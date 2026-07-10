#!/usr/bin/env python3
"""Verify the finite data in "A finite counterexample to the Beneš conjecture".

Paper: Raoui, A. A. (2026). A finite counterexample to the Beneš conjecture. 
Zenodo. https://doi.org/10.5281/zenodo.21287099

This script is intentionally dependency-free.  It checks the finite arithmetic
used in the manuscript:

* the partitions U_i and V_i;
* the permutation phi and the meet condition U wedge phi(U) = 0;
* the directed graph Gamma, its adjacency matrix A, and the displayed A^5;
* the permutation sigma and its demand matrix M;
* the 10-layer cut certificate used to prove sigma not in P^11;
* the equivalent graph-theoretic form: L^5 externally connected but L^10 not rearrangeable;
* an optional max-flow check for the same cut demand.

The script is not needed for the proof: the paper gives the cut certificate
explicitly.  Its role is to make the finite verification reproducible.
"""

from collections import deque
from typing import Deque, Dict, Iterable, List, Sequence, Set, Tuple

Point = int
Vertex = int
Edge = Tuple[Vertex, Vertex]
Matrix = List[List[int]]

E: List[Point] = list(range(1, 15))
VERTICES: List[Vertex] = list(range(7))

U: Dict[Vertex, Set[Point]] = {
    0: {1, 9},
    1: {2, 3},
    2: {10, 11},
    3: {5, 13},
    4: {4, 7},
    5: {6, 14},
    6: {8, 12},
}

V_EXPECTED: Dict[Vertex, Set[Point]] = {
    0: {1, 2},
    1: {3, 4},
    2: {5, 6},
    3: {7, 8},
    4: {9, 10},
    5: {11, 12},
    6: {13, 14},
}

# Each point of E is a directed edge of Gamma.
EDGES: Dict[Point, Edge] = {
    1: (0, 0),
    2: (0, 1),
    3: (1, 1),
    4: (1, 4),
    5: (2, 3),
    6: (2, 5),
    7: (3, 4),
    8: (3, 6),
    9: (4, 0),
    10: (4, 2),
    11: (5, 2),
    12: (5, 6),
    13: (6, 3),
    14: (6, 5),
}

A_EXPECTED: Matrix = [
    [1, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 1, 0, 1],
    [1, 0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 1],
    [0, 0, 0, 1, 0, 1, 0],
]

A5_EXPECTED: Matrix = [
    [7, 9, 4, 2, 6, 2, 2],
    [6, 7, 4, 4, 5, 4, 2],
    [4, 2, 4, 9, 2, 9, 2],
    [4, 4, 5, 4, 5, 4, 6],
    [5, 6, 5, 2, 6, 2, 6],
    [2, 2, 6, 2, 6, 2, 12],
    [4, 2, 4, 9, 2, 9, 2],
]

M_EXPECTED: Matrix = [
    [0, 0, 2, 0, 0, 0, 0],
    [2, 0, 0, 0, 0, 0, 0],
    [0, 2, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 2],
    [0, 0, 0, 0, 0, 2, 0],
]

X: Dict[int, Set[Vertex]] = {
    0: {1, 2, 3, 6},
    1: {3, 4, 5},
    2: {2, 6},
    3: {3, 5},
    4: {2, 6},
    5: {3, 5},
    6: {2, 6},
    7: {3, 5},
    8: {2, 6},
    9: {3, 5},
    10: {2, 4, 6},
}

ENTERING_EXPECTED: Dict[int, List[Edge]] = {
    0: [],
    1: [],
    2: [],
    3: [(4, 2)],
    4: [],
    5: [(4, 2)],
    6: [],
    7: [(4, 2)],
    8: [],
    9: [(1, 4), (4, 2)],
}


def fail(message: str) -> None:
    raise AssertionError(message)


def cycle_perm(cycles: Iterable[Sequence[Point]]) -> Dict[Point, Point]:
    """Return the permutation represented by disjoint cycle notation."""
    permutation = {x: x for x in E}
    for cycle in cycles:
        cycle_list = list(cycle)
        if len(cycle_list) == 0:
            continue
        for a, b in zip(cycle_list, cycle_list[1:] + cycle_list[:1]):
            if a not in permutation or b not in permutation:
                fail(f"cycle contains a point outside E: {cycle_list}")
            permutation[a] = b
    if set(permutation.values()) != set(E):
        fail("cycle notation did not define a bijection of E")
    return permutation


def image(permutation: Dict[Point, Point], subset: Set[Point]) -> Set[Point]:
    return {permutation[x] for x in subset}


def matmul(left: Matrix, right: Matrix) -> Matrix:
    n = len(left)
    return [
        [sum(left[i][k] * right[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def matrix_power(matrix: Matrix, exponent: int) -> Matrix:
    if exponent < 1:
        fail("this script only needs positive matrix powers")
    result = matrix
    for _ in range(exponent - 1):
        result = matmul(result, matrix)
    return result


def row_sums(matrix: Matrix) -> List[int]:
    return [sum(row) for row in matrix]


def column_sums(matrix: Matrix) -> List[int]:
    return [sum(matrix[i][j] for i in range(len(matrix))) for j in range(len(matrix[0]))]


def adjacency_matrix() -> Matrix:
    matrix = [[0 for _ in VERTICES] for _ in VERTICES]
    for a, b in EDGES.values():
        matrix[a][b] += 1
    return matrix


def check_partitions_and_graph(phi: Dict[Point, Point]) -> None:
    all_u_points = set().union(*U.values())
    all_v_points = set().union(*V_EXPECTED.values())
    assert all_u_points == set(E), "U blocks do not cover E"
    assert all_v_points == set(E), "V blocks do not cover E"
    assert sum(len(block) for block in U.values()) == len(E), "U blocks overlap"
    assert sum(len(block) for block in V_EXPECTED.values()) == len(E), "V blocks overlap"
    assert {len(block) for block in U.values()} == {2}, "U is not uniform of block size 2"

    incoming = {v: {e for e, (_, b) in EDGES.items() if b == v} for v in VERTICES}
    outgoing = {v: {e for e, (a, _) in EDGES.items() if a == v} for v in VERTICES}
    assert incoming == U, f"incoming edge sets differ from U: {incoming}"
    assert outgoing == V_EXPECTED, f"outgoing edge sets differ from V: {outgoing}"

    indegrees = {v: 0 for v in VERTICES}
    outdegrees = {v: 0 for v in VERTICES}
    ordered_pairs = set()
    for edge in EDGES.values():
        a, b = edge
        outdegrees[a] += 1
        indegrees[b] += 1
        if edge in ordered_pairs:
            fail(f"Gamma has a repeated directed edge {edge}")
        ordered_pairs.add(edge)
    assert set(indegrees.values()) == {2}, f"indegrees are not all 2: {indegrees}"
    assert set(outdegrees.values()) == {2}, f"outdegrees are not all 2: {outdegrees}"

    phi_u = {v: image(phi, U[v]) for v in VERTICES}
    assert phi_u == V_EXPECTED, f"phi(U_i) differs from V_i: {phi_u}"

    meet_max = max(len(U[j] & V_EXPECTED[i]) for i in VERTICES for j in VERTICES)
    assert meet_max <= 1, f"meet condition fails; maximum intersection is {meet_max}"


def check_matrices(sigma: Dict[Point, Point]) -> Matrix:
    A = adjacency_matrix()
    assert A == A_EXPECTED, f"adjacency matrix differs from displayed A: {A}"
    A5 = matrix_power(A, 5)
    assert A5 == A5_EXPECTED, f"A^5 differs from displayed matrix: {A5}"
    assert min(min(row) for row in A5) > 0, "A^5 is not positive"

    # This is not used in the proof, but records that the exponent 6 is not an
    # artifact of using a longer walk than necessary for this graph.
    A4 = matrix_power(A, 4)
    assert any(A4[i][j] == 0 for i in VERTICES for j in VERTICES), "A^4 unexpectedly positive"

    M = [[len(image(sigma, U[i]) & V_EXPECTED[j]) for j in VERTICES] for i in VERTICES]
    assert M == M_EXPECTED, f"demand matrix differs from displayed M: {M}"
    assert row_sums(M) == [2] * 7, f"row sums of M are not all 2: {row_sums(M)}"
    assert column_sums(M) == [2] * 7, f"column sums of M are not all 2: {column_sums(M)}"
    return M


def entering_edges_by_layer() -> Dict[int, List[Edge]]:
    result: Dict[int, List[Edge]] = {}
    for t in range(10):
        outside = set(VERTICES) - X[t]
        entering = sorted({(a, b) for a, b in EDGES.values() if a in outside and b in X[t + 1]})
        result[t] = entering
    return result


def check_cut_certificate(M: Matrix) -> None:
    C = set(VERTICES) - X[0]
    T = X[10]
    assert C == {0, 4, 5}, f"wrong initial outside set C: {C}"
    assert T == {2, 4, 6}, f"wrong final target set T: {T}"

    demand_crossing_cut = sum(M[i][j] for i in C for j in T)
    assert demand_crossing_cut == 6, f"demand crossing cut should be 6, got {demand_crossing_cut}"

    entering = entering_edges_by_layer()
    assert entering == ENTERING_EXPECTED, f"entering layer edges differ: {entering}"
    cut_capacity = sum(len(edges) for edges in entering.values())
    assert cut_capacity == 5, f"cut capacity should be 5, got {cut_capacity}"
    assert demand_crossing_cut > cut_capacity, "cut does not obstruct the demand"


def check_graph_theoretic_form(M: Matrix) -> None:
    """Check the Open Problem Garden graph-theoretic formulation.

    The associated ordered 2-stage graph L has adjacency matrix A.  Thus L^m
    is externally connected exactly when every entry of A^m is positive.  The
    displayed mask has multiplicity matrix M.  The same 10-layer cut proves that
    this mask is not routable in L^10.
    """
    A = adjacency_matrix()
    A5 = matrix_power(A, 5)
    assert min(min(row) for row in A5) > 0, "L^5 is not externally connected"

    A4 = matrix_power(A, 4)
    assert any(A4[i][j] == 0 for i in VERTICES for j in VERTICES), (
        "L^4 unexpectedly externally connected"
    )

    # The mask for L^10 has source-target multiplicities M.  Since L is
    # 2-regular, every valid mask must have all source and target degrees 2.
    assert row_sums(M) == [2] * 7, "mask source degrees are not all 2"
    assert column_sums(M) == [2] * 7, "mask target degrees are not all 2"

    C = set(VERTICES) - X[0]
    T = X[10]
    mask_demand_across_cut = sum(M[i][j] for i in C for j in T)
    cut_capacity = sum(len(edges) for edges in entering_edges_by_layer().values())
    assert mask_demand_across_cut == 6, "mask demand across graph-form cut should be 6"
    assert cut_capacity == 5, "graph-form cut capacity should be 5"
    assert mask_demand_across_cut > cut_capacity, "mask is not obstructed by the cut"


class Dinic:
    """Small integer max-flow implementation for an independent cut check."""

    def __init__(self, number_of_nodes: int) -> None:
        self.graph: List[List[List[int]]] = [[] for _ in range(number_of_nodes)]

    def add_edge(self, u: int, v: int, capacity: int) -> None:
        forward = [v, capacity, len(self.graph[v])]
        backward = [u, 0, len(self.graph[u])]
        self.graph[u].append(forward)
        self.graph[v].append(backward)

    def max_flow(self, source: int, sink: int) -> int:
        flow = 0
        infinity = 10**9
        while True:
            level = [-1] * len(self.graph)
            level[source] = 0
            queue: Deque[int] = deque([source])
            while queue:
                u = queue.popleft()
                for v, capacity, _ in self.graph[u]:
                    if capacity > 0 and level[v] < 0:
                        level[v] = level[u] + 1
                        queue.append(v)
            if level[sink] < 0:
                return flow

            index = [0] * len(self.graph)

            def dfs(u: int, pushed: int) -> int:
                if u == sink:
                    return pushed
                while index[u] < len(self.graph[u]):
                    edge = self.graph[u][index[u]]
                    v, capacity, reverse_index = edge
                    if capacity > 0 and level[v] == level[u] + 1:
                        sent = dfs(v, min(pushed, capacity))
                        if sent > 0:
                            edge[1] -= sent
                            self.graph[v][reverse_index][1] += sent
                            return sent
                    index[u] += 1
                return 0

            while True:
                pushed = dfs(source, infinity)
                if pushed == 0:
                    break
                flow += pushed


def node(t: int, v: Vertex) -> int:
    return 7 * t + v


def check_max_flow_for_selected_demand(M: Matrix) -> None:
    """Independent max-flow check for the C-to-T sub-demand over 10 steps."""
    source = 77
    sink = 78
    network = Dinic(79)
    C = set(VERTICES) - X[0]
    T = X[10]

    for i in C:
        supply = sum(M[i][j] for j in T)
        network.add_edge(source, node(0, i), supply)
    for t in range(10):
        for a, b in EDGES.values():
            network.add_edge(node(t, a), node(t + 1, b), 1)
    for j in T:
        demand = sum(M[i][j] for i in C)
        network.add_edge(node(10, j), sink, demand)

    flow = network.max_flow(source, sink)
    assert flow == 5, f"selected C-to-T max flow should be 5, got {flow}"


def main() -> None:
    phi = cycle_perm([(2, 3, 4, 9), (5, 7, 10), (6, 11), (8, 13), (12, 14)])
    sigma = cycle_perm([(1, 5, 7, 10, 3, 2), (4, 9, 6, 13, 8, 11)])

    check_partitions_and_graph(phi)
    M = check_matrices(sigma)
    check_cut_certificate(M)
    check_graph_theoretic_form(M)
    check_max_flow_for_selected_demand(M)

    print("All checks passed.")
    print("Verified: phi(U_i)=V_i, meet condition, A^5 positivity, demand matrix M,")
    print("          partition-form obstruction sigma not in P^11,")
    print("          graph-form obstruction L^5 externally connected but L^10 not rearrangeable,")
    print("          10-layer cut capacity 5 < selected demand 6, and max-flow value 5.")


if __name__ == "__main__":
    main()
