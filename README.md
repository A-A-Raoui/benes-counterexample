## A Finite Counterexample to the Beneš Conjecture - Verification Code

This repository contains a verification script for the paper:

> Raoui, A. A. (2026). A finite counterexample to the Beneš conjecture. *Zenodo*. [https://doi.org/10.5281/zenodo.21287099](https://doi.org/10.5281/zenodo.21287099)

### About

The [Beneš conjecture](https://www.openproblemgarden.org/op/bene_conjecture) (1975) asks whether transitivity of a product of partition-stabilizer cosets implies that a larger power equals the full symmetric group. Equivalently, in the [graph-theoretic formulation](https://www.openproblemgarden.org/op/bene_conjecture_graph_theoretic_form_0), it asks whether external connectivity of a repeated 2-stage graph implies rearrangeability.

We give an explicit finite counterexample to this conjecture. The example has a ground set of size |E| = 14, a uniform partition with seven blocks of size 2, and a permutation φ satisfying U ∧ φ(U) = 0. The set (φ S(U))^6 is transitive, but (φ S(U))^11 ≠ Sym(E). 

Equivalently, in the graph-theoretic formulation, the associated simple 2-regular ordered 2-stage graph L has L^5 externally connected, while L^10 is not rearrangeable. Non-rearrangeability is certified by a cut of capacity 5 in a 10-step time-expanded network, while the corresponding demand requires 6 edge-disjoint crossings. The proof is self-contained and requires no computer assistance - only hand-checkable finite arithmetic.

### Verification Script

The script [`verify_benes_counterexample.py`](verify_benes_counterexample.py) independently recomputes every finite quantity displayed in the paper:

- The partitions U and phi(U)
- The meet condition U ∧ phi(U) = **0**
- The adjacency matrix A and its fifth power A^5 (all entries positive)
- The demand matrix M for the permutation sigma
- The 10-layer time-expanded cut certificate (capacity 5 < demand 6)
- The equivalent graph-theoretic form: L^5 externally connected, L^10 not rearrangeable
- An independent max-flow computation confirming the cut is tight

### Requirements

Python 3.6 or later. **No external packages required** - the script uses only the standard library.

### Usage

```bash
python3 verify_benes_counterexample.py
```

Expected output:

```
All checks passed.
Verified: phi(U_i)=V_i, meet condition, A^5 positivity, demand matrix M,
          partition-form obstruction sigma not in P^11,
          graph-form obstruction L^5 externally connected but L^10 not rearrangeable,
          10-layer cut capacity 5 < selected demand 6, and max-flow value 5.
```

## License

This verification code is released under the [MIT License](LICENSE).
