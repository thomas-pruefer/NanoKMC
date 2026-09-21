# Scientific provenance and literature

NanoKMC is a research implementation. The public code therefore distinguishes the specific NanoKMC Active-Filtered formulation from established algorithmic ingredients and neighbouring solver architectures.

The references below are **scientific precedents and comparison points**, not claims that the present source code was copied from the cited implementations.

## Core references

### Bortz, Kalos & Lebowitz (1975) — BKL / n-fold way and partial filtering

A. B. Bortz, M. H. Kalos, J. L. Lebowitz, "A new algorithm for Monte Carlo simulation of Ising spin systems," *Journal of Computational Physics* **17** (1975) 10–18. DOI: `10.1016/0021-9991(75)90060-1`.

Relevance to this repository:

- foundational rate-class / rejection-free Monte Carlo construction;
- dense class populations and reverse/cross-addressing are established prior art;
- Appendix A describes an intermediate binary-alloy procedure that maintains A sites having at least one B neighbour, then chooses a neighbour and retains residual rejection.

`KMCPartialFilterOptimized` is a controlled modern reference implementation of this **partial structural-filter design point** for the present FCC benchmark. It is not presented as a new NanoKMC method.

### Sadiq (1984) — active bonds with exact energetic classes

A. Sadiq, "A new algorithm for the Monte Carlo simulation of spin-exchange kinetics of Ising systems," *Journal of Computational Physics* **55** (1984) 387–396. DOI: `10.1016/0021-9991(84)90028-7`.

Relevance:

- maintains structurally active unlike exchange bonds;
- additionally assigns active bonds to exact energetic/rate classes;
- performs rate-weighted selection so the selected move executes without a subsequent Metropolis rejection;
- locally maintains both activity and energetic class information.

`KMCExactClassOptimized` occupies this neighbouring **active-bond + exact finite-rate-class** design point for the symmetric binary FCC benchmark. It is a controlled in-code reference architecture, not a claim of a new exact-class family.

### Shida & Henriques (1997) — filtered-pair normalization/correctness boundary

C. S. Shida, V. B. Henriques, "Kawasaki dynamics and equilibrium distributions in simulations of phase separating systems," arXiv:cond-mat/9703198 (1997).

Relevance:

Uniformly selecting only unlike nearest-neighbour pairs changes the state-dependent proposal normalization. A filtered-pair implementation therefore cannot simply count each selected pair as an ordinary classical proposal without compensating for the changing candidate population.

NanoKMC addresses this at the common-MCS level. For an FCC lattice with `B` current undirected unlike bonds, `KMCActiveFilteredGeneric` and `KMCActiveFilteredBinaryNN` advance by

```text
Delta MCS = 6 / B
```

for every selected active bond. This preserves the classical expected proposal frequency of each current undirected nearest-neighbour bond while leaving the selected-pair Metropolis test explicit.

### Schulze (2002) — exact classes and minimal searching

T. P. Schulze, "Kinetic Monte Carlo simulations with minimal searching," *Physical Review E* **65** (2002) 036704. DOI: `10.1103/PhysRevE.65.036704`.

Relevance:

- finite rate classes;
- list-based event storage;
- reverse-addressed local reclassification;
- avoiding global searches when the set of possible local rates is finite.

These are important prior-art boundaries for the data-structure and exact-class aspects of the repository.

### Schulze (2008) — inverted lists and rejection trade-offs

T. P. Schulze, "Efficient kinetic Monte Carlo simulation," *Journal of Computational Physics* **227** (2008) 2455–2462. DOI: `10.1016/j.jcp.2007.10.021`.

Relevance:

This work explicitly explores the trade between maintained preselection information and residual rejection using inverted-list approaches. That trade is central to interpreting the solver set in this repository.

### Saum, Schulze & Ratsch (2009) — coarse rate categories and residual rejection

M. A. Saum, T. P. Schulze, C. Ratsch, "Inverted List Kinetic Monte Carlo with Rejection Applied to Directed Self-Assembly of Epitaxial Growth," *Communications in Computational Physics* **6** (2009) 553–564.

Relevance:

- partitions accessible events into coarse rate categories with category upper bounds;
- selects a category using category population times its upper bound;
- selects a member uniformly within the category;
- retains residual rejection using exact rate divided by the category bound;
- varies category resolution and demonstrates the rejection-versus-bookkeeping trade-off directly.

`KMCRateCategoryOptimized` is a controlled same-physics realization of this **coarse rate-category rejection design point** for the present eight-class binary FCC benchmark. It is not a claim that the present source code reproduces the cited implementation line-for-line.

## How the six public solvers should be interpreted

| Solver | Scientific role |
|---|---|
| `KMCClassical` | Blind random site-neighbour Kawasaki/Metropolis baseline. |
| `KMCActiveFilteredGeneric` | NanoKMC Active-Filtered formulation with explicit local-environment evaluation after structural selection. |
| `KMCActiveFilteredBinaryNN` | Optimized binary-NN backend of the same Active-Filtered formulation. |
| `KMCPartialFilterOptimized` | In-code literature/reference architecture representing partial structural filtering with residual rejection. |
| `KMCRateCategoryOptimized` | In-code Schulze/Saum-style reference architecture representing coarse rate categories with residual rejection. |
| `KMCExactClassOptimized` | In-code literature/reference architecture representing active bonds with exact finite-rate classes and rejection-free selection. |

The repository does **not** claim that active lists, reverse addressing, swap-delete list maintenance, partial filtering, finite rate classes, or rejection-free selection are novel ideas. The NanoKMC contribution is the specific structural-only filtering formulation and its implementation/benchmark characterization.

## Accompanying NanoKMC paper

The repository's `CITATION.cff` intentionally leaves the paper DOI unset until the final bibliographic record exists. Once published, the release metadata and README should be updated with the permanent paper citation and archival software DOI if available.
