# Solver reference

## `KMCClassical`

**Family:** blind/random proposal Kawasaki KMC.

**Proposal population:** all oriented site-neighbour proposals. A lattice site is selected uniformly, followed by one of its 12 FCC nearest neighbours.

**Structural nulls:** same-species pairs are rejected immediately.

**Post-selection:** unlike pairs use the same binary symmetric NN Metropolis rule as the other benchmark solvers.

**Clock:** `1/N` common MCS per proposal, so one MCS is `N` random site-neighbour attempts.

The implementation evaluates the local binary NN energy class directly from neighbour counts. This is algorithmically classical blind proposal.

## `KMCActiveFilteredGeneric`

**Family:** NanoKMC Active-Filtered.

**Candidate population:** complete current set of unlike FCC nearest-neighbour bonds.

**Selection:** uniform over active bonds.

**Post-selection:** reconstruct ordered local NN environment, then Metropolis acceptance.

**Clock:** `6/B` common MCS per selected active bond.

**Bundled Hamiltonian:** symmetric homogeneous NN model, equal-species bond energy `-Ea/2`, unlike-species bond energy `0`. The evaluator and structural interface can represent multiple species, although the bundled homogeneous initializer currently populates only species 0/1. The ordered environment interface can be overridden for richer local Hamiltonians.

## `KMCActiveFilteredBinaryNN`

**Family:** same NanoKMC Active-Filtered method.

**Candidate population / selection / clock:** identical to Generic.

**Post-selection:** exact O(1) lookup from a cached binary NN coordination descriptor.

**Scope:** exactly two species, symmetric homogeneous NN model. The implementation uses eight precomputed Metropolis factors for the FCC exchange benchmark.

This is the performance-oriented binary specialization of the Active-Filtered method.

## `KMCPartialFilterOptimized`

**Family:** intermediate partial structural filtering.

**Literature anchor:** Bortz, Kalos & Lebowitz (1975), especially the Appendix-A binary-alloy procedure. See [`literature.md`](literature.md).

The solver maintains species-0 sites that have at least one species-1 NN. It selects an eligible species-0 site and then one of its 12 neighbors. Same-species neighbor picks are structural nulls; unlike picks retain a Metropolis acceptance test. Cached local counts accelerate eligibility and energy evaluation.

This solver is included as a controlled in-code comparison architecture, not as another NanoKMC Active-Filtered backend.

## `KMCRateCategoryOptimized`

**Family:** rate-informed intermediate / coarse rate-category rejection.

**Literature anchors:** Schulze's inverted-list rejection architecture and the Saum–Schulze–Ratsch rate-category study. See [`literature.md`](literature.md).

The bundled binary FCC NN benchmark has eight exact Metropolis classes. The release setting `RateCategoryCount=4` groups these as `{0,1}`, `{2,3}`, `{4,5}`, and `{6,7}`. Each category stores only active bonds whose exact rate lies in that group and is assigned the upper probability `p_hat_g` of the group.

A category is selected with weight `N_g p_hat_g`, a bond is selected uniformly inside it, and only then is the selected bond's exact class evaluated. The bond executes with residual probability `p_exact/p_hat_g`. Thus the solver maintains more energetic preselection information than NanoKMC Active-Filtered but less than `KMCExactClassOptimized`, and it deliberately retains rejection.

The common-MCS increment is `6/Qhat`, with `Qhat = sum_g N_g p_hat_g`. Supported `RateCategoryCount` values are `1`, `2`, `4`, and `8`; `4` is the default validated comparison setting. `M=1` is the single-majorant limit, while `M=8` becomes the exact-class rejection-free limit for the bundled eight-class model.

This is a controlled literature/reference architecture, not a NanoKMC method.

## `KMCExactClassOptimized`

**Family:** exact finite-rate-class / rejection-free event selection.

**Literature anchors:** BKL/n-fold, Sadiq's active-bond exact-class construction, and finite-rate-class/list implementations such as Schulze (2002). See [`literature.md`](literature.md).

Current unlike bonds are assigned to one of eight exact Metropolis-rate classes for the symmetric binary FCC benchmark. A class is selected with weight `N_c p_c`, then a bond uniformly within the class. The selected event executes without a subsequent Metropolis rejection. The mean common-MCS increment is `6/Q` with `Q = sum_c N_c p_c`.

Cached local counts and differential class refresh reduce the bookkeeping cost of the reference rate-class formulation.

## Scientific naming policy

The words `Optimized` in the literature/reference solvers describe implementation optimization inside this common benchmark code base. They do **not** imply that the underlying literature algorithm families were invented by NanoKMC. Likewise, `Generic` and `BinaryNN` are two acceptance backends of one Active-Filtered structural method, not two different KMC families.
