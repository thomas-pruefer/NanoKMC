# Architecture

## Public solver architecture

```text
SystemClass
  |
  +-- SystemKMCHomogenous
        |
        +-- SystemKMCClassical
        |
        +-- SystemKMCActiveFilteredBase
        |     |
        |     +-- SystemKMCActiveFilteredGeneric
        |     +-- SystemKMCActiveFilteredBinaryNN
        |
        +-- internal binary-rate base
        |     +-- SystemKMCRateCategoryOptimized
        |     +-- internal exact-class core
        |           +-- SystemKMCExactClassOptimized
        |
        +-- internal partial-filter core
              +-- SystemKMCPartialFilterOptimized
```

`SystemKMCHomogenous` supplies the shared homogeneous benchmark parameters and incremental common-MCS interface. The six public `SystemID` values are the leaf solvers listed in the README.

## Active-Filtered structural engine

`SystemKMCActiveFilteredBase` owns the part that defines the NanoKMC method:

1. Fixed FCC site topology is built once.
2. The complete current set of structurally active unlike NN bonds is held as a dense canonical list.
3. A bond is selected uniformly from that list.
4. A derived backend returns `PairAcceptanceProbability(...)`.
5. If accepted, the existing NanoKMC lattice transaction performs the species exchange.
6. Only the union of candidate bonds incident on the two exchanged sites is refreshed.
7. The common clock uses `6/B` for every selected active bond.

No derived acceptance backend can bias structural selection through this interface. The state-dependent `6/B` clock is part of the physical correspondence of filtered selection; see the Shida–Henriques correctness boundary summarized in [`literature.md`](literature.md).

## Precomputed FCC topology

For each fixed lattice site `q`, the engine stores the 12 NN site indices. Runtime neighbor lookup therefore becomes a contiguous integer lookup rather than repeated periodic-coordinate arithmetic.

The optimization principle is not specific to FCC. The current release simply instantiates the FCC topology used by the paper.

## Differential active-set refresh

For an FCC NN exchange only bonds incident on either endpoint can change structural activity. The union contains 23 unique undirected bonds (`2*z - 1` with `z=12`). The engine checks only those canonical keys and modifies list membership only when active/inactive state actually changes.

The number 23 is FCC-specific; differential local refresh is not.

## Acceptance backends

### Generic

`SystemKMCActiveFilteredGeneric` gathers an ordered `LocalPairEnvironment` after a bond is selected. The structure contains endpoint species plus all directional NN species around both endpoints. The bundled evaluator applies the paper's symmetric homogeneous NN Hamiltonian, but `AcceptanceProbabilityFromEnvironment(...)` is virtual so a future system can use direction-dependent, motif-dependent or otherwise richer local physics without changing structural filtering.

### BinaryNN

`SystemKMCActiveFilteredBinaryNN` assumes exactly two species and the symmetric homogeneous NN Hamiltonian used by the paper. It stores one byte per fixed site: the number of species-1 nearest neighbours. That count is updated locally after an accepted exchange. The selected pair's exact Metropolis factor is then obtained through integer arithmetic and an eight-entry precomputed probability table.

The BinaryNN descriptor is used only after structural selection. It is not an energetic preselection catalogue.

## Classical baseline

`SystemKMCClassical` implements blind Kawasaki proposals without an active-event catalogue: select a site uniformly, select one of its 12 FCC neighbours uniformly, reject equal-species pairs, then apply the same Metropolis rule to unlike pairs. Its common clock is `1/N` MCS per proposal. The implementation computes the binary nearest-neighbour energy class directly from local coordination counts.

## Coarse rate-category reference architecture

`SystemKMCRateCategoryOptimized` shares only the fixed-topology/eight-rate utility base with the exact-class solver and occupies the rate-informed intermediate region between structural-only Active-Filtered selection and exact finite-rate classes. It stores each active bond only in one **coarse** category. Categories are selected with weight `N_g * p_hat_g`, a member is selected uniformly, and the selected bond is accepted with `p_exact / p_hat_g`. The default four-category partition of the eight binary FCC Metropolis classes is `{0,1}`, `{2,3}`, `{4,5}`, `{6,7}`.

The common-MCS increment is `6/Qhat`, where `Qhat = sum_g N_g p_hat_g`. `RateCategoryCount=1,2,4,8` is supported for controlled validation; the release/paper setting is `4`.

## Provenance boundary

Dense lists, reverse addressing, local list maintenance, partial filtering and finite rate classes are established algorithmic ideas. The public architecture intentionally makes these boundaries visible: the NanoKMC method stores only the complete **structural** unlike-bond population and evaluates energetics after selection, while the literature/reference solvers move less or more information into the maintained preselection state. See [`literature.md`](literature.md).
