# Extending NanoKMC

NanoKMC `0.1.0` is not a general reaction-KMC framework, but the Active-Filtered implementation deliberately separates structural event management from local energetic evaluation.

## New local Hamiltonians

For a new local Hamiltonian on the current FCC exchange topology, the least invasive extension path is to derive from `SystemKMCActiveFilteredGeneric` and override:

```cpp
double AcceptanceProbabilityFromEnvironment(
    const LocalPairEnvironment& env) override;
```

`LocalPairEnvironment` contains the two selected sites/species, selected bond direction, and ordered nearest-neighbour site/species arrays around both endpoints. A derived evaluator can therefore distinguish arrangements having the same composition but different directional patterns.

The structural Active-Filtered list continues to contain all unlike nearest-neighbour pairs; energetic information remains post-selection unless a new method is intentionally designed.

## Multi-species models

The packed FCC lattice representation supports `2..16` species labels, and the Generic Active-Filtered engine carries those labels without reducing them to a binary descriptor. The bundled homogeneous initializer is intentionally simple and currently assigns only species `0` and `1` using `clvl`.

A true multi-species model can therefore be implemented by providing a multi-species initialization/distribution rule and an appropriate `AcceptanceProbabilityFromEnvironment(...)`. The Generic solver architecture does not need to be replaced. The optimized BinaryNN, Partial-Filter, Rate-Category and Exact-Class implementations are binary specializations and would need separate generalized constructions if equivalent optimized multi-species variants are desired.

## New optimized acceptance backends

If a model admits a compact sufficient local-state descriptor, a new optimized backend can derive directly from `SystemKMCActiveFilteredBase` and implement:

```cpp
double PairAcceptanceProbability(...);
void InitializeAcceptanceState();
void UpdateAcceptanceStateAfterExchange(...);
```

The BinaryNN backend is the example: its descriptor is one species-1 NN count per site.

## New lattices or event topologies

Future lattice support should replace/encapsulate the current FCC topology functions (`NeighborCoords`, canonical directions and affected-bond enumeration) rather than duplicating the Active-Filtered solver logic. A different topology also requires re-deriving the common-clock mapping; it is not merely an input-file switch in `0.1.0`.

Likewise, non-exchange reactions, vacancies with distinct event types, or arbitrary reaction networks would require a broader event abstraction. Those are natural future framework directions but are outside the current executable interface.
