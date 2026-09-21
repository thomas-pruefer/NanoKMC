# Extending NanoKMC

This release intentionally separates structural Active-Filtered event management from local energetic evaluation.

For a new local Hamiltonian, the least invasive extension path is to derive from `SystemKMCActiveFilteredGeneric` and override:

```cpp
double AcceptanceProbabilityFromEnvironment(
    const LocalPairEnvironment& env) override;
```

`LocalPairEnvironment` contains the two selected sites/species, selected bond direction, and ordered nearest-neighbour site/species arrays around both endpoints. A derived evaluator can therefore distinguish arrangements having the same composition but different directional patterns.

If a model admits a compact sufficient local-state descriptor, a new optimized backend can instead derive directly from `SystemKMCActiveFilteredBase` and implement:

```cpp
double PairAcceptanceProbability(...);
void InitializeAcceptanceState();
void UpdateAcceptanceStateAfterExchange(...);
```

The BinaryNN backend is the example: its descriptor is one species-1 NN count per site.

Future lattice support should replace/encapsulate the current FCC topology functions (`NeighborCoords`, canonical directions and affected-bond enumeration) rather than duplicating the Active-Filtered solver logic.
