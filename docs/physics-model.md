# Bundled physical model

The validated examples supplied with NanoKMC use a periodic FCC lattice with nearest-neighbour conserved exchange (Kawasaki-type dynamics).

- coordination number: `z = 12`
- conserved composition
- structural active event: unlike nearest-neighbour pair
- bundled examples: binary species A/B
- homogeneous symmetric NN interaction scale: `Ea`
- thermal parameter: `kT`
- Metropolis acceptance: `min(1, exp(-DeltaE/kT))`

For the bundled symmetric model, an equal-species nearest-neighbour bond contributes `-Ea/2` and an unlike bond contributes `0`.

For a selected binary A-B FCC pair, the four shared neighbors cancel in the exchange-energy difference and only seven non-canceling external positions remain at each endpoint. The unfavorable Metropolis exponents reduce to the finite set `1 ... 7`, giving the eight factors

```text
1, exp(-Ea/kT), ..., exp(-7*Ea/kT)
```

used by `KMCActiveFilteredBinaryNN`, `KMCRateCategoryOptimized`, and `KMCExactClassOptimized`.

The method-level Active-Filtered structural selection does not depend on this finite binary spectrum. `KMCActiveFilteredGeneric` receives the selected pair plus ordered local species environments and can be extended to other local Hamiltonians, including multi-species models, while retaining the FCC structural event engine and its clock construction.

`KMCRateCategoryOptimized` groups the eight exact binary factors into coarse contiguous categories. With `RateCategoryCount=4`, the groups are `{0,1}`, `{2,3}`, `{4,5}`, and `{6,7}`. Each group uses its largest probability as the category majorant; the exact selected-bond probability is applied only in the residual acceptance test.
