# Benchmark physical model

The bundled paper examples use a periodic FCC lattice with nearest-neighbour Kawasaki exchange.

- coordination number: `z = 12`
- conserved composition
- active structural event: unlike nearest-neighbour pair
- binary paper benchmark: species A/B
- homogeneous symmetric NN interaction scale: `Ea` in the source input
- thermal parameter: `kT`
- Metropolis acceptance: `min(1, exp(-DeltaE/kT))`

For the symmetric model used by the paper, an equal-species nearest-neighbour bond contributes `-Ea/2` and an unlike bond contributes `0`. This convention defines the bundled homogeneous benchmark Hamiltonian used throughout the release.

For a selected binary A-B FCC pair, the four shared neighbors cancel in the exchange-energy difference and only seven non-canceling external positions remain at each endpoint. The unfavorable Metropolis exponents reduce to the finite set `1 ... 7`, giving the eight factors

```text
1, exp(-Ea/kT), ..., exp(-7*Ea/kT)
```

used by `KMCActiveFilteredBinaryNN`, `KMCRateCategoryOptimized`, and `KMCExactClassOptimized`.

The method-level Active-Filtered selection does not depend on this finite spectrum; only the BinaryNN acceptance backend does.


`KMCRateCategoryOptimized` groups the eight exact factors into coarse contiguous categories. With the frozen `RateCategoryCount=4`, the groups are `{0,1}`, `{2,3}`, `{4,5}`, and `{6,7}`. Each group uses its largest probability as the category majorant; the exact selected-bond probability is applied only in the residual acceptance test.
