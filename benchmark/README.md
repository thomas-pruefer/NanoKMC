# Internal solver benchmark helper

`run_internal_comparison.py` runs the six public NanoKMC solver paths sequentially for one common binary FCC state and writes results under `benchmark/results/`.

Example:

```bash
python benchmark/run_internal_comparison.py \
  --exe build/nanokmc \
  --k 6 --mcs 30000 --xa 0.20 --kt 0.75 --seed 12345
```

The helper reports the clean solver-evolution timing from `Benchmark.csv` and checks same-seed non-timing equivalence of the Active-Filtered Generic and BinaryNN backends for the bundled symmetric binary Hamiltonian.

This is intentionally **not** the full manuscript cross-software benchmark harness for SPPARKS, kmcos, or KMC_Lattice. Final publication timing should follow the manuscript's controlled timing protocol.
