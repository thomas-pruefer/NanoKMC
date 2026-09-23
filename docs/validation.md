# Validation and reproducibility

## Smoke validation

```bash
python validation/run_smoke.py --exe build/nanokmc
```

The smoke campaign checks all six public solvers and requires:

- successful completion;
- internal structural unlike-bond count equals an independent FCC recount at checkpoints;
- species conservation;
- `KMCClassical` structural bond count equals the independent FCC recount;
- exact-class selected events equal executed events;
- `KMCActiveFilteredGeneric` and `KMCActiveFilteredBinaryNN` give identical same-seed non-timing benchmark fields for the bundled symmetric binary Hamiltonian.

The last check is stronger than the scientific requirement (statistical equivalence) and is used as a release regression test. Future harmless changes to internal list ordering could remove same-seed trajectory identity without changing the stochastic method.

## Input/parser validation

```bash
python validation/run_input_validation.py --exe build/nanokmc
```

This regression test checks release-hardening behavior that is easy to break without affecting normal examples:

- `Seed` preserves the full `uint64_t` range without 32-bit truncation;
- parameter names are matched exactly (`T` is distinct from `kT`);
- blank lines in `CalcData.csv` are ignored rather than becoming extra zero checkpoints;
- invalid packed-FCC dimensions such as `knx=1` fail safely with a clear diagnostic;
- checkpoint MCS values must be non-negative and nondecreasing.

## Output-feature smoke validation

```bash
python validation/run_output_smoke.py --exe build/nanokmc
```

This runs the small `examples/output_features` case and requires non-empty benchmark, cluster-distribution, axial/cylindrical/spherical composition-profile, XZ projection, RasMol/XYZ, BlenderSimple/XYZ, coordinate-CSV and packed snapshot outputs. It additionally checks that per-bin species counts reproduce the bin population and the composition fractions sum to one.

Axis configurability is independently regression-tested with:

```bash
python validation/run_profile_axis_smoke.py --exe build/nanokmc
```

This test runs the axial profile on `Y` and the cylindrical profile on `Z`.

## Manuscript compiler setting

For GCC/Clang, `NANOKMC_PAPER_BUILD=ON` applies:

```text
-O1 -DNDEBUG
```

Use this mode only when reproducing benchmark timings tied to that compiler policy. General NanoKMC use can use the normal Release build.

## Internal comparison benchmark

```bash
python benchmark/run_internal_comparison.py --exe build/nanokmc --k 5
python benchmark/run_internal_comparison.py --exe build/nanokmc --k 6
```

These compare only the six solvers contained in this repository; external-software comparisons are outside this release.

## Rate-category architecture limits

`validation/run_rate_category_limits.py` checks both controlled ends of the category-resolution spectrum. With `RateCategoryCount=8`, every exact class is its own category and the solver must match the exact-class trajectory apart from diagnostic probability-evaluation counters and timing. With `RateCategoryCount=1`, all structurally active bonds share one majorant category and the solver must retain residual energetic rejection without structural null proposals.
