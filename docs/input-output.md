# Input and output reference

NanoKMC uses a folder-based run interface. The executable changes into the selected run directory and reads `nanokmc.in` plus `output/CalcData.csv`.

## Minimal run-directory layout

```text
run/
  nanokmc.in
  output/
    CalcData.csv
    CalcData.template.csv      # optional convenience copy
    bit/                       # created/used for packed snapshots
  evaluation/                  # generated checkpoint analyses
```

The `examples/` tree contains ready-to-run cases.

## `nanokmc.in`

The public release supports the following fields.

| Field | Meaning | Typical value / constraint |
|---|---|---|
| `SystemID` | public solver selection | one of the six IDs below |
| `Seed` | `std::mt19937_64` random seed | integer `0 .. 18446744073709551615` (`uint64_t`) |
| `NSpecies` | number of lattice species represented by the packed lattice | `2..16`; supplied runs are binary; see species note below |
| `knx`, `kny`, `knz` | FCC half-grid exponents used by the lattice builder | each must be `>=2`; practical maximum is memory/platform dependent |
| `lc` | lattice-coordinate scale used by visualization/evaluation helpers | e.g. `0.4338` |
| `kT` | thermal-energy scale | must be `>0` |
| `Ea` | homogeneous NN interaction parameter used by the bundled model | validated examples use `1.0` |
| `clvl` | bundled homogeneous initializer's species-0 fraction | `0 <= clvl <= 1` |
| `dev` | retained compatibility parameter of the homogeneous initializer | current bundled initializer does not use it beyond logging |
| `SpeciesName(i)` | output label for species `i` | e.g. `"A"`, `"B"` |
| `SpeciesColor(i)` | RasMol color command | e.g. `"[255,0,0]"` |
| `RasmolSpeciesPlot(i)` | include species in visualization/CSV outputs | `0` or `1` |
| `EvalParam` | space-delimited checkpoint output tokens | e.g. `" Benchmark AxialCompositionProfile "` |
| `SysEvalParam` | solver-specific validation tokens | normally left as supplied by examples |
| `AxialCompositionProfileAxis` | direction used by `AxialCompositionProfile` | `X`, `Y`, or `Z`; default `X` |
| `CylindricalCompositionProfileAxis` | cylinder axis used by `CylindricalCompositionProfile` | `X`, `Y`, or `Z`; default `X` |
| `RateCategoryCount` | coarse category count for `KMCRateCategoryOptimized` | `1`, `2`, `4`, or `8`; default/validated comparison value `4` |
| `Index`, `T`, `EB`, `dE`, `E01`, `f01`, `Fluence` | compatibility metadata retained by the folder/output workflow | keep the example values unless adapting the metadata workflow |

Supported public `SystemID` values:

```text
KMCClassical
KMCActiveFilteredGeneric
KMCActiveFilteredBinaryNN
KMCPartialFilterOptimized
KMCRateCategoryOptimized
KMCExactClassOptimized
```

The public examples also contain `SeedStart`, `SeedEnd` and `SeedStep` fields inherited from older batch workflows. The current executable seeds one invocation from `Seed`; the other three fields do not control the RNG of that invocation.

### Species note

The lattice storage supports `2..16` species labels. `KMCActiveFilteredGeneric` carries arbitrary species labels through the structural active-bond list and ordered local-environment interface. The bundled `SystemKMCHomogenous` initializer currently creates only species `0` and `1` from `clvl`, so a true multi-species simulation requires an extended initializer (or another compatible state source) and an acceptance model appropriate to that physics. The other five public solver implementations in `0.1.0` require exactly two species.

## `output/CalcData.csv`

`CalcData.csv` is the checkpoint schedule and progress record. The public examples ship a matching `CalcData.template.csv`; copy the template back before repeating a completed run.

For `NSpecies=2`, the current row layout is:

```text
requested_mcs ; bond_number ; record_count ; legacy_time_ms ; jumps_in_record ; atom_count_0 ; atom_count_1 ;
```

Important notes:

- column 0 is the requested common-MCS checkpoint; values must be non-negative and nondecreasing;
- column 1 being non-zero marks a previously completed checkpoint in the resume workflow;
- blank lines are ignored rather than interpreted as extra zero-valued checkpoint rows;
- the timing value in column 3 is retained for compatibility and is **not** the clean benchmark timer;
- solver-independent benchmark timing is written to `evaluation/Benchmark.csv`.

## `evaluation/Benchmark.csv`

Enable with:

```text
EvalParam=" Benchmark ";
```

The file contains, at each checkpoint:

- requested checkpoint and internal common-MCS coordinate;
- internally maintained unlike-bond count;
- independent FCC unlike-bond recount;
- interface density;
- proposal/selected and accepted/executed event counters;
- structural-null and probability-evaluation counters;
- Active-Filtered/list-update counters where relevant;
- clean cumulative KMC evolution wall time;
- total lattice size and species populations.

The clean evolution timer excludes initialization and checkpoint evaluation/output work outside the KMC evolution loop.

## Packed snapshots and restart semantics

Every requested checkpoint is written as a packed lattice bit-state. At run completion the checkpoint files are consolidated into:

```text
output/bit/data.zip
```

The archive stores the **lattice state only**. It does not serialize the `std::mt19937_64` engine state. Loading a checkpoint in a new process therefore starts from the saved physical configuration but reinitializes the RNG from `Seed`; the subsequent stochastic trajectory is not expected to be identical to an uninterrupted run. This is the intended `0.1.0` checkpoint boundary.

## Optional output tokens

The release-validated output tokens are grouped as follows.

Simulation / morphology observables:

```text
Benchmark
ClusterDistribution
AxialCompositionProfile
CylindricalCompositionProfile
SphericalCompositionProfile
ProjectedCompositionXZ
```

State / visualization exports:

```text
Rasmol
OnefileRasmol
BlenderSimple
WriteCSV
```

For example:

```text
EvalParam=" Benchmark ClusterDistribution AxialCompositionProfile CylindricalCompositionProfile SphericalCompositionProfile ProjectedCompositionXZ Rasmol OnefileRasmol BlenderSimple WriteCSV ";
AxialCompositionProfileAxis="X";
CylindricalCompositionProfileAxis="X";
```

See [`output-features.md`](output-features.md) for exact files, conventions and validation status.

## Output timing

Output/evaluation work can be much more expensive than the solver step for small systems. For performance comparisons, use the benchmark timer and keep output policy controlled. The included `benchmark/` helper compares only the six in-code solvers.
