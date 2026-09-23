# Contributing to NanoKMC

NanoKMC `0.1.0` is a release-stable scientific baseline. Changes should preserve
scientific correspondence, documented interfaces, and reproducibility before
optimizing convenience or performance.

## Before submitting a change

For normal development, build and run the full registered test suite:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

When a change may affect manuscript benchmark timings, additionally test the
frozen benchmark compiler policy:

```bash
cmake -S . -B build-paper -DCMAKE_BUILD_TYPE=Release -DNANOKMC_PAPER_BUILD=ON
cmake --build build-paper
ctest --test-dir build-paper --output-on-failure
```

On Windows/MSYS2, `scripts\windows\build_and_test.bat` builds/tests the normal
release and `scripts\windows\build_and_test.bat paper` selects the frozen
benchmark compiler policy.

## Scientific invariants

Changes to the six public solvers should preserve their documented proposal/
selection semantics, acceptance rule, and common-MCS mapping. In particular:

- `KMCClassical` must remain a blind random site-neighbour proposal baseline;
- the two Active-Filtered backends must share the same structural candidate
  population and selection/clock semantics;
- Active-Filtered structural selection must not use energetic information;
- `KMCPartialFilterOptimized`, `KMCRateCategoryOptimized`, and
  `KMCExactClassOptimized` are literature/reference architectures and should
  not be reframed as NanoKMC inventions;
- output/evaluation code must not mutate the subsequent KMC trajectory.

If a change intentionally alters one of these invariants, document it as a
method/model change rather than a refactor.

## Performance changes

Do not infer solver speed rankings from parallel batch wall time. Parallel
launchers are useful for correctness and trajectory generation, but controlled
performance timing should use sequential single-core runs on a documented host.

A performance optimization should be accompanied by a correctness regression
check. When deterministic list/RNG ordering is unchanged, same-seed comparison
is useful; otherwise statistical/observable correspondence is the scientific
criterion.

## Literature and comments

Comments should explain non-obvious algorithmic invariants, data-structure
intent, extension boundaries, and scientific provenance rather than restating
C++ syntax. New claims about algorithmic novelty should be checked against
`docs/literature.md` and the associated literature.

## Public interface

The stable `0.1.0` interface is the six public `SystemID` values plus the input
and output features documented under `docs/`. Internal or compatibility hooks
visible in the source are not automatically part of the supported interface.
