# Reproducibility guide

This file defines the reproducibility boundary of the NanoKMC `0.1.0` software release.

## 1. What this repository reproduces

This repository contains six in-code solver implementations for the bundled binary FCC exchange model:

1. Classical blind proposal
2. Active-Filtered Generic
3. Active-Filtered BinaryNN
4. Partial-Filter reference
5. Rate-Category reference
6. Exact-Class reference

It contains the physical model, examples, native output definitions, validation utilities, and internal comparison helpers needed to inspect those implementations. The code base is also structured so the Generic Active-Filtered acceptance model and initializer can be extended beyond the bundled binary case.

The broader cross-software manuscript benchmark campaign is intentionally maintained outside this repository.

## 2. Build modes

### Normal software build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

This uses the toolchain's normal Release optimization policy.

### Frozen manuscript compiler policy

```bash
cmake -S . -B build-paper \
  -DCMAKE_BUILD_TYPE=Release \
  -DNANOKMC_PAPER_BUILD=ON
cmake --build build-paper
```

For GCC/Clang this explicitly selects `-O1` and `NDEBUG`. Performance values should be compared only with compiler/toolchain and host metadata recorded for the corresponding benchmark campaign.

## 3. Randomness and checkpoints

NanoKMC uses `std::mt19937_64` with an explicit `uint64_t` seed. The accepted input range is the full `0 .. 2^64-1` range. The code uses custom rejection sampling for integer indices and a fixed 53-bit conversion for unit-uniform doubles rather than implementation-defined standard-library distributions.

Within one process, re-created solver objects retain the process-wide engine stream when the seed is unchanged. Packed checkpoint archives, however, contain the lattice state only and do **not** serialize the engine state. Starting a new process from a saved checkpoint therefore restores the physical configuration but reinitializes the engine from `Seed`; exact uninterrupted-trajectory identity is not promised across process restarts.

Same-seed identity is used as a strong regression test when two code paths are expected to consume randomness in the same order. It is not a universal requirement for scientifically equivalent solver implementations.

## 4. Validation

After every build intended for scientific use, run:

```bash
ctest --test-dir build --output-on-failure
```

The registered tests cover:

- CLI availability;
- all six public solvers;
- structural unlike-bond recounting;
- species conservation;
- Active-Filtered Generic/BinaryNN same-seed non-timing equivalence;
- exact-class rejection-free execution;
- rate-category limiting cases;
- native output generation and profile normalization;
- configurable profile axes;
- full `uint64_t` seed parsing;
- exact input-key matching;
- blank-line-safe checkpoint parsing and nondecreasing checkpoint schedules;
- safe rejection of invalid packed-lattice dimensions.

The release is additionally checked with AddressSanitizer and UndefinedBehaviorSanitizer in CI.

## 5. Timing boundary

NanoKMC's benchmark diagnostic reports `EvolutionWallSecondsCumulative`, measured with `std::chrono::steady_clock` around the state-evolution loop. Initialization, checkpoint writing, `RunEval`, morphology analysis, and rendering are outside that timer.

The process-time field retained in `CalcData.csv` is compatibility metadata and is not the clean benchmark timer.

Parallel campaign launchers may be used for trajectory generation or throughput, but performance timing should use single-thread solver instances under controlled execution conditions.

## 6. Scientific invariants

A refactor must not silently change:

- solver proposal/selection semantics;
- Active-Filtered structural eligibility;
- the distinction between structural preselection and post-selection energetics;
- common-MCS mappings;
- the bundled Hamiltonian when running the validated examples;
- output/evaluation non-interference with subsequent KMC evolution.

If a change intentionally alters one of these, it is a scientific method/model change and should be documented and revalidated as such.

## 7. Archiving a release

For a publication-associated release:

1. run the full CTest suite on the intended compiler/platform;
2. record compiler, version, build flags, operating system, and CPU for performance data;
3. ensure `CMakeLists.txt`, `CHANGELOG.md`, and `CITATION.cff` match the intended tag/version/date;
4. create an immutable Git tag and GitHub Release;
5. archive that release with a persistent research-software repository such as Zenodo;
6. cite the **version-specific** archived software DOI in work that depends on that exact implementation;
7. after archival, add the DOI/badge to the default branch's README/CITATION metadata as a follow-up commit if desired.
