# Reproducibility guide

This file defines the reproducibility boundary of the NanoKMC `0.1.0` scientific software release.

## 1. What this repository reproduces

This repository contains the six in-code solver implementations used to study the binary FCC Kawasaki benchmark:

1. Classical blind proposal
2. Active-Filtered Generic
3. Active-Filtered BinaryNN
4. Partial-Filter reference
5. Rate-Category reference
6. Exact-Class reference

It contains the benchmark physical model, examples, native output definitions, and validation utilities needed to inspect those implementations.

The broader cross-software benchmark campaign (SPPARKS, kmcos, KMC_Lattice, campaign orchestration, and manuscript figure generation) is intentionally maintained outside this focused public code release.

## 2. Build modes

### Normal software build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

This uses the toolchain's normal Release optimization policy.

### Frozen paper compiler policy

```bash
cmake -S . -B build-paper \
  -DCMAKE_BUILD_TYPE=Release \
  -DNANOKMC_PAPER_BUILD=ON
cmake --build build-paper
```

For GCC/Clang this explicitly selects `-O1` and `NDEBUG`. Performance values should be compared with the compiler/toolchain and host metadata recorded for the corresponding benchmark campaign.

## 3. Randomness

NanoKMC uses `std::mt19937_64` with an explicit integer seed. The code uses custom rejection sampling for integer indices and a fixed 53-bit conversion for unit-uniform doubles rather than implementation-defined standard-library distributions.

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
- configurable profile axes.

The release is additionally validated with AddressSanitizer, UndefinedBehaviorSanitizer, and leak detection.

## 5. Timing boundary

NanoKMC's benchmark diagnostic reports `EvolutionWallSecondsCumulative`, measured with `std::chrono::steady_clock` around the state-evolution loop. Initialization, checkpoint writing, `RunEval`, morphology analysis, and rendering are outside that timer.

The process-time field retained in `CalcData.csv` is compatibility metadata and is not the clean benchmark timer.

Parallel campaign launchers may be used for trajectory generation or throughput, but paper-style performance timing should use single-thread solver instances under controlled execution conditions.

## 6. Scientific invariants

A refactor must not silently change:

- solver proposal/selection semantics;
- Active-Filtered structural eligibility;
- the distinction between structural preselection and post-selection energetics;
- common-MCS mappings;
- the bundled Hamiltonian;
- output/evaluation non-interference with subsequent KMC evolution.

If a change intentionally alters one of these, it is a scientific method change and should be documented and revalidated as such.

## 7. Archiving a release

For a paper-associated release:

1. run the full CTest suite on the intended compiler/platform;
2. record compiler, version, build flags, operating system, and CPU in the benchmark metadata;
3. ensure `CITATION.cff` matches the tag;
4. create an immutable Git tag;
5. archive that tag with a persistent research-software repository (for example Zenodo) and add the resulting DOI to the citation metadata and manuscript.
