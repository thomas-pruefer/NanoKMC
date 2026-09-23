# NanoKMC quick start

This guide builds NanoKMC, runs the release tests, and starts a small Active-Filtered FCC exchange example.

## 1. Build

### Windows / MSYS2 UCRT64

From a normal Windows terminal in the repository root:

```text
scripts\windows\build_and_test.bat
```

This creates `build\nanokmc.exe` and runs the CTest suite.

For the frozen compiler policy used in the accompanying manuscript benchmark:

```text
scripts\windows\build_and_test.bat paper
```

which creates `build-paper\nanokmc.exe` with `NANOKMC_PAPER_BUILD=ON`.

### Linux / macOS-style toolchain

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Manuscript benchmark-policy build:

```bash
cmake -S . -B build-paper -DCMAKE_BUILD_TYPE=Release -DNANOKMC_PAPER_BUILD=ON
cmake --build build-paper
ctest --test-dir build-paper --output-on-failure
```

## 2. Inspect the executable

```bash
./build/nanokmc --version
./build/nanokmc --list-solvers
./build/nanokmc --help
```

## 3. Run the primary Active-Filtered example

Windows:

```text
build\nanokmc.exe examples\active_filtered_binary_nn
```

Linux/macOS-style shell:

```bash
./build/nanokmc examples/active_filtered_binary_nn
```

The example uses `KMCActiveFilteredBinaryNN` with the bundled binary FCC nearest-neighbour model.

## 4. Inspect and reset the example

Important outputs include:

```text
examples/active_filtered_binary_nn/
  output/CalcData.csv
  output/bit/data.zip
  evaluation/Benchmark.csv
```

`CalcData.csv` is updated as checkpoints complete. To reset an example, restore `output/CalcData.csv` from `output/CalcData.template.csv` (or restore the tracked file with Git).

## 5. Other solvers and outputs

Ready-to-run examples for all six public solver paths are under `examples/`. The `examples/output_features` case demonstrates the documented morphology, profile, visualization, and checkpoint outputs.

See:

- [Solver reference](solvers.md)
- [Physics model](physics-model.md)
- [Input/output reference](input-output.md)
- [Output features](output-features.md)
- [Validation](validation.md)
- [Reproducibility](reproducibility.md)
