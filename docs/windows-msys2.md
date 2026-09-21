# Windows / MSYS2

NanoKMC supports a native Windows workflow through MSYS2 UCRT64.

The optional convenience helper searches for MSYS2 under `C:\msys64` and `D:\msys64`, then configures, builds, and runs the CTest suite with CMake/Ninja.

Standard Release build:

```text
scripts\windows\build_and_test.bat
```

Paper compiler-policy build (`NANOKMC_PAPER_BUILD=ON`):

```text
scripts\windows\build_and_test.bat paper
```

The same operations can always be run directly from an MSYS2 UCRT64 shell:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Benchmark and validation helpers are Python scripts and are documented in `benchmark/README.md` and [`validation.md`](validation.md); separate batch wrappers are intentionally not maintained.
