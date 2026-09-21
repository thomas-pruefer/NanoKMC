# Known limitations and deliberate scope

NanoKMC `0.1.0` is a focused scientific-software release rather than a general-purpose KMC engine. The following limitations are explicit so users can distinguish validated scope from future development directions.

## Scientific scope

- The bundled topology is the periodic FCC lattice used by the paper benchmark.
- The optimized BinaryNN, Partial-Filter, Rate-Category, and Exact-Class paths assume the bundled binary nearest-neighbour benchmark model.
- The Generic Active-Filtered backend exposes a richer local-environment acceptance interface, but adding a different event topology or clock mapping requires a new scientific construction; it is not a configuration-only change.
- Solver instances are single-threaded by design in this release.

## Runtime / packaging dependencies

- Packed checkpoint archives currently use command-line `zip` / `unzip` tools. They must be available on `PATH`.
- The executable is the supported interface; NanoKMC `0.1.0` does not install a stable C++ library/API target.
- The key-value `nanokmc.in` parser is intentionally compatible with the existing research input format. It is documented and validated for the public parameters but is not a general schema-driven configuration system.

## Implementation debt retained deliberately

The lattice core uses compact raw-array storage and some C-style I/O/path construction. This code is covered by the release regression and sanitizer tests, but it has not been converted wholesale to higher-level containers because such a rewrite would expand the scientific regression surface without changing the published method. Some shell-assisted checkpoint/run-folder handling also remains and is part of the current packaging boundary.

Internal invariant failures in several solver implementations terminate with a non-zero process exit rather than propagating C++ exceptions. These checks are intended to make corrupted internal state fail loudly during research/validation runs.

Future releases may introduce stronger RAII containers, structured error types, an input schema, and a reusable library target. Such changes should be regression-tested against the `0.1.0` scientific reference behavior.

## Reproducibility boundary

The complete cross-software benchmark campaign and manuscript figure pipeline are maintained separately. This repository contains the six in-code solver paths and their internal comparison helpers, not the external packages or their full campaign orchestration.
