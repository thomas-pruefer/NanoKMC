# Known limitations and deliberate scope

NanoKMC `0.1.0` is reusable FCC exchange-KMC research software, but it is not yet a general-purpose KMC framework. The following boundaries are explicit so users can distinguish implemented/validated capability from extension directions.

## Scientific scope

- The bundled topology is a periodic FCC lattice with nearest-neighbour exchange.
- The supplied homogeneous initializer and release validation cases are binary.
- The packed lattice core and Generic Active-Filtered local-environment interface are designed for `2..16` species, but true multi-species initialization/physics must be supplied by an extension.
- The optimized BinaryNN, Partial-Filter, Rate-Category, and Exact-Class paths assume the bundled binary nearest-neighbour model; `KMCClassical` is also binary in this release.
- The Generic Active-Filtered backend exposes a richer local-environment acceptance interface, but adding a different event topology or clock mapping requires a new scientific construction; it is not a configuration-only change.
- Solver instances are single-threaded by design in this release.

## Checkpoint boundary

Packed checkpoints serialize the lattice configuration but not the `std::mt19937_64` engine state. A checkpoint reopened in a new process is therefore a restart from the same physical state with the RNG initialized from `Seed`, not a bit-for-bit continuation of the uninterrupted random stream.

## Runtime / packaging dependencies

- Packed checkpoint archives currently use command-line `zip` / `unzip` tools. They must be available on `PATH`.
- The executable is the supported interface; NanoKMC `0.1.0` does not install a stable C++ library/API target.
- The key-value `nanokmc.in` parser remains compatible with the existing research input format. It now matches parameter names exactly, but it is not a general schema-driven configuration system.
- The packed FCC representation requires `knx`, `kny`, `knz >= 2`; very large lattices are additionally limited by memory and platform integer ranges.

## Implementation debt retained deliberately

The lattice core uses compact raw-array storage and some C-style I/O/path construction. This code is covered by release regression and sanitizer tests, but it has not been converted wholesale to higher-level containers because such a rewrite would expand the scientific regression surface without changing the algorithms.

Some shell-assisted checkpoint/run-folder handling also remains and is part of the current packaging boundary. Internal invariant failures in several solver implementations terminate with a non-zero process exit rather than propagating structured C++ exceptions.

Future releases may introduce stronger RAII containers, structured error types, an input schema, a reusable library target, and a more general topology/event abstraction. Such changes should be regression-tested against the `0.1.0` reference behavior.

## Reproducibility boundary

The full cross-software benchmark campaign and manuscript figure pipeline are maintained separately. This repository contains the six in-code solver paths and their internal comparison helpers, not external KMC packages or the complete manuscript campaign orchestration.
