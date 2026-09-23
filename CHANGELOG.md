# Changelog

## 0.1.0 — 2026-09-23

First public research-software release of NanoKMC, including the validated
configuration used for the accompanying manuscript study.

Highlights:

- Six public FCC exchange-KMC solver architectures: Classical, Active-Filtered
  Generic, Active-Filtered BinaryNN, Partial-Filter, Rate-Category, and
  Exact-Class.
- Extensible Generic Active-Filtered local-environment interface for richer
  local Hamiltonians and multi-species extensions.
- Full `uint64_t` random-seed parsing with deterministic custom sampling from
  `std::mt19937_64`.
- Hardened input validation for packed-FCC dimensions, species count, thermal
  parameter, composition fraction, exact parameter-name matching, and
  checkpoint schedules.
- Blank lines in `CalcData.csv` are ignored safely instead of creating implicit
  zero-valued checkpoints.
- Reproducible CMake/C++17 build with a separate frozen manuscript-benchmark
  compiler-policy option.
- Regression validation for solver invariants, rate-category limits, input
  parsing, and documented output features.
- Example inputs, internal benchmark helper, scientific provenance, citation
  metadata, GitHub CI, cross-compiler validation, and sanitizer validation.
- Cleanup of demonstrably unused legacy members/debug fragments while retaining
  the validated lattice/storage and output machinery.
