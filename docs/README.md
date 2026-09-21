# NanoKMC documentation

Recommended reading order:

1. [`quickstart.md`](quickstart.md) — build, validate, and run a first simulation.
2. [`reproducibility.md`](reproducibility.md) — reproducibility boundary and paper-build policy.
3. [`physics-model.md`](physics-model.md) — benchmark lattice and Hamiltonian.
4. [`solvers.md`](solvers.md) — six public solver architectures.
5. [`architecture.md`](architecture.md) — implementation structure.
6. [`input-output.md`](input-output.md) — run-directory format and parameters.
7. [`output-features.md`](output-features.md) — observables and state/visualization exports.
8. [`literature.md`](literature.md) — scientific provenance and prior-art boundaries.
9. [`validation.md`](validation.md) — regression and scientific checks.
10. [`extending.md`](extending.md) — extending the Generic backend.
11. [`known-limitations.md`](known-limitations.md) — deliberate scope and technical limitations.
12. [`windows-msys2.md`](windows-msys2.md) — Windows/MSYS2 workflow.

## Documentation figures

The gallery in `docs/assets/` is generated from native NanoKMC output:

```bash
python -m pip install -r docs/requirements-plotting.txt
python validation/run_output_smoke.py --exe build/nanokmc
python docs/examples/plot_output_features.py validation/_output_run --outdir docs/assets
```

See [`output-features.md`](output-features.md) for definitions and limitations.
