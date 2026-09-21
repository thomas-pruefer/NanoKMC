# Output-feature example

This case uses `KMCActiveFilteredBinaryNN` and enables every native output feature that is release-smoke-tested in NanoKMC `0.1.0`:

```text
Benchmark ClusterDistribution AxialCompositionProfile CylindricalCompositionProfile SphericalCompositionProfile ProjectedCompositionXZ Rasmol OnefileRasmol BlenderSimple WriteCSV
```

Run from the repository root:

```text
build\nanokmc.exe examples\output_features
```

Expected outputs include:

- `evaluation/Benchmark.csv`
- `evaluation/clusters_S0.csv`, `clusters_S1.csv`
- `evaluation/surfaceatoms_S*.csv` and `clusterbonds_S*.csv`
- `evaluation/AxialCompositionProfile/axis_X.csv`
- `evaluation/CylindricalCompositionProfile/axis_X.csv`
- `evaluation/SphericalCompositionProfile/profile.csv`
- `evaluation/ProjectedCompositionXZ/projection.csv`
- `evaluation/Rasmol/*.xyz` and `*.rsm`
- `evaluation/OnefileRasmol/*.xyz` and `*.rsm`
- `evaluation/BlenderSimple/*.xyz`
- `evaluation/CSV/*.csv`
- `output/bit/data.zip`

See `docs/output-features.md` for conventions and scope.


The example uses the X axis for both axial and cylindrical profiles:

```text
AxialCompositionProfileAxis="X";
CylindricalCompositionProfileAxis="X";
```

The profile outputs are tidy semicolon-delimited CSVs with headers and per-bin species fractions.

## Plot the generated outputs

After this example completes, regenerate the documentation-style plots with:

```bash
python -m pip install -r docs/requirements-plotting.txt
python docs/examples/plot_output_features.py examples/output_features --outdir docs/assets
```

The plotting helper is optional and is not required by NanoKMC itself.
