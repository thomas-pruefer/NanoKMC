#!/usr/bin/env python3
"""Generate the NanoKMC documentation output gallery from one completed run.

Usage
-----
Run the output-feature example first, then generate the gallery::

    build/nanokmc examples/output_features
    python docs/examples/plot_output_features.py examples/output_features

The script reads only native NanoKMC outputs and writes PNG documentation
figures.  It does not rerun or alter the simulation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUN = ROOT / "examples" / "output_features"
DEFAULT_OUT = ROOT / "docs" / "assets"


def read_semicolon_dict(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=";"))


def latest_checkpoint(rows: list[dict[str, str]]) -> float:
    vals = [float(r["checkpoint_mcs"]) for r in rows]
    return max(vals)


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_distribution(run: Path, outdir: Path) -> None:
    path = run / "evaluation" / "clusters_S0.csv"
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = [x for x in line.rstrip().split(";") if x != ""]
            if not fields:
                continue
            checkpoint = float(fields[0])
            counts = np.asarray([int(v) for v in fields[1:]], dtype=float)
            # Native bins: 0..998 correspond to upper bounds
            # 10, 20, ..., 9990 sites; the final bin is >9990 sites.
            bounds = 10.0 * (np.arange(len(counts)) + 1)
            mask = counts > 0
            rows.append((checkpoint, bounds[mask], counts[mask]))

    fig, ax = plt.subplots(figsize=(7.3, 4.6))
    for checkpoint, x, y in rows:
        ax.plot(x, y, marker="o", linewidth=1.4, label=f"{checkpoint:g} MCS")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cluster-size bin upper bound [sites]")
    ax.set_ylabel("Number of A clusters")
    ax.set_title("Native ClusterDistribution output")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(title="Checkpoint")
    save(fig, outdir / "cluster_distribution.png")


def plot_axial(run: Path, outdir: Path) -> None:
    path = next((run / "evaluation" / "AxialCompositionProfile").glob("axis_*.csv"))
    rows = read_semicolon_dict(path)
    checkpoint = latest_checkpoint(rows)
    fig, ax = plt.subplots(figsize=(7.3, 4.6))
    species = sorted({(int(r["species"]), r["species_name"]) for r in rows})
    for idx, name in species:
        sub = [r for r in rows if float(r["checkpoint_mcs"]) == checkpoint and int(r["species"]) == idx]
        sub.sort(key=lambda r: int(r["coordinate_index"]))
        ax.plot([float(r["coordinate"]) for r in sub], [float(r["fraction"]) for r in sub],
                marker="o", linewidth=1.5, label=name)
    axis = rows[0]["axis"]
    ax.set_xlabel(f"{axis} coordinate")
    ax.set_ylabel("Species fraction")
    ax.set_ylim(0, 1)
    ax.set_title(f"AxialCompositionProfile at {checkpoint:g} MCS")
    ax.grid(True, alpha=0.2)
    ax.legend(title="Species")
    save(fig, outdir / "axial_composition_profile.png")


def plot_cylindrical(run: Path, outdir: Path) -> None:
    path = next((run / "evaluation" / "CylindricalCompositionProfile").glob("axis_*.csv"))
    rows = read_semicolon_dict(path)
    checkpoint = latest_checkpoint(rows)
    chosen = min(int(r["species"]) for r in rows)
    sub = [r for r in rows if float(r["checkpoint_mcs"]) == checkpoint and int(r["species"]) == chosen]
    axial = sorted({(int(r["axial_index"]), float(r["axial_coordinate"])) for r in sub})
    radial = sorted({(int(r["radial_bin"]), float(r["radial_coordinate"])) for r in sub})
    ai = {key: i for i, (key, _) in enumerate(axial)}
    ri = {key: i for i, (key, _) in enumerate(radial)}
    arr = np.full((len(radial), len(axial)), np.nan)
    for r in sub:
        arr[ri[int(r["radial_bin"])], ai[int(r["axial_index"])]] = float(r["fraction"])
    fig, ax = plt.subplots(figsize=(7.3, 4.9))
    im = ax.imshow(arr, origin="lower", aspect="auto", vmin=0, vmax=1,
                   extent=[axial[0][1], axial[-1][1], radial[0][1], radial[-1][1]])
    axis = rows[0]["axis"]
    species_name = sub[0]["species_name"]
    ax.set_xlabel(f"{axis}-axis coordinate")
    ax.set_ylabel("Radial coordinate")
    ax.set_title(f"CylindricalCompositionProfile: fraction {species_name}, {checkpoint:g} MCS")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(f"Fraction {species_name}")
    save(fig, outdir / "cylindrical_composition_profile.png")


def plot_spherical(run: Path, outdir: Path) -> None:
    path = run / "evaluation" / "SphericalCompositionProfile" / "profile.csv"
    rows = read_semicolon_dict(path)
    checkpoint = latest_checkpoint(rows)
    fig, ax = plt.subplots(figsize=(7.3, 4.6))
    species = sorted({(int(r["species"]), r["species_name"]) for r in rows})
    for idx, name in species:
        sub = [r for r in rows if float(r["checkpoint_mcs"]) == checkpoint and int(r["species"]) == idx]
        sub.sort(key=lambda r: int(r["radial_bin"]))
        ax.plot([float(r["radial_coordinate"]) for r in sub], [float(r["fraction"]) for r in sub],
                marker="o", linewidth=1.5, label=name)
    ax.set_xlabel("Radius from box centre")
    ax.set_ylabel("Species fraction")
    ax.set_ylim(0, 1)
    ax.set_title(f"SphericalCompositionProfile at {checkpoint:g} MCS")
    ax.grid(True, alpha=0.2)
    ax.legend(title="Species")
    save(fig, outdir / "spherical_composition_profile.png")


def plot_projection(run: Path, outdir: Path) -> None:
    path = run / "evaluation" / "ProjectedCompositionXZ" / "projection.csv"
    rows = read_semicolon_dict(path)
    checkpoint = latest_checkpoint(rows)
    chosen = min(int(r["species"]) for r in rows)
    sub = [r for r in rows if float(r["checkpoint_mcs"]) == checkpoint and int(r["species"]) == chosen]
    xs = sorted({(int(r["x_index"]), float(r["x_coordinate"])) for r in sub})
    zs = sorted({(int(r["z_index"]), float(r["z_coordinate"])) for r in sub})
    xi = {key: i for i, (key, _) in enumerate(xs)}
    zi = {key: i for i, (key, _) in enumerate(zs)}
    arr = np.full((len(zs), len(xs)), np.nan)
    for r in sub:
        arr[zi[int(r["z_index"])], xi[int(r["x_index"])]] = float(r["fraction"])
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    im = ax.imshow(arr, origin="lower", aspect="equal", vmin=0, vmax=1,
                   extent=[xs[0][1], xs[-1][1], zs[0][1], zs[-1][1]])
    species_name = sub[0]["species_name"]
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Z coordinate")
    ax.set_title(f"ProjectedCompositionXZ: fraction {species_name}, {checkpoint:g} MCS")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(f"Fraction {species_name}")
    save(fig, outdir / "projected_composition_xz.png")


def parse_xyz(path: Path) -> tuple[np.ndarray, np.ndarray]:
    species: list[str] = []
    xyz: list[tuple[float, float, float]] = []
    with path.open(encoding="utf-8") as handle:
        lines = handle.readlines()[2:]
    for line in lines:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            species.append(fields[0])
            xyz.append((float(fields[1]), float(fields[2]), float(fields[3])))
        except ValueError:
            continue
    return np.asarray(species), np.asarray(xyz, dtype=float)


def plot_morphology(run: Path, outdir: Path) -> None:
    xyz_files = sorted((run / "evaluation" / "OnefileRasmol").glob("*.xyz"))
    if not xyz_files:
        raise FileNotFoundError("No OnefileRasmol XYZ files found")
    path = xyz_files[-1]
    species, xyz = parse_xyz(path)
    fig = plt.figure(figsize=(6.3, 5.8))
    ax = fig.add_subplot(111, projection="3d")
    for name in sorted(set(species.tolist())):
        pts = xyz[species == name]
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=7, alpha=0.6, label=name)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(f"OnefileRasmol / BlenderSimple coordinate export ({int(path.stem)} MCS)")
    ax.legend(title="Species", loc="upper right")
    save(fig, outdir / "morphology_coordinate_export.png")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", nargs="?", type=Path, default=DEFAULT_RUN,
                        help="Completed NanoKMC run directory (default: examples/output_features)")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT,
                        help="Output directory (default: docs/assets)")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    outdir = args.outdir.resolve()

    required = run / "evaluation" / "Benchmark.csv"
    if not required.exists():
        raise SystemExit(
            f"No completed output-feature run found at {run}.\n"
            "Run `build/nanokmc examples/output_features` first, or pass a completed run directory."
        )

    plot_cluster_distribution(run, outdir)
    plot_axial(run, outdir)
    plot_cylindrical(run, outdir)
    plot_spherical(run, outdir)
    plot_projection(run, outdir)
    plot_morphology(run, outdir)
    print(f"Wrote NanoKMC documentation gallery to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
