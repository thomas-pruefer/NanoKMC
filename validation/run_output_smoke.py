#!/usr/bin/env python3
"""Validate the documented NanoKMC 0.1.0 native output surface."""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "output_features"
RUN_DIR = Path(__file__).resolve().parent / "_output_run"


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0



def validate_fraction_table(path: Path, group_fields: tuple[str, ...], errors: list[str]) -> None:
    if not nonempty(path):
        errors.append(f"missing or empty: {path.relative_to(RUN_DIR)}")
        return
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter=";"))
    if not rows:
        errors.append(f"no data rows: {path.relative_to(RUN_DIR)}")
        return
    required = set(group_fields) | {"species", "count", "site_count", "fraction"}
    missing = required - set(rows[0])
    if missing:
        errors.append(f"missing columns in {path.relative_to(RUN_DIR)}: {sorted(missing)}")
        return
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in rows:
        key = tuple(row[field] for field in group_fields)
        grouped.setdefault(key, []).append(row)
    for key, group in grouped.items():
        try:
            frac_sum = sum(float(row["fraction"]) for row in group)
            site_counts = {int(row["site_count"]) for row in group}
            count_sum = sum(int(row["count"]) for row in group)
        except (ValueError, TypeError) as exc:
            errors.append(f"invalid numeric data in {path.relative_to(RUN_DIR)} group {key}: {exc}")
            return
        if len(site_counts) != 1 or count_sum != next(iter(site_counts)):
            errors.append(f"count/site-count mismatch in {path.relative_to(RUN_DIR)} group {key}")
            return
        if abs(frac_sum - 1.0) > 1e-10:
            errors.append(f"fractions do not sum to 1 in {path.relative_to(RUN_DIR)} group {key}: {frac_sum}")
            return

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, default=ROOT / "build" / "nanokmc")
    args = parser.parse_args()
    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"NanoKMC executable not found: {exe}")

    shutil.rmtree(RUN_DIR, ignore_errors=True)
    shutil.copytree(EXAMPLE, RUN_DIR)
    (RUN_DIR / "evaluation").mkdir(exist_ok=True)
    (RUN_DIR / "output" / "bit").mkdir(parents=True, exist_ok=True)
    shutil.copy2(RUN_DIR / "output" / "CalcData.template.csv",
                 RUN_DIR / "output" / "CalcData.csv")

    print("[output-smoke] output_features")
    cp = subprocess.run([str(exe), str(RUN_DIR)], cwd=ROOT)
    if cp.returncode != 0:
        print(f"OUTPUT VALIDATION: FAIL - process exited with {cp.returncode}")
        return 1

    required = [
        RUN_DIR / "evaluation" / "Benchmark.csv",
        RUN_DIR / "evaluation" / "clusters_S0.csv",
        RUN_DIR / "evaluation" / "clusters_S1.csv",
        RUN_DIR / "evaluation" / "surfaceatoms_S0.csv",
        RUN_DIR / "evaluation" / "surfaceatoms_S1.csv",
        RUN_DIR / "evaluation" / "clusterbonds_S0.csv",
        RUN_DIR / "evaluation" / "clusterbonds_S1.csv",
        RUN_DIR / "output" / "bit" / "data.zip",
    ]

    errors: list[str] = []
    for path in required:
        if not nonempty(path):
            errors.append(f"missing or empty: {path.relative_to(RUN_DIR)}")

    profile_tables = [
        (RUN_DIR / "evaluation" / "AxialCompositionProfile" / "axis_X.csv",
         ("checkpoint_mcs", "axis", "coordinate_index")),
        (RUN_DIR / "evaluation" / "CylindricalCompositionProfile" / "axis_X.csv",
         ("checkpoint_mcs", "axis", "axial_index", "radial_bin")),
        (RUN_DIR / "evaluation" / "SphericalCompositionProfile" / "profile.csv",
         ("checkpoint_mcs", "radial_bin")),
        (RUN_DIR / "evaluation" / "ProjectedCompositionXZ" / "projection.csv",
         ("checkpoint_mcs", "x_index", "z_index")),
    ]
    for table, groups in profile_tables:
        validate_fraction_table(table, groups, errors)

    patterns = {
        "RasMol per-species XYZ": "evaluation/Rasmol/*_S0.xyz",
        "RasMol script": "evaluation/Rasmol/*.rsm",
        "OnefileRasmol XYZ": "evaluation/OnefileRasmol/*.xyz",
        "BlenderSimple XYZ": "evaluation/BlenderSimple/*.xyz",
        "coordinate CSV": "evaluation/CSV/*.csv",
    }
    for label, pattern in patterns.items():
        hits = [p for p in RUN_DIR.glob(pattern) if nonempty(p)]
        if not hits:
            errors.append(f"no non-empty {label} output ({pattern})")

    bench = RUN_DIR / "evaluation" / "Benchmark.csv"
    if nonempty(bench):
        with bench.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter=";"))
        if not rows:
            errors.append("Benchmark.csv has no data rows")
        else:
            for i, row in enumerate(rows):
                if row.get("BondNumberInternal") != row.get("UnequalBondCount"):
                    errors.append(f"Benchmark.csv bond recount mismatch at row {i}")
                    break

    if errors:
        print("\nOUTPUT VALIDATION: FAIL")
        for error in errors:
            print(" -", error)
        return 1

    print("\nOUTPUT VALIDATION: PASS")
    print(" - Benchmark.csv")
    print(" - packed checkpoint archive")
    print(" - cluster / surface / cluster-bond distributions")
    print(" - axial / cylindrical / spherical composition profiles")
    print(" - projected XZ composition map")
    print(" - per-species RasMol XYZ + scripts")
    print(" - combined RasMol XYZ")
    print(" - BlenderSimple XYZ")
    print(" - plain coordinate CSV")
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
