#!/usr/bin/env python3
"""Validate configurable axes for the public spatial-composition profiles."""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "output_features"
RUN_DIR = Path(__file__).resolve().parent / "_profile_axis_run"


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

    input_path = RUN_DIR / "nanokmc.in"
    text = input_path.read_text(encoding="utf-8")
    text = text.replace('AxialCompositionProfileAxis="X";',
                        'AxialCompositionProfileAxis="Y";')
    text = text.replace('CylindricalCompositionProfileAxis="X";',
                        'CylindricalCompositionProfileAxis="Z";')
    input_path.write_text(text, encoding="utf-8")

    print("[profile-axis-smoke] axial=Y cylindrical=Z")
    cp = subprocess.run([str(exe), str(RUN_DIR)], cwd=ROOT)
    if cp.returncode != 0:
        print(f"PROFILE AXIS VALIDATION: FAIL - process exited with {cp.returncode}")
        return 1

    checks = [
        (RUN_DIR / "evaluation" / "AxialCompositionProfile" / "axis_Y.csv", "Y"),
        (RUN_DIR / "evaluation" / "CylindricalCompositionProfile" / "axis_Z.csv", "Z"),
    ]
    errors: list[str] = []
    for path, expected_axis in checks:
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"missing or empty: {path.relative_to(RUN_DIR)}")
            continue
        with path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter=";"))
        if not rows:
            errors.append(f"no data rows: {path.relative_to(RUN_DIR)}")
            continue
        axes = {row.get("axis") for row in rows}
        if axes != {expected_axis}:
            errors.append(
                f"wrong axis values in {path.relative_to(RUN_DIR)}: {sorted(axes)}"
            )

    if errors:
        print("\nPROFILE AXIS VALIDATION: FAIL")
        for error in errors:
            print(" -", error)
        return 1

    print("\nPROFILE AXIS VALIDATION: PASS")
    print(" - AxialCompositionProfileAxis=Y")
    print(" - CylindricalCompositionProfileAxis=Z")
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
