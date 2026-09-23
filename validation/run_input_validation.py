#!/usr/bin/env python3
"""Regression checks for public input parsing and validation."""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "classical"
RUN_ROOT = Path(__file__).resolve().parent / "_input_validation"
MAX_UINT64 = "18446744073709551615"


def prepare_case(name: str) -> Path:
    dst = RUN_ROOT / name
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(EXAMPLE, dst)
    (dst / "evaluation").mkdir(exist_ok=True)
    (dst / "output" / "bit").mkdir(parents=True, exist_ok=True)
    shutil.copy2(dst / "output" / "CalcData.template.csv",
                 dst / "output" / "CalcData.csv")
    return dst


def replace_parameter(path: Path, key: str, value: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    prefix = f"{key}="
    replaced = False
    out: list[str] = []
    for line in lines:
        if line.startswith(prefix):
            out.append(f"{key}={value};")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        raise RuntimeError(f"parameter {key} not found in {path}")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def run(exe: Path, case: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(exe), str(case)],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def benchmark_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as f:
        return len(list(csv.DictReader(f, delimiter=";")))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, default=ROOT / "build" / "nanokmc")
    args = parser.parse_args()
    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"NanoKMC executable not found: {exe}")

    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    RUN_ROOT.mkdir(parents=True)
    errors: list[str] = []

    # 1. Seed is a true uint64_t, including the maximum representable value.
    case = prepare_case("uint64_seed")
    replace_parameter(case / "nanokmc.in", "Seed", MAX_UINT64)
    cp = run(exe, case)
    if cp.returncode != 0:
        errors.append(f"uint64 seed run failed with {cp.returncode}:\n{cp.stdout}")
    else:
        log = (case / "log.dat").read_text(encoding="utf-8", errors="replace")
        if f"Seed: {MAX_UINT64}" not in log:
            errors.append("maximum uint64_t seed was not preserved in the runtime log")
        eval_sce = (case / "evaluation" / "eval.sce").read_text(
            encoding="utf-8", errors="replace")
        # Exact key matching: T must resolve to T=100, not the earlier kT=1.0.
        if 'setIntMeta("CalcT",100);' not in eval_sce:
            errors.append("exact key parsing failed: metadata T was confused with kT")

    # 2. Blank lines in CalcData.csv are ignored rather than becoming zero rows.
    case = prepare_case("blank_calcdata_line")
    calc = case / "output" / "CalcData.csv"
    original = calc.read_text(encoding="utf-8").splitlines()
    calc.write_text(original[0] + "\n\n" + "\n".join(original[1:]) + "\n",
                    encoding="utf-8")
    cp = run(exe, case)
    if cp.returncode != 0:
        errors.append(f"blank-line CalcData run failed with {cp.returncode}:\n{cp.stdout}")
    else:
        bench = case / "evaluation" / "Benchmark.csv"
        if not bench.exists() or benchmark_rows(bench) != len(original):
            errors.append("blank CalcData line created an unintended checkpoint record")

    # 3. Packed FCC dimensions must be at least 2; reject safely before shifting.
    case = prepare_case("invalid_dimension")
    replace_parameter(case / "nanokmc.in", "knx", "1")
    cp = run(exe, case)
    if cp.returncode == 0:
        errors.append("knx=1 was accepted; expected a non-zero validation exit")
    elif "knx, kny and knz must each be >= 2" not in cp.stdout:
        errors.append("knx=1 failed, but without the expected clear validation message")

    # 4. Checkpoint schedules must be non-negative and nondecreasing.
    case = prepare_case("decreasing_checkpoint")
    calc = case / "output" / "CalcData.csv"
    lines = calc.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        errors.append("classical example unexpectedly has fewer than two checkpoints")
    else:
        first_fields = lines[0].split(";")
        second_fields = lines[1].split(";")
        first_fields[0] = "10"
        second_fields[0] = "5"
        lines[0] = ";".join(first_fields)
        lines[1] = ";".join(second_fields)
        calc.write_text("\n".join(lines) + "\n", encoding="utf-8")
        cp = run(exe, case)
        if cp.returncode == 0:
            errors.append("decreasing CalcData checkpoints were accepted")
        elif "checkpoints must be nondecreasing" not in cp.stdout:
            errors.append("decreasing checkpoints failed without the expected diagnostic")

    if errors:
        print("INPUT VALIDATION: FAIL")
        for err in errors:
            print(" -", err)
        return 1

    print("INPUT VALIDATION: PASS")
    print(" - full uint64_t seed range is parsed without truncation")
    print(" - parameter names are matched exactly (T is distinct from kT)")
    print(" - blank CalcData.csv lines are ignored")
    print(" - invalid packed-FCC dimensions fail safely")
    print(" - checkpoint schedules reject decreasing MCS values")
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
