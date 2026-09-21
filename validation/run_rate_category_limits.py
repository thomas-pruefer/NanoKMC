#!/usr/bin/env python3
"""Validate the two controlled limiting cases of KMCRateCategoryOptimized.

M=8 must reproduce the Exact-Class solver trajectory for the bundled binary-NN
model (apart from diagnostics that count the residual-probability evaluation).
M=1 is checked as a structurally active single-majorant population with residual
rejection and exact structural bond accounting.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = Path(__file__).resolve().parent / "_rate_category_limits"
TIMING_FIELDS = {"TimeMsCumulative", "EvolutionWallSecondsCumulative"}
DIAGNOSTIC_FIELDS = {
    "ProbabilityEvaluationCountRecord",
    "ProbabilityEvaluationsPerAccepted",
    "ProbabilityEvaluationFraction",
}


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))


def prepare_from_example(example: str, name: str) -> Path:
    src = ROOT / "examples" / example
    dst = RUN_ROOT / name
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    (dst / "evaluation").mkdir(exist_ok=True)
    (dst / "output" / "bit").mkdir(parents=True, exist_ok=True)
    shutil.copy2(dst / "output" / "CalcData.template.csv", dst / "output" / "CalcData.csv")
    return dst


def run(exe: Path, case: Path) -> list[dict[str, str]]:
    cp = subprocess.run([str(exe), str(case)], cwd=ROOT)
    if cp.returncode != 0:
        raise RuntimeError(f"{case.name}: executable returned {cp.returncode}")
    return read_rows(case / "evaluation" / "Benchmark.csv")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", type=Path, default=ROOT / "build" / "nanokmc")
    args = ap.parse_args()
    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"NanoKMC executable not found: {exe}")

    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    RUN_ROOT.mkdir(parents=True)

    exact_case = prepare_from_example("exact_class_optimized", "exact")
    m8_case = prepare_from_example("rate_category_optimized", "rate_M8")
    text = (m8_case / "nanokmc.in").read_text(encoding="utf-8")
    text = text.replace("RateCategoryCount=4;", "RateCategoryCount=8;")
    (m8_case / "nanokmc.in").write_text(text, encoding="utf-8")

    exact = run(exe, exact_case)
    m8 = run(exe, m8_case)
    if len(exact) != len(m8):
        print("RATE-CATEGORY LIMIT VALIDATION: FAIL")
        print(" - M=8 / Exact-Class row-count mismatch")
        return 1

    ignored = TIMING_FIELDS | DIAGNOSTIC_FIELDS
    for i, (a, b) in enumerate(zip(exact, m8)):
        for key in set(a) | set(b):
            if key in ignored:
                continue
            if a.get(key) != b.get(key):
                print("RATE-CATEGORY LIMIT VALIDATION: FAIL")
                print(f" - M=8 differs from Exact-Class at row {i}, {key}: {b.get(key)} != {a.get(key)}")
                return 1

    m1_case = prepare_from_example("rate_category_optimized", "rate_M1")
    text = (m1_case / "nanokmc.in").read_text(encoding="utf-8")
    text = text.replace("RateCategoryCount=4;", "RateCategoryCount=1;")
    (m1_case / "nanokmc.in").write_text(text, encoding="utf-8")
    m1 = run(exe, m1_case)
    saw_rejection = False
    for i, row in enumerate(m1):
        if row["BondNumberInternal"] != row["UnequalBondCount"]:
            print("RATE-CATEGORY LIMIT VALIDATION: FAIL")
            print(f" - M=1 bond recount mismatch at row {i}")
            return 1
        jumps = int(float(row["NJumpsRecord"]))
        accepted = int(float(row["NAcceptedRecord"]))
        inactive = int(float(row["InactiveRejectedCountRecord"]))
        if inactive != 0:
            print("RATE-CATEGORY LIMIT VALIDATION: FAIL")
            print(f" - M=1 produced structural null proposals at row {i}")
            return 1
        saw_rejection |= jumps > accepted
    if not saw_rejection:
        print("RATE-CATEGORY LIMIT VALIDATION: FAIL")
        print(" - M=1 showed no residual energetic rejection")
        return 1

    print("RATE-CATEGORY LIMIT VALIDATION: PASS")
    print(" - M=8 reproduced Exact-Class non-timing/non-diagnostic outputs")
    print(" - M=1 retained one active population with residual rejection")
    print(" - structural unlike-bond accounting remained exact")
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
