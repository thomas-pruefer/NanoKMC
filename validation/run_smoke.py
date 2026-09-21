#!/usr/bin/env python3
"""Run and validate the six public NanoKMC solvers on small examples."""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
RUN_ROOT = Path(__file__).resolve().parent / "_runs"
CASES = [
    "classical",
    "active_filtered_generic",
    "active_filtered_binary_nn",
    "partial_filter_optimized",
    "rate_category_optimized",
    "exact_class_optimized",
]
TIMING_FIELDS = {"TimeMsCumulative", "EvolutionWallSecondsCumulative"}


def read_semicolon_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))


def prepare_case(name: str) -> Path:
    src = EXAMPLES / name
    dst = RUN_ROOT / name
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    (dst / "evaluation").mkdir(exist_ok=True)
    (dst / "output" / "bit").mkdir(parents=True, exist_ok=True)
    shutil.copy2(dst / "output" / "CalcData.template.csv",
                 dst / "output" / "CalcData.csv")
    return dst


def check_benchmark(name: str, rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    if not rows:
        return [f"{name}: empty Benchmark.csv"]
    a0 = rows[0]["AtomNumber0"]
    b0 = rows[0]["AtomNumber1"]
    for i, row in enumerate(rows):
        if row["BondNumberInternal"] != row["UnequalBondCount"]:
            errors.append(f"{name}: bond recount mismatch at row {i}")
        if row["AtomNumber0"] != a0 or row["AtomNumber1"] != b0:
            errors.append(f"{name}: species conservation failure at row {i}")
    if name == "exact_class_optimized":
        for i, row in enumerate(rows):
            if row["NJumpsRecord"] != row["NAcceptedRecord"]:
                errors.append(f"{name}: exact-class rejection-free check failed at row {i}")
    if name == "rate_category_optimized":
        saw_residual_rejection = False
        for i, row in enumerate(rows):
            jumps = int(float(row["NJumpsRecord"]))
            accepted = int(float(row["NAcceptedRecord"]))
            inactive = int(float(row.get("InactiveRejectedCountRecord", "0")))
            prob = int(float(row.get("ProbabilityEvaluationCountRecord", "0")))
            if inactive != 0:
                errors.append(f"{name}: structurally inactive proposal at row {i}")
            if jumps != prob:
                errors.append(f"{name}: exact selected-candidate probability not evaluated once per proposal at row {i}")
            if jumps > accepted:
                saw_residual_rejection = True
        if not saw_residual_rejection:
            errors.append(f"{name}: default four-category solver showed no residual rejection")
    if name == "classical":
        for i, row in enumerate(rows):
            if int(float(row.get("ActiveTableUpdatesRecord", "0"))) != 0:
                errors.append(f"{name}: classical solver unexpectedly updated an active table at row {i}")
    return errors


def compare_active_backends(generic: list[dict[str, str]],
                            binary: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    if len(generic) != len(binary):
        return ["Active-Filtered Generic/BinaryNN row-count mismatch"]
    for i, (g, b) in enumerate(zip(generic, binary)):
        keys = set(g) | set(b)
        for key in keys - TIMING_FIELDS:
            if g.get(key) != b.get(key):
                errors.append(
                    f"Active-Filtered same-seed mismatch row {i}, field {key}: "
                    f"Generic={g.get(key)} BinaryNN={b.get(key)}"
                )
                return errors
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, default=ROOT / "build" / "nanokmc")
    args = parser.parse_args()
    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"NanoKMC executable not found: {exe}")

    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    RUN_ROOT.mkdir(parents=True)
    results: dict[str, list[dict[str, str]]] = {}
    errors: list[str] = []

    for name in CASES:
        case = prepare_case(name)
        print(f"[smoke] {name}")
        cp = subprocess.run([str(exe), str(case)], cwd=ROOT)
        if cp.returncode != 0:
            errors.append(f"{name}: process exited with {cp.returncode}")
            continue
        bench = case / "evaluation" / "Benchmark.csv"
        if not bench.exists():
            errors.append(f"{name}: Benchmark.csv not produced")
            continue
        rows = read_semicolon_csv(bench)
        results[name] = rows
        errors.extend(check_benchmark(name, rows))

    if "active_filtered_generic" in results and "active_filtered_binary_nn" in results:
        errors.extend(compare_active_backends(
            results["active_filtered_generic"],
            results["active_filtered_binary_nn"],
        ))

    if errors:
        print("\nSMOKE VALIDATION: FAIL")
        for err in errors:
            print(" -", err)
        return 1

    print("\nSMOKE VALIDATION: PASS")
    print(" - all six public solvers completed")
    print(" - internal structural unlike-bond count matched independent recount")
    print(" - species counts were conserved")
    print(" - rate-category solver retained residual rejection without structural nulls")
    print(" - exact-class solver was rejection-free")
    print(" - Generic and BinaryNN Active-Filtered backends matched same-seed non-timing outputs")
    shutil.rmtree(RUN_ROOT, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
