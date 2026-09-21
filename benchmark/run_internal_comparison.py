#!/usr/bin/env python3
"""Sequential six-solver NanoKMC benchmark used for release verification.

This script intentionally benchmarks only the six public solvers in this
repository. It does not replace the cross-software benchmark campaign used by
the paper.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "classical": "KMCClassical",
    "active_filtered_binary_nn": "KMCActiveFilteredBinaryNN",
    "active_filtered_generic": "KMCActiveFilteredGeneric",
    "partial_filter_optimized": "KMCPartialFilterOptimized",
    "rate_category_optimized": "KMCRateCategoryOptimized",
    "exact_class_optimized": "KMCExactClassOptimized",
}
DEFAULT_CHECKPOINTS = [0, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000]
TIMING_FIELDS = {"TimeMsCumulative", "EvolutionWallSecondsCumulative"}


def input_text(system_id: str, k: int, mcs: int, seed: int,
               x_a: float, kt: float) -> str:
    return f'''// NanoKMC public internal benchmark\nSystemID="{system_id}";\nSeed={seed};\nSeedStart={seed};\nSeedEnd={seed};\nSeedStep=1;\nNSpecies=2;\nknx={k};\nkny={k};\nknz={k};\nlc=0.4338;\nkT={kt};\nSpeciesName(1)="A";\nSpeciesName(2)="B";\nSpeciesColor(1)="[255,0,0]";\nSpeciesColor(2)="[0,0,255]";\nRasmolSpeciesPlot(1)=1;\nRasmolSpeciesPlot(2)=1;\nEvalParam=" Benchmark ";\nSysEvalParam=" ";\nIndex="StrCalcSystemID_IntCalcT_IntSysEB_IntSysdE_IntCalcSeed";\nT={mcs};\nEB=1.0;\ndE=0.0;\nE01=0.0;\nf01=1.0;\nEa=1.0;\ndev=0.0;\nclvl={x_a};\nFluence=0;\n'''


def prepare_case(base: Path, name: str, sid: str, args) -> Path:
    d = base / name
    shutil.rmtree(d, ignore_errors=True)
    (d / "output" / "bit").mkdir(parents=True)
    (d / "evaluation").mkdir(parents=True)
    (d / "nanokmc.in").write_text(
        input_text(sid, args.k, args.mcs, args.seed, args.xa, args.kt),
        encoding="utf-8",
    )
    checkpoints = [x for x in DEFAULT_CHECKPOINTS if x <= args.mcs]
    if checkpoints[-1] != args.mcs:
        checkpoints.append(args.mcs)
    calc = "".join(f"{x};0;0;0;0;0;0;\n" for x in checkpoints)
    (d / "output" / "CalcData.csv").write_text(calc, encoding="utf-8")
    (d / "output" / "CalcData.template.csv").write_text(calc, encoding="utf-8")
    return d


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))


def compare_active(generic, binary):
    for i, (g, b) in enumerate(zip(generic, binary)):
        for key in set(g) | set(b):
            if key in TIMING_FIELDS:
                continue
            if g.get(key) != b.get(key):
                return False, f"row {i}, field {key}: {g.get(key)} != {b.get(key)}"
    return True, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", type=Path, default=ROOT / "build" / "nanokmc")
    ap.add_argument("--k", type=int, choices=[5, 6], default=5)
    ap.add_argument("--mcs", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--xa", type=float, default=0.10)
    ap.add_argument("--kt", type=float, default=0.75)
    args = ap.parse_args()
    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"Executable not found: {exe}")

    tag = f"k{args.k}_xA{int(round(args.xa*100)):02d}_T{str(args.kt).replace('.', 'p')}_{args.mcs}mcs"
    base = ROOT / "benchmark" / "results" / tag
    shutil.rmtree(base, ignore_errors=True)
    base.mkdir(parents=True)

    summary = []
    all_rows = {}
    for name, sid in CASES.items():
        case = prepare_case(base, name, sid, args)
        print(f"[benchmark] {name}")
        cp = subprocess.run([str(exe), str(case)], cwd=ROOT)
        if cp.returncode != 0:
            return cp.returncode
        rows = read_rows(case / "evaluation" / "Benchmark.csv")
        all_rows[name] = rows
        last = rows[-1]
        accepted = sum(int(float(r["NAcceptedRecord"])) for r in rows)
        selected = sum(int(float(r["NJumpsRecord"])) for r in rows)
        sec = float(last["EvolutionWallSecondsCumulative"])
        summary.append({
            "case": name,
            "system_id": sid,
            "evolution_wall_seconds": sec,
            "selected_or_proposed_events": selected,
            "executed_exchanges": accepted,
            "microseconds_per_executed_exchange": sec * 1e6 / accepted if accepted else 0.0,
            "final_interface_bond_fraction": float(last["InterfaceDensityFCC"]),
            "final_common_mcs_exact": float(last["CommonMCSExact"]),
            "total_atoms": int(float(last["TotalAtoms"])),
        })

    ok, detail = compare_active(
        all_rows["active_filtered_generic"],
        all_rows["active_filtered_binary_nn"],
    )

    out_csv = base / "performance_summary.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]))
        w.writeheader(); w.writerows(summary)

    by = {r["case"]: r for r in summary}
    ratios = {
        "classical_over_binary_runtime": by["classical"]["evolution_wall_seconds"] /
                                         by["active_filtered_binary_nn"]["evolution_wall_seconds"],
        "generic_over_binary_runtime": by["active_filtered_generic"]["evolution_wall_seconds"] /
                                       by["active_filtered_binary_nn"]["evolution_wall_seconds"],
        "partial_over_binary_runtime": by["partial_filter_optimized"]["evolution_wall_seconds"] /
                                       by["active_filtered_binary_nn"]["evolution_wall_seconds"],
        "rate_category_over_binary_runtime": by["rate_category_optimized"]["evolution_wall_seconds"] /
                                               by["active_filtered_binary_nn"]["evolution_wall_seconds"],
        "exact_over_binary_runtime": by["exact_class_optimized"]["evolution_wall_seconds"] /
                                     by["active_filtered_binary_nn"]["evolution_wall_seconds"],
        "active_backend_same_seed_equivalent": ok,
        "active_backend_mismatch_detail": detail,
    }
    (base / "comparison.json").write_text(json.dumps(ratios, indent=2), encoding="utf-8")
    (base / "metadata.json").write_text(json.dumps({
        "k": args.k, "x_A": args.xa, "kT": args.kt, "target_mcs": args.mcs,
        "seed": args.seed, "run_policy": "sequential", "paper_build": "-O1 -DNDEBUG",
    }, indent=2), encoding="utf-8")

    print("\nInternal comparison complete")
    for r in summary:
        print(f"  {r['system_id']}: {r['evolution_wall_seconds']:.6g} s")
    print(f"  Generic/BinaryNN same-seed non-timing equivalence: {'PASS' if ok else 'FAIL'}")
    print(f"Results: {base}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
