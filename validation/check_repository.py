#!/usr/bin/env python3
"""Check release-tree invariants that should remain true for the public repository."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_SOLVERS = {
    "KMCClassical",
    "KMCActiveFilteredGeneric",
    "KMCActiveFilteredBinaryNN",
    "KMCPartialFilterOptimized",
    "KMCRateCategoryOptimized",
    "KMCExactClassOptimized",
}
EXAMPLE_BY_SOLVER = {
    "KMCClassical": "classical",
    "KMCActiveFilteredGeneric": "active_filtered_generic",
    "KMCActiveFilteredBinaryNN": "active_filtered_binary_nn",
    "KMCPartialFilterOptimized": "partial_filter_optimized",
    "KMCRateCategoryOptimized": "rate_category_optimized",
    "KMCExactClassOptimized": "exact_class_optimized",
}


def fail(message: str) -> None:
    raise SystemExit(f"REPOSITORY CHECK: FAIL - {message}")


def main() -> int:
    required = [
        "README.md", "LICENSE", "CITATION.cff", "CONTRIBUTING.md",
        "CHANGELOG.md", "CMakeLists.txt", "docs/quickstart.md",
        "docs/reproducibility.md", "docs/validation.md",
    ]
    for rel in required:
        if not (ROOT / rel).is_file():
            fail(f"missing required file: {rel}")

    # Implementation files must never again be included from headers.
    for header in list((ROOT / "src").glob("*.h")) + list((ROOT / "src/internal").glob("*.h")):
        text = header.read_text(encoding="utf-8")
        if re.search(r'#include\s+["<].*\.cpp[">]', text):
            fail(f"implementation file included from header: {header.relative_to(ROOT)}")

    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    version_match = re.search(r'project\(NanoKMC VERSION ([0-9.]+)', cmake)
    if not version_match:
        fail("could not read project version from CMakeLists.txt")
    version = version_match.group(1)
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if f"version: {version}" not in cff:
        fail(f"CITATION.cff version does not match CMake version {version}")

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for solver in sorted(PUBLIC_SOLVERS):
        if solver not in readme:
            fail(f"README does not mention public solver {solver}")
        example = ROOT / "examples" / EXAMPLE_BY_SOLVER[solver]
        if not (example / "nanokmc.in").is_file():
            fail(f"missing example nanokmc.in for {solver}")
        if not (example / "output/CalcData.template.csv").is_file():
            fail(f"missing CalcData template for {solver}")

    root_batch = list(ROOT.glob("*.bat"))
    if root_batch:
        fail("batch helpers must live under scripts/windows, not repository root")

    generated_validation_dirs = [
        "validation/_runs",
        "validation/_output_run",
        "validation/_profile_axis_run",
        "validation/_rate_category_limits",
        "validation/_input_validation",
    ]
    for rel in generated_validation_dirs:
        if (ROOT / rel).exists():
            fail(f"generated validation output must not ship in public release: {rel}")

    maintenance_docs = {
        "docs/release-hardening-audit.md",
        "docs/release-validation.md",
        "docs/source-cleanup.md",
    }
    for rel in maintenance_docs:
        if (ROOT / rel).exists():
            fail(f"internal maintenance note should not ship in public release: {rel}")

    print("REPOSITORY CHECK: PASS")
    print(f" - version metadata consistent: {version}")
    print(" - six public solver examples present")
    print(" - no .cpp implementation files included from headers")
    print(" - required release/documentation files present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
