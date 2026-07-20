#!/usr/bin/env python3
"""
Compare OpenPIC grid diagnostic files from two runs.

Usage:
    python tools/compare_grid_outputs.py RUN_A_DIR RUN_B_DIR
    python tools/compare_grid_outputs.py RUN_A_DIR RUN_B_DIR --atol 1e-10 --rtol 1e-8

The script matches files by relative path and filename pattern *grd*.dat.
It expects tab/space separated OpenPIC diagnostic tables with a header row.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable


def find_grid_files(root: Path) -> dict[str, Path]:
    return {
        str(path.relative_to(root)).replace("\\", "/"): path
        for path in root.rglob("*grd*.dat")
        if path.is_file()
    }


def parse_table(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    if not lines:
        raise ValueError(f"{path}: empty file")

    header = lines[0].split()
    rows: list[list[float]] = []
    for line_no, line in enumerate(lines[1:], start=2):
        fields = line.split()
        if len(fields) != len(header):
            raise ValueError(
                f"{path}:{line_no}: expected {len(header)} columns, got {len(fields)}"
            )
        try:
            rows.append([float(value) for value in fields])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: non-numeric value") from exc

    return header, rows


def rel_error(diff: float, a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1.0)
    return diff / scale


def compare_rows(
    rel_name: str,
    path_a: Path,
    path_b: Path,
    atol: float,
    rtol: float,
) -> tuple[bool, str]:
    header_a, rows_a = parse_table(path_a)
    header_b, rows_b = parse_table(path_b)

    if header_a != header_b:
        return False, f"{rel_name}: header mismatch"

    if len(rows_a) != len(rows_b):
        return False, f"{rel_name}: row count mismatch {len(rows_a)} != {len(rows_b)}"

    max_abs = 0.0
    max_rel = 0.0
    max_location = ""

    for row_idx, (row_a, row_b) in enumerate(zip(rows_a, rows_b), start=1):
        for col_idx, (value_a, value_b) in enumerate(zip(row_a, row_b)):
            if not (math.isfinite(value_a) and math.isfinite(value_b)):
                if value_a != value_b:
                    return False, (
                        f"{rel_name}: non-finite mismatch at row {row_idx}, "
                        f"column {header_a[col_idx]}"
                    )
                continue

            diff = abs(value_a - value_b)
            rel = rel_error(diff, value_a, value_b)
            if diff > max_abs or rel > max_rel:
                max_abs = max(max_abs, diff)
                max_rel = max(max_rel, rel)
                max_location = f"row {row_idx}, column {header_a[col_idx]}"

            if diff > atol and rel > rtol:
                return False, (
                    f"{rel_name}: max tolerance exceeded at row {row_idx}, "
                    f"column {header_a[col_idx]}: "
                    f"a={value_a:.17g}, b={value_b:.17g}, "
                    f"abs={diff:.3e}, rel={rel:.3e}"
                )

    detail = f"max_abs={max_abs:.3e}, max_rel={max_rel:.3e}"
    if max_location:
        detail += f" at {max_location}"
    return True, f"{rel_name}: OK ({detail})"


def sorted_names(names: Iterable[str]) -> list[str]:
    return sorted(names, key=lambda value: value.lower())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_a", type=Path)
    parser.add_argument("run_b", type=Path)
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    files_a = find_grid_files(args.run_a)
    files_b = find_grid_files(args.run_b)

    names_a = set(files_a)
    names_b = set(files_b)
    missing_in_b = sorted_names(names_a - names_b)
    missing_in_a = sorted_names(names_b - names_a)
    common = sorted_names(names_a & names_b)

    ok = True

    for name in missing_in_b:
        ok = False
        print(f"Missing in B: {name}")
    for name in missing_in_a:
        ok = False
        print(f"Missing in A: {name}")

    compared = 0
    for name in common:
        same, message = compare_rows(name, files_a[name], files_b[name], args.atol, args.rtol)
        compared += 1
        if not same:
            ok = False
        if not args.quiet or not same:
            print(message)

    print(
        f"\nCompared {compared} common grid files; "
        f"missing_in_a={len(missing_in_a)}, missing_in_b={len(missing_in_b)}."
    )

    if ok:
        print("PASS")
        return 0

    print("FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
