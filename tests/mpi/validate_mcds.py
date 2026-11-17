#!/usr/bin/env python3
"""Compare two PhysiCell MultiCellDS snapshots by labels, IDs, and coordinates."""

import argparse
from pathlib import Path
import sys

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from tools.compare_simulations import read_snapshot  # noqa: E402


def maximum_error(left, right):
    if left.size == 0:
        return 0.0
    return float(np.max(np.abs(left - right)))


def compare(reference, candidate, cell_atol, field_atol, rtol):
    failures = []

    if reference["fields"].keys() != candidate["fields"].keys():
        failures.append(
            f"substrate mismatch: {sorted(reference['fields'])} != {sorted(candidate['fields'])}"
        )
    if reference["labels"].keys() != candidate["labels"].keys():
        failures.append("cell-label sets differ")

    for key in ("bounds", "xyz", "weights"):
        if reference[key].shape != candidate[key].shape or not np.allclose(
            reference[key], candidate[key], atol=field_atol, rtol=rtol
        ):
            failures.append(f"microenvironment {key} differs")

    worst_field = ("none", 0.0)
    for name in reference["fields"].keys() & candidate["fields"].keys():
        left = reference["fields"][name]
        right = candidate["fields"][name]
        error = maximum_error(left, right)
        if error > worst_field[1]:
            worst_field = (name, error)
        if left.shape != right.shape or not np.allclose(left, right, atol=field_atol, rtol=rtol):
            failures.append(f"substrate {name!r} differs (max absolute error {error:.6e})")

    if "ID" not in reference["labels"] or "ID" not in candidate["labels"]:
        failures.append("saved cell data has no ID label")
        reference_order = candidate_order = np.array([], dtype=int)
    else:
        reference_ids = reference["labels"]["ID"][0]
        candidate_ids = candidate["labels"]["ID"][0]
        reference_order = np.argsort(reference_ids)
        candidate_order = np.argsort(candidate_ids)
        if not np.array_equal(reference_ids[reference_order], candidate_ids[candidate_order]):
            failures.append("cell ID sets differ")

    worst_cell = ("none", 0.0)
    if len(reference_order) == len(candidate_order):
        for name in reference["labels"].keys() & candidate["labels"].keys():
            left = reference["labels"][name][:, reference_order]
            right = candidate["labels"][name][:, candidate_order]
            error = maximum_error(left, right)
            if error > worst_cell[1]:
                worst_cell = (name, error)
            if left.shape != right.shape or not np.allclose(left, right, atol=cell_atol, rtol=rtol):
                failures.append(f"cell label {name!r} differs (max absolute error {error:.6e})")

    summary = (
        f"cells={reference['n']}, substrates={len(reference['fields'])}, "
        f"worst cell error={worst_cell[1]:.6e} ({worst_cell[0]}), "
        f"worst field error={worst_field[1]:.6e} ({worst_field[0]})"
    )
    return failures, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--cell-atol", type=float, default=1e-6)
    parser.add_argument("--field-atol", type=float, default=1e-10)
    parser.add_argument("--rtol", type=float, default=1e-10)
    args = parser.parse_args()

    reference = read_snapshot(args.reference)
    candidate = read_snapshot(args.candidate)
    failures, summary = compare(
        reference, candidate, args.cell_atol, args.field_atol, args.rtol
    )
    if failures:
        print("FAIL:", summary)
        for failure in failures[:20]:
            print(" -", failure)
        if len(failures) > 20:
            print(f" - ... {len(failures) - 20} additional differences")
        return 1

    print("PASS:", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
