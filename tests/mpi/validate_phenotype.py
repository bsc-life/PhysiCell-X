#!/usr/bin/env python3
"""Compare deterministic one-rank and multi-rank phenotype integrations."""

import argparse
import sys
from pathlib import Path

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from compare_simulations import read_snapshot  # noqa: E402


def relative_l2(reference, candidate):
    denominator = np.linalg.norm(reference)
    return float(np.linalg.norm(candidate - reference) / denominator) if denominator else 0.0


def relative_max_error(reference, candidate):
    scale = max(float(np.max(np.abs(reference))), float(np.max(np.abs(candidate))), 1e-30)
    return float(np.max(np.abs(candidate - reference)) / scale)


def sorted_label(snapshot, name, order):
    if name not in snapshot["labels"]:
        raise ValueError(f"{snapshot['path']}: missing required cell label '{name}'")
    return snapshot["labels"][name][:, order]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--position-max", type=float, default=1.0)
    parser.add_argument("--position-rmse", type=float, default=0.25)
    parser.add_argument("--field-relative-l2", type=float, default=0.02)
    parser.add_argument("--field-max", type=float, default=6.0)
    parser.add_argument("--internalized-relative-max", type=float, default=0.01)
    parser.add_argument("--pressure-max", type=float, default=0.1)
    args = parser.parse_args()

    reference = read_snapshot(args.reference)
    candidate = read_snapshot(args.candidate)
    failures = []

    if reference["n"] != candidate["n"]:
        failures.append(f"cell count differs: {reference['n']} versus {candidate['n']}")
    if reference["fields"].keys() != candidate["fields"].keys():
        failures.append("substrate names differ")

    reference_ids = reference["labels"]["ID"][0]
    candidate_ids = candidate["labels"]["ID"][0]
    reference_order = np.argsort(reference_ids)
    candidate_order = np.argsort(candidate_ids)
    if not np.array_equal(reference_ids[reference_order], candidate_ids[candidate_order]):
        failures.append("cell IDs differ")

    exact_labels = ("cell_type", "dead", "oncoprotein", "total_volume")
    for name in exact_labels:
        left = sorted_label(reference, name, reference_order)
        right = sorted_label(candidate, name, candidate_order)
        if not np.allclose(left, right, rtol=1e-12, atol=1e-10):
            failures.append(f"cell invariant '{name}' differs")

    reference_position = sorted_label(reference, "position", reference_order)
    candidate_position = sorted_label(candidate, "position", candidate_order)
    position_delta = candidate_position - reference_position
    position_max = float(np.max(np.abs(position_delta)))
    position_rmse = float(np.sqrt(np.mean(position_delta**2)))
    if position_max > args.position_max or position_rmse > args.position_rmse:
        failures.append("cell position error exceeds tolerance")

    reference_pressure = sorted_label(reference, "pressure", reference_order)
    candidate_pressure = sorted_label(candidate, "pressure", candidate_order)
    pressure_max = float(np.max(np.abs(candidate_pressure - reference_pressure)))
    if pressure_max > args.pressure_max:
        failures.append("cell pressure error exceeds tolerance")

    reference_internalized = sorted_label(
        reference, "internalized_total_substrates", reference_order
    )
    candidate_internalized = sorted_label(
        candidate, "internalized_total_substrates", candidate_order
    )
    internalized_relative_max = relative_max_error(
        reference_internalized, candidate_internalized
    )
    if internalized_relative_max > args.internalized_relative_max:
        failures.append("internalized substrate error exceeds tolerance")

    field_relative_l2 = 0.0
    field_max = 0.0
    for name in reference["fields"]:
        left = reference["fields"][name]
        right = candidate["fields"][name]
        field_relative_l2 = max(field_relative_l2, relative_l2(left, right))
        field_max = max(field_max, float(np.max(np.abs(right - left))))
    if field_relative_l2 > args.field_relative_l2 or field_max > args.field_max:
        failures.append("diffusive field error exceeds tolerance")

    status = "FAIL" if failures else "PASS"
    print(
        f"{status}: cells={reference['n']}, position max={position_max:.6g} um, "
        f"position RMSE={position_rmse:.6g} um, pressure max={pressure_max:.6g}, "
        f"internalized relative max={internalized_relative_max:.6%}, "
        f"field relative L2={field_relative_l2:.6%}, field max={field_max:.6g}"
    )
    for failure in failures:
        print(f" - {failure}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
