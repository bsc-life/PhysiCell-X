#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path

POSITION_TOL = 1e-3
VELOCITY_TOL = 1e-5
PRESSURE_TOL = 1e-10
INTERFACE_HALF_WIDTH = 10.0


def parse_neighbor_ids(raw_value: str):
    if not raw_value:
        return ()
    return tuple(sorted(int(value) for value in raw_value.split(";") if value))


def load_snapshots(path: Path):
    snapshots = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            iteration = int(row["iteration"])
            cell_id = int(row["cell_id"])
            snapshots.setdefault(iteration, {})[cell_id] = {
                "time": float(row["time"]),
                "position": (
                    float(row["x"]),
                    float(row["y"]),
                    float(row["z"]),
                ),
                "velocity": (
                    float(row["vx"]),
                    float(row["vy"]),
                    float(row["vz"]),
                ),
                "pressure": float(row["pressure"]),
                "mechanics_neighbor_candidates": int(row["mechanics_neighbor_candidates"]),
                "mechanics_neighbor_interactions": int(row["mechanics_neighbor_interactions"]),
                "neighbors": parse_neighbor_ids(row["neighbor_ids"]),
            }
    return snapshots


def vector_error(left, right):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def compare(reference, mpi):
    report_lines = []
    ref_iterations = sorted(reference)
    mpi_iterations = sorted(mpi)

    if ref_iterations != mpi_iterations:
        report_lines.append("FAIL")
        report_lines.append("First failing iteration: n/a")
        report_lines.append("Failing cell IDs: n/a")
        report_lines.append(f"Iteration mismatch: reference={ref_iterations}, mpi={mpi_iterations}")
        return False, "\n".join(report_lines)

    worst_position_error = 0.0
    worst_velocity_error = 0.0
    worst_pressure_error = 0.0
    first_failure_iteration = None
    failing_cells = []

    for iteration in ref_iterations:
        ref_cells = reference[iteration]
        mpi_cells = mpi[iteration]

        ref_ids = sorted(ref_cells)
        mpi_ids = sorted(mpi_cells)

        if len(ref_ids) != len(mpi_ids):
            report_lines.append("FAIL")
            report_lines.append(f"First failing iteration: {iteration}")
            report_lines.append("Failing cell IDs: n/a")
            report_lines.append(
                f"Iteration {iteration}: cell count mismatch reference={len(ref_ids)} mpi={len(mpi_ids)}"
            )
            return False, "\n".join(report_lines)

        if ref_ids != mpi_ids:
            report_lines.append("FAIL")
            report_lines.append(f"First failing iteration: {iteration}")
            report_lines.append("Failing cell IDs: n/a")
            report_lines.append(f"Iteration {iteration}: cell ID mismatch reference={ref_ids} mpi={mpi_ids}")
            return False, "\n".join(report_lines)

        for cell_id in ref_ids:
            ref_row = ref_cells[cell_id]
            mpi_row = mpi_cells[cell_id]

            if ref_row["neighbors"] != mpi_row["neighbors"]:
                if first_failure_iteration is None:
                    first_failure_iteration = iteration
                failing_cells.append(cell_id)
                report_lines.append("FAIL")
                report_lines.append(f"First failing iteration: {first_failure_iteration}")
                report_lines.append(f"Failing cell IDs: {sorted(set(failing_cells))}")
                report_lines.append(
                    f"Iteration {iteration}, cell {cell_id}: neighbor mismatch "
                    f"reference={ref_row['neighbors']} mpi={mpi_row['neighbors']}"
                )
                return False, "\n".join(report_lines)

            position_error = vector_error(ref_row["position"], mpi_row["position"])
            velocity_error = vector_error(ref_row["velocity"], mpi_row["velocity"])
            pressure_error = abs(ref_row["pressure"] - mpi_row["pressure"])
            worst_position_error = max(worst_position_error, position_error)
            worst_velocity_error = max(worst_velocity_error, velocity_error)
            worst_pressure_error = max(worst_pressure_error, pressure_error)

            count_fields = (
                "mechanics_neighbor_candidates",
                "mechanics_neighbor_interactions",
            )
            count_mismatches = {
                field: (ref_row[field], mpi_row[field])
                for field in count_fields
                if ref_row[field] != mpi_row[field]
            }

            if (
                position_error > POSITION_TOL
                or velocity_error > VELOCITY_TOL
                or pressure_error > PRESSURE_TOL
                or count_mismatches
            ):
                if first_failure_iteration is None:
                    first_failure_iteration = iteration
                failing_cells.append(cell_id)
                report_lines.append("FAIL")
                report_lines.append(f"First failing iteration: {first_failure_iteration}")
                report_lines.append(f"Failing cell IDs: {sorted(set(failing_cells))}")
                report_lines.append(
                    f"Iteration {iteration}, cell {cell_id}: "
                    f"position_error={position_error:.6e}, velocity_error={velocity_error:.6e}, "
                    f"pressure_error={pressure_error:.6e}, count_mismatches={count_mismatches}"
                )
                report_lines.append(
                    f"Worst errors so far: position={worst_position_error:.6e}, "
                    f"velocity={worst_velocity_error:.6e}, pressure={worst_pressure_error:.6e}"
                )
                return False, "\n".join(report_lines)

    mechanics_iterations = [iteration for iteration in mpi_iterations if iteration > 0]
    if mechanics_iterations:
        first_mechanics_iteration = mechanics_iterations[0]
        interface_cells = {
            cell_id: row
            for cell_id, row in mpi[first_mechanics_iteration].items()
            if abs(row["position"][0]) < INTERFACE_HALF_WIDTH
        }
        has_left_cell = any(row["position"][0] < 0.0 for row in interface_cells.values())
        has_right_cell = any(row["position"][0] > 0.0 for row in interface_cells.values())
        inactive_cells = [
            cell_id
            for cell_id, row in interface_cells.items()
            if row["pressure"] <= 0.0 or row["mechanics_neighbor_interactions"] == 0
        ]
        if not has_left_cell or not has_right_cell or inactive_cells:
            report_lines.append("FAIL")
            report_lines.append(f"First failing iteration: {first_mechanics_iteration}")
            report_lines.append(f"Failing cell IDs: {inactive_cells}")
            report_lines.append(
                "Cross-interface validation requires active, positive-pressure cells "
                f"on both sides of x=0; interface_cells={sorted(interface_cells)}, "
                f"inactive_cells={inactive_cells}"
            )
            return False, "\n".join(report_lines)

    report_lines.append("PASS")
    report_lines.append(f"Iterations compared: {len(ref_iterations)}")
    report_lines.append(f"Cells per iteration: {len(reference[ref_iterations[0]]) if ref_iterations else 0}")
    report_lines.append(f"Worst position error: {worst_position_error:.6e}")
    report_lines.append(f"Worst velocity error: {worst_velocity_error:.6e}")
    report_lines.append(f"Worst pressure error: {worst_pressure_error:.6e}")
    return True, "\n".join(report_lines)


def main(argv):
    if len(argv) not in (3, 4):
        print("Usage: compare_snapshots.py reference.csv mpi.csv [report.txt]", file=sys.stderr)
        return 2

    reference_path = Path(argv[1])
    mpi_path = Path(argv[2])
    report_path = Path(argv[3]) if len(argv) == 4 else None

    reference = load_snapshots(reference_path)
    mpi = load_snapshots(mpi_path)
    success, report = compare(reference, mpi)

    if report_path is not None:
        report_path.write_text(report + "\n")

    print(report)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
