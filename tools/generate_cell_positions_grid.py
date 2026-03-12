#!/usr/bin/env python3
import argparse
import math
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate evenly spaced cell positions inside a cube and "
            "output a PhysiCell CSV (x,y,z,typeID)."
        )
    )
    parser.add_argument("num_cells", type=int, help="Number of cells to generate")
    parser.add_argument("side_um", type=float, help="Cube side length in microns (um)")
    parser.add_argument(
        "num_types",
        type=int,
        help="Total number of cell types in the simulation",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="-",
        help="Output CSV path (default: stdout)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for type assignment (default: none)",
    )
    return parser.parse_args()


def generate_positions(num_cells: int, side_um: float):
    if num_cells <= 0:
        raise ValueError("num_cells must be positive")
    if side_um <= 0:
        raise ValueError("side_um must be positive")

    n_per_axis = math.ceil(num_cells ** (1.0 / 3.0))
    spacing = side_um / n_per_axis
    start = -side_um / 2.0 + spacing / 2.0

    for idx in range(num_cells):
        i = idx % n_per_axis
        j = (idx // n_per_axis) % n_per_axis
        k = idx // (n_per_axis * n_per_axis)
        x = start + i * spacing
        y = start + j * spacing
        z = start + k * spacing
        yield x, y, z


def main() -> int:
    args = parse_args()
    try:
        positions = generate_positions(args.num_cells, args.side_um)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    if args.num_types <= 0:
        print("Error: num_types must be positive", file=sys.stderr)
        return 2

    if args.seed is not None:
        import random
        rng = random.Random(args.seed)
    else:
        import random
        rng = random

    out_fh = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
    try:
        for x, y, z in positions:
            type_id = rng.randrange(args.num_types)
            out_fh.write(f"{x:.6f},{y:.6f},{z:.6f},{type_id}\n")
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
