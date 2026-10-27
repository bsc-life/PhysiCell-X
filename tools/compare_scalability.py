#!/usr/bin/env python3
"""Compare a PhysiCell OMP baseline with PhysiCell-X MPI+OMP CSV results.

Usage:
    python3 tools/compare_scalability.py PHYSICELL_CSV_OR_DIR \
        PHYSICELL_X_RESULTS_DIR OUTPUT_DIR

The PhysiCell-X directory is searched recursively for ``*_nodes.csv`` files.
Timing summaries use the mean of non-zero-duration rows.  The output contains
individual plots, a combined figure, and a CSV with the plotted summaries.

Requires matplotlib.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
import sys
from typing import Iterable, Sequence
import warnings


REQUIRED_COLUMNS = (
    "simulation_time",
    "elapsed_wall_seconds",
    "diffusion_seconds",
    "mechanics_seconds",
    "phenotype_seconds",
    "intracellular_seconds",
    "diffusion_percent",
    "mechanics_percent",
    "phenotype_percent",
    "intracellular_percent",
    "memory_gb",
    "mpi_ranks",
    "omp_threads_per_rank",
    "total_cells",
)

TIME_COLUMNS = (
    "elapsed_wall_seconds",
    "diffusion_seconds",
    "mechanics_seconds",
    "phenotype_seconds",
    "intracellular_seconds",
)

PERCENT_COLUMNS = (
    "diffusion_percent",
    "mechanics_percent",
    "phenotype_percent",
    "intracellular_percent",
)

COMPONENT_LABELS = ("Diffusion", "Mechanics", "Phenotype", "Intracellular")
NODE_FILE_PATTERN = re.compile(r"^(\d+)_nodes\.csv$")


class InputError(ValueError):
    """An actionable error in a CSV input or CLI path."""


@dataclass(frozen=True)
class Run:
    name: str
    path: Path
    nodes: int
    rows: tuple[dict[str, float], ...]

    @property
    def timing_rows(self) -> tuple[dict[str, float], ...]:
        return tuple(row for row in self.rows if row["elapsed_wall_seconds"] > 0)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare a single-node PhysiCell OMP CSV baseline with all "
            "*_nodes.csv PhysiCell-X results in a directory."
        )
    )
    parser.add_argument(
        "physicell_omp_results",
        type=Path,
        help="PhysiCell baseline CSV, or a directory containing exactly one CSV",
    )
    parser.add_argument(
        "physicellx_mpiomp_results_dir",
        type=Path,
        help="directory containing PhysiCell-X *_nodes.csv result files",
    )
    parser.add_argument(
        "output_dir", type=Path, help="directory in which plots and summary CSV are written"
    )
    return parser.parse_args(argv)


def resolve_baseline(path: Path) -> Path:
    if path.is_file():
        if path.suffix.lower() != ".csv":
            raise InputError(f"PhysiCell baseline is not a CSV file: {path}")
        return path
    if not path.exists():
        raise InputError(f"PhysiCell baseline path does not exist: {path}")
    if not path.is_dir():
        raise InputError(f"PhysiCell baseline path is neither a file nor a directory: {path}")

    candidates = sorted(candidate for candidate in path.glob("*.csv") if candidate.is_file())
    if len(candidates) != 1:
        listing = ", ".join(candidate.name for candidate in candidates) or "none"
        raise InputError(
            f"Expected exactly one CSV directly in baseline directory {path}; "
            f"found {len(candidates)} ({listing})"
        )
    return candidates[0]


def discover_physicellx(path: Path) -> list[tuple[int, Path]]:
    if not path.exists():
        raise InputError(f"PhysiCell-X results directory does not exist: {path}")
    if not path.is_dir():
        raise InputError(f"PhysiCell-X results path is not a directory: {path}")

    discovered: list[tuple[int, Path]] = []
    for candidate in sorted(path.rglob("*_nodes.csv")):
        match = NODE_FILE_PATTERN.fullmatch(candidate.name)
        if match:
            nodes = int(match.group(1))
            if nodes < 1:
                raise InputError(f"Node count must be positive in filename: {candidate}")
            discovered.append((nodes, candidate))

    if not discovered:
        raise InputError(f"No *_nodes.csv result files found under {path}")

    by_node: dict[int, list[Path]] = {}
    for nodes, candidate in discovered:
        by_node.setdefault(nodes, []).append(candidate)
    duplicates = {nodes: files for nodes, files in by_node.items() if len(files) > 1}
    if duplicates:
        details = "; ".join(
            f"{nodes} nodes: {', '.join(str(file) for file in files)}"
            for nodes, files in sorted(duplicates.items())
        )
        raise InputError(
            "Found multiple PhysiCell-X files for the same node count. "
            f"Pass one scenario directory or remove the ambiguity ({details})"
        )
    return sorted(discovered)


def read_run(path: Path, name: str, nodes: int) -> Run:
    try:
        handle = path.open(newline="", encoding="utf-8-sig")
    except OSError as error:
        raise InputError(f"Cannot open {path}: {error}") from error

    with handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise InputError(f"CSV has no header: {path}")
        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise InputError(f"{path}: missing required columns: {', '.join(missing)}")

        rows: list[dict[str, float]] = []
        seen_times: set[float] = set()
        for line_number, raw in enumerate(reader, start=2):
            row: dict[str, float] = {}
            for column in REQUIRED_COLUMNS:
                value = raw.get(column, "")
                try:
                    number = float(value)
                except (TypeError, ValueError) as error:
                    raise InputError(
                        f"{path}:{line_number}: {column} must be numeric, got {value!r}"
                    ) from error
                if not math.isfinite(number):
                    raise InputError(
                        f"{path}:{line_number}: {column} must be finite, got {value!r}"
                    )
                row[column] = number

            simulation_time = row["simulation_time"]
            if simulation_time in seen_times:
                raise InputError(
                    f"{path}:{line_number}: duplicate simulation_time {simulation_time:g}"
                )
            seen_times.add(simulation_time)

            nonnegative = TIME_COLUMNS + PERCENT_COLUMNS + ("memory_gb", "total_cells")
            negative = [column for column in nonnegative if row[column] < 0]
            if negative:
                raise InputError(
                    f"{path}:{line_number}: negative value in {', '.join(negative)}"
                )
            if row["mpi_ranks"] < 1 or row["omp_threads_per_rank"] < 1:
                raise InputError(
                    f"{path}:{line_number}: mpi_ranks and omp_threads_per_rank must be positive"
                )
            rows.append(row)

    if not rows:
        raise InputError(f"CSV has no data rows: {path}")
    if not any(row["elapsed_wall_seconds"] > 0 for row in rows):
        raise InputError(f"CSV has no positive elapsed_wall_seconds rows: {path}")
    rows.sort(key=lambda row: row["simulation_time"])
    return Run(name=name, path=path, nodes=nodes, rows=tuple(rows))


def mean(rows: Iterable[dict[str, float]], column: str) -> float:
    values = [row[column] for row in rows]
    if not values:
        raise InputError(f"Cannot calculate mean {column}: no applicable rows")
    return math.fsum(values) / len(values)


def common_time_speedups(baseline: Run, run: Run) -> list[tuple[float, float]]:
    baseline_by_time = {
        row["simulation_time"]: row["elapsed_wall_seconds"] for row in baseline.timing_rows
    }
    values = []
    for row in run.timing_rows:
        elapsed = baseline_by_time.get(row["simulation_time"])
        if elapsed is not None:
            values.append((row["simulation_time"], elapsed / row["elapsed_wall_seconds"]))
    return values


def summary_row(run: Run, baseline_mean_elapsed: float, implementation: str) -> dict[str, object]:
    timing_rows = run.timing_rows
    mean_elapsed = mean(timing_rows, "elapsed_wall_seconds")
    return {
        "implementation": implementation,
        "nodes": run.nodes,
        "source_csv": str(run.path.resolve()),
        "samples": len(run.rows),
        "timing_samples": len(timing_rows),
        "mean_elapsed_wall_seconds": mean_elapsed,
        "speedup_vs_physicell": baseline_mean_elapsed / mean_elapsed,
        "mean_diffusion_seconds": mean(timing_rows, "diffusion_seconds"),
        "mean_mechanics_seconds": mean(timing_rows, "mechanics_seconds"),
        "mean_phenotype_seconds": mean(timing_rows, "phenotype_seconds"),
        "mean_intracellular_seconds": mean(timing_rows, "intracellular_seconds"),
        "mean_diffusion_percent": mean(timing_rows, "diffusion_percent"),
        "mean_mechanics_percent": mean(timing_rows, "mechanics_percent"),
        "mean_phenotype_percent": mean(timing_rows, "phenotype_percent"),
        "mean_intracellular_percent": mean(timing_rows, "intracellular_percent"),
        "mean_memory_gb": mean(run.rows, "memory_gb"),
        "mean_total_cells": mean(run.rows, "total_cells"),
        "mpi_ranks": mean(run.rows, "mpi_ranks"),
        "omp_threads_per_rank": mean(run.rows, "omp_threads_per_rank"),
    }


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, output_dir: Path, stem: str) -> None:
    for suffix in ("png", "svg"):
        fig.savefig(output_dir / f"{stem}.{suffix}", dpi=180, bbox_inches="tight")


def plot_speedup_nodes(plt, output_dir: Path, x_rows: list[dict[str, object]]):
    fig, ax = plt.subplots(figsize=(7, 4.5), layout="constrained")
    nodes = [int(row["nodes"]) for row in x_rows]
    speedups = [float(row["speedup_vs_physicell"]) for row in x_rows]
    ax.plot(nodes, speedups, "o-", linewidth=2)
    ax.set(xlabel="PhysiCell-X nodes", ylabel="Speedup vs PhysiCell OMP", title="Scalability")
    ax.set_xticks(nodes)
    ax.grid(alpha=0.3)
    save_figure(fig, output_dir, "speedup_vs_nodes")
    return fig


def plot_speedup_time(plt, output_dir: Path, baseline: Run, x_runs: list[Run]):
    fig, ax = plt.subplots(figsize=(7, 4.5), layout="constrained")
    matched_counts: dict[int, int] = {}
    for run in x_runs:
        values = common_time_speedups(baseline, run)
        matched_counts[run.nodes] = len(values)
        if values:
            times, speedups = zip(*values)
            ax.plot(times, speedups, marker=".", label=f"{run.nodes} node(s)")
    if not any(matched_counts.values()):
        raise InputError("No shared positive-duration simulation_time values across the inputs")
    ax.set(
        xlabel="Simulation time",
        ylabel="Speedup vs PhysiCell OMP",
        title="Speedup over simulation time",
    )
    ax.grid(alpha=0.3)
    ax.legend()
    save_figure(fig, output_dir, "speedup_vs_time")
    return fig, matched_counts


def plot_breakdown(plt, output_dir: Path, all_runs: list[Run]):
    columns = 2
    rows = math.ceil(len(all_runs) / columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(9, 4 * rows), squeeze=False, layout="constrained"
    )
    colors = plt.get_cmap("tab20c").colors[::4][:4]
    for ax, run in zip(axes.flat, all_runs):
        percentages = [mean(run.timing_rows, column) for column in PERCENT_COLUMNS]
        ax.pie(
            percentages,
            autopct=lambda value: f"{value:.1f}%" if value >= 1 else "",
            startangle=90,
            colors=colors,
        )
        ax.set_title(run.name)
    for ax in list(axes.flat)[len(all_runs) :]:
        ax.set_visible(False)
    fig.legend(COMPONENT_LABELS, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncols=4)
    fig.suptitle("Mean runtime breakdown")
    save_figure(fig, output_dir, "runtime_breakdown")
    return fig


def plot_bars(
    plt,
    output_dir: Path,
    all_rows: list[dict[str, object]],
    column: str,
    ylabel: str,
    title: str,
    stem: str,
):
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    labels = [
        "PhysiCell\nbaseline"
        if row["implementation"] == "PhysiCell OMP"
        else f"PhysiCell-X\n{int(row['nodes'])} node(s)"
        for row in all_rows
    ]
    values = [float(row[column]) for row in all_rows]
    bars = ax.bar(labels, values, color=["#777777"] + ["#4472c4"] * (len(values) - 1))
    ax.bar_label(bars, fmt="%.3g", padding=3)
    ax.set(ylabel=ylabel, title=title)
    ax.grid(axis="y", alpha=0.3)
    save_figure(fig, output_dir, stem)
    return fig


def combined_figure(
    plt,
    output_dir: Path,
    baseline: Run,
    x_runs: list[Run],
    summaries: list[dict[str, object]],
    figures: list,
) -> None:
    """Create a vector-capable overview containing all requested chart types."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 15), layout="constrained")

    nodes = [run.nodes for run in x_runs]
    x_summaries = summaries[1:]
    axes[0, 0].plot(
        nodes, [float(row["speedup_vs_physicell"]) for row in x_summaries], "o-", linewidth=2
    )
    axes[0, 0].set(
        xlabel="PhysiCell-X nodes", ylabel="Speedup vs PhysiCell OMP", title="Scalability"
    )
    axes[0, 0].set_xticks(nodes)
    axes[0, 0].grid(alpha=0.3)

    for run in x_runs:
        values = common_time_speedups(baseline, run)
        if values:
            times, speedups = zip(*values)
            axes[0, 1].plot(times, speedups, marker=".", label=f"{run.nodes} node(s)")
    axes[0, 1].set(
        xlabel="Simulation time",
        ylabel="Speedup vs PhysiCell OMP",
        title="Speedup over simulation time",
    )
    axes[0, 1].grid(alpha=0.3)
    axes[0, 1].legend()

    breakdown_run = x_runs[-1]
    percentages = [mean(breakdown_run.timing_rows, column) for column in PERCENT_COLUMNS]
    colors = plt.get_cmap("tab20c").colors[::4][:4]
    axes[1, 0].pie(
        percentages,
        autopct=lambda value: f"{value:.1f}%" if value >= 1 else "",
        startangle=90,
        colors=colors,
    )
    axes[1, 0].set_title(f"Mean runtime breakdown: {breakdown_run.nodes} nodes")
    axes[1, 0].legend(
        COMPONENT_LABELS, loc="lower center", bbox_to_anchor=(0.5, -0.12), ncols=2
    )

    labels = ["PhysiCell"] + [f"X: {run.nodes} node(s)" for run in x_runs]
    bar_colors = ["#777777"] + ["#4472c4"] * len(x_runs)
    memory_bars = axes[1, 1].bar(
        labels, [float(row["mean_memory_gb"]) for row in summaries], color=bar_colors
    )
    axes[1, 1].bar_label(memory_bars, fmt="%.3g", padding=3)
    axes[1, 1].set(ylabel="Mean memory (GB)", title="Memory usage")
    axes[1, 1].grid(axis="y", alpha=0.3)

    wall_bars = axes[2, 0].bar(
        labels,
        [float(row["mean_elapsed_wall_seconds"]) for row in summaries],
        color=bar_colors,
    )
    axes[2, 0].bar_label(wall_bars, fmt="%.3g", padding=3)
    axes[2, 0].set(ylabel="Mean elapsed wall time (s)", title="Wall-clock time")
    axes[2, 0].grid(axis="y", alpha=0.3)

    axes.flat[-1].set_visible(False)
    fig.suptitle("PhysiCell / PhysiCell-X scalability comparison", fontsize=16)
    save_figure(fig, output_dir, "scalability_comparison")
    plt.close(fig)
    for source in figures:
        plt.close(source)


def run(args: argparse.Namespace) -> list[Path]:
    baseline_path = resolve_baseline(args.physicell_omp_results)
    x_files = discover_physicellx(args.physicellx_mpiomp_results_dir)
    baseline = read_run(baseline_path, "PhysiCell OMP baseline", 1)
    x_runs = [
        read_run(path, f"PhysiCell-X, {nodes} node(s)", nodes) for nodes, path in x_files
    ]

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    if not output_dir.is_dir():
        raise InputError(f"Output path is not a directory: {output_dir}")

    baseline_mean_elapsed = mean(baseline.timing_rows, "elapsed_wall_seconds")
    baseline_summary = summary_row(baseline, baseline_mean_elapsed, "PhysiCell OMP")
    x_summaries = [
        summary_row(run, baseline_mean_elapsed, "PhysiCell-X MPI+OMP") for run in x_runs
    ]
    summaries = [baseline_summary, *x_summaries]

    # Import plotting only after inputs have been validated, and use a writable
    # cache location on clusters where the user's config directory is read-only.
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/physicell-scalability-matplotlib")
    try:
        warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise InputError("matplotlib is required to generate scalability plots") from error

    figures = []
    figures.append(plot_speedup_nodes(plt, output_dir, x_summaries))
    speedup_time_figure, matched_counts = plot_speedup_time(plt, output_dir, baseline, x_runs)
    figures.append(speedup_time_figure)
    figures.append(plot_breakdown(plt, output_dir, [baseline, *x_runs]))
    figures.append(
        plot_bars(
            plt,
            output_dir,
            summaries,
            "mean_memory_gb",
            "Mean memory (GB)",
            "Memory usage",
            "memory_usage",
        )
    )
    figures.append(
        plot_bars(
            plt,
            output_dir,
            summaries,
            "mean_elapsed_wall_seconds",
            "Mean elapsed wall time (s)",
            "Wall-clock time",
            "wall_time_vs_nodes",
        )
    )

    for row in x_summaries:
        row["matched_time_samples"] = matched_counts[int(row["nodes"])]
    baseline_summary["matched_time_samples"] = len(baseline.timing_rows)
    summary_path = output_dir / "scalability_summary.csv"
    # Ensure every row has the final field before DictWriter fixes its schema.
    ordered_summaries = [baseline_summary, *x_summaries]
    write_summary(summary_path, ordered_summaries)

    plot_names = [
        "speedup_vs_nodes",
        "speedup_vs_time",
        "runtime_breakdown",
        "memory_usage",
        "wall_time_vs_nodes",
    ]
    combined_figure(plt, output_dir, baseline, x_runs, ordered_summaries, figures)

    return sorted(output_dir / f"{name}.{suffix}" for name in [*plot_names, "scalability_comparison"] for suffix in ("png", "svg")) + [summary_path]


def main(argv: Sequence[str] | None = None) -> int:
    try:
        outputs = run(parse_args(argv))
    except (InputError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print("Wrote scalability comparison artifacts:")
    for output in outputs:
        print(f"  {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
