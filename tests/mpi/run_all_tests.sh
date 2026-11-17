#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

JOBS=${JOBS:-2}
MPI_RANKS=${MPI_RANKS:-2}
CELL_SND_COUNT=${CELL_SND_COUNT:-128}
DIFFUSION_STEPS=${DIFFUSION_STEPS:-20}
MPIEXEC=${MPIEXEC:-mpiexec}
MODE=${1:-all}

if [[ "$MODE" != "all" && "$MODE" != "build" && "$MODE" != "run" ]]; then
    echo "Usage: $0 [all|build|run]" >&2
    exit 2
fi
if (( MPI_RANKS < 2 )); then
    echo "MPI_RANKS must be at least 2" >&2
    exit 2
fi

RESULTS_DIR=${RESULTS_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/physicell-x-mpi-tests.XXXXXX")}
mkdir -p "$RESULTS_DIR"
LOG_DIR="$RESULTS_DIR/logs"
mkdir -p "$LOG_DIR"

MPI_ARGS=()
if "$MPIEXEC" --version 2>&1 | grep -qi "Open MPI"; then
    MPI_ARGS+=(--oversubscribe)
    if [[ $(id -u) -eq 0 ]]; then
        MPI_ARGS+=(--allow-run-as-root)
    fi
fi

run_mpi()
{
    local ranks=$1
    shift
    "$MPIEXEC" "${MPI_ARGS[@]}" -n "$ranks" "$@"
}

run_logged()
{
    local phase=$1
    local name=$2
    local slug=$3
    shift 3

    local log_file="$LOG_DIR/$slug.log"
    local status=0

    printf "[%-5s] %s ... " "$phase" "$name"
    set +e
    ( set -e; "$@" ) >"$log_file" 2>&1
    status=$?
    set -e

    if (( status == 0 )); then
        printf "PASS\n"
        return 0
    fi

    printf "FAIL\n"
    printf "        Log: %s\n" "$log_file"
    return "$status"
}

build_directory()
{
    make -C "$SCRIPT_DIR/$1" -j "$JOBS" all
}

build_cell_snd_validation()
{
    make -C "$SCRIPT_DIR/cell_snd" -j "$JOBS" validate
}

build_suite()
{
    local directories=(MultiCellDS cell_snd cell_xml diffusion mechanics phenotype)
    for directory in "${directories[@]}"; do
        run_logged BUILD "$directory" "build-$directory" build_directory "$directory"
    done
    run_logged BUILD "cell_snd validation executable" "build-cell-snd-validation" \
        build_cell_snd_validation
}

test_xml_definitions()
{
    (
        cd "$SCRIPT_DIR/cell_xml"
        run_mpi "$MPI_RANKS" ./cellxml ./config/Asymmetric_division.xml --validate
    )
}

test_cell_serialization()
{
    (
        cd "$SCRIPT_DIR/cell_snd"
        OMP_NUM_THREADS=2 run_mpi 2 ./cell_snd_val "$CELL_SND_COUNT"
    )
}

test_diffusion()
{
    (
        cd "$SCRIPT_DIR/diffusion"
        run_mpi "$MPI_RANKS" ./diffusion ./config/PhysiCell_settings-serial.xml \
            --validate "$DIFFUSION_STEPS"
    )
}

test_mechanics()
{
    local mechanics_dir="$RESULTS_DIR/mechanics"
    mkdir -p "$mechanics_dir"
    (
        cd "$SCRIPT_DIR/mechanics"
        ./mechanics_ref ./config/validation.xml "$mechanics_dir/reference.csv"
        run_mpi "$MPI_RANKS" ./mechanics_val ./config/validation.xml "$mechanics_dir/mpi.csv"
    )
    python3 "$SCRIPT_DIR/mechanics/compare_snapshots.py" \
        "$mechanics_dir/reference.csv" "$mechanics_dir/mpi.csv" "$mechanics_dir/report.txt"
}

test_multicellds_output()
{
    local multids_ref="$RESULTS_DIR/multicellds/reference"
    local multids_mpi="$RESULTS_DIR/multicellds/mpi"
    mkdir -p "$multids_ref" "$multids_mpi"
    (
        cd "$SCRIPT_DIR/MultiCellDS"
        run_mpi 1 ./multids ./config/PhysiCell_settings.xml "$multids_ref" --validate
        run_mpi "$MPI_RANKS" ./multids ./config/PhysiCell_settings.xml "$multids_mpi" --validate
    )
    python3 "$SCRIPT_DIR/validate_mcds.py" \
        "$multids_ref/initial.xml" "$multids_mpi/initial.xml"
}

test_phenotype_integration()
{
    local phenotype_ref="$RESULTS_DIR/phenotype/reference"
    local phenotype_mpi="$RESULTS_DIR/phenotype/mpi"
    mkdir -p "$phenotype_ref" "$phenotype_mpi"
    (
        cd "$SCRIPT_DIR/phenotype"
        run_mpi 1 ./mechanics ./config/PhysiCell_settings.xml "$phenotype_ref" --validate
        run_mpi "$MPI_RANKS" ./mechanics ./config/PhysiCell_settings.xml "$phenotype_mpi" --validate
    )
    python3 "$SCRIPT_DIR/validate_phenotype.py" \
        "$phenotype_ref/final.xml" "$phenotype_mpi/final.xml"
}

run_suite()
{
    run_logged TEST "XML cell definitions" "test-xml-definitions" test_xml_definitions
    run_logged TEST "cell serialization semantics" "test-cell-serialization" \
        test_cell_serialization
    run_logged TEST "analytical MPI diffusion" "test-diffusion" test_diffusion
    run_logged TEST "serial versus MPI mechanics" "test-mechanics" test_mechanics
    run_logged TEST "one-rank versus MPI MultiCellDS output" "test-multicellds" \
        test_multicellds_output
    run_logged TEST "one-rank versus MPI phenotype integration" "test-phenotype" \
        test_phenotype_integration
}

if [[ "$MODE" == "all" || "$MODE" == "build" ]]; then
    build_suite
fi
if [[ "$MODE" == "all" || "$MODE" == "run" ]]; then
    run_suite
fi

if [[ "$MODE" == "build" ]]; then
    echo "All MPI test executables built successfully."
else
    echo "All requested MPI tests passed. Results: $RESULTS_DIR"
fi
