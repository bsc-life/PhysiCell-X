# PhysiCell-X MPI validation suite

Build and run every correctness test from the repository root with:

```bash
tests/mpi/run_all_tests.sh
```

The runner creates an isolated results directory under `/tmp` and prints its
location. It stops at the first failure. The `template` directory is not run
because it is development scaffolding rather than a test.

Console output contains only test metadata: test name and pass/fail status.
Complete output from each executable is stored under `logs/` in the
results directory. If a test fails, the runner prints the corresponding log
path.

Useful overrides:

```bash
JOBS=8 MPI_RANKS=2 tests/mpi/run_all_tests.sh
tests/mpi/run_all_tests.sh build
RESULTS_DIR=/path/to/results tests/mpi/run_all_tests.sh run
```

The suite checks:

- XML cell definitions and asymmetric-division probabilities on every rank;
- preservation of cell identity, kinematics, phenotype state, and custom data
  across `Cell::pack` / `unpack` communication;
- an analytical constant steady-state diffusion invariant across MPI boundaries;
- deterministic serial-versus-MPI mechanics trajectories;
- one-rank-versus-multi-rank MultiCellDS cell and field data;
- one-rank-versus-multi-rank deterministic phenotype integration.

The phenotype comparison reports paper-ready numerical equivalence metrics. It
requires identical cell IDs, types, viability, oncoprotein, and volume; a cell
position maximum error below 1 micrometre and RMSE below 0.25 micrometre; field
relative L2 error below 2%; internalized-substrate relative maximum error below
1%; and pressure maximum error below 0.1. Random cell orientation is excluded
because motility is disabled in this deterministic validation scenario.

The launch command defaults to `mpiexec`. Set `MPIEXEC` to another compatible
launcher when needed. Existing Slurm launch scripts remain available for large
benchmarks but are not used by this short correctness suite.
