Run the local voxel-index regression from the repository root:

```sh
python3 tests/unit/test_voxel_local_index.py
```

The runner requires `mpic++`. It compiles the production mesh and vector sources
into a temporary executable, then calls `nearest_voxel_local_index` directly.
It does not initialize MPI, launch simulations, or replace the simulator's object
files or executable.

Coverage includes every voxel centre in the 250 × 250 × 4 run geometry, one-,
two-, and three-rank layouts (including uneven partitions), singleton local
axes, exact voxel faces, local domain boundaries, and lower/upper clamping on
all three axes. Lower-bound undershoots reproduce the signed/unsigned comparison
bug. Failures return a nonzero exit code, independently of assertion settings.
