# Comparing simulation outputs

```bash
python3 tools/compare_simulations.py output ../PhysiCell/output \
  --out analysis/comparison_30um \
  --labels PhysiCell-X PhysiCell
```

Dependencies: Python 3.9+, NumPy, SciPy, Matplotlib. No GUI is required. Input
folders are read only; the report directory must be outside them and empty or
new. The script also works with older XML files whose `simplified_data` source is
`BioFVM`. MAT row labels, substrate names and voxel coordinates determine the
mapping; array order and cell IDs are never assumed to correspond between runs.

## Generated files

- `timeline_mechanics.png`: radius, vertical extent, volume, pressure, nearest
  neighbor distance, and overlap over time.
- `timeline_phenotype.png`: total/live/dead-state populations, cells with positive
  cycle or growth rate, arrested cells, mean live rate, dead volume, oncoprotein.
- `timeline_diffusion_<substrate>.png`: domain/core/cell-position concentrations,
  minimum, maximum and negative-voxel counts, for every common substrate.
- `timeline_errors_<substrate>.png`: relative field L2 error, maximum local error,
  and relative gradient L2 error.
- `spatial_diffusion_<substrate>_{xy,xz}.png`: paired concentration slices and
  signed difference maps, with common color scales across runs and times.
- `profiles_diffusion_<substrate>.png`: radial and x profiles at selected times.
  X profiles average over a cylinder around the specified center's x axis.
- `spatial_cells_*.png`: pressure, death and growth/cycle state in a thin z slice.
- `radial_evolution_*.png`: time-versus-radius maps of pressure, dead fraction,
  and live growth/cycle rate. Empty bins remain blank.
- CSV files: all timeline measurements, absolute/relative cell metric differences,
  and radial cell profiles. Missing values and undefined relative errors are `nan`.
- `report.md`: final comparison table and limitations.
- `provenance.json`: matched XML files/times, input sizes/modification times,
  configuration hashes, recorded seeds, initial-state checks, and warnings.

## Spatial controls

```bash
python3 tools/compare_simulations.py output ../PhysiCell/output \
  --out analysis/comparison_spatial \
  --center 0 0 0 --core-radius 100 --spatial-radius 300 \
  --interface-x 0 --slice-z 0 --slice-width 20 \
  --spatial-times 0 60 480 1000.01
```

Default spatial snapshots are the first, second, middle and last matched times.
Requested spatial times select the nearest matched snapshot; actual times are
printed on the graphs and recorded in provenance. Every matched time contributes
to the timelines and radial evolution maps.

The default center is `(0,0,0)`, held fixed between runs. The default spatial radius
is 1.2 times the first matched reference population's largest center-to-cell
distance; the core radius defaults to one third of that plot radius. For an
off-center colony, provide `--center`. Radius uses the 3D distance of cell centers,
not the physical outer cell surface. `--interface-x` only marks a diagnostic
position; it does not assert that an MPI boundary exists there.

For a metabolic model with a saved scalar `growth_rate` label:

```bash
python3 tools/compare_simulations.py output reference_output \
  --out analysis/metabolism --growth-field growth_rate --growth-threshold 1e-8
```

The threshold is in the selected field's XML units. By default, "growing" means
a live cell with positive `current_cycle_phase_exit_rate`; this is a proliferation
proxy, not a measured biomass flux or proof that division occurred. "Arrested"
means live and at/below the threshold. Alternative substrate consumption cannot
be inferred from concentrations alone. Missing model-specific labels produce
unavailable metrics rather than fabricated zeros.

## Definitions and limits

- Counts use the complete MAT population, while cell spatial maps explicitly use
  a slice of thickness `--slice-width`. Apoptosis and necrosis use the standard
  PhysiCell phase codes 100 and 101–103 together with the saved dead flag.
- Death counts are cells currently present, not cumulative death events. Total
  dead-cell volume helps assess swelling, but does not establish its cause.
- Nearest-neighbor overlap is `max(0, radius_i + radius_neighbor - distance)` in
  microns, averaged once per cell. It is not an all-pair overlap statistic.
- Domain means and L2 errors are voxel-volume weighted. Radial averages are
  voxel/cell means. Cell exposure samples the nearest geometric diffusion voxel
  using saved cell positions; it cannot reconstruct incorrect runtime ownership
  or stale cached indices. It is not necessarily the concentration the code used.
- Gradients use finite differences on the saved Cartesian grid, skipping singleton
  axes (so planar/2D simulations are supported). No interpolated mesh comparison
  or resampling in time is performed.
- `initial.xml`, numbered snapshots and `final.xml` are matched by XML simulated
  time, with `--time-tolerance` (default 0.0001 min). Same-time aliases must have
  identical labelled data. Snapshots older than the current initial snapshot are
  excluded with a warning. Reversed timestamps and conflicting duplicate times
  are rejected. This catches common stale-output cases but cannot prove a run's
  source revision or distinguish every possible mixed run without a run manifest.
- This is descriptive analysis of one run per simulator. It does not claim
  statistical equivalence, infer code bugs automatically, reconstruct cumulative
  divisions/deaths, or compute actual metabolic exchange fluxes from rate constants.
  Use matched initial conditions, replicate seeds, and explicit biological
  tolerances to establish equivalence.

Tests: `python3 -m unittest discover -s tests/analysis -p 'test_compare_simulations.py'`.
