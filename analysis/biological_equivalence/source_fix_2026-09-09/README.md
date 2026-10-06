# Deterministic substrate application reassessment

Work started 9 September 2026; the final coupled experiment started 10 September after an interrupted session. Production changes are applied to both PhysiCell-X and the sibling PhysiCelldFBA working tree. Original simulation settings, objects, executables and output directories were preserved. All executables and runs used here are isolated in this directory.

## Changes

`BioFVM/BioFVM_basic_agent.cpp` rebuilds source constants on every active secretion call using its actual `dt`. Changes in rates and timestep therefore take effect even when volume is unchanged. This corrects a cached 0.01-minute FBA increment being reused at a 0.001-minute secretion timestep.

`core/PhysiCell_cell_sources.h` groups eligible cells by microenvironment and substrate voxel, sorts each group by cell ID, and parallelizes across groups. Each voxel has one writer within its owning process. The complete implicit secretion/uptake and explicit net-export operation runs sequentially within a group. Both Cell_Container source paths in PhysiCell-X and the reference source path use this helper. The Makefiles include the new header dependency.

FBA optimization, update triggers, update ordering, diffusion algorithms and Dirichlet rules were not changed. The existing FBA-versus-secretion ordering differs between projects. This experiment tests source application without treating FBA solve frequency as the established cause of the historical mismatch.

`x.patch` and `reference.patch` describe this task's changes relative to the pre-task working trees. Previously present changes are preserved. The matching `*_final_changes.json` manifests hash the changed files. `patch_checks.json` records successful byte-identical reversal and reapplication in temporary copies. Initial patches/manifests are retained for the audit trail; the final staged headers include the corrected volume accessor.

## Validation and run conditions

`test_cell_sources.cpp` uses 64 cells sharing four voxels. Both fresh project builds must reproduce a sequential source update bit for bit with 1, 2 and 8 OpenMP threads, reversed input ordering and 20 repetitions per thread count. It tests a deliberately stale timestep cache, a changed net-export rate, an actual timestep change, implicit uptake/secretion, and internal/external accounting. Results: `test_results.json`.

The initial reference fixture failure was in uninitialized RNG setup during Cell construction, confirmed by a debugger backtrace. The final fixture calls SeedRandom before constructing cells. It was not a failure of the production colony or the source helper. The first absolute-path MPI smoke run omitted cell MAT files because the existing writer prepends `./`; the final smoke/coupled configurations use `..` relative to each output/work directory. This output-path workaround does not change simulation code. `smoke_final_comparison` contains the complete paired smoke comparison.

The main executables were built from 41 PhysiCell-X and 33 reference objects with matching optimization flags and each project's bundled Clp/CoinUtils/SBML libraries. The interrupted final test link was recovered with `build_case.py --recover-tests`. Linking now uses a temporary local file before copying the complete executable into the analysis build folder. `compile_initial.json`, `compile.json`, `provenance_initial.json` and `provenance.json` retain compilation commands and source/binary hashes. Successful final smoke runs validate the main executables whose final link preceded that interruption.

Both coupled runs use two OpenMP threads per process, seed 50, the same 250 × 250 × 4 substrate mesh, matching boundaries, a 0.001-minute diffusion/secretion step and 0.01-minute FBA interval. PhysiCell-X uses two MPI ranks; the reference uses one process. Thus OpenMP counts match per process, not total allocated cores. Initial cell positions, volumes and substrate fields match; cell orientations differ because a shared seed does not imply identical RNG streams in the two implementations. The two metabolic-model XML files have identical bytes.

The experiment runs through 180 minutes with full snapshots every 30 minutes. `runs.json`, the generated XML files and each run's `run.json` record settings, commands, times and completion checks. The runner checks for complete referenced cell and substrate data files, not just a successful exit code.

## Source accounting

Set `PHYSICELL_SUBSTRATE_BUDGET_PREFIX=/path/prefix` to enable per-rank `prefix.rankN.csv` logs. Logging is disabled by default. `PHYSICELL_SUBSTRATE_BUDGET_EVERY` controls the number of source steps per flush; coupled runs use 1000, or one simulated minute.

- `requested_net_export`: sum of actual secretion `dt × net_export_rate` over eligible cells.
- `requested_total`: predicted extracellular change, including implicit secretion/uptake.
- `applied_total`: measured concentration change across the source call, multiplied by voxel volume.
- `max_abs_step_error`: largest absolute difference between the rank's summed requested-total and applied-total amounts on an individual step in the interval.

These are native concentration × voxel-volume amount units, without an additional SI conversion. In this colony conventional secretion and uptake coefficients are zero, so the independent net-export sum must equal the measured extracellular source change. This is an accounting check before the next diffusion step; it is not a balance of the final field after Dirichlet boundary exchange. The main loop performs one extra tick after the last numbered snapshot, which the budget summarizer excludes from the requested 0–180-minute totals. Spatial cancellation is possible in aggregate residuals; deterministic per-voxel tests provide the complementary check.

## Commands used

From the PhysiCell-X project root:

```bash
python3 analysis/biological_equivalence/source_fix_2026-09-09/build_case.py x
python3 analysis/biological_equivalence/source_fix_2026-09-09/build_case.py reference
python3 analysis/biological_equivalence/source_fix_2026-09-09/run_source_tests.py
python3 analysis/biological_equivalence/source_fix_2026-09-09/prepare_runs.py
python3 analysis/biological_equivalence/source_fix_2026-09-09/run_cases.py smoke_final
python3 analysis/biological_equivalence/source_fix_2026-09-09/run_cases.py coupled
python3 analysis/biological_equivalence/source_fix_2026-09-09/analyse_budgets.py coupled
OPENBLAS_NUM_THREADS=1 python3 tools/compare_biological_outputs.py \
  analysis/biological_equivalence/source_fix_2026-09-09/runs/x_coupled \
  analysis/biological_equivalence/source_fix_2026-09-09/runs/reference_coupled \
  --out analysis/biological_equivalence/source_fix_2026-09-09/comparison
```

MPI needs permission to open local sockets. A new experiment should use fresh output directories; `run_cases.py` refuses to overwrite a completed run. `check_progress.py` reads flushed intervals while a coupled run is active.

Historical context and the unchanged-transport control are in `../reassessment_2026-09-09/report.md` and `../reassessment_2026-09-09/transport_diagnosis/report.md`. The latter found bitwise-identical serial and two-rank diffusion for matched input and boundaries. The new combined intervention also matches thread counts; it cannot separately assign the historical mismatch to the race, stale timestep or thread-count change. A single pair of stochastic colony runs does not establish statistical biological equivalence.

`replay_late_transport.py prepare` and `replay_late_transport.py run` add a matched-input control from the fresh reference 120-minute snapshot, using the coupled experiment's production object files. The second command needs MPI sockets. Results, commands, fixtures and hashes are in `transport_replay_120`. All four comparisons (X serial or two ranks versus reference serial, with or without fixed sources) were bitwise identical after 100 steps. No production source or original output is changed by this replay.

`replay_source_layout.py` uses the reference solver and a common initial field to compare the two saved 120-minute source layouts. One pair of fixtures retains the saved rates; another uses the same net-export vector per unit cell volume for every cell in both colonies. The latter reproduces the coupled spatial difference with 0.09–0.35% residual L2 error, normalized by the observed field difference itself. This isolates the dominant effect of differing source positions/volume distributions without introducing differing FBA rates or transport operators. Results and exact calls are saved as `transport_replay_120/layout_summary.json` and `layout_calls.json`.

The layout script accepts a snapshot minute, for example `replay_source_layout.py 150`, and stores each follow-up in its own `transport_replay_<minute>` folder. Positions outside the saved mesh are excluded from these fixtures, and exclusion counts are recorded. `--reuse-results` regenerates the analysis/figure from existing replay outputs without repeating transport. The follow-up distinguishes finite-domain source eligibility from substrate Dirichlet clamping; those are separate conditions.
