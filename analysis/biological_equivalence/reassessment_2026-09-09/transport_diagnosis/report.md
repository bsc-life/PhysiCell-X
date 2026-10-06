**The controlled tests point away from the diffusion algorithm and mismatched Dirichlet conditions. They identify defects in applying cellular substrate fluxes, rather than a difference in the FBA optimization itself. However, they do not yet establish which defect accounts for the entire historical 18–24% discrepancy.**

The similar curves in `substrates.png` are real: the implementations predict similar mean concentrations near the colony. The previous report's larger percentages describe depletion or accumulation relative to the initial field, not the percentage difference in oxygen or glucose concentration.

| Mean concentration within the plotted 25 µm radius | Time | PhysiCell-X | Reference | Difference relative to reference concentration |
|---|---:|---:|---:|---:|
| Oxygen | 180 min | 23.39884 | 23.65496 | −1.08% |
| Glucose | 180 min | 25.08416 | 25.46564 | −1.50% |
| Acetate | 180 min | 0.31256 | 0.26288 | +18.90% |
| CO₂ | 180 min | 1.95784 | 1.64363 | +19.12% |
| Oxygen | 330 min | 14.81678 | 15.13509 | −2.10% |
| Glucose | 330 min | 11.74561 | 12.22715 | −3.94% |
| Acetate | 330 min | 2.02558 | 1.96379 | +3.15% |
| CO₂ | 330 min | 12.52188 | 12.13304 | +3.20% |

At 30 minutes, oxygen/glucose concentration differences in this region are just −0.150%/−0.210%, even though whole-domain depletion differs by 23–24%. The large initial backgrounds, averaging over space, and plotting the entire growth trajectory make these discrepancies look small. For products, the approximately 19% difference at 180 minutes occupies a small part of the vertical range needed for the later values. [Relative-difference plots](substrate_relative_differences.png) show both concentration and depletion/accumulation denominators explicitly. The small effects on recorded growth in this nutrient-rich experiment are consistent with the similar concentration trajectories; that does not establish equal substrate accounting.

| Candidate explanation | Evidence from this investigation | Assessment |
|---|---|---|
| Diffusion discretization or MPI diffusion exchange | Serial reference, PhysiCell-X serial, and PhysiCell-X two-rank blocking solver give bit-for-bit identical fields after 100 steps from identical inputs, both without cells and with identical fixed cell sources. | No discrepancy in this controlled configuration. |
| Different configured Dirichlet faces/values | Saved settings agree; measured outer x/y faces equal their prescribed values in both 30-minute snapshots, and z faces are not clamped. Matching boundaries also produce identical diffusion results. | No evidence that a boundary-condition mismatch causes these outputs. |
| Different FBA optimization or exchange conversion | Current optimizer implementation differs only in blank lines; model XML and exchange-conversion formulas match. At 30 minutes, active mean stored export-rate magnitudes differ by at most 0.00384%. | Not supported as the explanation for the 23–24% early field difference. |
| Concurrent cell-to-voxel source updates | Both implementations write shared voxel densities without synchronization. Production-method tests lose flux with multiple threads; one-thread and locked controls conserve the prescribed flux. | Confirmed defect and a candidate contributor to the difference between runs. Its historical contribution is not quantified. |
| FBA-to-BioFVM cached timestep | Both cache a 0.01-minute increment which can be reused by a 0.001-minute secretion call. A production-method probe reproduces a tenfold increment. | Confirmed shared coupling defect; does not by itself explain the difference between projects. |

The diffusion experiments use the actual 250 × 250 × 4 mesh, 2 µm spacing, five diffusion coefficients, zero decay, and 0.001-minute step. All implementations start from the **same reference 30-minute concentration field**. The source replay also uses the same 211 reference cells, positions, volumes and saved export rates, with growth, mechanics and FBA disabled. Sources are applied with a single thread and the correct 0.001-minute increment, so this control deliberately removes both source-update defects. It runs 100 steps, or 0.1 minute, and compares every voxel after coordinate alignment, not just spatial means.

All four comparisons—PhysiCell-X serial/two-rank, with/without fixed sources—have **zero maximum absolute concentration difference for all five substrates**. The replay checks the cell-to-voxel mapping across x = 0 as well as the diffusion step. Results and raw-array hashes are in [diffusion_summary.json](diffusion_summary.json); the raw replay fields are stored in `/tmp` and can be regenerated. This is strong evidence against a diffusion implementation difference for the tested mesh and state, not certification of every solver variant, geometry or long trajectory.

The tested MPI solver is `diffusion_decay_solver__constant_coefficients_LOD_3D_BLOCKING`: [initialize_microenvironment](../../../../BioFVM/BioFVM_microenvironment.cpp) selects it for this 3-D colony, overriding the constructor's initial AVX solver pointer. The AVX alternative was not used as the basis of the conclusion. Dirichlet masks in the probe explicitly reproduce the saved face settings: x/y boundaries fixed at `[25, 27.5, 0, 0, 0]`, no imposed z-face values and no internal Dirichlet nodes. This tests boundary application under matching masks; it does not recover a historical executable's runtime masks. Independent [snapshot face measurements](boundary30.json) and source inspection support that the saved runs used those faces.

The concurrent source defect is in [PhysiCell-X Basic_Agent::simulate_secretion_and_uptake](../../../../BioFVM/BioFVM_basic_agent.cpp) and its [reference counterpart](../../../../../PhysiCelldFBA/BioFVM/BioFVM_basic_agent.cpp). The cell-container loops invoke it with `#pragma omp parallel for`. Multiple cells can perform read–modify–write operations on the same voxel's density. There is no atomic operation, lock, reduction or ownership by voxel. Even the zero-secretion and unit-denominator operations can overwrite another thread's export update; protecting only the final net-export addition is not a complete general repair.

In the actual 30-minute fixture, **211 cells occupy 120 substrate voxels; 64 voxels contain multiple cells, with up to three cells per voxel**. I ran 100 independent one-step trials with the production source routine, fixed rates, no diffusion, and 1, 2, 8 and 112 OpenMP threads. Examples from this machine:

| Test | Mean oxygen flux lost | Largest one-step oxygen flux loss |
|---|---:|---:|
| Reference, 1 thread | Roundoff only | Roundoff only |
| Reference, 2 threads | 0.682% | 4.739% |
| Reference, 8 threads | 0.370% | 2.370% |
| Reference, 112 threads | 0.019% | 0.474% |
| PhysiCell-X, 2 threads, all fixture cells | 0.180% | 3.791% |
| PhysiCell-X, 2 threads, left/right subsets tested separately | 0.028% / 0.078% | 0.926% / 0.971% |
| Either implementation, 112 threads with a lock around each complete cell update | Roundoff only | Roundoff only |

See [all substrate summaries](source_summary.json), [commands and environment](source_runs.json), and individual trial CSVs. **These measurements do not reproduce the full 18–24% historical gap.** This machine exposes only 8 CPUs, so 112 threads are oversubscribed; the tests use passive waiting and freshly compiled code. Race losses depend on scheduling, affinity, compiler, ordering and simultaneous execution, and need not increase monotonically with requested thread count. The historical reference used 112 threads versus 2 per rank for PhysiCell-X. It would be unjustified to extrapolate a numerical correction from these local losses, or to call the entire mismatch proven to be a reference-side race. The reference itself is not a reliable numerical ground truth until this defect is controlled.

The second defect is a **timestep mismatch in the shared dFBA/BioFVM interface**. In both projects, `update_dfba_outputs` calls `set_internal_uptake_constants(dfba_time_step)` with 0.01 minutes. That function stores the already integrated increment `net_export_rate × cached_dt / voxel_volume` and clears `volume_is_changed`. Later, `simulate_secretion_and_uptake(..., 0.001)` rebuilds the constants only if that flag is set; otherwise it adds the existing increment without rescaling it to the passed timestep. The colony's null volume callbacks do not guarantee a rebuild before secretion.

A probe using freshly compiled production `Basic_Agent` methods confirms the consequence in **both** implementations:

| Cached timestep | Secretion call timestep | Expected concentration increment | Actual increment |
|---|---:|---:|---:|
| 0.01 min | 0.001 min | 0.001 | 0.01 |
| 0.001 min | 0.001 min | 0.001 | 0.001 |

The fixture prescribes export rate 8 into a voxel of volume 8. It reproduces the method sequence without running the FBA optimizer; supporting BioFVM objects come from the existing builds. [Cache probe commands/results](cache_probe.json) and [source/object provenance](source_provenance.json) document that scope. This demonstrates a tenfold per-call error **when that cache state persists**, not that every historical concentration is exactly ten times too large. Division, volume changes or other cache rebuilds can alter when it occurs. Since both projects contain the issue, it can corrupt absolute accounting while their concentration curves remain similar.

There are also scheduling differences: PhysiCell-X updates intracellular state before secretion, whereas the reference secretes first; its FBA update trigger and next-update scheduling differ from the reference. The current FBA header also lacks the reference's optional custom-optimization dispatch, but the compared colony does not install such a callback. Near-identical pre-division volume and stored fluxes argue against a large difference in FBA update counts as the early cause. These controls do not reconstruct all historical source increments, cache states or MPI migration histories.

The most useful next step is to make substrate application deterministic and timestep-consistent in **both** projects: assign one writer per substrate voxel (or serialize updates as a diagnostic), and rebuild or rescale source increments using the actual secretion timestep. Then repeat the coupled colony with matching thread counts, logging the summed requested and applied substrate amounts before diffusion. If a discrepancy remains, a matched-input diffusion replay already provides a working control for separating that discrepancy from transport. No production simulation source, settings or original outputs were changed during this diagnosis.

Reproduce the focused experiments from the project root:

```bash
python3 analysis/biological_equivalence/reassessment_2026-09-09/reproduce_transport_probes.py
```

The runner builds temporary executables, reconstructs the common field/cell fixture, runs source tests and 100-step diffusion replays including two MPI ranks, and regenerates the comparison summaries and residual plot. MPI must be permitted to open local sockets. The cache probe additionally needs the existing supporting BioFVM objects in both projects. Source tests are intentionally nondeterministic and their loss percentages can change on repetition. These are focused diagnostic experiments, not a fresh 330-minute colony run or a statistical equivalence test.
