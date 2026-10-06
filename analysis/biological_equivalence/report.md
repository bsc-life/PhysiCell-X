The available outputs do **not support overall biological equivalence** between PhysiCell-X and PhysiCelldFBA. Growth and the recorded metabolic phenotype agree closely over the shared interval, but cell packing and pressure differ substantially, and the coupled substrate fields show systematic differences. These are descriptive results from one run per project, not a statistical equivalence test.

Compared on 2026-09-08: `PhysiCell-X/output` against `../PhysiCelldFBA/output` (the actual directory spelling). Seven paired full snapshots cover 0, 30, 60, 90, 120, 150 and 180 simulated minutes. PhysiCelldFBA has additional numbered snapshots through 510 minutes; PhysiCell-X has no matching full snapshots beyond 180 minutes in this analysis. The reference `final.*` files are dated July 9 and belong to an older run; they were excluded. Neither compared series establishes the configured 600-minute outcome.

| Component | Assessment from these outputs |
|---|---|
| Diffusion and substrate fields | Similar qualitative consumption/production, but appreciable differences in depletion, accumulation and spatial gradients. The diffusion operator alone has not been isolated. |
| Cell mechanics | Strong disagreement in relaxation, overlap and pressure, already present before division. The spatial pattern suggests a problem in interactions across the MPI partition. |
| Phenotype | Close agreement in population growth, total cell volume and recorded dFBA growth/exchange fluxes. Death, motility and alternative metabolic states are not meaningfully exercised. |

The saved settings have only two textual differences: OpenMP threads are 2 versus 112, and `domain/use_2D` is empty versus `false`. The local XML reader uses `as_bool()` with a false default, and both saved meshes have four z layers, so the empty element does not explain this as a 2D-versus-3D comparison. Both runs use seed 50, a 500 × 500 × 8 µm domain, a 250 × 250 × 4 voxel mesh, 2 µm spacing, and 0.001-minute diffusion/mechanics/phenotype steps. The dFBA step is 0.01 minutes. Diffusion coefficients, decay rates, boundary settings and cell parameters match in the saved XML. The current `config/Ecoli_core.xml` files are byte-identical (SHA-256 `5e0e011badcb1135d8a61e73beb2521cf5ac570ea75ced37010c7d8416f4aedf`); this does not prove historical binary or input provenance.

At time zero, both runs have exactly 211 cells. After sorting by cell ID, positions, volumes and all other stored fields match except orientation. All five initial concentration fields and voxel coordinates/volumes match exactly. Different orientations show that a shared seed did not produce an identical complete stochastic state. Later cell IDs should not be treated as matched lineages after division. The script therefore compares population distributions and geometry independently of cell ordering.

The following measurements describe the **coupled substrate fields at 180 minutes**. Concentrations retain the output's numerical units: the snapshot calls them dimensionless, while intracellular transport parameters use mM. No additional unit conversion was imposed.

| Quantity | PhysiCell-X | PhysiCelldFBA | Relative difference |
|---|---:|---:|---:|
| Minimum oxygen | 22.09692 | 22.60090 | −2.230% |
| Minimum glucose | 24.13557 | 24.72299 | −2.376% |
| Mean acetate | 0.055079 | 0.047700 | +15.469% |
| Mean CO₂ | 0.371320 | 0.320573 | +15.830% |
| Lactate, every voxel | 0 | 0 | Exact agreement in this inactive field |

Integrated oxygen depletion is 15.812% greater in PhysiCell-X, and integrated glucose depletion is 15.281% greater. These quantities are the voxel-volume-weighted deficit relative to the common initial field; **they are not cumulative uptake or a mass-balance measurement**, because boundary exchange is not recorded. Product accumulation uses the corresponding increase relative to time zero.

Whole-field relative L2 errors are only 0.230% for oxygen and 0.286% for glucose, because their large background concentrations dominate the denominator. Normalizing instead by the reference's change from the initial field gives 16.720% and 15.611%. Acetate and CO₂ relative L2 errors are 15.870% and 16.567%. Finite-difference gradient relative L2 errors are larger: oxygen 113.6%, glucose 87.1%, acetate 98.2%, CO₂ 111.6%. These are spatial derivative diagnostics of saved fields, not the solver's face fluxes. They establish spatial disagreement but cannot distinguish transport implementation errors from different cell positions, source placement or boundary handling.

The **mechanics disagreement is stronger than the phenotype disagreement**:

| Measurement | Time (min) | PhysiCell-X | PhysiCelldFBA |
|---|---:|---:|---:|
| Mean reported pressure (dimensionless) | 30 | 0.381316 | 0.015039 |
| Maximum pair overlap | 30 | 78.382% | 1.075% |
| Pairs exceeding 20% overlap | 30 | 16 | 0 |
| Mean reported pressure | 180 | 0.309726 | 0.035455 |
| Maximum reported pressure | 180 | 32.8156 | 0.172820 |
| Maximum pair overlap | 180 | 87.819% | 35.162% |
| Pairs exceeding 20% overlap | 180 | 90 | 1 |
| Radius of gyration (µm) | 180 | 9.99797 | 10.59738 |
| Mean nearest-neighbor distance (µm) | 180 | 1.43682 | 1.53219 |

Pair overlap is `max(0, 1 − distance/(radius_i + radius_j))`, calculated from the full MAT cell population using a spatial tree. Twenty percent is a descriptive cutoff for locating large overlaps, not a validated biological acceptance margin. The 30-minute comparison is especially informative: neither run has divided, total cell volume differs by only 0.000122%, but mean pressure differs by 25.36-fold. Matching IDs have a mean position difference of 0.8783 µm and a maximum of 3.3290 µm. All 16 large-overlap pairs straddle x = 0. At 60 minutes, all 12 such pairs also straddle x = 0. At 180 minutes, 77 of 90 do. This localization makes cross-partition neighbor/ghost interactions a priority for investigation; the precise cause is not established by these snapshots.

Independent pressure reconstruction from all geometric overlaps also reveals a discrepancy: at 30 minutes its mean is 1.60972 in PhysiCell-X versus 0.015024 in the reference, whereas the recorded means are 0.381316 and 0.015039. The reconstruction uses the local standard cell-cell pressure formula, `sum(overlap²)/0.027288820670331`. It excludes wall contributions and uses post-update snapshot positions, so exact equality is not expected. Nevertheless, the large discrepancy in PhysiCell-X reinforces the need to inspect interaction coverage. At 180 minutes, eight PhysiCell-X cell centers lie outside the saved z bounds versus one reference cell, another boundary-related diagnostic rather than proof of its cause.

Two output limitations matter here. Saved velocities are zero in all compared snapshots: the position-update code clears velocity after moving cells, so this cannot establish equal dynamics. Also, PhysiCell-X's neighbor graph has only 107 rows for 211 cells at 30 minutes and 887 rows for 1671 cells at 180 minutes. The reference has a row for every cell. The MPI output path calls a graph writer that iterates local `all_cells` without gathering a global graph. Consequently, graph files were not used to infer global contact counts; all overlaps above come from the full MAT positions and radii.

The **recorded growth phenotype remains close**:

| Time (min) | PhysiCell-X cells | Reference cells | Count difference |
|---|---:|---:|---:|
| 0 | 211 | 211 | 0% |
| 30 | 211 | 211 | 0% |
| 60 | 379 | 386 | −1.813% |
| 90 | 460 | 468 | −1.709% |
| 120 | 825 | 827 | −0.242% |
| 150 | 1289 | 1323 | −2.570% |
| 180 | 1671 | 1673 | −0.120% |

At 180 minutes, total cell volume is 3324.7769 versus 3329.8639 µm³ (−0.1528%), and mean growth rate is 0.832233761 versus 0.832232997 h⁻¹ (+0.0000918%). Across nonzero paired times, the relative difference in mean growth rate never exceeds 0.000356%. At 180 minutes, mean recorded glucose, oxygen and CO₂ specific exchange fluxes differ by less than 0.0005%; acetate differs by 0.00753%. Mean net export rates, which also depend on cell size/state, differ by about 0.17–0.18%.

All recorded cells remain alive, have cycle model 5 and phase 14, and are nonmotile. Saved adhesion, repulsion and maximum adhesion-distance parameters agree. No death-response or motility equivalence can be inferred from inactive processes. Oxygen and glucose remain far above the configured transport Km of 0.01 in both runs; this is consistent with the nearly identical growth rate despite differing extracellular fields. A nutrient-rich colony does not test whether the implementations agree during starvation or metabolic transitions. The zero lactate field likewise does not test active lactate transport.

For certification, the next experiments should first isolate mechanics with identical deterministic positions and no growth, then isolate diffusion with identical fixed sources and boundaries, and finally repeat the coupled model across seeds with biological equivalence margins chosen in advance. The existing evidence is sufficient to flag mechanics and substrate-field discrepancies, but not to assign a unique code defect or certify untested regimes. No simulation source, settings or original outputs were modified for this comparison, and no fresh simulations were run.

Reproduce from the project root with:

```bash
python3 tools/compare_biological_outputs.py output ../PhysiCelldFBA/output
```

The script requires NumPy, SciPy and Matplotlib. It verifies finite arrays, unique cell IDs and voxel coordinates, matching schemas/units and matching meshes/physical parameters; aligns voxels by coordinates; compares every stored phenotype component except IDs, positions, orientations and attack targets; and records input hashes in [provenance.json](provenance.json). Distribution distances are descriptive Wasserstein distances, without treating cells or voxels as independent biological replicates.

Validation passed for all seven paired summaries, the known initial 5% packing overlap, and a same-input comparison of the final substrate fields giving zero concentration and gradient errors. The reported zero death counts, constant phase 14 and zero saved velocities were also checked across every paired snapshot.

Detailed data: [cell metrics](cell_metrics.csv), [substrate metrics](substrate_metrics.csv), [phenotype distributions](phenotype_distributions.csv). Figures: [growth and mechanics](comparison.png), [substrate trajectories](substrates.png), [spatial pressure](spatial_pressure.png). Relevant current source: [XML boolean reader](../../modules/PhysiCell_pugixml.cpp), [cell mechanics and velocity reset](../../core/PhysiCell_cell.cpp), [MPI velocity/neighbor calculations](../../core/PhysiCell_standard_models.cpp), [graph output](../../modules/PhysiCell_MultiCellDS.cpp), [dFBA implementation](../../addons/dFBA/src/dfba_intracellular.cpp). Current source inspection supports the interpretations above but is not a reconstruction of the historical executable.
