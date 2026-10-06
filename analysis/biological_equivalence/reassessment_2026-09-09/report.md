**The latest changes substantially repaired the original MPI-boundary mechanics failure, but they did not resolve overall biological equivalence.** Early mechanics now agrees almost exactly with the reference. Growth remains close. Substrate depletion and product accumulation still differ systematically and are slightly worse at the original 180-minute endpoint. Later mechanics also retains differences in overlap and pressure distributions.

Follow-up: [isolated substrate-cause investigation](transport_diagnosis/report.md) finds identical diffusion/fixed-source replay results, matching boundary evidence, and shared defects in concurrent source application and cached source timesteps. It also explains why the concentration curves look closer than the depletion percentages suggest.

This reassessment was conducted on 2026-09-09 using the newer `PhysiCell-X/output` and the same `../PhysiCelldFBA/output` reference. It reruns the saved-output analysis; it does not launch another colony simulation. The [previous report](../report.md), [code audit](../code_audit.md), and original CSVs remain unchanged.

| Component | Assessment |
|---|---|
| Original cross-partition mechanics defect | Strong evidence of repair: production mesh probe passes and the pre-division packing failure disappears. |
| Mechanics over the complete trajectory | Much improved, but pressure and overlap distributions still differ after division. |
| Growth and recorded metabolism | Close agreement through 330 minutes in this nutrient-rich regime. |
| Coupled substrate fields | Mismatch persists; not explained by the saved colony geometry alone. |
| Overall biological equivalence | Not established by these single-run comparisons. |

There are now **12 paired full snapshots, 0–330 minutes at 30-minute intervals**. Before/after claims use the original shared 0–180-minute interval; there are no old PhysiCell-X results at 210–330 minutes in the previous analysis. The newer run log, `test-45570163.log`, ends with cancellation due to the job time limit, and neither comparison establishes the configured 600-minute outcome. SVGs beyond 330 minutes are not substitutes for full cell and concentration snapshots. Reference `final.*` files remain excluded.

Input checks confirm that all 44 reference files covered by the old provenance manifest are byte-identical to the previous analysis. Of the corresponding 44 PhysiCell-X files, 38 have changed. The earlier manifest's absolute checkout paths were relocated by project and filename for this check. Saved settings still differ only in OpenMP threads (2 versus 112) and the empty versus `false` `use_2D` element; both outputs have the same three-dimensional mesh. Current metabolic model XML files are byte-identical. At time zero, both runs have 211 cells with identical positions, volumes, and other stored fields except orientation, and exactly identical initial concentration fields. A shared seed therefore still does not establish identical complete stochastic states.

The current source contains both changes suggested by the earlier audit. [Mesh neighbor indexing](../../../BioFVM/BioFVM_mesh.cpp) now uses `ceil(1e-16 + extent/spacing)`, consistent with mesh allocation. The destructive radius-cache reset is commented out in the MPI overload of [update_all_cells](../../../core/PhysiCell_cell_container.cpp). This follows the reference's approach of retaining the cache; it is not an explicit rebuild before each mechanics step. The serial overload still contains its reset, so this finding should not be generalized to every execution path.

The existing production mesh probe was freshly compiled and run against the current source. It found **0 missing links out of 2,384 expected directed cross-partition links** in the actual 8 µm slab, versus **2,384/2,384 missing previously**. The divisible 10 µm control also passed. See [probe commands, results and source hashes](mesh_probe.json). The old probe runner expects the defective code and applies a temporary repair, so it was not reused unchanged. This probe verifies mesh construction, not a complete MPI communication exchange or a runtime cache invariant. Object and executable timestamps precede the newer output and are consistent with a rebuild; current binary/source hashes are recorded, but timestamps cannot prove historical build provenance or separate the effects of the two fixes.

The strongest mechanical evidence comes from **30 minutes, before either population has divided**:

| Measurement | Before | After | Reference |
|---|---:|---:|---:|
| Mean reported pressure | 0.381316 | 0.015039071 | 0.015039011 |
| Maximum fractional overlap | 78.382% | 1.075264% | 1.075284% |
| Pairs exceeding 20% overlap | 16 | 0 | 0 |
| Mean ID-matched position difference from reference (µm) | 0.8783 | 0.000003146 | — |

Mean pressure now differs by only **0.000399%**. Mean reconstructed geometric pressure is 0.015023642 versus 0.015023582. The maximum ID-matched position difference is 0.000004569 µm. These measurements strongly support repair of the original missing-contact failure rather than merely a change in population growth.

At the original **180-minute endpoint**, the improvement remains large:

| Measurement | Before | After | Reference |
|---|---:|---:|---:|
| Mean reported pressure | 0.309726 | 0.0355464 | 0.0354551 |
| Maximum reported pressure | 32.8156 | 0.152588 | 0.172820 |
| Maximum fractional overlap | 87.819% | 34.703% | 35.162% |
| Pairs exceeding 20% overlap | 90 | 3 | 1 |
| Such pairs crossing x = 0 | 77 | 1 | 0 |
| Radius of gyration (µm) | 9.99797 | 10.60506 | 10.59738 |

Large-overlap pairs fall by **96.7%**, and mean pressure is within **0.257%** of the reference. However, this endpoint does not establish agreement throughout the trajectory. At 60 and 90 minutes, current mean pressure is respectively 66.5% and 55.3% below the reference, although absolute pressures are small (0.000890 versus 0.002652 and 0.002541 versus 0.005680). Division timing and different stochastic trajectories can influence these sparse snapshots; their contribution has not been isolated.

At **330 minutes**, mean pressure is 0.328452 versus 0.342040 (−3.97%), radius of gyration differs by −0.269%, and mean nearest-neighbor distance by about −0.422%. Nevertheless, there are **243 versus 136 pairs exceeding 20% overlap**, involving 416 versus 216 cells. Only 3 of the 243 PhysiCell-X pairs cross x = 0, so the old midplane concentration of large overlaps is absent. Maximum overlap is actually higher in the reference (91.2% versus 73.2%), illustrating why a single maximum is insufficient to compare packing.

Geometric-pressure means remain 23.9% above the reference at 180 minutes and 16.4% above it at 330 minutes. This reconstruction uses all saved pair overlaps and the standard `sum(overlap²)/0.027288820670331` formula. It excludes wall forces and uses post-update positions, whereas recorded pressure comes from the mechanics update; differences cannot alone identify a remaining force defect. At 330 minutes, 195 versus 212 cell centers lie outside the saved z bounds, so boundary excursions occur in both implementations.

**Substrate-field agreement has not been repaired.** On the directly comparable 180-minute snapshot:

| Voxel-volume-weighted change relative to the common initial field | Previous excess versus reference | Current excess versus reference |
|---|---:|---:|
| Oxygen depletion | +15.81% | +18.50% |
| Glucose depletion | +15.28% | +18.07% |
| Acetate accumulation | +15.47% | +18.22% |
| CO₂ accumulation | +15.83% | +18.52% |

These differences worsen by approximately 2.7–2.8 percentage points. They measure the saved field's deficit or increase, **not cumulative cellular uptake or a mass balance**, because boundary exchange is not available. Concentrations retain their saved numerical units, without imposing a unit conversion.

At 180 minutes, field L2 errors normalized by the reference change from the initial field are 19.10% for oxygen, 18.28% for glucose, 18.48% for acetate and 19.00% for CO₂. Finite-difference gradient errors remain 83.7–109.4%. These spatial diagnostics mix transport, source placement and the differing cell trajectories; they do not isolate the diffusion operator.

Crucially, the mismatch already exists at **30 minutes**: integrated oxygen/glucose depletion is 23.02%/24.13% greater, and acetate/CO₂ accumulation is 23.45%/22.64% greater. At that time positions agree to a few millionths of a micrometre, total cell volume differs by only 0.000118%, and mean stored net export rates for active substrates differ in magnitude by at most 0.00384%. Thus the saved geometry and recorded exchange-rate differences are insufficient to explain the field mismatch. This points toward the extracellular source/transport coupling, boundary handling or update scheduling as the next investigation target; it does not identify a unique code defect or prove identical histories between snapshots.

At 330 minutes, integrated field differences decrease to 3.01–3.08%, with perturbation-normalized L2 errors of 3.07–3.36%. Gradient errors remain 36.2–53.7%. This later reduction cannot be credited as a before/after improvement because the original analysis has no matching old 330-minute output. Lactate remains exactly zero in both projects, and no negative concentrations were found.

**Growth remains close across all 12 paired times.** Population-count differences range from −2.137% to +1.799%; total-volume differences stay within 0.188%. At 180 minutes there are 1,670 versus 1,673 cells, with total volume −0.0297% relative to the reference. At 330 minutes there are 13,215 versus 13,227 cells (−0.0907%), and total volume differs by −0.1655%. The maximum relative difference in mean growth rate across the paired times is 0.000431%. At 330 minutes, mean recorded glucose, oxygen and CO₂ specific exchange fluxes differ by less than 0.0015%; acetate differs by 0.0246%.

All cells remain alive and nonmotile, with cycle model 5 and phase 14. Oxygen and glucose remain well above the configured Km of 0.01, including at 330 minutes. Agreement in this nutrient-rich case does not exercise starvation, death, motility or active lactate transport. Saved velocities are zero after position updates and cannot establish equal dynamics. MPI neighbor graph output remains incomplete (6,639 rows for 13,215 cells at 330 minutes); contact measurements here use the full MAT cell population instead.

The next discriminating experiment is a **fixed-cell, no-growth source/transport comparison**, starting from identical fields and positions. Compare serial/reference, one MPI rank and two MPI ranks, recording per-substrate field changes immediately before and after source application and diffusion, plus total owned-cell source contributions and boundary exchange. This would distinguish duplicate/missing source application, scaling or scheduling differences from diffusion and boundary effects. A repeated-seed coupled comparison with predetermined biological acceptance margins is still required for an equivalence claim. The present evidence supports a successful repair of the dominant mechanics defect, with extracellular coupling and later packing still requiring validation.

Reproduce from the project root:

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/compare_biological_outputs.py output ../PhysiCelldFBA/output --out analysis/biological_equivalence/reassessment_2026-09-09
python3 analysis/biological_equivalence/reassessment_2026-09-09/supplement.py
```

The comparison checks finite arrays, unique IDs and coordinates, matching schemas/units/meshes and snapshot physical parameters, and records input hashes. Supplementary validation checks unchanged reference data/metrics, inactive phenotypes, nonnegative fields, matching current model files, and zero concentration/gradient errors for a same-input final-field comparison. These are descriptive diagnostics from one run per implementation, without treating cells or voxels as independent biological replicates; 20% overlap is a diagnostic cutoff, not a validated biological acceptance margin. Simulation source, settings and original outputs were not modified.

Artifacts: [before/after figure](before_after.png), [before/after measurements](before_after.csv), [full trajectory figure](comparison.png), [substrate trajectories](substrates.png), [full-colony pressure at 330 minutes](spatial_pressure_full.png), [cell metrics](cell_metrics.csv), [substrate metrics](substrate_metrics.csv), [phenotype distributions](phenotype_distributions.csv), [input provenance](provenance.json), and [validation](validation.json). Use `spatial_pressure_full.png` for the later colony; the original comparison script's `spatial_pressure.png` retains its old ±22 µm crop and 0.2 pressure saturation.
