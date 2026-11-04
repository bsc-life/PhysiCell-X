# PhysiCell-X intracellular dFBA sample validation

Validation date: 2026-11-04

## Scope and provenance

- PhysiCell-X: `8b454333b2f3ad8ca0785213f9549355045bb497` (`dev_biofvmb`), plus the fixes in this working tree.
- PhysiCelldFBA reference: `7f2b1d0a52782918ea4319c2ef94563bf719f253` (`development`).
- Compiler: Open MPI `mpic++`, GCC 11.4.0.
- Comparison procedure: `tools/agent/biological_equivalence.md` and `tools/compare_simulations.py`.
- All comparisons used one MPI rank, one OpenMP thread, XML seed 50, identical paired XML inputs, identical SBML model hashes, and labelled XML/MAT output rather than file-index assumptions.

Project mapping:

| PhysiCell-X sample | PhysiCelldFBA reference |
|---|---|
| `bacterial_colony` | `bacterial_colony` |
| `cancer_metabolism` | `human_template` (renamed reference sample) |
| `dfba_unit_test` | `dfba_unit_test` |

## Repairs required before validation

- Added a root selector for `dfba_unit_test` and made it idempotent.
- Corrected the cancer and unit-test Makefiles so plain `make` builds the executable, resolves the downloaded Coin-OR archive names, and installs/checks both dFBA dependencies.
- Added the missing cancer SBML model and cell CSV, and the missing unit-test config/model/script tree.
- Synchronized sample XML parameters and observable reaction IDs with the reference projects.
- Made custom-data writes safe and exported realised exchange fluxes, applied lower bounds, and growth rate under declared XML labels.
- Made the libSBML setup tolerate the known shared-filesystem dangling-symlink metadata failure while still requiring the header and static archive.

## Build results

Each sample was selected at the repository root, cleaned, and built with `make -j4` in an isolated copy.

| PhysiCell-X target | Result | Wall time |
|---|---:|---:|
| `bacterial-colony` | PASS | 29.82 s |
| `cancer-metabolism-sample` | PASS | 17.73 s |
| `dfba_unit_test` | PASS | 45.71 s |

## Paired simulations

Every executable was protected by `timeout 180s`. The requested configured maxima were not exceeded; PhysiCell also emits its normal final snapshot one diffusion step after the configured maximum.

| Pair | Configured simulated time | Last matched snapshot | PhysiCell-X wall time | Reference wall time | Matched snapshots |
|---|---:|---:|---:|---:|---:|
| bacterial colony | 0.005 min | 0.006 min | 1.46 s | 1.05 s | 7 |
| cancer/human | 6.0 min | 6.01 min | 7.91 s | 3.16 s | 4 |
| dFBA unit test | 0.01 min | 0.011 min | 0.85 s | 0.16 s | 4 |

The complete initial labelled cell state and all initial substrate fields were bitwise equal for every pair.

## Biological comparison

Acceptance for this short deterministic smoke test was defined on biological observables: identical population fate/counts, less than 1% relative difference in growth and non-negligible substrate-field L2 norms, and absolute-error interpretation for fields whose reference norm is close to zero. Exact floating-point identity was not expected because PhysiCell-X uses its distributed BioFVM/cell-update path while the reference uses the serial path.

| Final metric | bacterial colony | cancer/human | dFBA unit test |
|---|---:|---:|---:|
| cells (X / reference) | 211 / 211 | 1 / 1 | 1 / 1 |
| live / dead (both) | 211 / 0 | 1 / 0 | 1 / 0 |
| mean growth difference | 0.000636% | 9.10e-11% | 0.000112% |
| largest non-negligible field relative L2 | 0.00250% (glucose) | 0.483% (oxygen) | 0.000110% (glucose) |
| mechanics/position outcome | equivalent; pressure mean differs 3.99e-7% | identical | identical |
| total-volume difference | 3.03e-8% | 0% | 1.23e-8% |

Near-zero product fields have misleading percentage errors but small absolute errors:

- Bacterial CO2: 10.77% relative L2, maximum absolute error 0.00291; acetate: 12.61%, maximum absolute error 0.000466.
- Unit-test CO2: 3.37% relative L2, maximum absolute error 3.78e-5.
- Human/cancer CO2 remained exactly zero. Lactate relative L2 was 0.472%.

The final mean realised fluxes also agree closely except where a low local substrate amplifies a relative percentage. The largest such case is human/cancer oxygen exchange: `-0.277814` versus `-0.258646` mmol/gDW/h (7.41%, absolute difference 0.0192), while growth is equal to displayed precision and population/volume outcomes are identical.

## Earliest divergence and interpretation

| Pair | Last exact matched state | First numerical divergence |
|---|---:|---:|
| bacterial colony | 0.001 min | fields at 0.002 min; cell flux/volume/mechanics labels at 0.003 min |
| cancer/human | 0 min | fields and flux labels at 3 min |
| dFBA unit test | 0 min | fields and flux/volume labels at 0.005 min |

The divergence pattern starts in diffusion/export coupling and remains small in absolute terms. This is consistent with the different distributed and serial execution/update paths, not an input, seed, model, mesh, field-label, or initial-cell mismatch. No divergence in fate counts was observed.

## Conclusion and limits

All three repaired PhysiCell-X intracellular dFBA samples pass clean-build and short, one-rank biological-equivalence smoke tests against their PhysiCelldFBA references.

This is one deterministic short run per sample. It does not establish statistical equivalence, long-horizon equivalence, or multi-rank equivalence; those require additional replicates and MPI decomposition tests.
