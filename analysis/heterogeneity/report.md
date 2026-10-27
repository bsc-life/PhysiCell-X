# Simulation comparison

Compared 18 matched times, 0–1000.01 min.

| Final metric | Simulation | Reference | Δ | Δ (%) |
|---|---:|---:|---:|---:|
| apoptotic_count | 499 | 517 | -18 | -3.482 |
| arrested_count | 1 | 0 | +1 | +nan |
| centroid_x | -0.273311 | -0.360702 | +0.087391 | +24.23 |
| centroid_y | -0.185461 | -0.73535 | +0.549889 | +74.78 |
| centroid_z | -0.37482 | -0.294331 | -0.0804888 | -27.35 |
| dead_count | 499 | 517 | -18 | -3.482 |
| dead_volume | 77104.6 | 83934 | -6829.44 | -8.137 |
| growing_count | 20203 | 20102 | +101 | +0.5024 |
| growth_mean | 0.000140451 | 0.000139583 | +8.68311e-07 | +0.6221 |
| live_count | 20204 | 20102 | +102 | +0.5074 |
| n_cells | 20703 | 20619 | +84 | +0.4074 |
| nearest_neighbor_mean | 14.8176 | 14.8494 | -0.031806 | -0.2142 |
| nearest_neighbor_overlap_mean | 1.09891 | 1.08082 | +0.0180815 | +1.673 |
| necrotic_count | 0 | 0 | +0 | +nan |
| oncoprotein_mean | 1.00947 | 1.01014 | -0.000673962 | -0.06672 |
| oncoprotein_p95 | 1.4135 | 1.41708 | -0.00357424 | -0.2522 |
| pressure_mean | 0.944006 | 0.924024 | +0.0199824 | +2.163 |
| pressure_p95 | 1.69496 | 1.66031 | +0.0346506 | +2.087 |
| radial_max | 267.165 | 267.799 | -0.634289 | -0.2369 |
| radial_p50 | 199.674 | 199.458 | +0.215844 | +0.1082 |
| radial_p95 | 246.164 | 245.906 | +0.257534 | +0.1047 |
| spread_x | 112.633 | 112.392 | +0.241051 | +0.2145 |
| spread_y | 112.341 | 112.44 | -0.0996459 | -0.08862 |
| spread_z | 112.439 | 112.359 | +0.0795344 | +0.07079 |
| total_volume | 4.85683e+07 | 4.83864e+07 | +181914 | +0.376 |
| total_volume_mean | 2345.95 | 2346.69 | -0.734587 | -0.0313 |
| total_volume_p95 | 2494 | 2494 | +0 | +0 |
| vertical_p95 | 203.506 | 202.912 | +0.594408 | +0.2929 |

## Interpretation limits

- simulation: no archived PhysiCell_settings.xml
- Recorded seeds differ: 798451835 versus 0
- One run per simulator: descriptive agreement only, not statistical equivalence. Source-code causes are not inferred by this script.

All population metrics use the full MAT population. Cell maps are thin slices. The dotted x line is a user-selected diagnostic location, not an automatically detected MPI boundary.

Diffusion summaries and physical cell-position exposure can conceal incorrect runtime cell ownership; use the spatial difference maps and x profiles. Actual fluxes, cumulative division/death events, solver state and code provenance cannot be recovered from these snapshots alone.
