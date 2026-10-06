Two concrete defects in PhysiCell-X's mechanics path can suppress repulsion between overlapping cells. The first directly explains the concentration of large overlaps at x = 0; the second independently rejects valid interactions and remains relevant after repairing the first. Both mechanisms were reproduced using production C++ functions. This audit assessed the current source; it did not apply fixes to simulation source or rerun the full colony.

**1. MPI neighbor IDs use a truncated mesh dimension — primary defect.**

The colony selects a 5 µm mechanics voxel size in [main.cpp](../../main.cpp#L154). With z bounds −4 to +4 µm, [Cartesian_Mesh::resize](../../BioFVM/BioFVM_mesh.cpp#L1471) allocates `ceil(8 / 5) = 2` z layers. It assigns global voxel IDs using that count. But [create_moore_neighborhood](../../BioFVM/BioFVM_mesh.cpp#L897) recomputes node counts by conversion to `int`, truncating `8 / 5` to 1:

```cpp
int z_nodes = (bounding_box[5]-bounding_box[2])/dz;
int global_y_jump = z_nodes;
int global_x_jump = z_nodes * y_nodes;
```

For this domain, `global_y_jump` becomes 1 instead of 2 and `global_x_jump` becomes 100 instead of 200. Consequently the left/right Moore lists contain incorrect global voxel IDs.

For example, on rank 0 the boundary voxel centered at (−2.5, 2.5, −1.5) has ID 9900. Its face neighbor on rank 1 is ID 10100, centered at (+2.5, 2.5, −1.5). The code instead requests ID 10000, which is centered at (+2.5, −247.5, −1.5). Some requested IDs are absent from the received face entirely; others identify the wrong voxel.

The failure continues through the stack:

1. [ghost_exchange / serialize_ghost_voxel](../../core/PhysiCell_cell_container_mechanics.cpp#L110) transmits the allocated voxels' actual IDs and centers.
2. [standard_update_cell_velocity](../../core/PhysiCell_standard_models.cpp#L755) searches those received voxels using the incorrect Moore-list IDs.
3. A failed lookup is skipped. A wrongly located voxel can also be rejected by the spatial filter.
4. `Cell::add_potentials` never receives the actual contacting cell across the boundary, so its repulsion and pressure contribution are omitted.

The reference PhysiCelldFBA uses the same 5 µm mechanics voxel size, but its serial local-neighborhood path does not construct these MPI global-neighbor lists. This explains why identical mechanical coefficients do not produce identical packing. The observed boundary localization is consistent with this defect: at both 30 and 60 minutes every pair exceeding 20% overlap straddles x = 0.

The mesh probe instantiated the two x partitions and called the actual production mesh constructor. It checked every geometrically expected adjacent voxel across the partition, in both directions:

| Test | Missing expected directed voxel links |
|---|---:|
| Original code, actual 8 µm slab | 2384 / 2384 |
| Original code, divisible 10 µm slab control | 0 / 2384 |
| Temporary rounding correction, original 8 µm slab | 0 / 2384 |

The temporary correction only made the three integer node-count calculations use the same `ceil(1e-16 + extent/spacing)` expression as `resize()`. It preserved the original 8 µm domain. This establishes the indexing cause without a change to the biological geometry. A production repair should use a single consistent source of global mesh dimensions for allocation, global IDs and neighbor strides, then check neighbor-list validity and uniqueness.

**2. The interaction-radius cache is erased and not rebuilt for the dFBA model.**

In [Cell_Container::update_all_cells](../../core/PhysiCell_cell_container.cpp#L481), PhysiCell-X sets every `max_cell_interactive_distance_in_voxel` entry to zero during the phenotype update. This happens after the intracellular update. The dFBA growth update calls [set_total_volume](../../addons/dFBA/src/dfba_intracellular.cpp#L621), which can populate the cache, but the subsequent reset erases those values.

The code expects the bundled phenotype/volume update to restore the cache. However, this colony explicitly sets both `update_phenotype` and `volume_update_function` to null in [custom.cpp](../../custom_modules/custom.cpp#L91). In [advance_bundled_phenotype_functions](../../core/PhysiCell_cell.cpp#L451), the relevant `set_total_volume` calls depend on those callbacks. Geometry update alone does not rebuild the voxel maximum. Thus ordinary non-dividing cells leave the cache empty before mechanics. Division or other explicit volume/position operations can populate some entries, making the failure state-dependent.

The corresponding reset is deliberately commented out in [PhysiCelldFBA's update_all_cells](../../../PhysiCelldFBA/core/PhysiCell_cell_container.cpp#L162).

Both local and MPI [is_neighbor_voxel](../../core/PhysiCell_cell.cpp#L2086) filters use the querying cell's interaction radius plus the other voxel's cached maximum. With a zero cache, the effective search reach is too short. Ghost serialization also transmits that zero value, so correcting MPI IDs alone does not repair this second problem.

A direct call to the freshly compiled production MPI filter reproduced the false rejection with two equal-radius cells:

| Quantity | Value |
|---|---:|
| Radius | 0.777748 µm |
| Cell centers | (−1, 0, 0) and (+0.4, 0, 0) µm |
| Physical overlap | 9.997% |
| Other-voxel cache = 0 | Voxel rejected |
| Other-voxel cache = 0.777748 | Voxel accepted |

The robust repair is to rebuild each voxel's maximum interaction reach from its current `agent_grid` cells immediately before halo exchange and mechanics, after phenotype/division updates. A parallel loop over voxels can give each maximum one writer and avoid a race between cells sharing a voxel. The reach must cover both repulsion and adhesion. It should not depend on whether a biological model happens to install a volume callback. Removing the reset would restore the reference's approach, but an explicit rebuild makes the invariant clearer and avoids stale maxima.

The existing mechanics validation setup is poorly positioned to reveal either defect: [validation.xml](../../tests/mpi/mechanics/config/validation.xml) has extents divisible by its 20 µm mechanics spacing, and [its cell setup](../../tests/mpi/mechanics/custom_modules/custom.cpp#L129) installs `standard_volume_update_function`, which refills the radius cache. Regression coverage should include a nondivisible thin slab and null phenotype/volume callbacks, with overlapping pairs both across the MPI split and across local voxel faces.

The local and ghost repulsion formulas both use `(1 − distance/(r_i+r_j))²` times the geometric mean of the two repulsion coefficients. The defects above prevent that calculation from being reached for real contacts. The incomplete neighbor-graph text output identified earlier is an additional diagnostic-output issue; it does not explain the independently measured geometric overlaps.

Validation is limited to mesh construction and the production geometric filter. The mesh probe supplies rank/topology values without MPI communication, because the constructor performs no communication. The filter probe recompiles `core/PhysiCell_cell.cpp` and links the existing project objects for its supporting definitions. A full simulation after both fixes is still needed to quantify the improvement in colony overlap, pressure and substrate fields; no claim is made that these are the only defects in the complete simulator.

Reproduce from the project root:

```bash
python3 analysis/biological_equivalence/run_code_probes.py
```

The runner compiles into a temporary directory, tests the original code and a temporary corrected mesh source, and leaves production code unchanged. [Probe results](code_probe_results.txt) and [commands/source hashes](code_probe_provenance.json) document the execution. Probe sources: [mesh-neighbor check](probe_mesh_neighbors.cpp), [interaction-radius filter check](probe_mechanics_filter.cpp).
