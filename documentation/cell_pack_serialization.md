# `Cell::pack()` Serialization Map

This document follows the current `PhysiCell::Cell::pack()` order in `core/PhysiCell_cell_mpi.cpp` and classifies each packed field by how it can be reconstructed on the receiving MPI rank.

## Legend

- `Definition`: copied from `Cell_Definition` when a cell is created. The receiver can rebuild it from the local `Cell_Definition` if both ranks share the same definitions.
- `State`: true per-cell mutable state. This must be serialized to preserve the migrated cell.
- `Mixed`: starts from the definition, but can become cell-specific during the run. A future state-only MPI format would need to split this object.
- `Alias/Cache`: runtime cache or pointer-backed duplicate of data already stored elsewhere. It can usually be rebound or recomputed locally.

## Where definition data enters a `Cell`

- `Cell::Cell()` and `Cell::Cell(int)` start from `cell_defaults`.
- `create_cell(Cell_Definition&)` copies `type`, `type_name`, `custom_data`, `parameters`, `functions`, `phenotype`, and `is_movable` from the chosen `Cell_Definition`.
- `Cell::convert_to_cell_definition()` shows the intended conserved state during a type change:
  - preserved: `phenotype.volume`, `phenotype.geometry`, `phenotype.molecular.internalized_total_substrates`, and custom variables marked `conserved_quantity`
  - replaced from the new definition: most of `custom_data`, `parameters`, `functions`, and `phenotype`

That means a receiver can rebuild a large part of a migrated cell from the local `Cell_Definition`, then apply only the true per-cell state.

## Ordered Map Of `Cell::pack()`

### 1. `ID`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `ID` | `State` | Not from `Cell_Definition` | Must serialize. This is the cell identity. |

### 2. `position`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `position` | `State` | Not from `Cell_Definition` | Must serialize. Spatial location is rank-local runtime state. |

### 3. `type_name`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `type_name` | `Definition` | `Cell_Definition.name` | Redundant if `type` is already serialized and the receiver has the same cell-definition table. |

### 4. `custom_data`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `custom_data.name_to_index_map.size()` | `Definition` | From `Cell_Definition.custom_data` schema | Does not need per-cell serialization if schemas are identical on all ranks. |
| `custom_data.name_to_index_map[key]` | `Definition` | From `Cell_Definition.custom_data` schema | Schema only. |
| `custom_data.variables.size()` | `Definition` | From definition schema | Schema only. |
| `custom_data.variables[i].name` | `Definition` | From definition schema | Schema only. |
| `custom_data.variables[i].value` | `Mixed` | Default comes from definition | Serialize if this custom scalar can differ per cell. |
| `custom_data.variables[i].units` | `Definition` | From definition schema | Schema only. |
| `custom_data.variables[i].conserved_quantity` | `Definition` | From definition schema | Metadata only. |
| `custom_data.vector_variables.size()` | `Definition` | From definition schema | Schema only. |
| `custom_data.vector_variables[i].name` | `Definition` | From definition schema | Schema only. |
| `custom_data.vector_variables[i].value` | `Mixed` | Default comes from definition | Serialize if this custom vector can differ per cell. |
| `custom_data.vector_variables[i].units` | `Definition` | From definition schema | Schema only. |
| `custom_data.vector_variables[i].conserved_quantity` | `Definition` | From definition schema | Metadata only. |

### 5. `parameters`

Overall mode: `Definition`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `parameters.o2_hypoxic_threshold` | `Definition` | From `Cell_Definition.parameters` | Can be rebuilt from local definition unless user code mutates it per cell. |
| `parameters.o2_hypoxic_response` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.o2_hypoxic_saturation` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.o2_proliferation_saturation` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.o2_proliferation_threshold` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.o2_reference` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.o2_necrosis_threshold` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.o2_necrosis_max` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.max_necrosis_rate` | `Definition` | From `Cell_Definition.parameters` | Same as above. |
| `parameters.necrosis_type` | `Definition` | From `Cell_Definition.parameters` | Same as above. |

Not packed here: `parameters.pReference_live_phenotype` is a pointer and is intentionally not serialized.

### 6. `functions`

Overall mode: `Definition`

Current code serializes only `functions.cycle_model`. All other function pointers are not serialized and must be supplied by the local definition.

#### 6.1 `functions.cycle_model`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `functions.cycle_model.inverse_index_maps.size()` | `Definition` | Rebuild from cycle topology in local definition | Topology metadata. |
| `functions.cycle_model.inverse_index_maps[i][end_phase]` | `Definition` | Rebuild from cycle topology in local definition | Topology metadata. |
| `functions.cycle_model.name` | `Definition` | From local definition | Redundant if definitions are synchronized. |
| `functions.cycle_model.code` | `Definition` | From local definition | Redundant if definitions are synchronized. |
| `functions.cycle_model.phases.size()` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phases[i].index` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phases[i].code` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phases[i].name` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phases[i].division_at_phase_exit` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phases[i].removal_at_phase_exit` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phase_links.size()` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phase_links[i].size()` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phase_links[i][j].start_phase_index` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phase_links[i][j].end_phase_index` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.phase_links[i][j].fixed_duration` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.default_phase_index` | `Definition` | From local definition | Topology metadata. |
| `functions.cycle_model.data.inverse_index_maps` | `Definition` | Rebuild from local cycle model | Default model metadata, not active cell cycle state. |
| `functions.cycle_model.data.time_units` | `Definition` | From local definition | Default model metadata. |
| `functions.cycle_model.data.transition_rates` | `Definition` | From local definition | Default model rates. Serialize only if these are intentionally changed per cell. |
| `functions.cycle_model.data.current_phase_index` | `Definition` | From local definition | This is model-default data here, not the active per-cell cycle phase. |
| `functions.cycle_model.data.elapsed_time_in_phase` | `Definition` | From local definition | This is model-default data here, not the active per-cell cycle phase. |

Not packed inside `functions`:

- `instantiate_cell`
- `volume_update_function`
- `update_migration_bias`
- `custom_cell_rule`
- `update_phenotype`
- `update_phenotype_parallel`
- `pre_update_intracellular`
- `post_update_intracellular`
- `update_velocity`
- `update_velocity_parallel`
- `add_cell_basement_membrane_interactions`
- `calculate_distance_to_membrane`
- `set_orientation`
- `contact_function`
- `cell_division_function`
- `plot_agent_SVG`
- `plot_agent_legend`
- `Phase::entry_function`
- `Phase_Link::arrest_function`
- `Phase_Link::exit_function`

### 7. `state`

Overall mode: `State`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `state.orientation` | `State` | Constructor gives a default or random orientation | Serialize to preserve the exact migrated cell state. |
| `state.simple_pressure` | `State` | Not from definition | Serialize. |
| `state.number_of_nuclei` | `State` | Constructor default is `1` | Serialize. |
| `state.total_attack_time` | `State` | Constructor default is `0` | Serialize. |
| `state.contact_with_basement_membrane` | `State` | Constructor default is `false` | Serialize. |

Not packed in `Cell_State`:

- `attached_cells`
- `spring_attachments`
- `neighbors`
- `mechanics_neighbor_candidates`
- `mechanics_neighbor_interactions`

These are local graph or neighbor-cache data and should be rebuilt after migration.

### 8. `phenotype`

Overall mode: `Mixed`

#### 8.1 Top-level phenotype flags

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.flagged_for_division` | `State` | Not from definition | Serialize. |
| `phenotype.flagged_for_removal` | `State` | Not from definition | Serialize. |

#### 8.2 `phenotype.cycle`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.cycle.data.inverse_index_maps` | `Definition` | Rebuild from local cycle model | Topology copy. |
| `phenotype.cycle.data.time_units` | `Definition` | From local cycle model | Definition metadata. |
| `phenotype.cycle.data.transition_rates` | `Mixed` | Defaults from local cycle model | Serialize if cycle rates can vary per cell. |
| `phenotype.cycle.data.current_phase_index` | `State` | Not from definition | Must serialize. This is active cycle state. |
| `phenotype.cycle.data.elapsed_time_in_phase` | `State` | Not from definition | Must serialize. |
| `phenotype.cycle.asymmetric_division.asymmetric_division_probabilities` | `Mixed` | Defaults from local definition | Serialize if asymmetric-division probabilities can vary per cell. |

Not packed in `phenotype.cycle`:

- `phenotype.cycle.pCycle_Model`

After `Cell::unpack()`, the current code rebinds this pointer with `find_cell_definition(this->type)` and `phenotype.cycle.sync_to_cycle_model(...)`.

#### 8.3 `phenotype.death`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.death.rates` | `Mixed` | Defaults from local definition | Serialize if death rates can vary per cell. |
| `phenotype.death.parameters[i].time_units` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.parameters[i].unlysed_fluid_change_rate` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.parameters[i].lysed_fluid_change_rate` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.parameters[i].cytoplasmic_biomass_change_rate` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.parameters[i].nuclear_biomass_change_rate` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.parameters[i].calcification_rate` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.parameters[i].relative_rupture_volume` | `Definition` | From local definition | Death-model parameter. |
| `phenotype.death.dead` | `State` | Not from definition | Must serialize. |
| `phenotype.death.current_death_model_index` | `State` | Not from definition | Must serialize. |

Not packed in `phenotype.death`:

- `phenotype.death.models`

These are pointers to local death models and must be rebound from the local definition.

#### 8.4 `phenotype.volume`

Overall mode: `Mixed`

The header already splits this class into state variables and user-set parameters.

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.volume.total` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.solid` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.fluid` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.fluid_fraction` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.nuclear` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.nuclear_fluid` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.nuclear_solid` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.cytoplasmic` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.cytoplasmic_fluid` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.cytoplasmic_solid` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.calcified_fraction` | `State` | Definition gives only the initial value | Must serialize. |
| `phenotype.volume.cytoplasmic_to_nuclear_ratio` | `State` | Can be recomputed, but current code preserves it explicitly | Serialize for exact continuity. |
| `phenotype.volume.rupture_volume` | `State` | Can be recomputed, but current code preserves it explicitly | Serialize for exact continuity. |
| `phenotype.volume.cytoplasmic_biomass_change_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.nuclear_biomass_change_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.fluid_change_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.calcification_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.target_solid_cytoplasmic` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.target_solid_nuclear` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.target_fluid_fraction` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.target_cytoplasmic_to_nuclear_ratio` | `Definition` | From local definition | Parameter. |
| `phenotype.volume.relative_rupture_volume` | `Definition` | From local definition | Parameter. |

#### 8.5 `phenotype.geometry`

Overall mode: `State`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.geometry.radius` | `State` | Can be recomputed from `volume.total` | Serialize for exact continuity, or recompute after unpack. |
| `phenotype.geometry.nuclear_radius` | `State` | Can be recomputed from `volume.nuclear` | Serialize for exact continuity, or recompute after unpack. |
| `phenotype.geometry.surface_area` | `State` | Can be recomputed from volume and radius | Serialize for exact continuity, or recompute after unpack. |
| `phenotype.geometry.polarity` | `State` | Not from definition in general | Serialize. |

#### 8.6 `phenotype.mechanics`

Overall mode: `Definition`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.mechanics.cell_cell_adhesion_strength` | `Definition` | From local definition | Serialize only if mechanics parameters can vary per cell. |
| `phenotype.mechanics.cell_BM_adhesion_strength` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.cell_cell_repulsion_strength` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.cell_BM_repulsion_strength` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.relative_maximum_adhesion_distance` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.attachment_rate` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.detachment_rate` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.relative_maximum_attachment_distance` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.relative_detachment_distance` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.attachment_elastic_constant` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.maximum_attachment_rate` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.maximum_number_of_attachments` | `Definition` | From local definition | Same as above. |
| `phenotype.mechanics.cell_adhesion_affinities` | `Definition` | From local definition | Definition-sized affinity table. |

#### 8.7 `phenotype.motility`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.motility.is_motile` | `Definition` | From local definition | Serialize only if motility enablement can vary per cell. |
| `phenotype.motility.persistence_time` | `Definition` | From local definition | Motility parameter. |
| `phenotype.motility.migration_speed` | `Definition` | From local definition | Motility parameter. |
| `phenotype.motility.migration_bias_direction` | `State` | Not reliably reconstructible from `Cell_Definition` | Serialize if you want the same biased direction immediately after migration. |
| `phenotype.motility.restrict_to_2D` | `Definition` | From local definition | Motility parameter. |
| `phenotype.motility.motility_vector` | `State` | Not from definition | Serialize. |
| `phenotype.motility.chemotaxis_index` | `Definition` | From local definition | Motility parameter. |
| `phenotype.motility.chemotaxis_direction` | `Definition` | From local definition | Motility parameter. |
| `phenotype.motility.chemotactic_sensitivities` | `Definition` | From local definition and local microenvironment sizing | Motility parameter vector. |

Important omission:

- `phenotype.motility.migration_bias` exists in the class, but `Motility::pack()` does not serialize it.

#### 8.8 `phenotype.secretion`

Overall mode: `Definition`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.secretion.secretion_rates` | `Definition` | From local definition and local microenvironment sizing | Serialize only if secretion is changed per cell at runtime. |
| `phenotype.secretion.uptake_rates` | `Definition` | From local definition and local microenvironment sizing | Same as above. |
| `phenotype.secretion.saturation_densities` | `Definition` | From local definition and local microenvironment sizing | Same as above. |
| `phenotype.secretion.net_export_rates` | `Definition` | From local definition and local microenvironment sizing | Same as above. |

Not packed in `phenotype.secretion`:

- `phenotype.secretion.pMicroenvironment`

This must be rebound locally.

#### 8.9 `phenotype.molecular`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.molecular.internalized_total_substrates` | `State` | Not from definition | Must serialize. This is conserved intracellular substrate state. |
| `phenotype.molecular.fraction_released_at_death` | `Definition` | From local definition and local microenvironment sizing | Parameter vector. |
| `phenotype.molecular.fraction_transferred_when_ingested` | `Definition` | From local definition and local microenvironment sizing | Parameter vector. |

Not packed in `phenotype.molecular`:

- `phenotype.molecular.pMicroenvironment`

This must be rebound locally.

#### 8.10 `phenotype.cell_integrity`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.cell_integrity.damage` | `State` | Not from definition | Must serialize. |
| `phenotype.cell_integrity.damage_rate` | `Definition` | From local definition | Serialize only if damage kinetics vary per cell. |
| `phenotype.cell_integrity.damage_repair_rate` | `Definition` | From local definition | Serialize only if damage kinetics vary per cell. |

#### 8.11 `phenotype.intracellular`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.intracellular` type tag (`"maboss"` / `"not-maboss"`) | `Definition` | The receiver can allocate the correct local intracellular implementation from this tag and the local build | Must serialize the discriminator if intracellular models are enabled. |
| `phenotype.intracellular` model payload | `State` | Not from `Cell_Definition` alone | Must serialize for stateful intracellular models such as MaBoSS. |

#### 8.12 `phenotype.cell_interactions`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.cell_interactions.apoptotic_phagocytosis_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.cell_interactions.necrotic_phagocytosis_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.cell_interactions.other_dead_phagocytosis_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.cell_interactions.live_phagocytosis_rates` | `Definition` | From local definition | Definition-sized parameter vector. |
| `phenotype.cell_interactions.attack_rates` | `Definition` | From local definition | Definition-sized parameter vector. |
| `phenotype.cell_interactions.immunogenicities` | `Definition` | From local definition | Definition-sized parameter vector. |
| `phenotype.cell_interactions.attack_damage_rate` | `Definition` | From local definition | Parameter. |
| `phenotype.cell_interactions.total_damage_delivered` | `State` | Not from definition | Must serialize. |
| `phenotype.cell_interactions.attack_duration` | `Definition` | From local definition | Parameter. |
| `phenotype.cell_interactions.fusion_rates` | `Definition` | From local definition | Definition-sized parameter vector. |

Not packed in `phenotype.cell_interactions`:

- `phenotype.cell_interactions.pAttackTarget`

This pointer must not cross ranks as raw data.

#### 8.13 `phenotype.cell_transformations`

Overall mode: `Definition`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `phenotype.cell_transformations.transformation_rates` | `Definition` | From local definition | Definition-sized parameter vector. Serialize only if transformation rates can vary per cell. |

### 9. `is_out_of_domain`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `is_out_of_domain` | `State` | Constructor default is `false` | Must serialize. |

### 10. `is_movable`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `is_movable` | `Definition` | From `Cell_Definition.is_movable` | Redundant if the receiver instantiates from the correct local definition. |

### 11. `displacement`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `displacement` | `State` | Constructor default is zero | Serialize to preserve the in-flight mechanical update state. |

### 12. `Basic_Agent` state and caches

#### 12.1 `get_total_volume()`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `Basic_Agent::volume` via `get_total_volume()` | `Alias/Cache` | Can be set from `phenotype.volume.total` with `set_total_volume()` | Duplicates volume information already carried in `phenotype.volume.total`. |

#### 12.2 `get_is_volume_changed()`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `Basic_Agent::volume_is_changed` via `get_is_volume_changed()` | `Alias/Cache` | Can be derived from whether total volume changed during local setup | Runtime cache flag. |

#### 12.3 BioFVM source/sink solver temp vectors

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `cell_source_sink_solver_temp1` | `Alias/Cache` | Can be recomputed locally | Solver cache. |
| `cell_source_sink_solver_temp2` | `Alias/Cache` | Can be recomputed locally | Solver cache. |
| `cell_source_sink_solver_temp_export1` | `Alias/Cache` | Can be recomputed locally | Solver cache. |
| `cell_source_sink_solver_temp_export2` | `Alias/Cache` | Can be recomputed locally | Solver cache. |

#### 12.4 Velocity history and extracellular bookkeeping

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `previous_velocity` | `State` | Constructor default depends on `Basic_Agent` setup | Serialize if the integrator relies on velocity history. |
| `total_extracellular_substrate_change` | `Alias/Cache` | Can be recomputed locally | Runtime bookkeeping vector. |

#### 12.5 Pointer-backed BioFVM vectors duplicated from `phenotype.secretion`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `*secretion_rates` | `Alias/Cache` | Rebind to `phenotype.secretion.secretion_rates` | Duplicate of `phenotype.secretion`. |
| `*saturation_densities` | `Alias/Cache` | Rebind to `phenotype.secretion.saturation_densities` | Duplicate of `phenotype.secretion`. |
| `*uptake_rates` | `Alias/Cache` | Rebind to `phenotype.secretion.uptake_rates` | Duplicate of `phenotype.secretion`. |
| `*net_export_rates` | `Alias/Cache` | Rebind to `phenotype.secretion.net_export_rates` | Duplicate of `phenotype.secretion`. |

#### 12.6 Pointer-backed BioFVM vectors duplicated from `phenotype.molecular`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `*internalized_substrates` | `Alias/Cache` | Rebind to `phenotype.molecular.internalized_total_substrates` | Duplicate of `phenotype.molecular`. |
| `*fraction_released_at_death` | `Alias/Cache` | Rebind to `phenotype.molecular.fraction_released_at_death` | Duplicate of `phenotype.molecular`. |
| `*fraction_transferred_when_ingested` | `Alias/Cache` | Rebind to `phenotype.molecular.fraction_transferred_when_ingested` | Duplicate of `phenotype.molecular`. |

### 13. `type`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `type` | `Definition` | From `Cell_Definition.type` | This is the key field that should identify the local definition on the receiving rank. It is currently packed very late. |

### 14. `velocity`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `velocity` | `State` | Not from definition | Must serialize. |

## Important Fields That Exist But Are Not Packed By `Cell::pack()`

### `Cell`

- `container`
- `current_mechanics_voxel_index`
- `updated_current_mechanics_voxel_index`
- `crossed_to_left_subdomain`
- `crossed_to_right_subdomain`

These are local container or migration-management details and should be rebuilt on the receiving rank.

### `Basic_Agent`

- `microenvironment`
- `selected_microenvironment`
- `current_microenvironment_voxel_index`
- `current_voxel_index`
- `current_voxel_global_index`
- `is_active`
- `index`

These are local environment/container bindings or bookkeeping.

### Additional pointer fields rebuilt locally

- `parameters.pReference_live_phenotype`
- `phenotype.cycle.pCycle_Model`
- `phenotype.death.models`
- `phenotype.secretion.pMicroenvironment`
- `phenotype.molecular.pMicroenvironment`
- `phenotype.cell_interactions.pAttackTarget`

## Practical Split For A Future State-Only MPI Payload

If the goal is to minimize rank-to-rank traffic, the current pack order naturally splits into four groups:

1. Keep serialized:
   - `ID`, `position`, `velocity`
   - `state`
   - phenotype flags
   - active cycle state
   - death state
   - volume state
   - geometry or post-unpack geometry recomputation inputs
   - `phenotype.molecular.internalized_total_substrates`
   - `phenotype.cell_integrity.damage`
   - `phenotype.cell_interactions.total_damage_delivered`
   - intracellular dynamic state
   - `is_out_of_domain`, `displacement`, `previous_velocity`

2. Rebuild from local `Cell_Definition` after unpacking `type`:
   - `type_name`
   - most of `custom_data` schema
   - `parameters`
   - `functions`
   - most of `phenotype.mechanics`
   - most of `phenotype.motility`
   - most of `phenotype.secretion`
   - most of `phenotype.cell_interactions`
   - `phenotype.cell_transformations`
   - `is_movable`

3. Serialize only if the model allows per-cell mutation of what was originally definition data:
   - `custom_data` values
   - cycle transition rates
   - death rates
   - secretion vectors
   - mechanics parameters
   - motility enablement and sensitivities
   - molecular release/transfer fractions
   - damage kinetics

4. Rebind or recompute locally:
   - `Basic_Agent::volume`
   - `volume_is_changed`
   - BioFVM temporary solver vectors
   - duplicated BioFVM pointer vectors mirrored from `phenotype.secretion` and `phenotype.molecular`

The largest current redundancies are:

- `type_name` in addition to `type`
- full `functions.cycle_model`
- duplicated BioFVM vectors after `phenotype.secretion` and `phenotype.molecular`
- `Basic_Agent::volume` in addition to `phenotype.volume.total`
