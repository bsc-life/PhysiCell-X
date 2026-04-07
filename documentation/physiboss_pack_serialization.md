# PhysiBoSS Serialization Map

This document is the PhysiBoSS companion to [cell_pack_serialization.md](/home/jose/GPFS/development/PhysiCell-X/documentation/cell_pack_serialization.md). It follows the current PhysiBoSS serialization order used by the intracellular MaBoSS addon.

The entry path is:

1. `Phenotype::pack()` writes the intracellular type tag.
2. If the tag is `"maboss"`, `MaBoSSIntracellular::pack()` writes the addon payload.
3. `MaBoSSIntracellular::pack()` then delegates to nested serializers:
   - `MaBoSSInput::pack()`
   - `MaBoSSOutput::pack()`
   - `MaBoSSNetwork::pack()`

## Legend

- `Definition`: copied from the XML / cell-definition intracellular configuration. The receiver can rebuild it from the local `Cell_Definition` if ranks share the same definitions and MaBoSS files.
- `State`: true per-cell mutable state. This must be serialized to preserve a migrated PhysiBoSS cell.
- `Mixed`: starts from definition/config data but can become cell-specific during the run.
- `Alias/Cache`: runtime cache or local lookup data that can be rebound or recomputed from local registries and the already available definition/state.

## Where PhysiBoSS Data Enters A Cell

- `MaBoSSIntracellular::initialize_intracellular_from_pugixml()` loads the intracellular XML block into the `Cell_Definition` intracellular object.
- `Phenotype::operator=` clones `intracellular` through `Intracellular::clone()`.
- `MaBoSSIntracellular::clone()` uses the copy constructor:
  - copied: configuration fields such as filenames, timings, maps, and mappings
  - not copied: the live `maboss` engine/network state and `next_physiboss_run`
- `create_cell(Cell_Definition&)` calls `phenotype.intracellular->start()` on the new cell.
- `MaBoSSIntracellular::start()` initializes the MaBoSS engine, restarts node values, and sets `next_physiboss_run`.

That means a fresh cell created from a `Cell_Definition` gets the PhysiBoSS configuration from the definition, then generates its own runtime MaBoSS state locally. A migrated cell needs the runtime MaBoSS state serialized if exact continuity is required.

## Ordered Map Of PhysiBoSS Serialization

### 0. `Phenotype` intracellular type tag

Before `MaBoSSIntracellular::pack()` runs, `Phenotype::pack()` serializes a type discriminator.

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| intracellular type tag (`"maboss"` / `"not-maboss"`) | `Definition` | The receiver uses it to allocate the correct intracellular subclass | Must serialize the discriminator before the addon payload. |

## `MaBoSSIntracellular::pack()` In Order

### 1. `bnd_filename`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `bnd_filename` | `Definition` | From intracellular XML in the local `Cell_Definition` | Redundant if the receiver instantiates from the correct local definition. |

### 2. `cfg_filename`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `cfg_filename` | `Definition` | From intracellular XML in the local `Cell_Definition` | Redundant if the receiver instantiates from the correct local definition. |

### 3. Timing and execution settings

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `time_step` | `Definition` | From intracellular XML | Definition-level config. |
| `discrete_time` | `Definition` | From intracellular XML | Definition-level config. |
| `time_tick` | `Definition` | From intracellular XML | Definition-level config. |
| `scaling` | `Definition` | From intracellular XML | Definition-level config. |
| `time_stochasticity` | `Definition` | From intracellular XML | Definition-level config. |

These are later re-applied to the live `MaBoSSNetwork` after unpacking.

### 4. Inheritance settings

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `inherit_state` | `Definition` | From intracellular XML | Division/inheritance policy, not dynamic migrated-cell state. |
| `inherit_nodes.size()` | `Definition` | From intracellular XML | Definition metadata. |
| `inherit_nodes[key]` | `Definition` | From intracellular XML | Definition metadata. |

### 5. `start_time`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `start_time` | `Definition` | From intracellular XML | Used by `start()` to schedule the first run. |

### 6. `initial_values`

Overall mode: `Definition`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `initial_values.size()` | `Definition` | From intracellular XML | Definition metadata. |
| `initial_values[node_name]` | `Definition` | From intracellular XML | Initial node-value override probabilities. |

These are used when a cell is started or restarted. They are not the current live node states of a migrated cell.

### 7. `mutations`

Overall mode: `Definition`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `mutations.size()` | `Definition` | From intracellular XML | Definition metadata. |
| `mutations[node_name]` | `Definition` | From intracellular XML | Definition-level mutation setup. |

These are re-applied to the local `MaBoSSNetwork` after unpacking.

### 8. `parameters`

Overall mode: `Mixed`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `parameters.size()` | `Definition` | From intracellular XML | Definition metadata. |
| `parameters[param_name]` | `Mixed` | Default comes from intracellular XML | This outer map is configuration-level parameter override data. Serialize if code may mutate it per cell. |

Important distinction:

- `MaBoSSIntracellular::parameters` is the config map loaded from XML.
- The live parameter values actually used by the simulation are serialized later inside `maboss` through the symbol table payload.

## `listOfInputs`

Overall mode: `Mixed`

The map structure comes from XML, but one field inside each `MaBoSSInput` is true runtime state.

### 9. `listOfInputs` container

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `listOfInputs.size()` | `Definition` | From intracellular XML | Definition metadata. |
| input map key | `Definition` | From intracellular XML | Duplicates the intracellular target name already represented in the `MaBoSSInput` object. |

### 10. Each `MaBoSSInput`

Pack order inside `MaBoSSInput::pack()`:

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `physicell_name` | `Definition` | From intracellular XML mapping | Mapping metadata. |
| `type` | `Definition` | From XML mapping constructor | Distinguishes node-input vs parameter-input mapping. |
| `intracellular_name` | `Definition` | From intracellular XML mapping | Mapping metadata. |
| `intracellular_parameter` | `Definition` | From intracellular XML mapping | Mapping metadata. |
| `action` | `Definition` | From intracellular XML mapping | Mapping metadata. |
| `threshold` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `inact_threshold` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `scaling` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `smoothing` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `smoothed_value` | `State` | Starts at `0` on a new cell | Runtime filtered signal state. Must serialize for exact continuity if smoothing is in use. |
| `use_for_dead` | `Definition` | From intracellular XML mapping | Mapping parameter. |

`smoothed_value` is the key field that makes `MaBoSSInput` a mixed object.

### 11. `indicesOfInputs`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `indicesOfInputs` | `Alias/Cache` | Can be rebuilt from `listOfInputs` with `find_signal_index(...)` | Runtime lookup cache. Current code serializes it, but conceptually it should be recomputed locally. |

## `listOfOutputs`

Overall mode: `Mixed`

The mapping definition comes from XML, but `probability` and `initialized` hold live smoothed output state.

### 12. `listOfOutputs` container

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `listOfOutputs.size()` | `Definition` | From intracellular XML | Definition metadata. |
| output map key | `Definition` | From intracellular XML | Duplicates the PhysiCell behavior name already present in the `MaBoSSOutput` object. |

### 13. Each `MaBoSSOutput`

Pack order inside `MaBoSSOutput::pack()`:

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `physicell_name` | `Definition` | From intracellular XML mapping | Mapping metadata. |
| `intracellular_name` | `Definition` | From intracellular XML mapping | Mapping metadata. |
| `action` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `value` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `base_value` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `smoothing` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `probability` | `State` | Starts at `0.5` for a new mapping | Runtime smoothed activation probability. Serialize for exact continuity if smoothing is in use. |
| `initialized` | `State` | Starts at `false` for a new mapping | Runtime initialization flag for the smoothing state. |
| `steepness` | `Definition` | From intracellular XML mapping | Mapping parameter. |
| `use_for_dead` | `Definition` | From intracellular XML mapping | Mapping parameter. |

### 14. `indicesOfOutputs`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `indicesOfOutputs` | `Alias/Cache` | Can be rebuilt from `listOfOutputs` with `find_behavior_index(...)` | Runtime lookup cache. Current code serializes it, but conceptually it should be recomputed locally. |

## `maboss.pack()`

This is the live MaBoSS simulation payload. It is the part that most directly represents the intracellular state of a migrated cell.

### 15. `time_to_update`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `maboss.time_to_update` | `State` | Not from definition | Must serialize to preserve the stochastic interval until the next MaBoSS run. |

### 16. Network node-state block

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| node count | `Definition` | Implied by the local BND/CFG model | Serialized as a consistency check. The receiver assumes the same network topology and node order. |
| per-node current Boolean state | `State` | A fresh cell would regenerate these from initial conditions, not from the migrated cell | Must serialize for exact continuity. |

The sender and receiver rely on identical node ordering from the local MaBoSS network built with the same `bnd_filename` and `cfg_filename`.

### 17. Symbol-table block

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| symbol-table count | `Definition` | Implied by the local BND/CFG model | Serialized as a consistency check. |
| symbol-table values in `getSymbolsNames()` order | `Mixed` | Default symbol values come from the local config files and outer parameter overrides | Must serialize if any parameter or symbol can be changed at runtime. |

This block is important because it captures live MaBoSS parameter values even if they no longer match the outer `parameters` map.

## 18. `next_physiboss_run`

| Packed field | Mode | Receiver init | MPI note |
|---|---|---|---|
| `next_physiboss_run` | `State` | For a fresh cell it is set in `start()` from `start_time` and current simulation time | Must serialize to preserve the exact future scheduling of the next PhysiBoSS update. |

## Important Fields That Exist But Are Not Packed

### `MaBoSSIntracellular`

- `intracellular_type`
  - This is serialized one level higher by `Phenotype::pack()`, not inside `MaBoSSIntracellular::pack()`.
- `counter`
  - Static bookkeeping, not per-cell state.

### `MaBoSSNetwork`

- `network`
- `config`
- `engine`
- `output_mask`
- `nodesByName`
- `parametersByName`
- `update_time_step`
- `initial_values`
- `mutations`
- `time_stochasticity`
- `scaling`

These are rebuilt or reapplied locally through:

- `maboss.init_maboss(bnd_filename, cfg_filename)`
- `maboss.unpack(...)`
- `maboss.mutate(mutations)`
- `maboss.set_time_stochasticity(time_stochasticity)`
- `maboss.set_scaling(scaling)`
- `maboss.set_discrete_time(discrete_time, time_tick)`
- `maboss.set_update_time_step(time_step)`

Some of these are duplicated conceptually between `MaBoSSIntracellular` and `MaBoSSNetwork`. The outer object carries config; the inner object carries live runtime state.

## Practical Split For A Future State-Only PhysiBoSS Payload

If the goal is to minimize MPI traffic while preserving the current addon behavior, the PhysiBoSS data naturally splits like this:

1. Keep serialized:
   - `maboss.time_to_update`
   - current MaBoSS node states
   - live MaBoSS symbol-table values
   - `next_physiboss_run`
   - `MaBoSSInput.smoothed_value` when smoothing is used
   - `MaBoSSOutput.probability`
   - `MaBoSSOutput.initialized`

2. Rebuild from local intracellular definition:
   - `bnd_filename`, `cfg_filename`
   - time settings
   - inheritance settings
   - `start_time`
   - `initial_values`
   - `mutations`
   - most of `parameters`
   - input and output mapping definitions

3. Recompute locally:
   - `indicesOfInputs`
   - `indicesOfOutputs`
   - MaBoSS pointers and lookup tables
   - node count and symbol-table count

## Current Redundancies

The largest addon redundancies in the current serialization are:

- outer `parameters` map plus inner live symbol-table values
- XML-derived input/output mapping definitions plus cached `indicesOfInputs` / `indicesOfOutputs`
- `bnd_filename` / `cfg_filename` and other XML config that the receiver could rebuild from the local `Cell_Definition`

The genuinely essential migrated PhysiBoSS state is much smaller:

- current Boolean node states
- live symbol-table values
- update scheduling
- smoothed input/output runtime memory
