# LLM Guide: Convert PhysiCell / PhysiCelldFBA OpenMP Sample Projects to PhysiCell-X MPI + OpenMP

This document is an instruction set for an LLM or coding agent that must transform an OpenMP-only sample project from `PhysiCell` or `PhysiCelldFBA` into a `PhysiCell-X` hybrid MPI + OpenMP project.

The target implementation style is the one used by `PhysiCell-X/sample_projects/*_mpi` and `PhysiCell-X/sample_projects_intracellular/*/*_mpi`:

- MPI process-level decomposition through `DistPhy::mpi`.
- OpenMP threads inside each MPI process through `omp_set_num_threads(...)`.
- MPI-aware microenvironment setup, cell container creation, diffusion, mechanics/phenotype update, output, and rank-gated logging.
- Project-specific biology preserved from the source OpenMP sample.

## Scope

Use these instructions when converting:

- `PhysiCell/sample_projects/<project>`
- `PhysiCelldFBA/sample_projects/<project>`
- `PhysiCelldFBA/sample_projects_intracellular/fba/<project>`
- Other OpenMP-only PhysiCell-style sample projects with `main.cpp`, `Makefile`, `config/`, and `custom_modules/`

into:

- `PhysiCell-X/sample_projects/<project>_mpi`, or
- `PhysiCell-X/sample_projects_intracellular/fba/<project>_mpi` for dFBA projects.

Do not treat this as a text-only rename. The MPI build changes runtime ownership of cells, voxel locality, output functions, and initialization semantics.

## Required Inputs

Before editing, inspect these files:

- Source project:
  - `main.cpp`
  - `Makefile`
  - `custom_modules/*.h`
  - `custom_modules/*.cpp`
  - `config/PhysiCell_settings.xml`
  - any `config/*.csv`, `config/*.xml`, SBML/FBA models, rule files, or scripts referenced by XML/custom code
- Target examples in `PhysiCell-X`:
  - `sample_projects/template_mpi/main.cpp`
  - `sample_projects/template_mpi/custom_modules/custom.h`
  - `sample_projects/template_mpi/custom_modules/custom.cpp`
  - `sample_projects/template_mpi/Makefile`
  - for dFBA: `sample_projects_intracellular/fba/dfba_unit_test/main.cpp`
  - for dFBA: `sample_projects_intracellular/fba/dfba_unit_test/Makefile`
  - for a project with an existing MPI analogue, prefer that analogue over the generic template.

## High-Level Migration Algorithm

1. Create a new target project directory under `PhysiCell-X`.
2. Copy source biology files from the OpenMP project: `config/`, `custom_modules/`, scripts, SBML/FBA files, and project-specific `main` logic.
3. Use a `PhysiCell-X` MPI sample as the structural base, not the OpenMP source main.
4. Port source custom biology into MPI-compatible custom hooks.
5. Replace serial PhysiCell APIs with MPI-aware PhysiCell-X APIs.
6. Gate filesystem writes and human-readable logging to the MPI I/O rank.
7. Ensure cell placement creates only cells owned by the local MPI rank.
8. Keep OpenMP enabled and configured through `PhysiCell_settings.omp_num_threads`.
9. Build with `mpic++` and include `DistPhy` and MPI-specific PhysiCell objects.
10. Validate with one MPI rank, then multiple MPI ranks, then compare biological outputs to the original OpenMP run.

## Directory Strategy

For a normal sample project:

```text
PhysiCell-X/sample_projects/<project>_mpi/
  Makefile
  main.cpp
  main-serial.cpp            optional, useful reference copy of source main
  config/
  custom_modules/
  scripts/
```

For a dFBA sample project:

```text
PhysiCell-X/sample_projects_intracellular/fba/<project>_mpi/
  Makefile
  main.cpp
  main-serial.cpp            optional
  config/
  custom_modules/
  scripts/
```

If the converted project is meant to be installed at the `PhysiCell-X` repository root, still use the same MPI API substitutions.

## Main Program Conversion

### Add MPI Headers and Namespace

Add:

```cpp
#include <mpi.h>
```

Keep:

```cpp
#include <omp.h>
#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"
#include "./custom_modules/custom.h"
```

Add:

```cpp
using namespace DistPhy::mpi;
```

### Initialize MPI Before Reading Configuration

At the start of `main`:

```cpp
mpi_Environment world;
world.Initialize();

mpi_Cartesian cart_topo;
cart_topo.Build_Cartesian_Topology(world);
cart_topo.Find_Cartesian_Coordinates(world);
cart_topo.Find_Left_Right_Neighbours(world);
```

Do this before `load_PhysiCell_config_file(...)`, microenvironment setup, or any code that may depend on rank count or domain topology.

### Load Configuration With the MPI-Aware API

Replace:

```cpp
XML_status = load_PhysiCell_config_file(argv[1]);
XML_status = load_PhysiCell_config_file("./config/PhysiCell_settings.xml");
```

with:

```cpp
XML_status = load_PhysiCell_config_file(argv[1], world);
XML_status = load_PhysiCell_config_file("./config/PhysiCell_settings.xml", world);
```

If config files are copied to output, only the I/O rank should run `system(...)`:

```cpp
if (IOProcessor(world))
{
    system(copy_command);
}
```

Do not allow every rank to write the same file simultaneously.

### Preserve OpenMP Setup

Keep:

```cpp
omp_set_num_threads(PhysiCell_settings.omp_num_threads);
```

This controls OpenMP threads per MPI rank. If running with `N` MPI ranks and `T` OpenMP threads, total execution slots are roughly `N * T`.

### Use MPI-Aware Microenvironment and Cell Container Setup

Replace:

```cpp
setup_microenvironment();
Cell_Container* cell_container =
    create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size);
```

with:

```cpp
setup_microenvironment(world, cart_topo);
Cell_Container* cell_container =
    create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size, world, cart_topo);
```

### Use MPI-Aware Custom Setup Hooks

Replace:

```cpp
create_cell_types();
setup_tissue();
```

with one of these patterns, depending on the target custom signatures:

```cpp
create_cell_types(world, cart_topo);
setup_tissue(microenvironment, world, cart_topo);
```

or, for dFBA samples that use the shorter current local pattern:

```cpp
create_cell_types(world, cart_topo);
setup_tissue(world, cart_topo);
```

Choose one pattern and keep `custom.h`, `custom.cpp`, and `main.cpp` consistent.

### Use MPI-Aware Output

Replace serial MultiCellDS output:

```cpp
save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);
```

with the MPI overload available in the target branch:

```cpp
save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time, world, cart_topo);
```

or, in projects using the pugi XML writer:

```cpp
save_PhysiCell_to_MultiCellDS_xml_pugi(filename, microenvironment, PhysiCell_globals.current_time, world, cart_topo);
```

Use whichever function is used by the closest working `PhysiCell-X` sample. Do not mix output APIs casually in one project.

Replace serial SVG output:

```cpp
SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time,
         cell_coloring_function, substrate_coloring_function);
```

with:

```cpp
SVG_plot_mpi(filename, microenvironment, 0.0, PhysiCell_globals.current_time,
             cell_coloring_function, substrate_coloring_function, world, cart_topo);
```

Only the I/O rank should create legends and print citations:

```cpp
if (IOProcessor(world))
{
    create_plot_legend(filename, cell_coloring_function);
    display_citations();
}
```

### Use MPI-Aware Simulation Status and Legacy Logs

Replace:

```cpp
display_simulation_status(std::cout);
log_output(time, index, microenvironment, report_file);
```

with:

```cpp
display_simulation_status(std::cout, world, cart_topo);
log_output(time, index, microenvironment, report_file, world, cart_topo);
```

If the project opens `simulation_report.txt`, only open and write it on `IOProcessor(world)`.

### Use MPI-Aware Diffusion and Cell Updates

Replace:

```cpp
microenvironment.simulate_diffusion_decay(diffusion_dt);
((Cell_Container*) microenvironment.agent_container)->update_all_cells(
    PhysiCell_globals.current_time);
```

with:

```cpp
microenvironment.simulate_diffusion_decay(diffusion_dt, world, cart_topo);
((Cell_Container*) microenvironment.agent_container)->update_all_cells(
    PhysiCell_globals.current_time, world, cart_topo);
```

### Finalize MPI

At every normal exit path:

```cpp
world.Finalize();
return 0;
```

If configuration loading fails after MPI initialization, finalize before exiting where practical:

```cpp
if (!XML_status)
{
    world.Finalize();
    return -1;
}
```

## Custom Module Conversion

### Include DistPhy Helpers

In MPI custom modules, include the MPI helper headers used by the local samples:

```cpp
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"
```

Add:

```cpp
using namespace DistPhy::mpi;
```

### Add MPI Prototypes to `custom.h`

For ordinary projects, keep the serial prototypes if the file still supports serial compilation, and add MPI overloads:

```cpp
void create_cell_types(void);
void setup_tissue(void);
void setup_microenvironment(void);

void create_cell_types(mpi_Environment& world, mpi_Cartesian& cart_topo);
void setup_tissue(Microenvironment& m, mpi_Environment& world, mpi_Cartesian& cart_topo);
void setup_microenvironment(mpi_Environment& world, mpi_Cartesian& cart_topo);
```

For dFBA projects that are MPI-only, the target may use:

```cpp
void create_cell_types(mpi_Environment& world, mpi_Cartesian& cart_topo);
void setup_tissue(mpi_Environment& world, mpi_Cartesian& cart_topo);
void setup_microenvironment(mpi_Environment& world, mpi_Cartesian& cart_topo);
```

Match `main.cpp` exactly.

### Convert `create_cell_types`

Replace serial XML initialization:

```cpp
initialize_cell_definitions_from_pugixml();
setup_signal_behavior_dictionaries();
setup_cell_rules();
display_cell_definitions(std::cout);
```

with MPI-aware calls:

```cpp
initialize_cell_definitions_from_pugixml(world, cart_topo);
setup_signal_behavior_dictionaries(world, cart_topo);
setup_cell_rules(&world, &cart_topo);       // if this signature exists in target branch
display_cell_definitions(std::cout, world, cart_topo);
```

Some `PhysiCell-X` samples use a shared helper that accepts nullable MPI pointers:

```cpp
setup_cell_rules(world_ptr, cart_topo_ptr);
if (world_ptr == nullptr || IOProcessor(*world_ptr))
{
    display_cell_definitions(std::cout);
}
```

Use the pattern already present in the nearest target sample.

Set MPI-capable velocity hooks if the source sets velocity behavior:

```cpp
cell_defaults.functions.update_velocity = standard_update_cell_velocity;
cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;
```

Preserve all source project biology:

- custom phenotype functions
- contact functions
- cell interactions
- rule setup
- intracellular pre/post hooks
- custom data initialization
- custom cell definition changes

Do not discard project-specific cell function assignments when replacing boilerplate.

### Convert `setup_microenvironment`

Replace:

```cpp
initialize_microenvironment();
```

with:

```cpp
initialize_microenvironment(world, cart_topo);
```

If the source customizes Dirichlet conditions, initial substrate gradients, or density vectors before initialization, preserve those edits. If it writes directly to voxel data after initialization, ensure it only writes local voxels or use MPI-safe helper patterns from an existing `PhysiCell-X` sample.

PhysiCell-X works as a 3-D MPI domain decomposition. For nominal 2-D samples, keep `default_microenvironment_options.simulate_2D = true` semantics if the source uses them, but do not assume a truly 2-D MPI mesh unless the target API explicitly supports it.

### Convert Cell Placement in `setup_tissue`

Serial code often does this:

```cpp
Cell* pC = create_cell(*pCD);
pC->assign_position(position);
```

MPI code must create cells on the rank that owns their position. There are two acceptable patterns.

#### Pattern A: Root Generates Global Positions, Then Distributes

Use this for deterministic equivalence with serial random placement or CSV-derived global positions.

```cpp
mpi_CellPositions cp;
mpi_MyCells mc;

if (world.rank == 0)
{
    std::vector<std::vector<double>> positions;
    // Fill positions exactly as the serial project would.

    int start_cell_ID = Basic_Agent::get_max_ID_in_parallel();
    cp.positions_to_rank_list(
        positions,
        m.mesh.bounding_box[0], m.mesh.bounding_box[3],
        m.mesh.bounding_box[1], m.mesh.bounding_box[4],
        m.mesh.bounding_box[2], m.mesh.bounding_box[5],
        m.mesh.dx, m.mesh.dy, m.mesh.dz,
        world, cart_topo, start_cell_ID);

    Basic_Agent::set_max_ID_in_parallel(start_cell_ID + positions.size());
}

distribute_cell_positions(cp, mc, world, cart_topo);

for (int i = 0; i < mc.my_no_of_cell_IDs; i++)
{
    Cell* pC = create_cell(*pCD, mc.my_cell_IDs[i]);
    pC->assign_position(
        mc.my_cell_coords[3*i],
        mc.my_cell_coords[3*i+1],
        mc.my_cell_coords[3*i+2],
        world, cart_topo);
}
```

This pattern is usually best for biological equivalence because random numbers are drawn once on rank 0, in serial order.

#### Pattern B: Each Rank Generates Local Cells in Its Local Bounding Box

Use this only when exact serial equivalence is not required and the source biology is naturally spatially uniform.

```cpp
double Xmin = microenvironment.mesh.local_bounding_box[0];
double Ymin = microenvironment.mesh.local_bounding_box[1];
double Zmin = microenvironment.mesh.local_bounding_box[2];
double Xmax = microenvironment.mesh.local_bounding_box[3];
double Ymax = microenvironment.mesh.local_bounding_box[4];
double Zmax = microenvironment.mesh.local_bounding_box[5];

int local_cells = parameters.ints("number_of_cells") / world.size
    + (world.rank < parameters.ints("number_of_cells") % world.size ? 1 : 0);

for (int n = 0; n < local_cells; n++)
{
    std::vector<double> position(3, 0.0);
    position[0] = Xmin + UniformRandom() * (Xmax - Xmin);
    position[1] = Ymin + UniformRandom() * (Ymax - Ymin);
    position[2] = Zmin + UniformRandom() * (Zmax - Zmin);

    Cell* pC = create_cell(*pCD);
    pC->assign_position(position, world, cart_topo);
}
```

This pattern can change stochastic trajectories because random draws are split across MPI ranks.

### Convert CSV Cell Loading

Replace:

```cpp
load_cells_from_pugixml();
```

with:

```cpp
load_cells_from_pugixml(world, cart_topo);
```

This is mandatory when `config/PhysiCell_settings.xml` references `config/cells.csv` or another initial cell file.

### Assign Position With MPI Context

Replace every serial placement:

```cpp
pC->assign_position(position);
pC->assign_position(x, y, z);
```

with:

```cpp
pC->assign_position(position, world, cart_topo);
pC->assign_position(x, y, z, world, cart_topo);
```

Only create the cell on the rank that owns the position. Do not create every cell on every rank.

### Gate Logging

Replace unconditional printing in initialization paths:

```cpp
std::cout << "Placing cells ..." << std::endl;
display_cell_definitions(std::cout);
```

with:

```cpp
if (IOProcessor(world))
{
    std::cout << "Placing cells ..." << std::endl;
}
```

Use MPI-aware display helpers when available:

```cpp
display_cell_definitions(std::cout, world, cart_topo);
```

Do not gate per-cell coloring functions unless they perform output. Coloring functions are called during rendering and must return valid colors on every rank participating in SVG generation. Remove debug prints from coloring functions unless rank-gated and explicitly needed.

## Makefile Conversion

### Use MPI Compiler

Replace:

```make
CC := g++
```

with:

```make
CC := mpic++
```

Keep support for `PHYSICELL_CPP` if present.

### Keep OpenMP

Ensure compile flags include:

```make
-fopenmp
```

For PhysiCell-X samples, a typical value is:

```make
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -m64 -std=c++11 -g -fopenmp
```

### Add DistPhy Objects

Add:

```make
DistPhy_OBJECTS := DistPhy_Environment.o DistPhy_Cartesian.o DistPhy_Utils.o DistPhy_Collective.o
```

Add these objects into `PhysiCell_OBJECTS`:

```make
PhysiCell_OBJECTS := $(BioFVM_OBJECTS) $(pugixml_OBJECTS) \
    $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS) $(DistPhy_OBJECTS)
```

Add compile rules:

```make
DistPhy_Environment.o: ./DistPhy/DistPhy_Environment.cpp
	$(COMPILE_COMMAND) -c ./DistPhy/DistPhy_Environment.cpp

DistPhy_Cartesian.o: ./DistPhy/DistPhy_Cartesian.cpp
	$(COMPILE_COMMAND) -c ./DistPhy/DistPhy_Cartesian.cpp

DistPhy_Utils.o: ./DistPhy/DistPhy_Utils.cpp
	$(COMPILE_COMMAND) -c ./DistPhy/DistPhy_Utils.cpp

DistPhy_Collective.o: ./DistPhy/DistPhy_Collective.cpp
	$(COMPILE_COMMAND) -c ./DistPhy/DistPhy_Collective.cpp
```

### Add MPI PhysiCell Core Objects

Ensure `PhysiCell_core_OBJECTS` includes:

```make
PhysiCell_cell_mpi.o
PhysiCell_standard_models_mpi.o
PhysiCell_cell_container_mechanics.o
```

Add compile rules:

```make
PhysiCell_cell_mpi.o: ./core/PhysiCell_cell_mpi.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_mpi.cpp

PhysiCell_standard_models_mpi.o: ./core/PhysiCell_standard_models_mpi.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models_mpi.cpp

PhysiCell_cell_container_mechanics.o: ./core/PhysiCell_cell_container_mechanics.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container_mechanics.cpp
```

### Add MPI-Aware Modules Expected by PhysiCell-X

Compare against `PhysiCell-X/sample_projects/template_mpi/Makefile`. MPI projects often include:

```make
PhysiCell_pathology_SVG.o
PhysiCell_geometry.o
```

Use the object list and rules from the closest working `PhysiCell-X` sample.

### dFBA Makefile Requirements

For dFBA projects, merge the `PhysiCell-X` MPI Makefile with the dFBA sections from the source `PhysiCelldFBA` Makefile.

Use:

```make
CC := mpic++
INC := -DADDON_PHYSIDFBA -I$(EXT_DIR)/libsbml/include -I$(EXT_DIR)/coin-or/include -I$(CUR_DIR)/addons/dFBA/src
PhysiCelldFBA_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o dfba_intracellular.o
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS) $(PhysiCelldFBA_OBJECTS)
```

Prefer static bzip2 detection if the target environment has runtime library problems:

```make
BZ2_STATIC_LIB ?= $(shell $(CC) -print-file-name=libbz2.a 2>/dev/null)

LIB := -L$(EXT_DIR)/libsbml/lib \
    $(EXT_DIR)/coin-or/lib/libClp-static.a \
    $(EXT_DIR)/coin-or/lib/libCoinUtils-static.a \
    -lsbml-static -lxml2 $(BZ2_STATIC_LIB) -lz -ggdb -llapack -no-pie
```

Add a check target if using static bzip2:

```make
.PHONY: check-static-bz2

check-static-bz2:
	@if [ "$(BZ2_STATIC_LIB)" = "libbz2.a" ] || [ ! -f "$(BZ2_STATIC_LIB)" ]; then \
		echo "ERROR: static libbz2.a was not found by $(CC)."; \
		echo "Load a bzip2 development module or run: make BZ2_STATIC_LIB=/path/to/libbz2.a"; \
		exit 1; \
	fi
```

Then make `all` depend on it:

```make
all: check-static-bz2 main.cpp $(ALL_OBJECTS) $(libFBA)
	$(COMPILE_COMMAND) $(INC) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp $(LIB)
	make name
```

Compile dFBA-dependent objects with `$(INC)`:

```make
PhysiCell_cell.o: ./core/PhysiCell_cell.cpp $(libFBA)
	$(COMPILE_COMMAND) $(INC) -c ./core/PhysiCell_cell.cpp

dfba_intracellular.o: addons/dFBA/src/dfba_intracellular.cpp $(libFBA)
	$(COMPILE_COMMAND) $(INC) -c addons/dFBA/src/dfba_intracellular.cpp

custom.o: ./custom_modules/custom.cpp
	$(COMPILE_COMMAND) $(INC) -c ./custom_modules/custom.cpp
```

## dFBA-Specific Code Conversion

### Preserve Intracellular Hooks

If the source assigns pre/post intracellular hooks, preserve them:

```cpp
cell_defaults.functions.pre_update_intracellular = NULL;
cell_defaults.functions.post_update_intracellular = post_update_intracellular;
cell_defaults.functions.update_phenotype = NULL;
cell_defaults.functions.volume_update_function = NULL;

for (auto* pCD : cell_definitions_by_index)
{
    pCD->functions.pre_update_intracellular = NULL;
    pCD->functions.post_update_intracellular = post_update_intracellular;
    pCD->functions.update_phenotype = NULL;
    pCD->functions.volume_update_function = NULL;
}
```

### Start Intracellular Models on Local Cells

After MPI-aware cell loading/placement:

```cpp
for (Cell* pCell : *all_cells)
{
    if (pCell->phenotype.intracellular != NULL)
    {
        pCell->phenotype.intracellular->start();
    }
}
```

This loop runs on each rank and only starts intracellular models for local cells.

### Keep FBA Flux Extraction Local

Custom per-cell post-update logic remains local:

```cpp
#ifdef ADDON_PHYSIDFBA
if (pCell->phenotype.intracellular != NULL &&
    pCell->phenotype.intracellular->intracellular_type == "dfba")
{
    PhysiCelldFBA::dFBAIntracellular* dfba =
        static_cast<PhysiCelldFBA::dFBAIntracellular*>(pCell->phenotype.intracellular);

    pCell->custom_data["growth_rate"] = dfba->get_growth_rate();
    pCell->custom_data["oxygen_flux"] = dfba->get_flux_value("R_EX_o2_e");
    pCell->custom_data["glucose_flux"] = dfba->get_flux_value("R_EX_glc__D_e");
    pCell->custom_data["acetate_flux"] = dfba->get_flux_value("R_EX_ac_e");
}
#endif
```

Do not gather FBA objects to root. Each rank owns and updates its local cells.

### Avoid Unranked Debug Output in FBA Coloring

Serial dFBA samples often print fluxes from color functions. Remove or gate those prints because `SVG_plot_mpi` may call coloring on many cells/ranks:

```cpp
if (IOProcessor(world))
{
    std::cout << "flux_value: " << flux_value << std::endl;
}
```

If `world` is not available inside the coloring function, prefer removing the debug print.

## Common API Replacement Table

| Serial OpenMP API | PhysiCell-X MPI + OpenMP API |
|---|---|
| `load_PhysiCell_config_file(path)` | `load_PhysiCell_config_file(path, world)` |
| `initialize_microenvironment()` | `initialize_microenvironment(world, cart_topo)` |
| `setup_microenvironment()` | `setup_microenvironment(world, cart_topo)` |
| `create_cell_container_for_microenvironment(m, dx)` | `create_cell_container_for_microenvironment(m, dx, world, cart_topo)` |
| `initialize_cell_definitions_from_pugixml()` | `initialize_cell_definitions_from_pugixml(world, cart_topo)` |
| `setup_signal_behavior_dictionaries()` | `setup_signal_behavior_dictionaries(world, cart_topo)` |
| `setup_cell_rules()` | `setup_cell_rules(&world, &cart_topo)` or target branch equivalent |
| `display_cell_definitions(std::cout)` | `display_cell_definitions(std::cout, world, cart_topo)` or gate with `IOProcessor(world)` |
| `create_cell(*pCD)` then `assign_position(pos)` | create only on owner rank, then `assign_position(pos, world, cart_topo)` |
| `load_cells_from_pugixml()` | `load_cells_from_pugixml(world, cart_topo)` |
| `save_PhysiCell_to_MultiCellDS_v2(..., time)` | `save_PhysiCell_to_MultiCellDS_v2(..., time, world, cart_topo)` |
| `save_PhysiCell_to_MultiCellDS_xml_pugi(..., time)` | `save_PhysiCell_to_MultiCellDS_xml_pugi(..., time, world, cart_topo)` |
| `SVG_plot(...)` | `SVG_plot_mpi(..., world, cart_topo)` |
| `display_simulation_status(std::cout)` | `display_simulation_status(std::cout, world, cart_topo)` |
| `log_output(..., report_file)` | `log_output(..., report_file, world, cart_topo)` |
| `microenvironment.simulate_diffusion_decay(dt)` | `microenvironment.simulate_diffusion_decay(dt, world, cart_topo)` |
| `update_all_cells(time)` | `update_all_cells(time, world, cart_topo)` |
| `std::cout` during global setup | `if (IOProcessor(world)) std::cout ...` |
| `system("cp ...")` | `if (IOProcessor(world)) system("cp ...")` |

## Filesystem and Output Rules

Only `IOProcessor(world)` should:

- copy input config files into output
- create `legend.svg`
- open/write `simulation_report.txt`
- print citations
- print simulation headers or project summaries
- perform cleanup or directory creation if the code does it manually

MPI-aware PhysiCell-X output functions handle distributed `.xml`, `.mat`, SVG, and graph outputs. Do not wrap `save_PhysiCell...` or `SVG_plot_mpi` in `IOProcessor(world)` unless the target sample does so; those functions need the MPI context to gather/write correctly.

## Biological Equivalence Rules

The conversion should preserve biological behavior, not byte-level output identity.

Preserve:

- XML model parameters
- cell definitions
- substrate names and ordering
- custom data names and units
- random seed semantics when possible
- initial positions from CSV exactly
- phenotype functions
- custom rules
- contact functions
- intracellular model paths and mappings
- save intervals and simulation time

Expected differences:

- cell IDs may differ across rank counts
- output file byte order may differ
- cell ordering in `.mat` files may differ
- voxel ordering may differ
- local floating-point trajectories may differ slightly
- wall time and runtime metadata will differ

For equivalence checks, compare after rematching:

- cells by physical position or by stable biological ID if available
- voxels by `(x, y, z)`
- cell counts over time
- live/dead counts
- total and per-cell volume distributions
- custom metabolic state such as `growth_rate`, `oxygen_flux`, `glucose_flux`, `acetate_flux`
- substrate domain means and biologically relevant local ranges

Do not fail the migration solely because `diff -qr outputA outputB` reports differences.

## Validation Procedure

### Build

From the target `PhysiCell-X` project directory or repository root, run:

```bash
make clean
make
```

For dFBA, ensure `addons/dFBA/ext` exists. If dependencies are missing, run the existing setup script:

```bash
python3 beta/setup_fba.py
```

Then rebuild.

### One-Rank Smoke Test

Run with one MPI rank first:

```bash
mpirun -np 1 ./<program_name> ./config/PhysiCell_settings.xml
```

Check:

- simulation reaches final time
- no duplicate output races
- initial/final XML are written
- SVG creation works if enabled
- cell count matches source OpenMP run

### Multi-Rank Test

Run with multiple ranks and controlled OpenMP threads:

```bash
OMP_NUM_THREADS=<threads_per_rank> mpirun -np <ranks> ./<program_name> ./config/PhysiCell_settings.xml
```

If the code calls `omp_set_num_threads(PhysiCell_settings.omp_num_threads)`, the XML setting overrides `OMP_NUM_THREADS`.

Check:

- all ranks finalize cleanly
- no rank hangs during diffusion, cell updates, or output
- output files are valid
- total agents are conserved or evolve as expected
- no obvious duplicate cell creation

### Biological Comparison

Compare source and target with scripts that:

1. Load each `output%08d_cells.mat`.
2. Sort or match cells by position or stable ID.
3. Compare relevant rows.
4. Load each `output%08d_microenvironment0.mat`.
5. Sort voxels by `(x, y, z)`.
6. Compare substrate rows.

For dFBA projects, explicitly track:

- `growth_rate`
- `oxygen_flux`
- `glucose_flux`
- `acetate_flux`
- substrate depletion/accumulation
- total cell count
- final volume distribution

Use tolerance-based assessment. A typical acceptable result is:

- exact or near-exact initial cell count and positions
- same qualitative growth/death/metabolic phases
- small final phenotype and substrate differences
- no duplicated cells
- no missing cells
- no negative concentrations introduced by the migration beyond behavior already present in the source or target baseline

## Common Failure Modes

### Every Rank Creates Every Cell

Symptom: cell count is multiplied by MPI rank count.

Cause: serial `setup_tissue()` was reused and each rank called `create_cell(...)`.

Fix: distribute positions to owner ranks or load cells with `load_cells_from_pugixml(world, cart_topo)`.

### Duplicate or Corrupt Output Files

Symptom: malformed XML, overwritten reports, strange config copies.

Cause: every MPI rank ran `system("cp ...")`, wrote `simulation_report.txt`, or wrote legend/citations.

Fix: gate non-MPI-aware filesystem writes with `IOProcessor(world)`.

### Link Errors for MPI Symbols

Symptom: undefined references to `mpi_Environment`, `mpi_Cartesian`, or MPI collectives.

Cause: Makefile still uses `g++`, lacks `DistPhy_OBJECTS`, or lacks MPI-specific object files.

Fix: use `mpic++`; add `DistPhy_*`, `PhysiCell_cell_mpi.o`, and `PhysiCell_standard_models_mpi.o`.

### Missing dFBA Symbols

Symptom: undefined references to `dFBAIntracellular`, `dfba_Model`, CLP, libSBML, or bzip2.

Cause: dFBA object/library sections were not merged into the MPI Makefile.

Fix: include `PhysiCelldFBA_OBJECTS`, `INC`, `LIB`, dFBA compile rules, and static/dynamic bzip2 handling.

### Non-Reproducible Random Initial Conditions

Symptom: MPI result has the right counts but different placement/distribution.

Cause: each rank draws random numbers independently.

Fix: rank 0 should generate the full serial-order position list and distribute it with `mpi_CellPositions`.

### Rank Hangs During Output or Diffusion

Symptom: run stalls at save step or diffusion step.

Cause: only some ranks call an MPI collective function.

Fix: all ranks must call MPI-aware collective functions such as `save_PhysiCell...(..., world, cart_topo)`, `SVG_plot_mpi(...)`, `simulate_diffusion_decay(...)`, and `update_all_cells(...)`. Only gate ordinary non-collective I/O.

## Minimal LLM Task Prompt

Use this prompt when asking another LLM to perform a conversion:

```text
Convert the OpenMP PhysiCell sample project at <source_path> into a PhysiCell-X MPI + OpenMP project at <target_path>.

Use the existing PhysiCell-X MPI samples as the structural template. Preserve all project-specific biology, XML parameters, config assets, custom modules, intracellular hooks, cell rules, and output settings.

Required transformations:
- initialize `mpi_Environment` and `mpi_Cartesian` at the start of `main`
- use `load_PhysiCell_config_file(path, world)`
- keep `omp_set_num_threads(PhysiCell_settings.omp_num_threads)`
- call `setup_microenvironment(world, cart_topo)`
- create the cell container with `create_cell_container_for_microenvironment(..., world, cart_topo)`
- call MPI-aware `create_cell_types` and `setup_tissue`
- use MPI-aware MultiCellDS and SVG output
- use `simulate_diffusion_decay(dt, world, cart_topo)`
- use `update_all_cells(time, world, cart_topo)`
- gate ordinary output/logging/system calls with `IOProcessor(world)`
- finalize MPI with `world.Finalize()`
- update the Makefile to use `mpic++`, `DistPhy` objects, MPI PhysiCell objects, OpenMP flags, and dFBA dependencies if present
- ensure cells are created only on their owner rank
- validate with 1 rank and multiple ranks

Return a summary of files changed, build/test commands run, and any biological equivalence caveats.
```

## Final Review Checklist

Before declaring the conversion complete:

- `main.cpp` includes MPI initialization and `world.Finalize()`.
- `main.cpp` has no remaining serial-only calls for config, output, diffusion, or cell updates.
- `custom.h` prototypes match `custom.cpp` and `main.cpp`.
- `custom.cpp` uses MPI-aware initialization and cell loading.
- Cell placement does not duplicate cells across ranks.
- Non-collective file writes are gated by `IOProcessor(world)`.
- Makefile uses `mpic++`.
- Makefile includes `DistPhy_OBJECTS`.
- Makefile includes `PhysiCell_cell_mpi.o` and `PhysiCell_standard_models_mpi.o`.
- dFBA projects include `-DADDON_PHYSIDFBA`, dFBA objects, CLP/CoinUtils/libSBML libraries, and required include paths.
- Project builds from clean state.
- One-rank run reaches final time.
- Multi-rank run reaches final time.
- Biological comparison is based on matched cells/voxels, not raw file byte identity.
