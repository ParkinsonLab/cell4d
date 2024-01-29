# Cell4D
---

Cell4D is a Linux-and-MacOS C++\-based graphical spatial stochatic cell simulator capable of simulating a wide variety of cellular pathways. Molecules are simulated as particles within a user-defined simulation space under a Smoluchowski-based reaction-diffusion system on a static time-step basis. At each timestep, particles will diffuse under Brownian-like motion and any potential reactions between molecules will be resolved. Simulation space is divided into cubic sub-partitions called c-voxels, groups of these c-voxels can be used to define spatial compartments that can have optional rules that govern particle permeability, and reactions can be compartment-specific as well.

The program can be run in a graphical mode where the user can view the simulation state as it is being modelled, or a non-graphical mode if only the model data output is desired. Data output includes a count table for all molecules in each compartment in a .tsv format, JSON particle position logs that record the position and state of all particles at regular time-step intervals, as well as intermediate JSON checkpoint files that can be used to resume simulations without needing to start over from the beginning.

## Installation
---

This program can be run on MacOS or Linux systems only. To install for Ubuntu, navigate to the directory containing the Cell4D files and type:

```bash
./install_ubuntu.sh 
./build_and_run.sh
```

To install for MacOS, use:

```bash
./install_macos.sh
./build_and_run.sh
```

This script will download all required dependencies of Cell4D and compile the source files into an executable named "Cell4D".


## Usage
---


#### Graphical mode
![menu](./readme_files/menu.png)
## Menu Options
---
- Glycolysis: launches a glycolysis reaction simulation demo
- Cascade: launches a cascade reaction simulation demo
- CEACAM: launches a CEACAM reaction simulation demo
- AB Fusion: launches an AB Fusion reaction simulation demo

- Workbench - Simulation: Prompts the user to supply a SBML-compliant model file in XML
- Open Model Editor: Launches a FireFox web browser instance to the Online Model Editor at https://compsysbio.org/cell4d_ui/build

- About: Displays information on where to get additional help
- License: Displays information on where to locate the license for this tool

## Usage options 
---
To run Cell4D, navigate to the repository folder, and enter:
```bash
./Cell4D <model> <sim-parameters> <flags>
```
The parameters are described below:
```bash
where <model> is a cell4D formatted XML file, and <sim-parameters> include:
	--check-in or --checkpoint: checkpointed molecule file from previous cell4D runs, formatted as json.
	--whole-particle-output: name of particle count .tsv at each time step. default: total_container_stats.tsv
	--compartment-particle-output: name of particle count file for each compartment. default: compart_container_stats.tsv
	--whole-bulk-output: name of bulk count .tsv at each time step. default: total_bulk_stats.tsv
	--compartment-bulk-output: name of bulk count file for each compartment. default: compart_bulk_stats.tsv
	--check-out or --checkpoint-output: name of checkpoint file created. default: checkpoint
	--checkpoint-timesteps: timestep interval for regular checkpointing. default: 25000
	--timestep: number of time cycles to run this simulation, overrides MAX_CYCLES field in model.
	--timescale: timestep length (seconds) to use in this simulation, overrides TIMESCALE field in model.
	--seed: simulation seed to use for RNG.
<flags> are:
	--no-graphics: turns off graphical display of simulation.
	--multi: compartmental molecule output will be separate files by compartment rather than long concatenated file
<sim-parameters> can be called as "--parameter_name" value OR "--parameter_name=value"
	
Experimental features:
	--poslog-out: name of molecule position log output
	--poslog-comparts: name of compartment of logged molecule positions
	--poslog-interval: interval of molecule log recordings, default at 500
```
Leaving the model, sim-parameters, and flags unfilled will launch Cell4D into GUI-mode.
To explicitly invoke No-Graphics mode, use:
```bash
./Cell4D <model> <sim-parameters> <flags> --no-graphics
``` 

To access the help display, and show the full options:
```bash
./Cell4D -h
```

Using either option to open the selected model file will cause the simulation UI to appear.

## Workbench launch simulation
---
![load_model_window](./readme_files/load_model_window.png)

If the user opts to load their own simulation, the user will be prompted to click the "Launch Simulation" button.

## Simulation window
---
![sim_window](./readme_files/sim_window.png)

When the simulation is loaded, the user will be see this simulation window
On the sidebar on the right are the runtime options of Cell4D graphical mode. 

![play](./readme_files/play.png)

Simulations can be freely started/paused/stopped by using the Control Panel buttons.  Clicking on the stop button will ask the user if they wish to quit the program entirely.  

![speed](./readme_files/speed.png)

The simulation speed input counter controls the simulation speed. 

![proteins](./readme_files/proteins.png)

Particles (proteins) and bulk molecules (metabolites) can be selectively hidden by the Display Proteins and Display Concentrations checkboxes respectively. The left radio button selects all.  The right button deselects all.  Each particle or bulk molecule can be selected individually.

![metabolites](./readme_files/metabolites.png)

Particles(proteins) and bulk molecules(metabolites) are displayed separately.

![compart](./readme_files/compart.png)

Compartments can be highlighted by Display Environment.  The left button selects all compartments to be highlighted. The right button deselects all compartments.  The compartments can be selected individually from the list.

![scale](./readme_files/scale.png)

The camera angle on the graphical simulation space can be freely rotated and zoomed, but not translated. Camera will always be centered on the simulation space. The camera angle is rotated by dragging the simulation.  The mouse scroll wheel is used to zoom in/out of the simulation window. Unchecking each box disables their function.  The last checkbox (Display Statistics will be explained below)

![close](./readme_files/close.png)

The close button exits the simulation without saving. 

![stats](./readme_files/stats.png)

Checking the "Display Statistics" box produces a time-series chart of the relative population of each particle and/or bulk molecule.
To choose which line to display on the graph, use the Display Proteins, and Display Metabolite checkboxes.
 


#### Non-graphical mode
To run in non-graphical mode, navigate to the Cell4D folder and type

```bash
./Cell4D <file-name> --no-graphics <additional-options>
```

where <file-name> is the path to the Cell4D XML model file. A comprehensive list of run options can be found with:

```bash
./Cell4D --help
```

In this mode, no graphics will be loaded and the simulation will simply start with events such as reactions being logged in the console output. Once the simulation is completed, the model data output will be automatically generated.

## Input files
---

Model files are written in XML, using a modified form of level 2 SBML. There are 6 main sections of a Cell4D model file: Environmental Variables, List of Compartments, List of Annotation Species, List of Species, List of Reactions, List of Events.

#### Environmental Variables

Environmental variables specify general simulation parameters such as the XYZ-dimensions of the simulation space (c-voxels in each cartesian direction), time and space-scale of the simulation, shape of the space, and total number of timesteps to simulate.

#### List of Compartments

The list of compartments describes all spatial compartments within the model as well as their locations. This is done by using 6 coordinates: x1 y1 z1, x2 y2 z2, to define rectangular prisms of arbitrary size within the bounds of the simulation space. Any number of these coordinate groups can be used in the definition of a compartment, meaning any arbirary shape for a compartment can be defined using these parameters.

#### List of Annotation Species

The list of annotation (or base) species describes the attributes of each unit that serves as a building block for a molecular species. For example, in a simple dimer reaction system where molecules ```A + B -> AB```, molecules ```A``` and ```B``` are comprised of a single unit of species ```A``` and ```B``` respectively. ```AB``` is a molecule that is formed by one unit of ```A``` and one unit of ```B```. This section describes the properties of these units, including their species type (discrete particle vs bulk molecule), diffusion constant, displayed color, a list of valid spatial compartments, and optional modification sites and binding sites.

There are two ways to model molecules in Cell4D: they can be simulated as discrete dimensionless point particles, or as local concentrations of molecules within a c-voxel. The particle approach should be used for most cases, as it is the most fine-grain method to simulate molecules and their interactions with one another. However, for molecules that are particularly high in concentration, such as ions, it would be more efficient to simulate them as local concentrations using the "bulk molecule" approach, which reduces the number of computations required to simulate the diffusion and reactions of that molecule type.

Species binding sites are additional optional parameters for the species definition. each site consists of a binding site name and optional binding state parameters

> Binding site name: string  

> Site states: Array\<string> | null

Sites with no site states will be considered "binding sites" with two possible states: ```bound``` and ```unbound```. These are flags indicating if a site is in use for molecule complexes with multiple species bound together.  
Site states can also be filled with multiple ```states```, which indicate the modification state of the species. For example, a species that can be phosphorylated can have the ```P``` and ```U``` site states, representing the phosphorylated and unphosphorylated states respectively. Binding sites can have as many sites as the user would like, but it is recommended that modification sites have at least 2 states.  

#### List of Species

Molecular species are comprised of one or more building block units of base species, described above. Using the previous example, species ```AB``` is a molecule that is formed by one unit of ```A``` and one unit of ```B```. This section lists a number of combinations of base species as well as their modification states and defines the molecule as a "species." By default, molecule properties are automatically inherited from the species that make up the molecule. 

#### List of Reactions

Reactions are patterns (reaction rules) that describe how molecules interact with one another within the simulation. 

> Reaction name: string  

> Reaction compartment: Array\<string> | null  

The reaction compartment fields allow the user to optionally select certain compartments where the reaction can occur in. If left blank, the reaction can occur anywhere within the simulation space.

Reactions are divided into interaction reactions and enzymatic reactions. Enzymatic reactions involve the conversion of a bulk molecule reactant into a bulk molecule product by a "modifier" enzyme, which is not changed by the reaction. The reaction kinetics are defined by the four reversible Michaelis-Menten parameters in the following equation:


![reversible_MM](./readme_files/reversible_MM_eq.png)

> Vfmax: number  
> Vrmax: number  
> Ksm: number  
> Kpm: number  

Interactions are the other type of reactions, and include unimolecular reactions and bimolecular reactions, which have one and two reactants respectively, and can produce one or more products. Interaction reactions can be set to reversible, which will require an additional reaction rate input.

> Reaction rate constant  
> Reverse reaction rate constant  

The units of the reaction rate constant depends on the number of reactants. For unimolecular reactions, the units of the rate constant are s<sup>-1</sup>. For bimolecular reactions it will be s<sup>-1</sup>M<sup>-1</sup>.

> List of reactants: molecule: base species  
> List of products: molecule: base species  

For interaction reactions, reactants are limited to a maximum of two molecules, but there is no limit on the number of product molecules. Each of these molecules must have at least 1 base species, but the binding and modification sites of the species are optional. When determining if a molecule is a valid reactant, the parameters listed in the list of reactants are compared to the molecule, so unspecified reactant sites will not be checked and are valid by default, regardless of its state. For product molecules, unfilled site parameters will be unchanged after the reaction.  

For enzymatic reactions, there may only be one reactant and one product, and they must both be bulk molecules. The modifier must be a particle.

#### List of Events

Events are custom behaviors that can be added to the simulation model that are unrelated to events. They can be set to happen from various triggers, and also set to execute probabilistically after the trigger occurs.

> Event name: string  

> Event probability: number  

> Event type: "Add molecules" | "Remove molecules" | "Molecule transport"  

Currently there are three types of possible events. A molecule, configured in the [List of Molecules](#List-of-Molecules), can be added to or deleted from a certain location, or be transported from one location to a different destination.

> Event location: string | X,Y,Z | X1,Y1,Z1, X2,Y2,Z2

Locations in an event definition can be formulated as a compartment string, a 3D coordinate representing a voxel, or two coordinates that form opposite corners of a rectangular prism, similar to a compartment lattice definition.

> Event trigger: "Time triggered" | "State triggered" | "Event triggered"  

There are three possible triggers for an event

> Time triggered: 
>> Repeat interval: integer | null  
>> Initial timestep: integer  
>> Final timestep: integer | null  

Time triggered events require an initial timestep that indicates the first time the event will occur, and can optionally be repeated after an interval of timesteps. The event can trigger repeatedly until the end of the simulation, or a final timestep can be specified which terminates further event activation after the selecte timestep.

> State triggered: 
>> Repeat interval: integer | null  
>> Initial timestep: integer  
>> Final timestep: integer | null  
>> Trigger location: string | X,Y,Z | X1,Y1,Z1, X2,Y2,Z2  
>> Trigger condition:  
>>> Molecule: string  
>>> Type: "more than" | "less than"  
>>> amount: integer  

State triggered reactions will have a trigger condition where the event will activate if a selected molecule, defined in the [List of Molecules](#List-of-Molecules), is greater or lesser than a specified amount within the trigger location when the state trigger condition is checked. The interval at which this condition is checked is specified by the repeat interval, so a repeat interval of 1 will cause the simulation to check for this trigger condition at every timestep.

> Event triggered: 
>> Triggering event: string  
>> Triggering delay: integer  

Finally, events can be triggered by other events. The event-triggered event can occur immediately after the triggering event occurs, or can be set to fire after a timestep delay.

## Authors
---

Donny Chan, Graham Cromar, Billy Taj & John Parkinson
