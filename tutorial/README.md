# MolPainter tutorial
## Building phase separated 2:1 x:y aspect ratio lipid bilayers of DPPC:DIPC:Chol 4:4:2 in the MARTINI 2 force field.

In this tutorial we’re going to construct a phase separated lipid bilayer of DPPC, DIPC, and Chol in a 2:1 x:y aspect ratio system using MolPainter, solvate this system using MolSolvator, and prepare this system for production MD simulation in GROMACS.

In the first part of this tutorial, we’ll make a system of 32x16=512 molecules per leaflet. Because CHOL (102 molecules per leaflet) has a 4:1 preference for DPPC over DIPC (see Pantelopulos et al. doi:10.1016/j.bpj.2018.10.011), we’ll place 82 CHOL alongside the DPPC on one side of the system and 20 CHOL alongside the DIPC on the other side of the system.

## Table of Contents
- [MolPainter tutorial](#molpainter-tutorial)
  * [Building phase separated 2:1 x:y aspect ratio lipid bilayers of DPPC:DIPC:Chol 4:4:2 in the MARTINI 2 force field.](#building-phase-separated-2-1-x-y-aspect-ratio-lipid-bilayers-of-dppc-dipc-chol-4-4-2-in-the-martini-2-force-field)
  * [Table of Contents](#table-of-contents)
    + [MolPainter](#molpainter)
      - [1: Gather Files needed to build a system in MolPainter.](#1--gather-files-needed-to-build-a-system-in-molpainter)
      - [2: Set the size of the MolPainter canvases.](#2--set-the-size-of-the-molpainter-canvases)
      - [3: Load the molecules into MolPainter.](#3--load-the-molecules-into-molpainter)
      - [4: Use Blends to get the right amount of CHOL on each side of the system](#4--use-blends-to-get-the-right-amount-of-chol-on-each-side-of-the-system)
      - [5: Paint the first layer! This is going to be the upper leaflet of the bilayer.](#5--paint-the-first-layer--this-is-going-to-be-the-upper-leaflet-of-the-bilayer)
      - [6: Make a new layer and edit the z-axis position of the first layer.](#6--make-a-new-layer-and-edit-the-z-axis-position-of-the-first-layer)
      - [7: Load new instances of DPPC, DIPC, CHOL, flipping them, and also new instances of the Lo and Ld phase blends. Paint this new leaflet.](#7--load-new-instances-of-dppc--dipc--chol--flipping-them--and-also-new-instances-of-the-lo-and-ld-phase-blends-paint-this-new-leaflet)
      - [8: Save your work for later use.](#8--save-your-work-for-later-use)
      - [9: Export the lipid bilayer on a hexagonal lattice.](#9--export-the-lipid-bilayer-on-a-hexagonal-lattice)
      - [10: Appreciate your fine art in VMD.](#10--appreciate-your-fine-art-in-vmd)
    + [MolSolvator](#molsolvator)
      - [1: Solvating the volume above the lipid bilayer](#1--solvating-the-volume-above-the-lipid-bilayer)
      - [2: Solvating the volume below the lipid bilayer](#2--solvating-the-volume-below-the-lipid-bilayer)
    + [Molecular Dynamics](#molecular-dynamics)
      - [1: Writing the topology file.](#1--writing-the-topology-file)
      - [2: Minimization of the system](#2--minimization-of-the-system)
      - [2: Thermalization and pre-equilibration of the system via annealing in NPT](#2--thermalization-and-pre-equilibration-of-the-system-via-annealing-in-npt)
      - [3: Test production simulation](#3--test-production-simulation)

### MolPainter
#### 1: Gather Files needed to build a system in MolPainter.
The PDB files needed to build this system are located in the directory "./PDBfiles" for DPPC, DIPC, CHOL, W, WF, NA+, and CL-.

MolPainter only needs one PDB-formatted file of each of the molecules you wish to paint onto the canvas layers of MolPainter. If you are building a lipid bilayer, the conformation of the molecule should be very straight, aligned along the z-axis. We must specify a group of one or more atoms to define both the z-axial and relative xy-plane positions of the molecule to be used for painting.

These groups are selected by setting the temperature/b-factor value of atoms defining the group in the molecule PDB file. If multiple groups are selected, the center of geometry will be used to determine where the molecule will be placed into each cell. As an example, let’s look at the DPPC pdb file "DPPC.pdb" in Figure 1.A. The selected PO4 atom in the PDB structure will define the center of geometry demonstrated in Figure 1.B.

For the purposes of this tutorial, I’ve set the selections for DPPC and DIPC both to PO4. For CHOL I’ve set the selection to ROH. This will make the center of the lipid bilayer relatively smooth, making for a decent initial condition.

*Figure 1. (A) PDB file of DPPC selecting PO4  by setting the b-factor to 1 to define the center of geometry used by MolPainter to define the center of each cell. (B) Cartoon of MARTINI DPPC showing the location of the center of geometry used by MolPainter defined by selecting PO4.*

#### 2: Set the size of the MolPainter canvases.
The canvases in MolPainter are defined by a square grid (which can be transformed into a hexagonal lattice when you export your painting). The number of cells in this grid and their representative area can be defined using the *Grid Size* button. Let’s set the width and height of the lattice to 32 cells x 16 cells, set the cell spacing to 8Å, and set the major gridlines period (for guiding the eye) to 4 cells.

*Figure 2. Modifying the MolPainter canvas grid using the pop-out window from Grid Size.*

#### 3: Load the molecules into MolPainter.
MolPainter molecules are objects with several features, including the molecule PDB file path (for exporting the painting), the molecule color (for painting), and options for transforming the molecule, like flipping it 180 degrees, turning it upside-down.

Let’s load the DPPC, DIPC, and CHOL pdb files, making new *DPPC*, *DIPC*, and *CHOL* molecules represented by the colors blue, red, and green, respectively. An example for DPPC is shown in Figure 3.

*Figure 3. A DPPC molecule in MolPainter added via New Molecule.*

#### 4: Use Blends to get the right amount of CHOL on each side of the system
MolPainter blends make it easy to paint areas of your system that are meant to achieve a particular molecular composition. In a blend each molecule has a weight that determines how likely it is to get drawn when you paint the blend. For example, to define the Blend to use for the left side of the system, which we will compose with 205 DPPC and 82 CHOL, we can use a blend where we simply set the weights of DPPC and CHOL to 205 and 82 (Figure 4). Then let’s make a second blend with weights of DIPC and CHOL at 205 and 20 for painting the right side.

*Figure 4. A DPPC+CHOL blend in MolPainter added via New Blend.*

#### 5: Paint the first layer! This is going to be the upper leaflet of the bilayer.
Paint using:
* The pen tool (hold left click and drag on canvas).
* The box tool (control+left click, hold, drag, and release).
* The spray tool, whose nozzle width can be adjusted with the *Spray radius* bar (shift+left click, hold, and drag on canvas).
Rather than an eraser tool, we have an "Empty" molecule. If you wish to have a cell that has no molecule, paint it with the "Empty" molecule.

Let’s paint the left-hand side of the system with DPPC+CHOL, which will form the liquid ordered (Lo) phase. We can make a quick initial painting using the box tool, then touch up the number and location of DPPC and CHOL using the pen tool. Let’s aim for something like Figure 5.A. After that, let’s paint the right-hand side of the system with DIPC+CHOL, which will form the liquid disordered (Ld) phase, first using the rectangle tool then touching up as needed with the pen, coming to a painting like Figure 5.B.

Figure 5. (A) Painting using the DPPC molecule, CHOL molecule, and the DPPC+CHOL "Lo Phase Blend" blend. (B) Painting using the DIPC molecule, CHOL molecule, and the DIPC+CHOL "Ld Phase Blend" blend.

#### 6: Make a new layer and edit the z-axis position of the first layer.
Great, our first layer is painted. Now let’s make a second layer to serve as the lower leaflet of the lipid bilayer. To do this click the *New Layer* button and let’s name it "Lower Leaflet" and set the z-axis position to be -17 Å (Figure 6). Let’s also edit Layer 1 by right clicking on it and selecting "edit" and rename it to "Upper Leaflet" and set the z-axis position to 17 Å.

*Figure 6. Editing the name and z position of a layer.*

#### 7: Load new instances of DPPC, DIPC, CHOL, flipping them, and also new instances of the Lo and Ld phase blends. Paint this new leaflet.
Now we need to load conformations of DPPC, DIPC, and CHOL and have them flipped 180 degrees over the x-axis. We’ll color these light blue, red, and green. And then we’ll also make the similar blends as in previous steps to achieve the 205+82 DPPC+CHOL and 205+20 DIPC+CHOL ratios phases on the left and right halves of the system. To load and flip DPPC, for example, you just check the "Flip 180 degrees over x-axis" box (Figure 7.A). After you’ve finished loading everything, completing our palette of molecules for painting, paint the lower leaflet, and arrive at something like Figure 7.B.

*Figure 7. (A) Loading and flipping a new DPPC molecule via New Molecule. (B) Painted lower leaflet and completed molecule palette.*

#### 8: Save your work for later use.
MolPainter has a JSON file format for saving and loading paintings and the corresponding palette of molecules. You can save using the *Save* button, giving a path to the file location. You can then load your work for later editing via the *Open* button. This makes it very easy to make small tweaks to the painting to accomplish tasks such as preparing new initial conditions with some randomization of the initial positions. The system built here is provided as "tutorial.json".

#### 9: Export the lipid bilayer on a hexagonal lattice.
Let’s export the painting to a PDB structure of the lipid bilayer. Use the *Export* button, name the file, let’s say, "solute.pdb", and hit the radio button "Hexagonal lattice" to export this structure on a hexagonal lattice to the active directory. You could also export it to some other location via writing the path or browsing to it using the "Browse" button. This hexagonal lattice is defined by shifting every other column up on the y-axis by one half a cell length.

#### 10: Appreciate your fine art in VMD.
Let’s give your painting a look (Figure 8). A tcl script is included in this tutorial to set up graphics in VMD to quickly visualize at the system you’ve made "visualize.tcl", which can be executed from the command line with
```
vmd phase_separation.pdb -e visualize.tcl
```

*Figure 8. Our painting exported from MolPainter visualized in VMD with DPPC in blue, DIPC in red, and CHOL in green.*

Nice.

### MolSolvator
Now that we’ve got our lipid bilayer, let’s solvate it. We'll do this using MolSolvator, which was specifically designed for use with MolPainter in mind.

Given there are 1024 molecules in the lipid bilayer, we’ll add 10,024 water particles. 10% of these must be the "anti-freeze" MARTINI particle, which employs very strong anti-freeze to normal water interactions to prevent anomalous freezing at ambient temperature. Thus we’ll insert 9022 normal "W" water particles and 1002 anti-freeze "WF" water particles. One straightforward way we can typically determine the number of ions is to divide the molar concentration of water, 55.5 M, by the molar concentration of ion, 0.150 M, to determine the waters per ion, 55.5 M / 0.150 M = 370. However, the MARTINI water effectively represents 4 waters, while the MARTINI ions still just represent one ion each. As such, we would have 92.5 waters per ion, 108 Na+ and 108 Cl-.

#### 1: Solvating the volume above the lipid bilayer


#### 2: Solvating the volume below the lipid bilayer


### Molecular Dynamics
Now let's get this system running in GROMACS.

#### 1: Writing the topology file.
GROMACS has a wonderfully easy way to input the force field parameters and any additional restraints via it’s topology file system. In short, rather than define the parameters for every individual copy of a molecule, GROMACS only needs to be told the parameters of each unique type molecule, defined by a "[ moleculetype ]".

There is a main top format file which can load parameters defined in additional itp files. The MARTINI developers have provided the parameters for all canonically-defined MARTINI 2 molecules on their website. We will load these in the top file (let’s call it "system.top") at the beginning via the GROMACS "#include" statement.

```
#include "martini_v2.2.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "martini_v2.0_ions.itp"
```

Next we should give the system a name.

```
[ system ]
MolPainter Tutorial
```

And now all we have to do is tell GROMACS how many of each molecule type are present in the system (these happen to also the residue names in the pdb files, though they are defined in the MARTINI force field itp files [ moleculetype ] entries). These must appear in the same order that they appear in the PDB file.

MolPainter writes out each molecule in the order that they were loaded into MolPainter, thus we will have in order DPPC, DIPC, CHOL, DPPC, DIPC, CHOL. Following this will be the two layers of solvent we added, thus NA+, CL-, WF, W, NA+ CL-, WF, W. With each of these, we just give the number of molecules appearing at that place in the order.

```
[ molecules ]
; name  number
DPPC    205
DIPC    205
CHOL    102
DPPC    205
DIPC    205
CHOL    102
NA+      54
CL-      54
WF      501
W      5012
NA+      54
CL-      54
WF      501
W      5012
```

#### 2: Minimization of the system
The *editconf* tool of GROMACS can convert our PDB file of the system to a GRO file. We need to set the PBC dimensions too. We’re going to add 2 Å to xedge and 1 Å to yedge and zedge to make sure there is a tiny bit of empty space to ensure minimization and then thermalization of the system goes smoothly.

```
gmx editconf -f solute+water+ion.pdb -o system.gro -box 24.5362 12.2681 8.9
```

For minimization and simulation we’ll be using the 1.1 Å nonbonded cutoff scheme "martini_v2.x_new-rf.mdp" [http://www.cgmartini.nl/index.php/force-field-parameters/input-parameters] tested for use in GPU-accelerated simulation. The exact scheme we’ll use is defined in "min.mdp".  Let’s minimize the system.

```
gmx grompp -f min.mdp -c system.gro -p system.top -o min.tpr
gmx mdrun -v -deffnm min
```

#### 2: Thermalization and pre-equilibration of the system via annealing in NPT
Now let’s prepare the system for production simulation in one step. We’ll anneal the system from 200 K to 295 K over 1-ns and equilibrate at 295 K for an additional 1-ns both at 1 atm pressure. We’ll use the aggressive Berendsen thermostat and barostat scheme during this process and use a conservative 10 fs time step.

The membrane molecules and solvent molecules will be separated into two groups for thermostatting. We’ll need to make an index file defining the atom numbers belonging to these groups. Normally when you use the *gmx make_ndx* program it will pull up an interactive prompt. Let’s just input the directions by piping a the corresponding string of text from *echo* to *gmx make_ndx*.

```
echo -e "2 | 3 | 4 \n name 8 memb \n 5 | 6 | 7 \n name 9 solv \n q \n" | gmx make_ndx -f system.gro -o system.ndx
```

Now let’s run the annealing. You’ll likely get a warning about verlet-buffer-tolerance – let’s just ignore that, we know the temperature of the system will not exactly be correct during annealing. We’re also going to try to avoid LINCS warnings during this by setting the environment variable "GMX_MAXCONSTRWARN" to -1 , as the system will be a bit testy when changing size during annealing. If the system is not stable, try making the annealing state and simulation time longer ("annealing-time" and "nsteps", respectively).

```
gmx grompp -f anneal.mdp -c min.gro -p system.top -n system.ndx -o anneal.tpr -maxwarn 1

export GMX_MAXCONSTRWARN=-1
gmx mdrun -v -deffnm anneal
```

OK! We can see our system should look something like Figure 9. We’re ready for production simulation.

*Figure 9. The system after annealing.*

#### 3: Test production simulation
For production simulation we’ll use the instructions in "run.mdp". Typically, I like to save only the coordinates MARTINI trajectories, and only in the compressed xtc format at a frequency of one frame per ns, and to just save thermodynamic and log information at this same interval. Just to give this a try, let’s run 5 ns of production simulation. We’ll now use a time step of 25 fs, the Parrinello-Rahman barostat, and the Bussi velocity rescaling thermostat. Another change here is that we won’t generate new initial velocities, but will use those present in the last 3 columns of "anneal.gro". Let’s give it a try.

```
gmx grompp -f run.mdp -c anneal.gro -r anneal.gro -p system.top -n system.ndx -o run.tpr
gmx mdrun -v -deffnm run
```

If this has run,  you should be clear for production simulation!