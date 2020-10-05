# MolPainter

## Table of Contents
- [Introducing MolPainter and MolSolvator](#introducing-molpainter-and-molsolvator)
- [Installation](#installation)
- [MolPainter Usage](#molpainter-usage)
  * [Buttons](#buttons)
  * [Molecule palette](#molecule-palette)
  * [Layers](#layers)
  * [PDB files for molecules](#pdb-files-for-molecules)
- [MolSolvator Usage](#molsolvator-usage)
- [Give it a try!](#give-it-a-try)

## Introducing MolPainter and MolSolvator
#### Tools for building and solvating complex, planar molecular systems of arbitrary molecular composition and placement via painting.

![Screenshot](https://raw.githubusercontent.com/gpantel/MolPainter/master/tutorial/images/TutorialFigure.png)

MolPainter is a novel graphical tool that enables users to specifically define the location of molecules in multi-layered, planar molecular systems. MolPainter achieves this by treating each plane of a hypothetical molecular system, defined by a z-axial position, as a two dimensional grids which serve as canvases. By associating molecular structures (in PDB format) to colors, these canvases can be painted to precisely define molecular environments.

MolSolvator, the sister program of MolPainter, is a command line tool that can rapidly solvate such planar systems within the context of the lattices of the "solute" systems produced by MolPainter.

## Installation

Install through pip:
```
pip install MolPainter
```

Run from the command line:
```
molpainter
```

```
molsolvator -i input.toml
```

## MolPainter Usage
### Buttons
* Topbar buttons
    - **New** resets everything in MolPainter.
    - **Open Painting** loads a saved MolPainter state so you can continue changing your painting.
    - **Save Painting** saves a MolPainter state for later! This is a JSON-format file.
    - **Export System** constructes your painted system in PDB format. You can construct your system on square lattices or hexagonal lattices.
    - **Import Solute** *NEW IN 1.1!!* inserts a solute (PDBfile) to the MolPainter canvases, obstructing cells occupied by the solute coordinates and outputting the solute with the painting when exported.
    - **New Molecule** loads a PDB file corresponding to a molecule you wish to paint with. You can also set all painted copies of the molecule to be randomly rotated around the z-axis or be flipped 180˚ over the x-axis (like for lipid bilayers).
    - **New Blend** makes a "blend" of molecules. Blends draw molecules onto layer canvases with a weighted probability, running a RNG to determine which of the molecules that compose the blend are painted given their weighted probability.
    - **New Layer** makes new layers to paint on, defined by a z-axial position.
    - **Grid Size** sets the lattice width - the number of cells in each layer along the x-axis, the lattice height - the number of cells in each layer along the y-axis, the lattice spacing - the length of *each* cell, and the major gridline spacing, used for visual guidance when painting. 

* Toolbox buttons
    - **Magnifying glass -** zooms out from the lattice view.
    - **Magnifying glass +** zooms in on the lattice view.
    - **Pencil** enables the pencil drawing tool. Click and hold to draw selected molecule/blend.
    - **Square** enables the square selection tool. Click position to start drawing a rectilinear shape, hold, move, and release at position to fill selected molecule/blend into the rectilinear shape. *Can also be used via a ctrl+click shortcut*
    - **Spray can** sprays selected molecule/blend into an area. *Can also be used via a shift+click shortcut*
    - **Spray radius** adjust the radius of the spray can.

### Molecule palette
The molecule palette, located to the lower left of MolPainter, contains a list of all molecules and blends loaded via the *New Molecule* and *New Blend* buttons, in order. The *Empty* molecule is always present, and is essentially an "eraser".

When you want to paint with a molecule or blend in the palette, click on the radio button. You can *edit*, *copy*, or *delete* any molecule via a drop-down menu that appears when right-clicked.

### Layers
The canvas layers to which molecules and blends are painted, created using *New Layer*, can be selected for painting by clicking on the corresponding layer button, which can be named. Layer names and z-axial positions can also be edited via a drop-down menu that appears when the layer button is right-clicked.

### PDB files for molecules
The conformation of the molecule in PDB files loaded to MolPainter are only translated or translated + rotated. Typically you would want to orient an input PDB file along the z-axis. The position that defines where the molecule will be translated into the layer grid (e.g. atoms defining the hydrophilic plane of a lipid bilayer) is defined by the tempfactor/b-factors column. Typically temperature/b-factors are set to 0, but if set to 1, MolPainter will use the center of geometry of all atoms set to 1 to define where the molecule will be placed once the painting is exported.

For example, in a lipid bilayer, a user might select lipid ester oxygen atoms and cholesterol hydroxyl atoms to define a lipid leaflet, as these atom groups should be at approximatly the same z-axial position in a lipid bilayer.

## MolSolvator Usage
MolSolvator is a simple tool which adds solvent molecules, also defined by PDB files, to 3D lattices whose bounds are defined *in terms* of the grids used to generate systems in MolPainter.

MolSolvator takes an *input* solute system (intended to be systems produced by MolPainter) defined within the approximate lower bounds of min(*x*), min(*y*) = 0, 0. It is to this solute system that solvent is *added*, such that no inserted solvent overlaps with with the solute or other solvent molecules.

MolSolvator inserts solvent molecules into a grid which has a separate lattice spacing from the solute system. This *solvent lattice spacing* would typically be chosen as something slightly larger than the typical Lennard-Jones interaction minimum energy distance of the model. Each solvent molecule is inserted to the lattice centered on its center of geometry. The minimum and maximum x, y, and z-dimensions of the molecule will mark neighboring cells as "occupied" if touching the center of other neighboring cells. This typically would not prevent enough space between solvent molecules, so in addition to this a *buffer* length is added to the solvent dimensions to make sure there is enough space for the solvent molecules. Typically, the *bufffer* length should be shorter than the *solvent lattice spacing*.

The input file that controls the i/o for MolSolvator is a TOML-format file. Each entry is explained here.

* Input/output files
    - solute = *string*, the path to the input solute system.
    - output = *string*, the path for the output solvated system.
* Solute lattice
    - shape = *string*, input solute is a "hex" or "square" lattice.
    - width = *integer*, solute lattice width
    - height = *integer*, solute lattice height
    - spacing = *float*, solute lattice cell spacing, in angstroms
* Solvent lattice
    - spacing = *float*, solvent lattice cell spacing, in angstroms
    - buffer = *float*, solvent buffer length, in angstroms
    - lower_z_position = *float*, lower bonds of the solvent grid
    - upper_z_position = *float*, upper bounds of the solvent grid
* Solvent molecules
    - paths = [*string*, ..., *string*], list of paths to files defining structure of solvents
    - numbers = [*integer*, ..., *integer*], list of numbers of each solvent molecule to be inserted. Must be of the same length as *paths*.
    - max_iterations = *integer*, maximum number of attempts to insert a solvent molecule to the lattice before aborting.
    - rotate = *lower case boolean*, boolean to enable random 3D rotations of solvent molecules. ("true", not "True". "false", not "False")

In addition to the i/o file, MolSolvator has additional command line actions

* Command line actions
    - *-zeroz*, translates the system such that min(*z*) = solvent lattice spacing / 2, following solvation.
    - *-centerc*, translates the system such that <x,y,z> = 0,0,0, following solvation.

## What's new in 1.1?

**Solute** systems can now be added during a MolPainter session! These are systems in PDB format to which you want to *paint* using MolPainter while avoiding clashes.

Let's say for example you wish to use MolPainter to paint lipids around one or multiple proteins. In 1.0 you would have to do this by *guessing* where the protein will be located and paint the system to try to avoid clasing with the protein. Not very convenient...

Now you can *Insert a Solute* into MolPainter canvases. A solute will "obstruct" cells in the layers of MolPainter,  making them unavailable for painting, based on the real coordinates of the Solute system described in the input PDB file. Just how far these obstructed cells extend from the coordinates of the Solute can be tweaked by adjusting the *Buffer space*. This is the same kind of "buffer" used in MolSolvator to avoid clashes when solvating.

When the system is ultimately exported, the solute will be written first, followed by the molecules that have been painted using MolPainter.


The new *Insert Solutes* button has the following fields

* **Path** path to a PDB file of the solute you want to add. There can only be one solute in the system at a time!
* **Buffer space (Å)** the length of the buffer added to obstruct more cells neighboring coordinates of the solute.
* **Center solute at z (Å)** centers the (x,y,z) coordinates of the solute to the center of the MolPainter grid at a particular plane in the z-dimension.
* **Expand grid to fit solute** expands the MolPainter grid to fit the (x,y,z) coordinates of the solute system. If the solute system has negative coordinates, these coordinates will still lay outside of the MolPainter canvases!



## Give it a try!

A complete tutorial showing how to use MolPainter and MolSolvator to construct a system and perform a simulation using GROMACS is available in the GitHub repo at: https://github.com/gpantel/MolPainter/tree/master/tutorial
