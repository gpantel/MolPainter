# MolPainter

## Table of Contents
- [Introducing MolPainter and MolSolvator](#introducing-molpainter-and-molsolvator)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Test Dependencies](#test-dependencies)
- [Documentation](#documentation)
- [Tutorial](#tutorial)
- [What's new in 1.1?](#whats-new-in-11)

## Introducing MolPainter and MolSolvator
#### Tools for building and solvating complex, planar molecular systems of arbitrary molecular composition and placement via painting.

![Screenshot](https://raw.githubusercontent.com/gpantel/MolPainter/master/docs/tutorial/images/TutorialFigure.png)

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

### Dependencies

MolPainter's GUI makes use of `tkinter` in order to provide a native interface on Linux, Mac, and Windows.

The additional package dependencies are `numpy`, `scipy`, `toml`, and `MDAnalysis`.

### Test Dependencies

Tests can be run with `pytest` from the root of the source tree:
```
python3 -m pytest
```

### Documentation

Descriptions of the objects and functions of MolPainter and MolSolvator are available [here](https://github.com/gpantel/MolPainter/blob/master/docs/README.md)

### Tutorial

A tutorial demonstrating the major functions of MolPainter and MolSolvator on a complex mixture lipid bilayer is available [here](https://github.com/gpantel/MolPainter/tree/master/docs/tutorial)

### What's new in 1.1?

**Solute** systems can now be added during a MolPainter session! These are systems in PDB format to which you want to *paint* using MolPainter while avoiding clashes.

Let's say for example you wish to use MolPainter to paint lipids around one or multiple proteins. In 1.0 you would have to do this by *guessing* where the protein will be located and paint the system to try to avoid clasing with the protein. Not very convenient...

Now you can *Insert a Solute* into MolPainter canvases. A solute will "obstruct" cells in the layers of MolPainter,  making them unavailable for painting, based on the real coordinates of the Solute system described in the input PDB file. Just how far these obstructed cells extend from the coordinates of the Solute can be tweaked by adjusting the *Buffer space*. This is the same kind of "buffer" used in MolSolvator to avoid clashes when solvating.

When the system is ultimately exported, the solute will be written first, followed by the molecules that have been painted using MolPainter.

The new *Insert Solutes* button has the following fields

* **Path** path to a PDB file of the solute you want to add. There can only be one solute in the system at a time!
* **Buffer space (Å)** the length of the buffer added to obstruct more cells neighboring coordinates of the solute.
* **Center solute at z (Å)** centers the (x,y,z) coordinates of the solute to the center of the MolPainter grid at a particular plane in the z-dimension.
* **Expand grid to fit solute** expands the MolPainter grid to fit the (x,y,z) coordinates of the solute system. If the solute system has negative coordinates, these coordinates will still lay outside of the MolPainter canvases!

