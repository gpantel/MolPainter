# Table of Contents
{:toc}

## Introducing MolPainter and MolSolvator
<h5>Tools for building and solvating complex, planar molecular systems of arbitrary molecular composition and placement via painting.</h5>

MolPainter is a novel graphical tool that enables users to specifically define the location of molecules in multi-layered, planar molecular systems. MolPainter achieves this by treating each plane of a hypothetical molecular system, defined by a z-axial position, as a two dimensional grids which serve as canvases. By associating molecular structures (in PDB format) to colors, these canvases can be painted to precisely define molecular environments.

MolSolvator, the sister program of MolPainter, is a command line tool that can rapidly solvate such planar systems within the context of the lattices of the "solute" systems produced by MolPainter.



## Usage
### Buttons
* Topbar buttons
    - **New** resets everything in MolPainter.
    - **Open** loads a saved MolPainter state so you can continue changing your painting.
    - **Save** saves a MolPainter state for later! This is a JSON-format file.
    - **Export** constructes your painted system in PDB format. You can construct your system on square lattices or hexagonal lattices.
    - **New Molecule** loads a PDB file corresponding to a molecule you wish to paint with. You can also set all painted copies of the molecule to be randomly rotated around the z-axis or be flipped 180˚ over the x-axis (like for lipid bilayers).
    - **New Blend** makes a "blend" of molecules. Blends draw molecules onto layer canvases with a weighted probability, running a RNG to determine which of the molecules that compose the blend are painted given their weighted probability.
    - **New Layer** makes new layers to paint on, defined by a z-axial position.
    - **Grid Size** sets the lattice width - the number of cells in each layer along the x-axis, the lattice height - the number of cells in each layer along the y-axis, the lattice spacing - the length of *each* cell, and the major gridline spacing, used for visual guidance when painting. 

* Toolbox buttons
    - **Magnifying glass -** zooms out from the lattice view.
    - **Magnifying glass -** zooms in on the lattice view.
    - **Pencil** enables the pencil drawing tool. Click and hold to draw selected molecule/blend..
    - **Square** enables the square selection tool. Click position to start drawing a rectilinear shape, hold, move, and release at position to fill selected molecule/blend into the rectilinear shape. *Can also be used via a ctrl+click shortcut*
    - **Spray can** sprays selected molecule/blend into an area. *Can also be used via a ctrl+click shortcut*
    - **Spray radius** adjust the radius of the spray can.

### Molecule palette
The molecule palette, located to the lower left of MolPainter, contains a list of all molecules and blends loaded via the *New Molecule* and *New Blend* buttons, in order. The *Empty* molecule is always present, and is essentially an "eraser".

When you want to paint with a molecule or blend in the palette, click on the radio button. You can *edit*, *copy*, or *delete* any molecule via a drop-down menu that appears when right-clicked.

### Layers
The canvas layers to which molecules and blends are painted, created using *New Layer*, can be selected for painting by clicking on the corresponding layer button, which can be named. Layer names and z-axial positions can also be edited

### PDB files for molecules
The conformation of the molecule in PDB files loaded to MolPainter are only translated or translated + rotated. Typically you would want to orient an input PDB file along the z-axis. The position that defines where the molecule will be translated into the layer grid (e.g. atoms defining the hydrophilic plane of a lipid bilayer) is defined by the tempfactor/b-factors column. Typically temperature/b-factors are set to 0, but if set to 1, MolPainter will use the center of geometry of all atoms set to 1 to define where the molecule will be placed once the painting is exported.

For example, in a lipid bilayer, a user might select lipid ester oxygen atoms and cholesterol hydroxyl atoms to define a lipid leaflet, as these atom groups should be at approximatly the same z-axial position in a lipid bilayer.

## Install

Install locally via pip:
```
python3 setup.py sdist bdist_wheel
pip install dist/MolPainter-0.9.0-py3-none-any.whl
```

Run from the command line:
```
molpainter
```

```
molsolvator -i input.tmol
```

## Usage




### See it in action

A complete demo with molecule source files will be added in the `example/` directory.


