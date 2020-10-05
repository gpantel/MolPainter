import json
import random
import numpy as np
import MDAnalysis as mda


class ZLayer():
    def __init__(self):
        self.lattice = None
        self.name = None
        self.zdepth = 0


class Molecule():
    def __init__(self):
        self.path = None
        self.name = "Empty"
        self.color = "#FFFFFF"
        self.index = 0
        self.flipbool = 0
        self.rotatebool = 0

    def get_paint_props(self):
        """
        Gives a tuple of properties needed to paint this molecule: ID and color.
        """
        return (self.index, self.color)


class Blender():
    def __init__(self):
        self.name = "Blender"
        self.index = -1
        self.molecules = {}

    def get_paint_props(self):
        """
        Gives a tuple of properties needed to paint a single molecule from this 
        blend: ID and color. The molecule to be painted is given by the distribution function.
        """
        mol = self.distribute_next()
        return (mol.index, mol.color)

    def distribute_next(self):
        """
        Distribution function that picks one molecule from this blend
        """
        return random.choices(population=list(self.molecules.keys()), weights=list(self.molecules.values()))[0]


class SolventLayer():
    def __init__(self):
        self.name = "Solvent layer"
        self.lattice_spacing = 4.0
        self.start_offset = 0
        self.end_offset = -10
        self.molecules = {}


class Project():
    def __init__(self):
        pass

    def init_defaults(self):
        """
        Set everything to default values for a new project
        """
        self.lattice_width = 25
        self.lattice_height = 25
        self.lattice_spacing = 8
        self.lattice_major_gridlines = 5
        self.layer_count = 1
        self.molecule_count = 0
        self.solvent_molecule_count = 0
        self.solvent_layer_count = 0
        self.blender_count = 0
        self.blender_offset = 1000
        self.import_solute = None
        self.solute_buffer_space = 3
        self.solute_z = None

        self.layers = []
        self.molecules = []
        self.blenders = []

        self.solvent_molecules = []
        self.solvent_layers = []

    def edit_lattice_params(self, width, height, spacing, lines):
        """
        Change the lattice size of this project. The layers must resize to match.
        """
        self.lattice_spacing = spacing
        self.lattice_width = width
        self.lattice_height = height
        self.lattice_major_gridlines = lines

        for layer in self.layers:
            new_lattice = np.zeros((self.lattice_height, self.lattice_width))
            for i in range(np.amin([layer.lattice.shape[0], new_lattice.shape[0]])):
                for j in  range(np.amin([layer.lattice.shape[1], new_lattice.shape[1]])):
                    new_lattice[i,j] = layer.lattice[i,j]
            layer.lattice = new_lattice

    def new_layer(self):
        """
        Come up with starting parameters for a potential new layer
        A name and Z depth may already be provided
        """
        layer = ZLayer()
        layer.lattice = np.zeros((self.lattice_height, self.lattice_width), dtype='int')
        layer.name = "Layer {}".format(self.layer_count)
        layer.zdepth = 0

        return layer

    def add_layer(self, layer):
        """
        Add this Z layer to the project
        """
        self.layer_count += 1
        self.layers.append(layer)

    def delete_layer(self, layer):
        """
        Delete this layer from the project
        """
        self.layers.remove(layer)
        layer = None

    def add_molecule(self, molecule):
        """
        Add this molecule to the project
        """
        self.molecules.append(molecule)
        self.molecule_count += 1

    def add_solvent_molecule(self, molecule):
        """
        Add this as a solvent molecule to the project
        """
        self.solvent_molecules.append(molecule)
        self.solvent_molecule_count += 1

    def delete_molecule(self, molecule):
        """
        Delete this molecule from the project
        """
        self.molecules.remove(molecule)
        molecule = None

    def delete_solvent_molecule(self, molecule):
        """
        Delete this solvent molecule from the project
        """
        self.solvent_molecules.remove(molecule)
        molecule = None

    def new_molecule(self):
        """
        Come up with starting parameters for a potential new molecule
        """
        default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        molecule = Molecule()
        molecule.index = self.molecule_count
        molecule.name = "Mol {}".format(self.molecule_count + self.solvent_molecule_count)
        molecule.color = default_colors[(self.molecule_count - 1) % 10]
        return molecule

    def new_solvent_molecule(self):
        """
        Come up with starting parameters for a potential new solvent molecule
        """
        default_colors = ['#8175aa', '#6fb899', '#31a1b3', '#ccb22b', '#a39fc9', '#94d0c0', '#959c9e', '#027b8e', '#9f8f12', '#767676']
        molecule = Molecule()
        molecule.index = self.solvent_molecule_count
        molecule.name = "Mol {}".format(self.molecule_count + self.solvent_molecule_count)
        molecule.color = default_colors[(self.solvent_molecule_count) % 10]
        return molecule

    def project_loaded(self):
        """
        Makes sure the counts match what was in the previous project
        """
        for mol in self.molecules:
            if mol.index >= self.molecule_count:
                self.molecule_count = mol.index + 1
        for mol in self.solvent_molecules:
            if mol.index >= self.solvent_molecule_count:
                self.solvent_molecule_count = mol.index + 1

    def new_blender(self):
        """
        Create an empty blender. 
        To make it 1-indexed, the default name is blender_count + 1, since
        there is no "default" or "empty" blender on a new project.
        """
        blender = Blender()
        blender.name = "Blend {}".format(self.blender_count + 1)
        blender.index = self.blender_count + self.blender_offset
        return blender

    def add_blender(self, blender):
        """
        Add this blender to the project
        """
        self.blenders.append(blender)
        self.blender_count += 1

    def delete_blender(self, blender):
        """
        Delete this blender from the project
        """
        self.blenders.remove(blender)
        blender = None

    def new_solvent_layer(self):
        """
        Come up with starting parameters for a new solvent layer
        """
        layer = SolventLayer()
        layer.name = "Solvent layer {}".format(self.solvent_layer_count + 1)
        return layer

    def add_solvent_layer(self, layer):
        """
        Add this solvent layer to the project
        """
        self.solvent_layers.append(layer)
        self.solvent_layer_count += 1

    def delete_solvent_layer(self, layer):
        """
        Remove this solvent layer from the project
        """
        self.solvent_layers.remove(layer)
        layer = None

    def edit_solute_settings(self, file, bufspace, center):
        """
        Apply solute configuration so it can be loaded
        """
        self.import_solute = file
        self.solute_buffer_space = float(bufspace)
        if center is not None:
            self.solute_z = float(center)
        else:
            self.solute_z = None

    def load_solute(self, should_expand):
        """
        Import a solute into the project. Expands the lattice if needed.
        """
        if self.import_solute is not None:
            self.solute = mda.Universe(self.import_solute)
        else:
            return

        if self.solute_z is not None:
            positions = self.solute.atoms.positions
            solute_indices = list(set(list(np.where(positions[:,2] > self.solute_z-self.lattice_spacing)[0])) & set(list(np.where(positions[:,2] < self.solute_z+self.lattice_spacing)[0])))
            solute_center = np.array([np.mean(positions[:,0]), np.mean(positions[:,1]), np.mean(positions[:,2])])
            lattice_center = np.array([(self.lattice_spacing*self.lattice_width)/2., (self.lattice_spacing*self.lattice_height)/2., self.solute_z])
            self.solute.atoms.positions = positions + (lattice_center - solute_center)
            pass

        if should_expand:
            xmax = np.amax(self.solute.atoms.positions[:,0])
            ymax = np.amax(self.solute.atoms.positions[:,1])
            newwidth = self.lattice_width
            newheight = self.lattice_height
            if xmax > (self.lattice_width*self.lattice_spacing):
                newwidth = int((xmax + self.solute_buffer_space + self.lattice_spacing) // self.lattice_spacing)
            if ymax > (self.lattice_height*self.lattice_spacing):
                newheight = int((ymax + self.solute_buffer_space + self.lattice_spacing) // self.lattice_spacing)
            self.edit_lattice_params(newwidth, newheight, self.lattice_spacing, self.lattice_major_gridlines)

        for layer in self.layers:
            self.overlay_solute(layer)

    def overlay_solute(self, layer):
        """
        Overlay the imported solute onto one layer so that the lattice regions
        occupied by the solute become obstructed.
        """
        positions = self.solute.atoms.positions
        solute_indices = list(set(list(np.where(positions[:,2] > layer.zdepth-self.lattice_spacing)[0])) & set(list(np.where(positions[:,2] < layer.zdepth+self.lattice_spacing)[0])))
        self.remove_overlay(layer)
        for ind in solute_indices:
            x = positions[ind][0]
            y = positions[ind][1]
            n_pos_x = int((x+self.solute_buffer_space) // self.lattice_spacing)
            n_neg_x = int((x-self.solute_buffer_space) // self.lattice_spacing)
            n_pos_y = int((y+self.solute_buffer_space) // self.lattice_spacing)
            n_neg_y = int((y-self.solute_buffer_space) // self.lattice_spacing)
            n_range_x = list(range(n_neg_x, n_pos_x+1))
            n_range_y = list(range(n_neg_y, n_pos_y+1))
            for row in n_range_y:
                for col in n_range_x:
                    if (row >= 0) and (row < layer.lattice.shape[0]) and (col >= 0) and (col < layer.lattice.shape[1]):
                        np.flip(layer.lattice.swapaxes(0,1), axis=1)[col][row] = -1

    def remove_overlay(self, layer):
        """
        Set obstructed lattice site IDs from -1 back to 0
        """
        for row in range(layer.lattice.shape[0]):
            for col in range(layer.lattice.shape[1]):
                if layer.lattice[row][col] == -1:
                    layer.lattice[row][col] = 0
    