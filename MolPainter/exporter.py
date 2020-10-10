import os
import _thread
import tkinter as tk
from tkinter import filedialog
import numpy as np
import MDAnalysis as mda


class Exporter(tk.Frame):
    """
    This defines a window for creating a PDB structure out of the project contents.
    It contains any post-processing options, such as creating a hexagonal lattice.
    """
    def __init__(self, project, file, projlocation, master=None):
        super().__init__(master)
        self.master = master
        self.project = project
        self.grid(column=0, row=0)
        self.exportfilevar = file
        self.projfilevar = projlocation

        self.lblf = tk.Label(self, text="Path")
        self.lblf.grid(row=0, column=0, padx=2, pady=4)
        self.file_entry = tk.Entry(self, textvariable=self.exportfilevar)
        self.file_entry.grid(row=0, column=1, columnspan=3, padx=2, pady=4, sticky=("E", "W"))
        self.file_button = tk.Button(self, text="Browse", command=self.file_select_action)
        self.file_button.grid(row=0, column=4, padx=2, pady=4)

        self.lattice_choice = tk.StringVar()
        self.lattice_choice.set("Square") 

        self.square_button    = tk.Radiobutton(self, 
                                text="Square lattice",
                                padx = 2, 
                                variable=self.lattice_choice, 
                                value="Square")
        self.square_button.grid(row=1, column=0, padx=40, pady=4)
        self.hexagonal_button = tk.Radiobutton(self, 
                                text="Hexagonal lattice",
                                padx = 2, 
                                variable=self.lattice_choice, 
                                value="Hexagon")
        self.hexagonal_button.grid(row=1, column=1, padx=5, pady=4)

        self.cancel_button = tk.Button(self, text="Close", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="Export", command=self.ok_action)
        self.cancel_button.grid(row=2, column=0, padx=2, pady=4)
        self.ok_button.grid(row=2, column=1, padx=2, pady=4)

        self.statuslabel = tk.Label(self, text="Export status: not started.")
        self.statuslabel.grid(row=3, column=0, columnspan=5, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)
        self.lattice_spacing = self.project.lattice_spacing

    def file_select_action(self):
        self.exportfilevar.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("PDB file","*.pdb"),("all files","*.*"))))

    def cancel_action(self):
        self.master.destroy()

    def ok_action(self):
        self.ok_button["state"] = "disabled"
        export_thread = _thread.start_new_thread(self.start_export, ())

    def return_event(self, event):
        self.ok_button.invoke()

    def start_export(self):
        self.statuslabel["text"] = "Export status: started..."
        try:
            self.export_system()
        except Exception as e:
            self.statuslabel["text"] += " failed\n{}".format(e)

        self.ok_button["state"] = "normal"
        self.statuslabel["text"] += "\nStopped."

    def export_system(self):    
        layers = []
        if self.project.import_solute != None:
            layers.append(self.project.solute.atoms)
        for i in range(len(self.project.layers)):
            layers.append(self.constructionate_layer(self.project.layers[i]))
        layer = mda.Merge(*layers)
        # assign x and y dimensions to the box. The z-dimension will equal the smaller of the two.
        self.lattice_spacing = self.project.lattice_spacing
        x = (self.project.lattice_width*self.lattice_spacing)
        y = (self.project.lattice_height*self.lattice_spacing)
        # add 1 more lattice_spacing to each dimensions -- if one dim. is longer, add more
        if self.project.lattice_width > self.project.lattice_height:
            x += (float(self.project.lattice_width)/float(self.project.lattice_height))*(self.lattice_spacing/2.)
            y += self.lattice_spacing/2.
        elif self.project.lattice_height < self.project.lattice_width:
            x += self.lattice_spacing/2.
            y += (float(self.project.lattice_width)/float(self.project.lattice_height))*(self.lattice_spacing/2.)
        else:
            x += self.lattice_spacing/2.
            y += self.lattice_spacing/2.

        layer.dimensions = np.array([x, y, np.amin([x,y]), 90., 90., 90.])
        # numer all residue IDs from 1 up to the number of molecules
        layer.residues.resids = np.arange(1,layer.residues.n_residues+1)
        destfile = self.exportfilevar.get()
        layer.atoms.write(destfile)
        self.statuslabel["text"] += "\nOutput saved as {}".format(destfile)
        return

    def constructionate_layer(self, layer):
        if self.lattice_choice.get() == "Square":
            lattice_xy_coordinates = self.map_lattice_to_xy_squares()
        elif self.lattice_choice.get() == "Hexagon":
            lattice_xy_coordinates = self.map_lattice_to_xy_hexagons()

        all_pdb_molecule_groups = []
        for i in range(len(self.project.molecules)):
            molecule_id = self.project.molecules[i].index
            # swap the layer.lattice x,y coordinates, then flip the y-coordinates s.t. x > and y ^
            pdb_lattice_sites = np.where(np.flip(layer.lattice.swapaxes(0,1), axis=1) == molecule_id)
            # if there are 0 instances of this molecule, just continue
            if len(pdb_lattice_sites[0]) == 0:
                continue
            path = self.project.molecules[i].path
            if path == None: 
                continue
            if os.path.dirname(path) == '':
                path = os.path.join(os.path.dirname(self.projfilevar.get()), path)
            if not os.path.exists(path):
                continue
            pdbu = mda.Universe(path)
            if len(np.where(pdbu.atoms.tempfactors == 1)[0]) == 0:
                head_indices = pdbu.atoms.indices
            else:
                head_indices = np.where(pdbu.atoms.tempfactors == 1)[0]
            # create an empty Universe to store system parameters
            n_mols = len(pdb_lattice_sites[0]) # number of molecules
            mol_n_residues = pdbu.atoms.n_residues
            n_residues = n_mols*mol_n_residues
            mol_n_atoms = pdbu.atoms.n_atoms
            n_atoms = mol_n_atoms*n_mols
            # determine residue index of each atom
            resindices = []
            for j in range(n_mols):
                for k in range(mol_n_atoms):
                    resindices.append(pdbu.atoms.resindices[k]+(j*mol_n_residues))
            segindices = [0]*n_residues

            # initialize an empty universe and populate it
            molecule_universe = mda.Universe.empty(n_atoms, n_residues=n_residues, atom_resindex=resindices, 
                                                            residue_segindex=segindices, trajectory=True)
            # assign topology attributes
            molecule_universe.add_TopologyAttr('name', list(pdbu.atoms.names)*n_mols)
            molecule_universe.add_TopologyAttr('type', list(pdbu.atoms.types)*n_mols)
            molecule_universe.add_TopologyAttr('resnames', list(pdbu.atoms.residues.resnames)*n_mols)
            molecule_universe.add_TopologyAttr('resid', list(range(1, (n_residues)+1)))
            molecule_universe.add_TopologyAttr('segid', [list(pdbu.atoms.segids)[0]])
            molecule_coordinates = []

            if self.project.molecules[i].flipbool == 1:
                pdbu.atoms.positions = self.flip_molecule(pdbu.atoms.positions)
            pdbuxyz = np.copy(pdbu.atoms.positions)
            if 1 in pdbu.atoms.tempfactors:
                head_indices = np.where(pdbu.atoms.tempfactors == 1)[0]
                head_positions = pdbu.atoms.positions[head_indices]
                head_cog = np.mean(pdbu.atoms.positions[head_indices], axis=0)
            else:
                head_cog = np.mean(pdbu.atoms.positions, axis=0)
            pdbuxyz = pdbuxyz - head_cog # center the coordinates around the head
            pdbuxyz[:,2] += layer.zdepth # move the coordinates to zdepth such that the head is at zdepth


            for j in range(len(pdb_lattice_sites[0])):
                # in units angstroms
                y, x = lattice_xy_coordinates[pdb_lattice_sites[1][j]][pdb_lattice_sites[0][j]]
                pdbxyz = np.copy(pdbuxyz)
                if self.project.molecules[i].rotatebool == 1:
                    pdbxyz = self.rotate_molecule(pdbxyz)
                pdbxyz[:,0] += x
                pdbxyz[:,1] += y
                molecule_coordinates.append(pdbxyz)
            molecule_coordinates = np.concatenate(molecule_coordinates)
            molecule_universe.atoms.positions = molecule_coordinates
            all_pdb_molecule_groups.append(molecule_universe.atoms)
        return mda.Merge(*all_pdb_molecule_groups).atoms # this is an atomgroup for the layer

    def flip_molecule(self, xyz):
        R = np.array([[1,  0,  0],
                      [0, -1,  0],
                      [0,  0,  -1]])
        return np.dot(xyz,R)

    def rotate_molecule(self, xyz):
        theta = np.random.random()*2*np.pi
        R = np.array([[np.cos(theta), -np.sin(theta),  0],
                      [np.sin(theta), np.cos(theta),  0],
                      [0,  0,  1]])
        return np.dot(xyz,R)

    def map_lattice_to_xy_squares(self):
        self.lattice_spacing = self.project.lattice_spacing
        lattice_xy_coordinates = []
        for i in range(self.project.lattice_height):
            lattice_xy_coordinates_row = []
            for j in range(self.project.lattice_width):
                x = (i*self.lattice_spacing)+(self.lattice_spacing/2)
                y = (j*self.lattice_spacing)+(self.lattice_spacing/2)
                lattice_xy_coordinates_column = (x, y)
                lattice_xy_coordinates_row.append(lattice_xy_coordinates_column)
            lattice_xy_coordinates.append(lattice_xy_coordinates_row)
        return lattice_xy_coordinates

    def map_lattice_to_xy_hexagons(self):
        lattice_xy_coordinates = []
        for i in range(self.project.lattice_height): # height or y-dimension
            lattice_xy_coordinates_row = []
            for j in range(self.project.lattice_width): # width or x-dimension
                if j%2 == 1: x = (self.lattice_spacing*i) + self.lattice_spacing
                elif j%2 == 0: x = (self.lattice_spacing*i) + (self.lattice_spacing/2)
                y = (self.lattice_spacing*j) + (self.lattice_spacing/2)
                lattice_xy_coordinates_column = (x, y)
                lattice_xy_coordinates_row.append(lattice_xy_coordinates_column)
            lattice_xy_coordinates.append(lattice_xy_coordinates_row)
        return lattice_xy_coordinates
