import os
import shutil
import json

import tkinter as tk
from tkinter import filedialog


class ProjSaver(tk.Frame):
    """
    This defines a window for saving the project information to a text file.
    The goal is to be able to copy the workspace from one machine to another.
    For this reason, molecule source files can be copied, and relative paths can be used.
    """
    def __init__(self, project, savfile, copybool, master=None):
        super().__init__(master)
        self.master = master
        self.project = project
        self.grid(column=0, row=0)
        self.savefilevar = savfile
        self.copyfilesbool = copybool

        self.lblf = tk.Label(self, text="Path")
        self.lblf.grid(row=0, column=0, padx=2, pady=4)
        self.file_entry = tk.Entry(self, textvariable=self.savefilevar)
        self.file_entry.grid(row=0, column=1, columnspan=3, padx=2, pady=4, sticky=("E", "W"))
        self.file_button = tk.Button(self, text="Browse", command=self.file_select_action)
        self.file_button.grid(row=0, column=4, padx=2, pady=4)

        self.copy_checkbox = tk.Checkbutton(self, text="Copy molecule source files", variable=self.copyfilesbool)
        self.copy_checkbox.grid(row=1, column=1, columnspan=2, padx=2, pady=4, sticky=("W"))

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="Save", command=self.ok_action)
        self.cancel_button.grid(row=2, column=0, padx=2, pady=4)
        self.ok_button.grid(row=2, column=4, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)

    def file_select_action(self):
        self.savefilevar.set(filedialog.asksaveasfilename(initialdir = "./",title = "Select file",filetypes = (("JSON file","*.json"),("all files","*.*"))))

    def cancel_action(self):
        self.master.destroy()

    def ok_action(self):
        destdir = self.get_destdir()
        if not destdir:
            return
        if self.copyfilesbool.get():
            if not self.copy_molecules(destdir):
                return
        if self.copyfilesbool.get():
            if self.project.import_solute is not None:
                self.copy_solute(destdir)
        self.write_project_file()
        self.master.destroy()

    def return_event(self, event):
        self.ok_button.invoke()

    def copy_molecules(self, destdir):
        """
        Make copies of molecule source files in the same location as the saved
        project. This is for sharing the project data among multiple users or
        workstations. The project workspace will then be able to load molecule
        data from a known location, rather than ask the user to browse to the
        correct files for each one.
        """
        for mol in self.project.molecules:
            if mol.index == 0:
                continue
            success = self.validate_mol(mol)
            if not success:
                return False

        for mol in self.project.molecules:
            if mol.index == 0:
                continue
            self.copy_mol(mol, destdir)

        return True

    def copy_solute(self, destdir):
        """
        Make a copy of the solute file along with the other project files.
        Check if the file exists and prompt to update it if it does not.
        Modifies the path to only refer to the filename.
        """
        exists = True
        try:
            exists = os.path.exists(self.project.import_solute)
        except:
            exists = False
        while not exists:
            newpath = None
            result = tk.messagebox.askyesno(
                message="The source path for solute doesn't exist.\n{}\nDo you want to choose a new path?".format(self.project.import_solute))
            if result:
                newpath = filedialog.askopenfilename(filetypes = (("PDB file","*.pdb"),("all files","*.*")))
            else:
                self.project.import_solute = None
                return
            try:
                exists = os.path.exists(self.project.import_solute)
            except:
                exists = False
        try:
            shutil.copy2(self.project.import_solute, destdir)
            self.project.import_solute = os.path.basename(self.project.import_solute)
        except OSError as e:
            tk.messagebox.showinfo(message=e.strerror)

    def get_destdir(self):
        """
        Get the destination directory based on the savefile location. If the
        location is empty or not valid, returns None.
        """
        try:
            dir = os.path.dirname(self.savefilevar.get())
            return dir
        except:
            return None

    def validate_mol(self, molecule):
        """
        Check if this molecule has a valid source file. If yes, it is OK to 
        copy. If no, ask the user if they need to modify this molecule.
        """
        exists = self.mol_exists(molecule)
        
        if not exists:
            result = tk.messagebox.askyesno(
                message="The source path for {} doesn't exist. Do you want to modify this molecule?".format(molecule.name))
            if result:
                self.master.master.cmds.molresumeaction = self.ok_action
                self.master.master.cmds.newmolecule = molecule
                self.master.master.cmds.edit_molecule_action()
                return False
            else:
                return True
        else:
            return True

    def mol_exists(self, molecule):
        """
        Check if a molecule's source file exists
        """
        exists = True
        try:
            exists = os.path.exists(molecule.path)
        except:
            exists = False
        return exists

    def copy_mol(self, molecule, destdir):
        """
        Copy the molecule source file to the destination directory. If the 
        source file does not exist, suppress any errors and do not copy it.
        This is called after validation, which gives the user the opportunity
        to modify any molecules whose source file does not exist.

        This will also modify the molecule to refer to only the filename with 
        no path component.
        """
        if self.mol_exists(molecule):
            try:
                shutil.copy2(molecule.path, destdir)
                molecule.path = os.path.basename(molecule.path)
            except OSError as e:
                tk.messagebox.showinfo(message=e.strerror)

    def write_project_file(self):
        """
        Write a file corresponding to the current project's contents.
        TODO: modify the JSON exporter so that the layers' lattices output more
        than one element per line.
        """
        proj = {
            'settings' : {
                'lattice_width' : self.project.lattice_width,
                'lattice_height' : self.project.lattice_height,
                'lattice_spacing' : self.project.lattice_spacing,
                'lattice_major_gridlines' : self.project.lattice_major_gridlines,
                'import_solute' : self.project.import_solute,
                'solute_buffer_space' : self.project.solute_buffer_space,
                'solute_z' : self.project.solute_z
            },
            'molecules' : [],
            'blenders' : [],
            'layers' : [],
        }

        for mol in self.project.molecules:
            if mol.index == 0:
                continue
            proj['molecules'].append({
                'path' : mol.path,
                'name' : mol.name,
                'color' : mol.color,
                'index' : mol.index,
                'flip' : mol.flipbool,
                'rotate' : mol.rotatebool,
            })

        for blend in self.project.blenders:
            weights = {}
            for key, val in blend.molecules.items():
                weights[key.index] = val
            proj['blenders'].append({
                'name' : blend.name,
                'index' : blend.index,
                'weights' : weights,
            })

        for layer in self.project.layers:
            proj['layers'].append({
                'name' : layer.name,
                'zpos' : layer.zdepth,
                'lattice' : layer.lattice.tolist(),
            })
        
        try:
            with open(self.savefilevar.get(), 'w') as savefile:
                json.dump(proj, savefile, indent=4)
        except Exception as e:
            tk.messagebox.showinfo(message=str(e))
