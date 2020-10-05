import json
import os
import tkinter as tk
import MDAnalysis as mda
import numpy as np

from .project import Molecule, Blender, ZLayer
from .layer_configurator import LayerConfigurator
from .molecule_configurator import MoleculeConfigurator
from .blend_configurator import BlendConfigurator
from .grid_configurator import GridConfigurator
from .proj_saver import ProjSaver
from .exporter import Exporter
from .solute_importer import SoluteImporter


class Commands:
    """
    This defines the actions that happen when buttons are clicked.
    """
    def __init__(self, gui, project):
        self.gui=gui
        self.project=project

        self.layername = tk.StringVar()
        self.layerdepth = tk.StringVar()
        self.layeraction = tk.StringVar()
        self.layeraction.trace_variable("w", self.layer_modified)
        self.layer_settings = None

        self.molname = tk.StringVar()
        self.molcolor = tk.StringVar()
        self.molpath = tk.StringVar()
        self.molflipbool = tk.IntVar()
        self.molrotatebool = tk.IntVar()
        self.molaction = tk.StringVar()
        self.molaction.trace_variable("w", self.molecule_modified)
        self.molecule_settings = None
        self.molresumeaction = None

        self.blendname = tk.StringVar()
        self.blendvals = tk.StringVar()
        self.blendaction = tk.StringVar()
        self.blendaction.trace_variable("w", self.blend_modified)
        self.blend_settings = None

        self.gridwidth = tk.StringVar()
        self.gridheight = tk.StringVar()
        self.gridspacing = tk.StringVar()
        self.gridaction = tk.StringVar()
        self.gridlines = tk.StringVar()
        self.gridaction.trace_variable("w", self.grid_modified)
        self.grid_settings = None

        self.exportfilevar = tk.StringVar()
        self.exportfilevar.set("output.pdb")
        self.importfilevar = tk.StringVar()

        self.projfilevar = tk.StringVar()
        self.projfilevar.set(os.path.join(os.path.expanduser('~'), "project.json"))
        self.copyfilesbool = tk.BooleanVar()
        self.copyfilesbool.set(False)


    def new_project_action(self):
        """
        Command to start a new project. Clears the existing settings and 
        creates a default project with one empty layer.
        """
        self.project.init_defaults()
        self.gui.mlist.molecules = []
        self.gui.mlist.blends = []
        self.gui.drawarea.layers = []
        self.gui.drawarea.forget_all_layers()
        newlayer = self.project.new_layer()
        self.gui.layer_created(newlayer)
        self.project.add_layer(newlayer)
        empty_molecule = Molecule()
        self.project.add_molecule(empty_molecule)
        self.gui.molecule_created(empty_molecule)
        self.gui.molecule_selected(empty_molecule)

    def save_project_action(self):
        """
        Command to save the working data for later use. Does not export the 
        actual molecules, for which 'export' should be used instead.
        """
        save_popup = tk.Toplevel(self.gui)
        ProjSaver(self.project, self.projfilevar, self.copyfilesbool, save_popup)

    def load_project_action(self):
        """
        Command to load project data. Attempts to validate molecule paths and
        prompts for update if a path is missing.
        """
        self.projfilevar.set(tk.filedialog.askopenfilename(initialdir = os.path.dirname(self.projfilevar.get()), title = "Select file",filetypes = (("JSON file","*.json"),("all files","*.*"))))
        path = self.projfilevar.get()
        if not self.file_exists(path):
            return
        with open(path, "r") as projfile:
            proj = json.load(projfile)
            self.project.init_defaults()
            self.gui.mlist.molecules = []
            self.gui.mlist.blends = []
            self.gui.drawarea.layers = []
            self.gui.drawarea.forget_all_layers()
            empty_molecule = Molecule()
            self.project.add_molecule(empty_molecule)
            self.gui.molecule_created(empty_molecule)
            self.project.edit_lattice_params(
                proj['settings']['lattice_width'],
                proj['settings']['lattice_height'],
                proj['settings']['lattice_spacing'],
                proj['settings']['lattice_major_gridlines']
                )
            for mol in proj['molecules']:
                newmol = Molecule()
                newmol.name = mol['name']
                newmol.color = mol['color']
                newmol.index = mol['index']
                newmol.flipbool = mol['flip']
                newmol.rotatebool = mol['rotate']
                newmol.path = mol['path']
                if newmol.path and not self.file_exists(newmol.path):
                    # Validation of paths - first check if a copy is included
                    # in the project location, then prompt the user to select
                    # a file if it could not be determined.
                    fname = os.path.basename(newmol.path)
                    copied_file = os.path.join(os.path.dirname(path), fname)
                    if not self.file_exists(copied_file):
                        confirm = tk.messagebox.askyesno(message="The path '{}' does not exist for molecule '{}'. Do you want to select a new path for '{}'?".format(newmol.path, newmol.name, newmol.name))
                        if confirm:
                            newmol.path = tk.filedialog.askopenfilename(initialdir = os.path.dirname(path), title = "Select file",filetypes = (("PDB file","*.pdb"),("all files","*.*")))
                self.project.add_molecule(newmol)
                self.gui.molecule_created(newmol)
            for blend in proj['blenders']:
                newblend = Blender()
                newblend.name = blend['name']
                newblend.index = blend['index']
                for key, val in blend['weights'].items():
                    for mol in self.project.molecules:
                        if mol.index == int(key):
                            newblend.molecules[mol] = val
                self.gui.blend_created(newblend)
                self.project.add_blender(newblend)
            for layer in proj['layers']:
                newlayer = ZLayer()
                newlayer.name = layer['name']
                newlayer.zdepth = layer['zpos']
                newlayer.lattice = np.array(layer['lattice'])
                self.gui.layer_created(newlayer)
                self.project.add_layer(newlayer)
            # New solute import feature. Suppress errors for backward compatibility.
            try:
                if proj['settings']['import_solute']:
                    solute_path = proj['settings']['import_solute']
                    if not self.file_exists(solute_path):
                        solute_path = os.path.join(os.path.dirname(path), os.path.basename(solute_path))
                    self.project.edit_solute_settings(solute_path, proj['settings']['solute_buffer_space'], proj['settings']['solute_z'])
                    try:
                        self.project.load_solute(False)
                    except Exception as e:
                        tk.messagebox.showinfo(message=e.strerror)
                        self.project.import_solute = None
            except KeyError:
                pass
            self.gui.drawarea.refresh_tabs()
            self.project.project_loaded()

    def file_exists(self, path):
        """
        Check if a file exists
        """
        exists = True
        try:
            exists = os.path.exists(path)
        except:
            exists = False
        return exists

    def grid_setup_action(self):
        """
        Command to change grid size and lattice spacing. Uses a configurator to
        collect values 
        """
        if self.grid_settings is not None:
            return
        self.gridwidth.set(self.project.lattice_width)
        self.gridheight.set(self.project.lattice_height)
        self.gridspacing.set(self.project.lattice_spacing)
        self.gridlines.set(self.project.lattice_major_gridlines)
        grid_popup = tk.Toplevel(self.gui)
        self.grid_settings = GridConfigurator(self.gridwidth, self.gridheight, self.gridspacing, self.gridlines, self.gridaction, "modify", grid_popup)

    def grid_modified(self, arg1, arg2, arg3):
        """
        Determine if the grid should be modified with the new settings or if the
        changes should be discarded.
        """
        if self.gridaction.get() == "modify":
            self.project.edit_lattice_params(int(self.gridwidth.get()), int(self.gridheight.get()), int(self.gridspacing.get()), int(self.gridlines.get()))
            if self.project.import_solute is not None:
                for layer in self.project.layers:
                    self.project.overlay_solute(layer)
            self.grid_settings = None
            self.gui.drawarea.refresh_tabs()
        elif self.gridaction.get() == "cancel":
            self.grid_settings = None

    def new_layer_action(self):
        """
        Command to create a new layer. Sets up a layer with default values,
        then hands it off to the configurator to be modified. If the layer
        was accepted, notify the GUI and add it to the project. If it was
        rejected, delete the layer.
        """
        if self.layer_settings is not None:
            return
        self.newlayer = self.project.new_layer()
        self.layername.set(self.newlayer.name)
        self.layerdepth.set(str(self.newlayer.zdepth))
        layer_popup = tk.Toplevel(self.gui)
        self.layer_settings = LayerConfigurator(self.layername, self.layerdepth, self.layeraction, "create", layer_popup)
    
    def edit_layer_action(self):
        """
        Command to edit a layer. Set self.newlayer prior to calling. This takes
        the layer and hands it off to the configurator to be modified.
        If the config was accepted, modify the layer in the project.
        """
        if self.layer_settings is not None:
            return
        self.layername.set(self.newlayer.name)
        self.layerdepth.set(str(self.newlayer.zdepth))
        layer_popup = tk.Toplevel(self.gui)
        self.layer_settings = LayerConfigurator(self.layername, self.layerdepth, self.layeraction, "modify", layer_popup)

    def delete_layer_action(self):
        """
        Command to delete a layer. Set self.newlayer before calling. Asks for 
        confirmation before deleting.
        """
        if self.layer_settings is not None:
            return
        confirm = tk.messagebox.askyesno(message='Really delete layer "{}"?'.format(self.newlayer.name))
        if confirm:
            self.gui.layer_removed(self.newlayer)
            self.project.delete_layer(self.newlayer)
        
        self.layer_settings = None
        self.newlayer = None

    def layer_modified(self, arg1, arg2, arg3):
        """
        Determine if a layer should be modified with the new settings or if the
        changes should be discarded.
        """
        if self.layeraction.get() == "create":
            self.newlayer.name = self.layername.get()
            self.newlayer.zdepth = int(self.layerdepth.get())
            self.project.add_layer(self.newlayer)
            self.gui.layer_created(self.newlayer)
            self.layer_settings = None
            self.newlayer = None
        elif self.layeraction.get() == "modify":
            self.newlayer.name = self.layername.get()
            self.newlayer.zdepth = int(self.layerdepth.get())
            self.gui.layer_modified(self.newlayer)
            self.layer_settings = None
            self.newlayer = None
        elif self.layeraction.get() == "cancel":
            self.layer_settings = None
            self.newlayer = None

    def new_molecule_action(self):
        """
        Command to create a new molecule type. Creates a template with default
        values, then hands it off to the configurator to be modified. If it was
        accepted, notify the GUI and add it to the project. If rejected, 
        simply delete it.
        """
        if self.molecule_settings is not None:
            return
        self.newmolecule = self.project.new_molecule()
        self.molname.set(self.newmolecule.name)
        self.molcolor.set(self.newmolecule.color)
        molecule_popup = tk.Toplevel(self.gui)
        self.molecule_settings = MoleculeConfigurator(self.molname, self.molcolor, self.molpath, self.molflipbool, self.molrotatebool, self.molaction, "create", molecule_popup)

    def copy_molecule_action(self):
        """
        Command to copy a molecule. Takes all properties of the source molecule,
        creates a new molecule with those properties, and opens the configurator
        for further modification.
        """
        if self.molecule_settings is not None:
            return
        self.molname.set(self.newmolecule.name + " copy")
        newcolorval = self.newmolecule.color.lstrip('#')
        red = int(newcolorval[0:2], base=16)
        green = int(newcolorval[2:4], base=16)
        blue = int(newcolorval[4:6], base=16)
        red = int((red + 0xFF) / 2)
        green = int((green + 0xFF) / 2)
        blue = int((blue + 0xFF) / 2)
        self.molcolor.set('#{:02x}{:02x}{:02x}'.format(red, green, blue))
        self.molpath.set(self.newmolecule.path)
        self.molflipbool.set(self.newmolecule.flipbool)
        self.molrotatebool.set(self.newmolecule.rotatebool)
        molecule_popup = tk.Toplevel(self.gui)
        self.newmolecule = self.project.new_molecule()
        self.molecule_settings = MoleculeConfigurator(self.molname, self.molcolor, self.molpath, self.molflipbool, self.molrotatebool, self.molaction, "create", molecule_popup)

    def edit_molecule_action(self):
        """
        Command to edit a molecule. Set self.newmolecule prior to calling. This
        takes the molecule and hands it off to the configurator to be modified.
        If the config was accepted, modify the molecule in the project.
        """
        if self.molecule_settings is not None:
            return
        self.molname.set(self.newmolecule.name)
        self.molcolor.set(self.newmolecule.color)
        self.molpath.set(self.newmolecule.path)
        self.molflipbool.set(self.newmolecule.flipbool)
        self.molrotatebool.set(self.newmolecule.rotatebool)
        molecule_popup = tk.Toplevel(self.gui)
        self.molecule_settings = MoleculeConfigurator(self.molname, self.molcolor, self.molpath, self.molflipbool, self.molrotatebool, self.molaction, "modify", molecule_popup)

    def delete_molecule_action(self):
        """
        Command to delete a molecule. Set self.newmolecule before calling. Asks
        for confirmation before deleting.
        """
        if self.molecule_settings is not None:
            return
        confirm = tk.messagebox.askyesno(message='Really delete molecule "{}"?'.format(self.newmolecule.name))
        if confirm:
            self.gui.molecule_removed(self.newmolecule)
            self.project.delete_molecule(self.newmolecule)
            self.gui.molecule_modified()
        
        self.molecule_settings = None
        self.newmolecule = None

    def molecule_resume_action(self):
        """
        Resume an operation that was interrupted because a molecule needed to
        be modified.
        """
        if self.molresumeaction:
            self.molresumeaction()
            self.molresumeaction = None

    def molecule_modified(self, arg1, arg2, arg3):
        """
        Determines if a molecule should be modified with the new settings or if
        the changes should be discarded.
        """
        if self.molaction.get() == "create":
            self.newmolecule.name = self.molname.get()
            self.newmolecule.path = self.molpath.get()
            self.newmolecule.color = self.molcolor.get()
            self.newmolecule.flipbool = self.molflipbool.get()
            self.newmolecule.rotatebool = self.molrotatebool.get()
            self.gui.molecule_created(self.newmolecule)
            self.project.add_molecule(self.newmolecule)
            self.molecule_settings = None
            self.newmolecule = None
        elif self.molaction.get() == "modify":
            self.newmolecule.name = self.molname.get()
            self.newmolecule.path = self.molpath.get()
            self.newmolecule.color = self.molcolor.get()
            self.newmolecule.flipbool = self.molflipbool.get()
            self.newmolecule.rotatebool = self.molrotatebool.get()
            self.gui.molecule_modified()
            self.molecule_settings = None
            self.newmolecule = None
            self.molecule_resume_action()
        elif self.molaction.get() == "cancel":
            self.molecule_settings = None
            self.newmolecule = None
            self.molecule_resume_action()

    def new_blend_action(self):
        """
        Command to create a new blend of molecules. Creates a template with default
        values, then hands it off to the configurator to be modified. If it was
        accepted, notify the GUI and add it to the project. If rejected, 
        simply delete it.
        """
        if self.blend_settings is not None:
            return
        self.newblend = self.project.new_blender()
        self.blendname.set(self.newblend.name)
        self.blendvals.set("")
        blend_popup = tk.Toplevel(self.gui)
        self.blend_settings = BlendConfigurator(self.blendname, self.blendvals, self.blendaction, "create", blend_popup)

    def copy_blend_action(self):
        """
        Command to copy a blend definition. Takes the definition of the source blend
        and passes it to the configurator to be modified.
        """
        if self.blend_settings is not None:
            return
        self.blendname.set(self.newblend.name + " copy")
        self.blendvals.set(self.get_val_string_from_blend(self.newblend))
        self.newblend = self.project.new_blender()
        blend_popup = tk.Toplevel(self.gui)
        self.blend_settings = BlendConfigurator(self.blendname, self.blendvals, self.blendaction, "create", blend_popup)

    def edit_blend_action(self):
        """
        Command to edit a molecule blender. Set self.newblend prior to calling. This
        takes the selected blend and hands it off to the configurator to be modified.
        If the config was accepted, modify the blender in the project.
        """
        if self.blend_settings is not None:
            return
        self.blendname.set(self.newblend.name)
        self.blendvals.set(self.get_val_string_from_blend(self.newblend))
        blend_popup = tk.Toplevel(self.gui)
        self.blend_settings = BlendConfigurator(self.blendname, self.blendvals, self.blendaction, "modify", blend_popup)

    def delete_blend_action(self):
        """
        Command to delete a blend. Set self.newblend before calling. Asks
        for confirmation before deleting.
        """
        if self.blend_settings is not None:
            return
        confirm = tk.messagebox.askyesno(message='Really delete the blend definition "{}"?\nThe molecules will not be deleted.'.format(self.newblend.name))
        if confirm:
            self.gui.blend_removed(self.newblend)
            self.project.delete_blender(self.newblend)
            self.gui.molecule_modified()
        
        self.blend_settings = None
        self.newblend = None

    def get_val_string_from_blend(self, blend):
        """
        Get a string representation of the molecules in this blend
        """
        weights = {}
        for key, val in blend.molecules.items():
            weights[key.index] = val
        return json.dumps(weights)

    def set_blend_vals_from_string(self, blend, valstr):
        """
        Update the molecules in this blend using a string representation
        """
        blend.molecules = {}
        try:
            srcdic = json.loads(valstr)
        except json.decoder.JSONDecodeError:
            return
        for key, val in srcdic.items():
            for mol in self.project.molecules:
                if mol.index == int(key):
                    blend.molecules[mol] = val

    def blend_modified(self, arg1, arg2, arg3):
        """
        Determines if a blend should be modified with the new settings or if
        the changes should be discarded.
        """
        if self.blendaction.get() == "create":
            self.newblend.name = self.blendname.get()
            self.set_blend_vals_from_string(self.newblend, self.blendvals.get())
            self.gui.blend_created(self.newblend)
            self.project.add_blender(self.newblend)
            self.blend_settings = None
            self.newblend = None
        elif self.blendaction.get() == "modify":
            self.newblend.name = self.blendname.get()
            self.set_blend_vals_from_string(self.newblend, self.blendvals.get())
            self.gui.molecule_modified()
            self.blend_settings = None
            self.newblend = None
        elif self.blendaction.get() == "cancel":
            self.blend_settings = None
            self.newblend = None

    def zoom_in_action(self):
        """
        Command to increase zoom level
        """
        self.gui.drawarea.zoom_in()

    def zoom_out_action(self):
        """
        Command to decrease zoom level
        """
        self.gui.drawarea.zoom_out()

    def no_action(self):
        """
        Placeholder for buttons that have no action
        """
        print("That button doesn't do anything")

    def export_action(self):
        """
        Command to save the output PDB file
        """
        export_popup = tk.Toplevel(self.gui)
        Exporter(self.project, self.exportfilevar, self.projfilevar, export_popup)

    def import_solute_action(self):
        """
        Command to import a PDB system that will be merged with the current project
        """
        if self.project.import_solute is not None:
            self.importfilevar.set(self.project.import_solute)
        import_popup = tk.Toplevel(self.gui)
        SoluteImporter(self.project, self.importfilevar, import_popup)
