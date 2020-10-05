import tkinter as tk
from tkinter import ttk

import numpy as np


class MoleculeLister(tk.Frame):
    """
    This defines the list of molecules and blends shown at the side of the window.
    """
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0, padx=5, pady=5)
        self.canvas = tk.Canvas(self, borderwidth=1, relief="solid", background="#FFFFFF")
        self.canvas.grid(column=0, row=0, padx=1, pady=1, sticky=("E", "S", "W", "N"))
        self.mlist = tk.Frame(self.canvas, background="#FFFFFF")
        self.mlist.grid(column=0, row=0)
        
        self.selected = tk.IntVar()
        self.selected.set(0)
        self.selected.trace_variable("w", self.activate_molecule)

        self.molecules = []
        self.blends = []

        self.vscroll = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vscroll.set)

        self.actions = tk.Menu(self)
        self.actions.add_command(label="Modify Molecule", command=self.master.cmds.edit_molecule_action)
        self.actions.add_command(label="Copy Molecule", command=self.master.cmds.copy_molecule_action)
        self.actions.add_command(label="Delete Molecule", command=self.master.cmds.delete_molecule_action)

    def show_actions(self, event, molecule):
        """
        Display a context menu for the molecule that was clicked. The empty
        molecule can not be modified or deleted. If a blend was clicked, set
        the menu actions to modify the blend.
        """
        state="normal"
        if molecule.index == 0:
            state="disabled"
            edit_cmd=None
            copy_cmd=None
            del_cmd=None
        elif molecule.index < self.master.project.blender_offset:
            self.master.cmds.newmolecule = molecule
            edit_cmd = self.master.cmds.edit_molecule_action
            copy_cmd = self.master.cmds.copy_molecule_action
            del_cmd = self.master.cmds.delete_molecule_action
        else:
            self.master.cmds.newblend = molecule
            edit_cmd = self.master.cmds.edit_blend_action
            copy_cmd = self.master.cmds.copy_blend_action
            del_cmd = self.master.cmds.delete_blend_action

        self.actions.entryconfigure(0, label='Edit "{}"'.format(molecule.name), state=state, command=edit_cmd)
        self.actions.entryconfigure(1, label='Copy "{}"'.format(molecule.name), state=state, command=copy_cmd)
        self.actions.entryconfigure(2, label='Delete "{}"'.format(molecule.name), state=state, command=del_cmd)
        
        self.actions.post(event.x_root, event.y_root)

    def add_molecule(self, molecule):
        """
        Adds this molecule to the list. Call this after creating a molecule.
        """
        number_tkvar = tk.StringVar()
        percent_tkvar = tk.StringVar()
        number_tkvar.set('0')
        percent_tkvar.set('0.0%')
        self.molecules.append({
            "color": tk.Canvas(self.mlist, width=20, height=20, background=molecule.color),
            "entry": tk.Radiobutton(self.mlist, text=str(molecule.name), variable=self.selected, value=molecule.index, background="#FFFFFF"),
            "number_tkvar": number_tkvar,
            "percent_tkvar": percent_tkvar,
            "number": tk.Label(self.mlist, textvariable=number_tkvar, background="#FFFFFF", width=5),
            "percent": tk.Label(self.mlist, textvariable=percent_tkvar, background="#FFFFFF", width=5),
            "molecule": molecule
            })
        self.display_molecules()
        self.molecules[-1]["entry"].select()

    def forget_molecule(self, molecule):
        """
        Remove this molecule from the list. Call this before deleting a molecule.
        """
        for mol in self.molecules:
            if mol["molecule"] == molecule:
                self.molecules.remove(mol)
                if self.selected.get() == molecule.index:
                    self.molecules[0]["entry"].select()
                return

    def add_blend(self, blend):
        """
        Adds this blend definition to the list.
        """
        newblend = {
            "color": tk.Canvas(self.mlist, width=20, height=20, background="#FFFFFF"),
            "image": tk.PhotoImage(master=self.mlist, width=20, height=20),
            "entry": tk.Radiobutton(self.mlist, text=str(blend.name), variable=self.selected, value=blend.index, background="#FFFFFF"),
            "blend": blend
        }
        newblend["color"].create_image("0", "0", anchor="nw", image=newblend["image"])
        self.draw_preview(blend, newblend["image"])
        self.blends.append(newblend)
        self.display_molecules()
        self.blends[-1]["entry"].select()

    def draw_preview(self, blend, image):
        """
        Updates the colorbox for a blend with a weighted blend of colors.
        """
        for row in range(int(image["height"])):
            for col in range(int(image["width"])):
                image.put(blend.get_paint_props()[1], to=(row, col))

    def forget_blend(self, blend):
        """
        Remove this blend from the list. Call this before deleting a blend.
        """
        for b in self.blends:
            if b["blend"] == blend:
                self.blends.remove(b)
                if self.selected.get() == blend.index:
                    self.molecules[0]["entry"].select()
                return

    def display_molecules(self):
        """
        Refreshes the molecule list to show current info on all molecules. Call
        this after modifying any molecules.
        """
        for element in self.mlist.grid_slaves():
            element.grid_forget()
        molcount = len(self.molecules)
        for i in range(molcount):
            colorbox = self.molecules[i]["color"]
            label = self.molecules[i]["entry"]
            number = self.molecules[i]["number"]
            percent = self.molecules[i]["percent"]
            molecule = self.molecules[i]["molecule"]
            colorbox["background"] = molecule.color
            label["text"] = str(molecule.name)
            self.display_element(colorbox, row=i, col=0, molecule=molecule)
            self.display_element(label, row=i, col=1, molecule=molecule)
            if i == 0:
                continue
            else:
                self.display_element(number, row=i, col=2, molecule=molecule)
                self.display_element(percent, row=i, col=3, molecule=molecule)
        for i in range(len(self.blends)):
            colorbox = self.blends[i]["color"]
            label = self.blends[i]["entry"]
            blend = self.blends[i]["blend"]
            label["text"] = str(blend.name)
            self.draw_preview(self.blends[i]["blend"], self.blends[i]["image"])
            self.display_element(colorbox, row=i+molcount, col=0, molecule=blend)
            self.display_element(label, row=i+molcount, col=1, molecule=blend)

    def display_element(self, element, row, col, molecule):
        """
        Puts the individual components of a molecule list item on the screen
        along with a context menu that is aware of the related molecule.
        """
        element.grid(column=col, row=row, padx=0, pady=4, sticky="W")
        element.bind('<Control-1>', lambda e: self.show_actions(e, molecule))
        element.bind('<2>', lambda e: self.show_actions(e, molecule))
        element.bind('<3>', lambda e: self.show_actions(e, molecule))

    def activate_molecule(self, arg1, arg2, arg3):
        """
        Sets the active molecule to be the one that was just clicked. If it is
        a blend, treat it as a molecule.
        """
        for mol in self.molecules:
            if mol["molecule"].index == self.selected.get():
                self.master.molecule_selected(mol["molecule"])
                return
        for b in self.blends:
            if b["blend"].index == self.selected.get():
                self.master.molecule_selected(b["blend"])
                return
