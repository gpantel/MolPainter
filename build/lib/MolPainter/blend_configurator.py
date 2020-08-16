import random
import json
import tkinter as tk
from tkinter import colorchooser, filedialog


class BlendConfigurator(tk.Frame):
    """
    This defines a window for defining a blend of molecules.
    """
    def __init__(self, name, vals, action, task="create", master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        self.action = action
        self.action.set(None)
        self.task = task
        self.molecules = []
        self.valsvar = vals
        self.resultdict = {}
        self.set_vals_from_string(self.valsvar.get())

        self.preview = tk.Canvas(self, width=20, height=20, background="#FFFFFF")
        self.preview.grid(row=0, column=0, padx=2, pady=4)
        self.preview_image = tk.PhotoImage(master=self, width=20, height=20)
        self.preview.create_image("0", "0", anchor="nw", image=self.preview_image)
        self.lbl = tk.Label(self, text="Name")
        self.lbl.grid(row=0, column=1, padx=2, pady=4)
        self.name_entry = tk.Entry(self, textvariable=name)
        self.name_entry.grid(row=0, column=2, columnspan=2, padx=2, pady=4)

        self.canvas = tk.Canvas(self, borderwidth=0)
        self.canvas.grid(row=1, column=0, columnspan=3, padx=4, pady=1, sticky=("E", "S", "W", "N"))
        self.mlist = tk.Frame(self.canvas)
        self.mlist.grid(column=0, row=0)
        self.molheader = tk.Label(self.mlist, text="Molecule")
        self.weightheader = tk.Label(self.mlist, text="Weight")
        self.molheader.grid(row=0, column=1)
        self.weightheader.grid(row=0, column=2)

        self.lbl3 = tk.Label(self)
        self.lbl3.grid(row=2, column=0, columnspan=4, padx=2, pady=4)

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="OK", command=self.ok_action)
        self.cancel_button.grid(row=3, column=0, padx=2, pady=4)
        self.ok_button.grid(row=3, column=4, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)
        self.master.bind('<Destroy>', self.destroy_event)
        self.list_molecules()
        self.draw_preview()

    def list_molecules(self):
        for mol in self.master.master.project.molecules:
            var=tk.StringVar()
            var.set(str(self.resultdict.get(mol.index,"0")))
            self.molecules.append({
                "color": tk.Canvas(self.mlist, width=20, height=20, background=mol.color),
                "label": tk.Label(self.mlist, text=str(mol.name)),
                "weight": var,
                "entry": tk.Entry(self.mlist, textvariable=var, validate="none", validatecommand=self.evaluate),
                "molecule": mol
            })
        self.display_molecules()

    def display_molecules(self):
        for i in range(len(self.molecules)):
            self.molecules[i]["color"].grid(row=i+1, column=0, padx=0, pady=4, sticky="W")
            self.molecules[i]["label"].grid(row=i+1, column=1, padx=0, pady=4, sticky="W")
            self.molecules[i]["entry"].grid(row=i+1, column=2, padx=0, pady=4, sticky="E")
            self.molecules[i]["entry"]["validate"] = "all"

    def evaluate(self):
        self.resultdict = {}
        for mol in self.molecules:
            try:
                index = mol["molecule"].index
                weight = int(mol["weight"].get())
                if weight > 0:
                    self.resultdict[index] = weight
            except ValueError as e:
                self.resultdict = {}
                return 1
        if not self.resultdict:
            self.lbl3["text"] = "The weights must be more than zero"
        else:
            self.lbl3["text"] = ""
        self.draw_preview()
        return 1

    def draw_preview(self):
        for row in range(int(self.preview_image["height"])):
            for col in range(int(self.preview_image["width"])):
                self.preview_image.put(self.get_color(), to=(row, col))

    def get_color(self):
        if not self.resultdict:
            return "#FFFFFF"
        choice = random.choices(population=list(self.resultdict.keys()), weights=list(self.resultdict.values()))[0]
        for mol in self.molecules:
            if mol["molecule"].index == choice:
                return mol["molecule"].color
        return "#FFFFFF"

    def set_vals_from_string(self, valstr):
        try:
            srcdic = json.loads(valstr)
        except json.decoder.JSONDecodeError:
            return
        for key, val in srcdic.items():
            self.resultdict[int(key)] = val

    def cancel_action(self):
        self.action.set("cancel")
        self.master.destroy()

    def ok_action(self):
        self.evaluate()
        if not self.resultdict:
            self.lbl3["text"] = "The weights must be positive integers"
            return
        self.valsvar.set(json.dumps(self.resultdict))
        self.action.set(self.task)
        self.master.destroy()

    def destroy_event(self, event):
        self.action.set("cancel")

    def return_event(self, event):
        self.ok_button.invoke()
        
