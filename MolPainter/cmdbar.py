import tkinter as tk
from tkinter import ttk, messagebox

from .layer_configurator import LayerConfigurator
from .molecule_configurator import MoleculeConfigurator


class CmdBar(tk.Frame):
    """
    This defines the row of buttons across the top of the window.
    """
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        self.colindex = 0

        self.btn_new = tk.Button(self, text="New", command=self.master.cmds.new_project_action)
        self.btn_open = tk.Button(self, text="Open Painting", command=self.master.cmds.load_project_action)
        self.btn_save = tk.Button(self, text="Save Painting", command=self.master.cmds.save_project_action)
        self.btn_export = tk.Button(self, text="Export System", command=self.master.cmds.export_action)
        self.btn_import = tk.Button(self, text="Import Solute", command=self.master.cmds.import_solute_action)
        self.btn_newmol = tk.Button(self, text="New Molecule", command=self.master.cmds.new_molecule_action)
        self.btn_newblend = tk.Button(self, text="New Blend", command=self.master.cmds.new_blend_action)
        self.btn_newlayer = tk.Button(self, text="New Layer", command=self.master.cmds.new_layer_action)
        self.btn_gridsize = tk.Button(self, text="Grid Settings", command=self.master.cmds.grid_setup_action)
        # self.btn_zoomout = tk.Button(self, text="Zoom Out", command=self.master.cmds.zoom_out_action)
        # self.btn_zoomin = tk.Button(self, text="Zoom In", command=self.master.cmds.zoom_in_action)
        #self.btn_solvents = tk.Button(self, text="Solvents", command=self.master.cmds.no_action)
        
        self.btn_new.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_open.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_save.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_export.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_import.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_newmol.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_newblend.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_newlayer.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        self.btn_gridsize.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        # self.btn_zoomout.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        # self.btn_zoomin.grid(column=self.col, row=0, columnspan=1, rowspan=1)
        # self.btn_solvents.grid(column=self.col, row=0, columnspan=1, rowspan=1)
    
    @property
    def col(self):
        self.colindex += 1
        return self.colindex
