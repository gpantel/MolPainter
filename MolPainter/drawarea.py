import tkinter as tk
from tkinter import ttk

from .layer_painter import LayerPainter


class LayerGrid(tk.Frame):
    """
    Represents the paintable canvas associated with one layer. One is contained
    inside each tab.
    """
    def __init__(self, zlayer, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        self.zlayer = zlayer
        self.paintarea = LayerPainter(self)
        self.paintarea.grid(column=0, row=0)


class Drawarea(tk.Frame):
    """
    The main area of the window where drawing takes place.
    """
    def __init__(self, master=None):
        super().__init__(master, relief="sunken", borderwidth=1)
        self.master = master
        self.grid(column=0, row=0)

        self.tabs = ttk.Notebook(self)
        self.tabs.grid(column=0, row=0)
        self.layers = []
        self.zoomvalues = [5, 7, 10, 15, 20]
        self.zoom = 2

        self.tool = 'pencil'

        self.molecule = None

        self.actions = tk.Menu(self)
        self.actions.add_command(label="Modify Layer", command=self.master.cmds.edit_layer_action)
        self.actions.add_command(label="Delete Layer", command=self.master.cmds.delete_layer_action)
        self.tabs.bind('<<NotebookTabChanged>>', self.count_molecules)
        self.tabs.bind('<Control-1>', self.show_actions)
        self.tabs.bind('<2>', self.show_actions)
        self.tabs.bind('<3>', self.show_actions)

    def count_molecules(self, event):
        cur_tab = self.tabs.nametowidget(self.tabs.select())
        cur_tab.paintarea.count_molecules()

    def show_actions(self, event):
        cur_layer = self.tabs.nametowidget(self.tabs.select()).zlayer
        self.master.cmds.newlayer = cur_layer
        self.actions.entryconfigure(0, label='Edit "{}"'.format(cur_layer.name))
        self.actions.entryconfigure(1, label='Delete "{}"'.format(cur_layer.name))
        self.actions.post(event.x_root, event.y_root)

    def add_layer(self, zlayer):
        """
        Add a tab for this layer.
        """
        newlayer = LayerGrid(zlayer=zlayer, master=self.tabs)
        self.layers.append(newlayer)
        self.tabs.add(newlayer, text=newlayer.zlayer.name)

    def forget_layer(self, zlayer):
        """
        Find the tab associated with this layer and remove it. Call this before
        deleting a layer.
        """
        for tabid in self.tabs.tabs():
            if self.tabs.nametowidget(tabid).zlayer == zlayer:
                self.tabs.forget(tabid)
        
    def forget_all_layers(self):
        """
        Delete all tabs - can be used when creating a new project
        """
        for tabid in self.tabs.tabs():
            self.tabs.forget(tabid)
        
    def refresh_tabs(self):
        """
        Makes every tab display the name of the associated layer. Call this after
        changing a layer name so that the tabs match.
        """
        for tabid in self.tabs.tabs():
            self.tabs.tab(tabid, text=self.tabs.nametowidget(tabid).zlayer.name)
            self.tabs.nametowidget(tabid).paintarea.count_molecules()

        for layer in self.layers:
            layer.paintarea.refresh_image()

    def zoom_in(self):
        """
        Increase zoom level, if possible
        """
        if self.zoom < len(self.zoomvalues)-1:
            self.zoom += 1
            self.refresh_tabs()

    def zoom_out(self):
        """
        Decrease zoom level, if possible
        """
        if self.zoom > 0:
            self.zoom -= 1
            self.refresh_tabs()

    def set_tool(self, tool):
        """
        Change what tool will be used to paint on the layer
        """
        self.tool = tool
        for layer in self.layers:
            if tool == 'pencil':
                layer.paintarea.canvas.config(cursor="pencil")
            elif tool == 'rect':
                layer.paintarea.canvas.config(cursor="cross")
            elif tool == 'spray':
                layer.paintarea.canvas.config(cursor="spraycan")
        
