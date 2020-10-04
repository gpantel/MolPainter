import tkinter as tk

from .project import Project
from .commands import Commands
from .cmdbar import CmdBar
from .toolbox import Toolbox
from .molecule_lister import MoleculeLister
from .drawarea import Drawarea


class Painter():
    def __init__(self):
        self.project = Project()

        self.root = tk.Tk()
        self.root.title("MolPainter")
        self.root.option_add('*tearOff', False)
        self.gui = PainterWindow(project=self.project, master=self.root)

        self.gui.cmds.new_project_action()

        self.gui.mainloop()


class PainterWindow(tk.Frame):
    def __init__(self, project, master=None):
        super().__init__(master)
        self.project = project
        self.master = master
        self.grid(column=0, row=0)

        self.cmds = Commands(gui=self, project=project)
        self.cmdbar = CmdBar(self)
        self.cmdbar.grid(column=0, row=0, columnspan=9, rowspan=1, sticky=("W", "N"))
        self.toolbox = Toolbox(self)
        self.toolbox.grid(column=0, row=1, columnspan=1, rowspan=7, sticky=("W", "N"))
        self.mlist = MoleculeLister(self)
        self.mlist.grid(column=0, row=8, columnspan=1, rowspan=9, sticky=("S", "W", "N"))
        self.drawarea = Drawarea(self)
        self.drawarea.grid(column=1, row=1, columnspan=8, rowspan=14, sticky=("E", "S", "W", "N"))

        # Start out with Pencil tool selected
        self.toolbox.btn_pencil.invoke()

    def layer_created(self, layer):
        """
        Used to notify the GUI that a new Z Layer was added
        """
        self.drawarea.add_layer(layer)
        if self.project.import_solute is not None:
            self.project.overlay_solute(layer)
            self.drawarea.refresh_tabs()

    def layer_modified(self, layer):
        """
        Used to notify the GUI that a Z Layer was modified
        """
        if self.project.import_solute is not None:
            self.project.overlay_solute(layer)
        self.drawarea.refresh_tabs()

    def layer_removed(self, layer):
        """ 
        Used to notify the GUI that a Z Layer will be deleted
        """
        self.drawarea.forget_layer(layer)

    def molecule_created(self, molecule):
        """
        Used to notify the GUI that a new molecule was added
        """
        self.mlist.add_molecule(molecule)

    def molecule_modified(self):
        """
        Used to notify the GUI that a molecule was modified
        """
        self.mlist.display_molecules()
        self.drawarea.refresh_tabs()

    def molecule_removed(self, molecule):
        """ 
        Used to notify the GUI that a molecule will be deleted
        """
        self.mlist.forget_molecule(molecule)

    def molecule_selected(self, molecule):
        """
        Used to notify the GUI that a molecule is selected
        """
        self.drawarea.molecule = molecule

    def blend_created(self, blend):
        """
        Used to notify the GUI that a new blender was added
        """
        self.mlist.add_blend(blend) 

    def blend_removed(self, blend):
        """ 
        Used to notify the GUI that a blend will be deleted
        """
        self.mlist.forget_blend(blend)

def main():
    painter = Painter()

if __name__ == "__main__":
    main()
