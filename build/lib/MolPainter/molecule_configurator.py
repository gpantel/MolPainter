import tkinter as tk
from tkinter import colorchooser, filedialog


class MoleculeConfigurator(tk.Frame):
    """
    This defines a window for setting properties of a molecule.
    """
    def __init__(self, name, color, file, flipbool, rotatebool, action, task="create", master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        self.action = action
        self.action.set(None)
        self.task = task
        self.filevar = file

        self.lblf = tk.Label(self, text="Path")
        self.lblf.grid(row=0, column=0, padx=2, pady=4)
        self.file_entry = tk.Entry(self, textvariable=self.filevar)
        self.file_entry.grid(row=0, column=1, columnspan=3, padx=2, pady=4, sticky=("E", "W"))
        self.file_button = tk.Button(self, text="Browse", command=self.file_select_action)
        self.file_button.grid(row=0, column=4, padx=2, pady=4)

        self.lbl = tk.Label(self, text="Name")
        self.lbl.grid(row=1, column=0, padx=2, pady=4)
        self.name_entry = tk.Entry(self, textvariable=name)
        self.name_entry.grid(row=1, column=1, columnspan=3, padx=2, pady=4, sticky=("E", "W"))

        self.colorvar = color
        self.lblc = tk.Label(self, text="Color")
        self.lblc.grid(row=2, column=0, padx=2, pady=4)
        self.color_entry = tk.Entry(self, textvariable=self.colorvar)
        self.canvas = tk.Canvas(self, width=20, height=20, background=self.colorvar.get())
        self.canvas.grid(row=2, column=1, padx=2, pady=4)
        self.color_entry.grid(row=2, column=2, columnspan=2, padx=2, pady=4, sticky=("E", "W"))
        self.color_button = tk.Button(self, text="Select", command=self.color_select_action)
        self.color_button.grid(row=2, column=4, padx=2, pady=4)

        self.flipbool = flipbool
        self.flip_checkbox = tk.Checkbutton(self, text="Flip 180"+u"\N{DEGREE SIGN}"+" over x-axis", variable=self.flipbool)
        self.flip_checkbox.grid(row=3, column=2, columnspan=2, padx=2, pady=4, sticky=("W"))
        self.rotatebool = rotatebool
        self.rotate_checkbox = tk.Checkbutton(self, text="Randomly rotate around z-axis", variable=self.rotatebool)
        self.rotate_checkbox.grid(row=4, column=2, columnspan=2, padx=2, pady=4, sticky=("W"))

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="OK", command=self.ok_action)
        self.cancel_button.grid(row=5, column=0, padx=2, pady=4)
        self.ok_button.grid(row=5, column=4, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)
        self.master.bind('<Destroy>', self.destroy_event)

    def file_select_action(self):
        self.filevar.set(filedialog.askopenfilename())

    def color_select_action(self):
        newcolor = self.colorvar.get()
        try:
            newcolor = colorchooser.askcolor(initialcolor=newcolor)
        except Exception:
            newcolor = colorchooser.askcolor()

        if newcolor[1] is not None:
            self.colorvar.set(newcolor[1])
            self.canvas["background"] = newcolor[1]

    def cancel_action(self):
        self.action.set("cancel")
        self.master.destroy()

    def ok_action(self):
        self.action.set(self.task)
        self.master.destroy()

    def destroy_event(self, event):
        self.action.set("cancel")

    def return_event(self, event):
        self.ok_button.invoke()
        
