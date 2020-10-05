import tkinter as tk


class GridConfigurator(tk.Frame):
    """
    This defines a window for setting dimensions of the project at a global level.
    """
    def __init__(self, gridwidth, gridheight, gridspacing, gridlines, gridaction, task="modify", master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        self.gridaction = gridaction
        self.gridaction.set(None)
        self.task = task
        
        self.widthvar = gridwidth
        self.lbl_width = tk.Label(self, text="Lattice width")
        self.lbl_width.grid(row=0, column=0, padx=2, pady=4)
        self.width_entry = tk.Entry(self, textvariable=self.widthvar)
        self.width_entry.grid(row=0, column=1, padx=2, pady=4)

        self.heightvar = gridheight
        self.lbl_height = tk.Label(self, text="Lattice height")
        self.lbl_height.grid(row=1, column=0, padx=2, pady=4)
        self.height_entry = tk.Entry(self, textvariable=self.heightvar)
        self.height_entry.grid(row=1, column=1, padx=2, pady=4)

        self.spacingvar = gridspacing
        self.lbl_area = tk.Label(self, text="Cell spacing " + u"\u212B")
        self.lbl_area.grid(row=2, column=0, padx=2, pady=4)
        self.area_entry = tk.Entry(self, textvariable=self.spacingvar)
        self.area_entry.grid(row=2, column=1, padx=2, pady=4)

        self.linesvar = gridlines
        self.lbl_lines = tk.Label(self, text="Major gridlines")
        self.lbl_lines.grid(row=3, column=0, padx=2, pady=4)
        self.lines_entry = tk.Entry(self, textvariable=self.linesvar)
        self.lines_entry.grid(row=3, column=1, padx=2, pady=4)

        self.lbl3 = tk.Label(self)
        self.lbl3.grid(row=4, column=0, columnspan=2, padx=2, pady=4)
        if self.master.master.project.import_solute is not None:
            self.lbl3["text"] = "Warning: a solute has been imported.\nChanging cell spacing may disrupt the painting."

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="OK", command=self.ok_action)
        self.cancel_button.grid(row=5, column=0, padx=2, pady=4)
        self.ok_button.grid(row=5, column=1, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)
        self.master.bind('<Destroy>', self.destroy_event)

    def cancel_action(self):
        self.gridaction.set("cancel")
        self.master.destroy()

    def ok_action(self):
        try:
            int(self.widthvar.get())
            int(self.heightvar.get())
            int(self.spacingvar.get())
        except ValueError:
            self.lbl3["text"] = "Only numeric values can be used"
            return
        self.gridaction.set(self.task)
        self.master.destroy()

    def destroy_event(self, event):
        self.gridaction.set("cancel")

    def return_event(self, event):
        self.ok_button.invoke()
        
