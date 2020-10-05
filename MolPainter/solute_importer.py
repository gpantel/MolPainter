import tkinter as tk
from tkinter import filedialog


class SoluteImporter(tk.Frame):
    """
    This defines a window for importing a PDB file that will be merged with the project.
    """
    def __init__(self, project, file, master=None):
        super().__init__(master)
        self.master = master
        self.project = project
        self.grid(column=0, row=0)
        self.filevar = file
        self.bufvar = tk.StringVar()
        self.bufvar.set(self.project.solute_buffer_space)

        self.lblf = tk.Label(self, text="Path")
        self.lblf.grid(row=0, column=0, padx=2, pady=4)
        self.file_entry = tk.Entry(self, textvariable=self.filevar)
        self.file_entry.grid(row=0, column=1, columnspan=3, padx=2, pady=4, sticky=("E", "W"))
        self.file_button = tk.Button(self, text="Browse", command=self.file_select_action)
        self.file_button.grid(row=0, column=4, padx=2, pady=4)

        self.lblb = tk.Label(self, text="Buffer space " + u"(\u212B)")
        self.lblb.grid(row=1, column=1, columnspan=2, padx=2, pady=4)
        self.buf_entry = tk.Entry(self, textvariable=self.bufvar)
        self.buf_entry.grid(row=1, column=3, padx=2, pady=4)

        self.centerbool = tk.BooleanVar()
        self.centervar = tk.StringVar()
        self.centervar.set("0")
        if self.project.solute_z is not None:
            self.centervar.set(str(self.project.solute_z))
            self.centerbool.set(True)
        self.center_checkbox = tk.Checkbutton(self, text="Center solute at z "+ u"(\u212B)", variable=self.centerbool)
        self.center_checkbox.grid(row=2, column=1, columnspan=2, padx=2, pady=4, sticky=("W"))
        self.center_entry = tk.Entry(self, textvariable=self.centervar)
        self.center_entry.grid(row=2, column=3, padx=2, pady=4)

        self.expandbool = tk.BooleanVar()
        self.expand_checkbox = tk.Checkbutton(self, text="Expand grid to fit solute", variable=self.expandbool)
        self.expand_checkbox.grid(row=3, column=1, columnspan=3, padx=2, pady=4, sticky=("W"))

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="OK", command=self.ok_action)
        self.cancel_button.grid(row=4, column=0, padx=2, pady=4)
        self.ok_button.grid(row=4, column=4, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)

    def file_select_action(self):
        self.filevar.set(filedialog.askopenfilename(filetypes = (("PDB file","*.pdb"),("all files","*.*"))))

    def cancel_action(self):
        self.master.destroy()

    def ok_action(self):
        # Apply settings before loading
        if self.centerbool.get():
            center = self.centervar.get()
        else:
            center = None
        self.project.edit_solute_settings(self.filevar.get(), self.bufvar.get(), center)
        # Try to load the solute into the project
        try:
            self.project.load_solute(self.expandbool.get())
        except Exception as e:
            tk.messagebox.showinfo(message=e.strerror)
            self.project.import_solute = None
            return
        self.master.master.drawarea.refresh_tabs()
        self.master.destroy()

    def return_event(self, event):
        self.ok_button.invoke()
        
