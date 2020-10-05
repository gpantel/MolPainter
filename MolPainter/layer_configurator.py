import tkinter as tk


class LayerConfigurator(tk.Frame):
    """
    This defines a window for setting properties of a single layer.
    """
    def __init__(self, layername, layerdepth, layeraction, task="create", master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        self.layeraction = layeraction
        self.layeraction.set(None)
        self.task = task
        
        self.lbl = tk.Label(self, text="Name")
        self.lbl.grid(row=0, column=0, padx=2, pady=4)
        self.name_entry = tk.Entry(self, textvariable=layername)
        self.name_entry.grid(row=0, column=1, padx=2, pady=4)

        self.depthvar = layerdepth
        self.lbl2 = tk.Label(self, text="z position " + u"(\u212B)")
        self.lbl2.grid(row=1, column=0, padx=2, pady=4)
        self.name_entry = tk.Entry(self, textvariable=self.depthvar)
        self.name_entry.grid(row=1, column=1, padx=2, pady=4)
        
        self.lbl3 = tk.Label(self)
        self.lbl3.grid(row=2, column=0, columnspan=2, padx=2, pady=4)

        self.cancel_button = tk.Button(self, text="Cancel", command=self.cancel_action)
        self.ok_button = tk.Button(self, text="OK", command=self.ok_action)
        self.cancel_button.grid(row=3, column=0, padx=2, pady=4)
        self.ok_button.grid(row=3, column=1, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)
        self.master.bind('<Destroy>', self.destroy_event)

    def cancel_action(self):
        self.layeraction.set("cancel")
        self.master.destroy()

    def ok_action(self):
        try:
            int(self.depthvar.get())
        except ValueError:
            self.lbl3["text"] = "The z position must be a number"
            return
        self.layeraction.set(self.task)
        self.master.destroy()

    def destroy_event(self, event):
        self.layeraction.set("cancel")

    def return_event(self, event):
        self.ok_button.invoke()
        