import tkinter as tk

import MolPainter

class AboutDialog(tk.Frame):
    """
    This defines a window that directs users to the documentation repository.
    """
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)

        text_header = "MolPainter {}\n\n".format(MolPainter.__version__)
        text_url = "For documentation, see\n{}\n\n".format(MolPainter.__url__)
        text_tutorial = "An end to end tutorial can be found at\n{}".format(MolPainter.__tutorial_url__)
        
        self.about_text = tk.Text(self, width=62, height=7)
        self.about_text.insert('end', text_header, ('header'))
        self.about_text.insert('end', text_url)
        self.about_text.insert('end', text_tutorial)
        self.about_text.grid(row=0, column=0, padx=2, pady=4)
        # self.about_text.tag_configure('header', font='bold')
        self.about_text['state'] = 'disabled'

        self.ok_button = tk.Button(self, text="OK", command=self.ok_action)
        self.ok_button.grid(row=5, column=0, padx=2, pady=4)

        self.master.bind('<Return>', self.return_event)

    def ok_action(self):
        self.master.destroy()

    def return_event(self, event):
        self.ok_button.invoke()
        
