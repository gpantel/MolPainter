import os
import tkinter as tk
from tkinter import ttk

class Toolbox(tk.Frame):
    """
    This defines the painting tool buttons at the side of the window.
    """
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)

        icondir = os.path.join(os.path.dirname(__file__), 'icons')

        self.pencil_icon = tk.PhotoImage(file=os.path.join(icondir, 'pencil.png'))
        self.rect_icon = tk.PhotoImage(file=os.path.join(icondir, 'square.png'))
        self.spray_icon = tk.PhotoImage(file=os.path.join(icondir, 'spraycan.png'))
        self.zoomin_icon = tk.PhotoImage(file=os.path.join(icondir, 'zoomin.png'))
        self.zoomout_icon = tk.PhotoImage(file=os.path.join(icondir, 'zoomout.png'))

        self.btn_zoomout = tk.Button(self, image=self.zoomout_icon, command=self.master.cmds.zoom_out_action, height=40, width=40)
        self.btn_zoomin = tk.Button(self, image=self.zoomin_icon, command=self.master.cmds.zoom_in_action, height=40, width=40)
        self.btn_pencil = tk.Button(self, image=self.pencil_icon, command=self.activate_pencil, height=40, width=40)
        self.btn_rect = tk.Button(self, image=self.rect_icon, command=self.activate_rect, height=40, width=40)
        self.btn_spray = tk.Button(self, image=self.spray_icon, command=self.activate_spray, height=40, width=40)

        self.spray_width_scale = tk.Scale(self, from_=1, to=6, orient="horizontal", label='Spray radius', length=130)
        
        self.btn_zoomout.grid(column=0, row=0, columnspan=1, rowspan=1, sticky=("W"))
        self.btn_zoomin.grid(column=1, row=0, columnspan=1, rowspan=1, sticky=("W"))
        self.btn_pencil.grid(column=0, row=1, columnspan=1, rowspan=1, sticky=("W"))
        self.btn_rect.grid(column=1, row=1, columnspan=1, rowspan=1, sticky=("W"))
        self.btn_spray.grid(column=2, row=1, columnspan=1, rowspan=1, sticky=("W"))
        self.spray_width_scale.grid(column=0, row=5, columnspan=3, rowspan=1, sticky=("N", "S", "W"))

    def activate_pencil(self):
        #self.btn_pencil.state(['pressed'])
        #self.btn_rect.state(['!pressed'])
        #self.btn_spray.state(['!pressed'])
        self.master.drawarea.set_tool('pencil')

    def activate_rect(self):
        #self.btn_pencil.state(['!pressed'])
        #self.btn_rect.state(['pressed'])
        #self.btn_spray.state(['!pressed'])
        self.master.drawarea.set_tool('rect')

    def activate_spray(self):
        #self.btn_pencil.state(['!pressed'])
        #self.btn_rect.state(['!pressed'])
        #self.btn_spray.state(['pressed'])
        self.master.drawarea.set_tool('spray')

    def open_noddle_popout(self):
        """
        Opens a popous window to the scale bar for the spray can
        """
        nozzle_popup = tk.Toplevel(self.gui)
