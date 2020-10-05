import tkinter as tk
import numpy as np


class LayerPainter(tk.Frame):
    """
    This handles the painting actions
    """
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0, row=0)
        
        self.virtual_spacing = self.master.master.master.zoomvalues[self.master.master.master.zoom] # number of pixels defining width and height of lattice points in pixels
        self.hardpad = 0 #int(self.virtual_spacing/2 )

        self.canvas = tk.Canvas(self, width=(self.lattice_width*self.virtual_spacing)+self.hardpad,\
                        height=(self.lattice_height*self.virtual_spacing)+self.hardpad, borderwidth=0,\
                        background="white", cursor="pencil")
        self.canvas.grid(column=0, row=0, padx=10, pady=10)

        # initialize the lattice as zeros -- this may be a list of lists forming a 2D array
        # will get these from elsewhere
        self.lattice = np.array([[0 for _ in range(self.lattice_width)] for _ in range(self.lattice_height)])
        # make the image we will draw on
        self.rectangles_image = tk.PhotoImage(master=self, width=(self.lattice_width*self.virtual_spacing)+self.hardpad,\
                                                    height=(self.lattice_height*self.virtual_spacing)+self.hardpad)
        # make the image all white, initially
        self.rectangles_image.put('white', to=(self.hardpad,self.hardpad,\
                                 int((self.lattice_width*self.virtual_spacing)+self.hardpad),\
                                 int((self.lattice_height*self.virtual_spacing)+self.hardpad)))
        # add the image to the canvas
        self.canvas.create_image(self.hardpad, self.hardpad, anchor="nw", image=self.rectangles_image)

        # make a grid of guiding lines, blocking out each 5x5 square
        #@TODO: Make this able to be turned on and off with a button
        self.fine_gridlines()
        self.guiding_gridlines()

        self.canvas.bind("<B1-Motion>", self.paint_action) # left click, dragged

        # bindings for selection rectangle, activated holding ctrl
        self.rectangle_active = False # if rectangle is being drawn, stop the pen!     
        self.canvas.bind("<Control-ButtonPress-1>", self.set_first_rectangle_coordinate)
        self.canvas.bind("<Control-ButtonRelease-1>", self.set_second_rectangle_coordinate_and_fill)
        # this is going to draw a rectangle to help guide the eyes during selection/drawing
        self.guiding_rectangle = None
        # spray can
        self.canvas.bind("<Shift-B1-Motion>", self.spraycan_on)
        # action completed
        self.canvas.bind("<ButtonRelease-1>", self.complete_action)

    @property
    def lattice_height(self):
        return self.master.master.master.master.project.lattice_height

    @property
    def lattice_width(self):
        return self.master.master.master.master.project.lattice_width

    @property
    def lattice_major_gridlines_index(self):
        return self.master.master.master.master.project.lattice_major_gridlines

    @property
    def spray_index(self):
        return self.master.master.master.master.toolbox.spray_width_scale.get()

    @property
    def active_tool(self):
        return self.master.master.master.tool

    def refresh_image(self):
        '''
        This function will redraw the layer image when called
        '''
        self.virtual_spacing = self.master.master.master.zoomvalues[self.master.master.master.zoom]
        self.canvas["width"] = (self.lattice_width*self.virtual_spacing)+self.hardpad
        self.canvas["height"] = (self.lattice_height*self.virtual_spacing)+self.hardpad
        self.rectangles_image["width"] = (self.lattice_width*self.virtual_spacing)+self.hardpad
        self.rectangles_image["height"] = (self.lattice_height*self.virtual_spacing)+self.hardpad
        self.canvas.delete('grid_line')
        self.rectangles_image.blank()
        self.rectangles_image.put('white', to=(self.hardpad,self.hardpad,\
                                 int((self.lattice_width*self.virtual_spacing)+self.hardpad),\
                                 int((self.lattice_height*self.virtual_spacing)+self.hardpad)))
        
        for row in range(self.master.zlayer.lattice.shape[0]):
            for col in range(self.master.zlayer.lattice.shape[1]):
                color = "#FFFFFF"
                molid = self.master.zlayer.lattice[row][col]
                if molid == -1:
                    # Obstructed cell
                    color = "#000000"
                else:
                    for mol in self.master.master.master.master.project.molecules:
                        if mol.index == molid:
                            color = mol.color
                            break
                self.rectangles_image.put(color, to=(int((col*self.virtual_spacing)),\
                      int((row*self.virtual_spacing)),\
                      int(((col+1)*self.virtual_spacing)),\
                      int(((row+1)*self.virtual_spacing))))
        self.fine_gridlines()
        self.guiding_gridlines()
        # cursor update
        self.set_cursor()

    def fine_gridlines(self):
        for i in range(0, self.lattice_width):
            # self.rectangles_image.put('black', to=((i*self.virtual_spacing), 0,\
            #                             (i*self.virtual_spacing) + 1, (self.lattice_height*self.virtual_spacing)))
            self.canvas.create_line([((i*self.virtual_spacing)+self.hardpad, 0),\
                ((i*self.virtual_spacing)+self.hardpad, (self.lattice_height*self.virtual_spacing)+self.hardpad)], tag='grid_line', width=1)
        for i in range(0, self.lattice_height):
            # self.rectangles_image.put('black', to=(0, (i*self.virtual_spacing),\
            #                             (self.virtual_spacing*self.lattice_width), (i*self.virtual_spacing) + 1))
            self.canvas.create_line([(0, (i*self.virtual_spacing)+self.hardpad),\
                ((self.virtual_spacing*self.lattice_width)+self.hardpad, (i*self.virtual_spacing)+self.hardpad)], tag='grid_line', width=1)
        return

    def guiding_gridlines(self):
        for i in range(0, self.lattice_width, self.lattice_major_gridlines_index):
            self.canvas.create_line([((i*self.virtual_spacing)+self.hardpad, 0),\
                ((i*self.virtual_spacing)+self.hardpad, (self.lattice_height*self.virtual_spacing)+self.hardpad)], tag='grid_line', width=2)
        for i in range(0, self.lattice_height, self.lattice_major_gridlines_index):
            self.canvas.create_line([(0, (i*self.virtual_spacing)+self.hardpad),\
                ((self.virtual_spacing*self.lattice_width)+self.hardpad, (i*self.virtual_spacing)+self.hardpad)], tag='grid_line', width=2)
        return

    def set_molecule(self, row, col):
        """
        Paint a single molecule onto this layer
        """
        if self.master.zlayer.lattice[row][col] == -1:
            # This is an obstructed site and can't be painted over
            return
        id, color = self.master.master.master.molecule.get_paint_props()
        if self.master.zlayer.lattice[row][col] != id:
            self.master.zlayer.lattice[row][col] = id
            self.rectangles_image.put(color, to=(int((col*self.virtual_spacing)),\
                                                int((row*self.virtual_spacing)),\
                                                int(((col+1)*self.virtual_spacing)),\
                                                int(((row+1)*self.virtual_spacing))))

    def count_molecules(self):
        number_of_molecules_in_layer = len(np.where(self.master.zlayer.lattice > 0)[0])
        molecules = self.master.master.master.master.mlist.molecules
        for i in range(1,len(molecules)):
            number_of_molecule = len(np.where(self.master.zlayer.lattice == molecules[i]["molecule"].index)[0])
            molecules[i]["number_tkvar"].set(str(number_of_molecule))
            if number_of_molecules_in_layer > 0:
                percent_of_molecule = str(np.around((number_of_molecule / number_of_molecules_in_layer)*100, decimals=1))+"%"
                molecules[i]["percent_tkvar"].set(percent_of_molecule)
            else:
                molecules[i]["percent_tkvar"].set("0.0%")
        return

    def set_cursor(self):
        """
        Make the cursor show the active tool
        """
        if self.active_tool == 'pencil':
            self.canvas.config(cursor="pencil")
        elif self.active_tool == 'rect':
            self.canvas.config(cursor="cross")
        elif self.active_tool == 'spray':
            self.canvas.config(cursor="spraycan")

    def paint_action(self, event):
        """
        Handle clicks on the canvas by painting with the currently selected tool.
        """
        if self.rectangle_active:
            self.draw_guiding_rectangle(event=event)
            return
        if self.active_tool == 'pencil':
            self.canvas.config(cursor="pencil")
            self.pencil(event=event)
        elif self.active_tool == 'rect':
            self.canvas.config(cursor="cross")
            if not self.rectangle_active:
                self.set_first_rectangle_coordinate(event=event)
        elif self.active_tool == 'spray':
            self.canvas.config(cursor="spraycan")
            self.spraycan_on(event=event)

    def spraycan_on(self, event):
        self.canvas.config(cursor="spraycan")
        spray_width = self.virtual_spacing*self.spray_index

        # within the square inscribed by minx, maxs...
        min_spray_x, max_spray_x = (event.x-self.hardpad)-spray_width, (event.x-self.hardpad)+spray_width
        min_spray_y, max_spray_y = (event.y-self.hardpad)-spray_width, (event.y-self.hardpad)+spray_width
        # spray paint all lattice points whose centers fall inside of spray_width from the click point
        rows = range(int(min_spray_y//self.virtual_spacing),int(max_spray_y//self.virtual_spacing)+1)
        cols = range(int(min_spray_x//self.virtual_spacing),int(max_spray_x//self.virtual_spacing)+1)

        # first compute and store polar radii (in units of lattice spacing) for points within the spray
        polar_radii = np.zeros((len(rows),len(cols)))
        for i in range(len(rows)):
            for j in range(len(cols)):
                polar_radii[i][j] = ( ((rows[i] + 0.5) - ((event.y-self.hardpad)/self.virtual_spacing))**2 + ((cols[j] + 0.5) - ((event.x-self.hardpad)/self.virtual_spacing))**2)

        # now paint the points from the spray that fall within the spray index
        for i in range(len(rows)):
            for j in range(len(cols)):
                row = rows[i]
                col = cols[j]
                if polar_radii[i][j] < self.spray_index**2:
                    # wrap around the PBC!
                    row = row % self.lattice_height
                    col = col % self.lattice_width
                    self.set_molecule(row, col)
        self.count_molecules()

    def complete_action(self, event):
        """
        Handle the release so the rectangle can complete
        """
        if self.rectangle_active:
            self.set_second_rectangle_coordinate_and_fill(event=event)
        self.set_cursor()

    def draw_guiding_rectangle(self, event):
        self.canvas.coords(self.guiding_rectangle, self.first_rectangle_x, self.first_rectangle_y, event.x, event.y)

    def set_first_rectangle_coordinate(self, event):
        self.rectangle_active = True # this will perclude the pen from activating at the same time
        self.canvas.config(cursor="cross")
        self.first_rectangle_x = event.x
        self.first_rectangle_y = event.y
        self.guiding_rectangle = self.canvas.create_rectangle(self.first_rectangle_x, self.first_rectangle_y,\
                         self.first_rectangle_x, self.first_rectangle_y, fill=None, outline='#00ff00', width=2)


    def set_second_rectangle_coordinate_and_fill(self, event):
        self.rectangle_active = False # this will allow the pen to be activated again
        self.second_rectangle_x = event.x
        self.second_rectangle_y = event.y
        self.canvas.delete(self.guiding_rectangle) # delete the guiding rectangle
        self.fill_selected_rectangle()
        self.set_cursor()

    def fill_selected_rectangle(self):
        # higher indices are toward the bottom right, lower indices are toward the top left
        col_first = int((self.first_rectangle_x-self.hardpad)//self.virtual_spacing)
        row_first = int((self.first_rectangle_y-self.hardpad)//self.virtual_spacing)
        col_second = int((self.second_rectangle_x-self.hardpad)//self.virtual_spacing)
        row_second = int((self.second_rectangle_y-self.hardpad)//self.virtual_spacing)

        if (col_first > col_second) and (row_first > row_second):
            rows = range(row_second, row_first+1)
            cols = range(col_second, col_first+1)
        elif (col_second > col_first) and (row_second > row_first):
            rows = range(row_first, row_second+1)
            cols = range(col_first, col_second+1)
        elif (col_second > col_first) and (row_first > row_second):
            rows = range(row_second, row_first+1)
            cols = range(col_first, col_second+1)
        elif (col_first > col_second) and (row_second > row_first):
            rows = range(row_first, row_second+1)
            cols = range(col_second, col_first+1)
        elif (col_first > col_second) and (row_first == row_second):
            cols = range(col_second, col_first+1)
            rows = [row_first]
        elif (col_second > col_first) and (row_first == row_second):
            cols = range(col_first, col_second+1)
            rows = [row_first]
        elif (col_first == col_second) and (row_first > row_second):
            cols = [col_first]
            rows = range(row_second, row_first+1)
        elif (col_first == col_second) and (row_second > row_first):
            cols = [col_first]
            rows = range(row_first, row_second+1)
        else: return # no boxes have been encapsulated, do nothing
        # let the selection wrap around the PBC
        rows = [row % self.lattice_height for row in rows]
        cols = [col % self.lattice_width for col in cols]
        for row in rows:
            for col in cols:
                self.set_molecule(row, col)
        self.count_molecules()

    def pencil(self, event):
        # Calculate column and row number
        col = int((event.x-self.hardpad)//self.virtual_spacing)
        row = int((event.y-self.hardpad)//self.virtual_spacing)
        # avoid attempting to under- or over-paint
        if (col < 0) or (row < 0): return
        if (col > self.lattice_width-1) or (row > self.lattice_height-1): return

        # If the tile is not filled, create a rectangle
        self.set_molecule(row, col)
        self.count_molecules()
