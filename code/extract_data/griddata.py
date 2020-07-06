from tqdm import tqdm
import numpy as np

class GridData:
    '''
    GridData is a class that implements the counting and iteration of grids based on the zoom scale, and map iteration counter
    '''
    def __init__(self, data:np.array, zoom:tuple, nd):
        '''
        Used to initialise the GridData class structure with the data array, and zoom_scale

        PARAMETERS
        ----------

            data_array: a NumPy array containing your data (Required)

            zoom: a tuple containing the x and y zoom factors in the form (x_zoom, y_zoom) (Required)

            nd: a nodata value (Required)

        EXAMPLES
        --------

        grids = GridData(data = np.array([1,2,3,4,5]), zoom_scale = 6)

        '''
        self.data = data
        self.x_zoom, self.y_zoom = zoom
        self.col_counter = 0
        self.total_counted_grids = 0
        self.nd = nd
        
        print("DATA SIZE", data.shape, "NODATA", self.nd)
        
    def __call__(self, counter:int):
        '''
        Function to traverse grids of the data based on the counter and count the number of filled grids within the search grid

        PARAMETERS
        ----------

            counter: a numerical counter to be used in offsetting the computed row index (Required)
            
        OUTPUT
        ------
        
            no_filled_grids: an integer representing the number of filled grids based on the zoom scale
            
        '''
        
        # count the next zoom_scale * zoom_scale grids
        no_filled_grids = 0

        # need to account for grid increments in row and column
        for row in range(self.x_zoom):
            
            for col in range(self.y_zoom):
            
                row_index = (counter * self.x_zoom) + row
                col_index = (self.col_counter * self.y_zoom) + col
                
                self.total_counted_grids += 1
                    
                if self.data[row_index][col_index] != self.nd:
                    no_filled_grids += 1
                    
        self.col_counter += 1

        return no_filled_grids

    def reset_col_counter(self):
        '''
        Function to reset the class variable col_counter that keeps track of how many columns have been traversed
        
        '''
        self.col_counter = 0
    
