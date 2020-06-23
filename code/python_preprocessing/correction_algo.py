import rasterio as rio, geopandas as gpd, numpy as np, os, fiona as fn

def set_zoom_scale(meta1:dict, meta2:dict) -> int:
    '''
    This function is used to the zoom scale by comparing the x-coordinate in the transforms of two sets of metadata. No order
    required as it automatically compares the width and height of both meta
    
    PARAMETERS
    ----------
    
        meta1: a dictionary for the first metadata (Required)
        
        meta2: a dictionary for the second metadata
        
    OUTPUT
    ------
    
        (x_zoom, y_zoom): a tuple giving the x and y zoom factors
        
    EXAMPLES
    --------
    
    zoom = set_zoom_scale({'driver': 'GTiff', 'dtype': 'float64', 'nodata': 0.0001, 'width': 5791, 'height': 2093, 'count': 1, 'crs': CRS.from_epsg(4326), 'transform': Affine(4.491576420597607e-05, 0.0, 106.71250031345062,
       0.0, -4.491576420597607e-05, -6.089141142605335), 'compress': 'lzw'}, \
       {'driver': 'GTiff', 'dtype': 'float32', 'nodata': 0.0001, 'width': 964, 'height': 348, 'count': 1, 'crs': CRS.from_epsg(4326), 'transform': Affine(0.00026994924954108627, 0.0, 106.5,
       0.0, -0.00026994924954108627, -6.0), 'compress': 'lzw'})
    
    '''
    print("META 1:", meta1)
    print("META 2:", meta2)
    # assert that the coordinate reference system is the same, or else it doesn't make sense to work out a zooms cale
    assert meta1['crs'] == meta2['crs']
    
    # TODO: probably add in a function in the future to convert the crs and make them the same. 
    
    # zoom_scale is a lambda function that takes in 3 arguments to get the floor of either the first or second index
    zoom_scale = lambda trans1, trans2, idx: int(trans1[idx] // trans2[idx])
    
    x_idx, y_idx = 0, 4
    
    x_zoom = zoom_scale(meta1['transform'], meta2['transform'], x_idx) if \
                meta1['height'] < meta2['height'] else zoom_scale(meta2['transform'], meta1['transform'], x_idx)
    
    y_zoom = zoom_scale(meta1['transform'], meta2['transform'], y_idx) if \
                meta1['width'] < meta2['width'] else zoom_scale(meta2['transform'], meta1['transform'], y_idx)
    
    return (x_zoom, y_zoom)


# 6*6 grid
def generate_grids(grid_size:int = 5) -> list:
    '''
    This function is used to generate a grid_size which corresponds from (-grid_size, grid_size) to (grid_size, grid_size)
    in single increments
    
    PARAMETERS
    ----------
    
        grid_size: an integer representing the grid size you wish to search in the correction algorithm. (Required)
        
    OUTPUT
    ------
    
        grid: a grid_size by grid_size list
        
    EXAMPLES
    --------
    
    grid = generate_grids()
    
    '''
    grid = [[(i,j) for j in range(-grid_size, grid_size)] for i in range(-grid_size, grid_size)]
    return grid


def generate_river_fill(data:np.array, zoom:tuple, meta:dict, grids, error_thres:float = 4e-1, algo:str = 'river'):
    '''
    This function is used to generate data points using a zoom factor, the metadata, and the data array
    
    PARAMETERS
    ----------
    
        data: a NumPy array with your data
        
        zoom: a tuple containing the x and y zoom factors in the form (x_zoom, y_zoom) (Required)

        meta: a dictionary containing metadata corresponding to the data array (Required)
        
        grids: a GridData class object (Required)

        error_thres: an error threshold for correcting the grid between 0 and 1. 
        If you give a higher value, less cells will be corrected and some may even be deleted because the grid count 
        that are filled do not meet this error_thres. If you give a lower value, more cells will be corrected. If no value given,
        defaults to 0.4 or 4e-1. (Optional)
        
        algo: whether to use the river correction or building correction algo. If you use the river correction algo,
        you will end up taking the fill to be 1 if at least error_thres % of the grids have a fill value of 1. Or else,
        you will take that the n * n percent of cells as the fill value. (Optional)
        
    OUTPUT
    ------
    
        geos: a list of shapely.ops.Points that contain the Point coordinate for each feature
        
        points: a list of points and containing fill values corresponding to an entry in geos
        
    EXAMPLES
    --------
    
    src_file = '../data/experiments/jakarta/fixed_raster/dem.tif'
    building_data, building_meta = open_rio(src_file)
    grids = GridData(data = hi_res_data.data, zoom = (6,6))
    generate_river_fill(data = building_data.data, zoom = (6, 6), grids = grids, meta = building_meta, error_thres = 2e-1)
    
    '''
    
    zoom_scale = zoom[0] if zoom[0] > zoom[1] else zoom[1]
    
    geos, points = list(), list()

    # get the nodata value
    nd = meta['nodata']
    map_trans = meta['transform']
    
    print("---- PROCESSING ----")
    print("NODATA", nd)
    print("DATA TRANSFORM", map_trans)
    
    no_grid_cells = zoom[0] * zoom[1]
    
    traverse_grid = generate_grids()
    

    def fill_algo(counter:int, item):
        '''
        Fill algo is an inner function of generate_river_fill. Created so that it can take the variables instantiated
        in generate_river_fill without having to redeclare any of the variables.
        
        PARAMETERS
        ----------

            counter: a numerical counter (Required)
            
            item: a data entry (Required)

        OUTPUT
        ------

            no_filled_grids: an integer number for the number of filled grids
        
        '''
        
        if algo == "river":
            # account for fill
            fill = 1 if grids(counter) / no_grid_cells > error_thres or item == 1 else 0 

            if fill == 0 and item != nd:

                # if a n * n grid around it has no filled values, then declare it as nothing 

                tmp_row, tmp_col = counter, edx
                tmp_filled_grids = 0 

                for ent1 in traverse_grid:

                    for ent2 in ent1:
                        tmp_col += ent2[1] 

                        if data[tmp_row][tmp_col] != nd:
                            tmp_filled_grids += 1

                    tmp_row += ent2[0]

                if tmp_filled_grids < (len(traverse_grid)*len(traverse_grid[0])):
                    fill = 0 
                else:
                    fill = 1
            
        else:
            
            fill = grids(counter) / no_grid_cells 
                
        return fill
    
    
    # set the required counter and fill variables, then count the number of grid cells based on the zoom
    counter, fill = 0, 0
    
    # need to go by using a 6*6 grid since dem data is 5x5. So as I iterate through the river grid, the corresponding map grid is a 
    # 6x6 grid where I calculate the percentage of cell that has river filling. 
    for idx, entry in tqdm(enumerate(data)):
        
        # get the top left x and y coordinates on every iteration
        top_left_x, top_left_y = map_trans[2], map_trans[5]+(counter * map_trans[4])
        temp_points, temp_geos = list(), list()

        for edx, item in enumerate(entry):
            
            fill = fill_algo(counter, item)

            temp_points.append(fill)
            temp_geos.append(Point([top_left_x, top_left_y]))

            top_left_x += map_trans[0] 

        points.append(temp_points)
        geos.append(temp_geos)
        counter += 1    
        grids.reset_col_counter()
        
        
    print("TOTAL COUNTED GRIDS:", grids.total_counted_grids)

    return points, geos


def get_percentage_corrected(original_points:list, corrected_points:list):
    def corrected_pts(points):
        total_corrected = 0 
        for item in points:
            data_arr = np.array(item) if type(item) != np.array else item
            total_corrected += len(data_arr[data_arr>0])
        return total_corrected
    
    original_num = corrected_pts(original_points)
    corrected_num = corrected_pts(corrected_points)
    perc_corrected = (corrected_num - original_num) / original_num
    
    print("ORIGINAL NUMBER OF POINTS:", original_num)
    print("CORRECTED NUMBER OF POINTS:", corrected_num)
    print("PERCENTAGE OF POINTS CORRECTED:", round(perc_corrected*100, 3), "%")
    
