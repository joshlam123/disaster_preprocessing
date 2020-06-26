import os, shutil, logging, sys
import earthpy.plot as ep, rasterio as rio, matplotlib.pyplot as plt, rioxarray as rx, numpy as np, geopandas as gpd, fiona as fn
from rasterio.plot import show, plotting_extent
from subprocess import Popen, run
from affine import Affine
from pyproj import Proj, transform
from rasterio import features
from osgeo import gdal, osr
from rasterio.crs import CRS
from distutils.version import LooseVersion
from pyproj.enums import WktVersion
from osgeo import osr
from tqdm import tqdm
from shapely.geometry import Point
from disaster_preprocessing.code.extract_data.griddata import GridData

def change_to_mercator(new_crs:int = 900913):
    '''
    This function is similar to gconvert from one CRS into a rasterio approved one. For more information about comptability,
    check # https://pyproj4.github.io/pyproj/stable/crs_compatibility.html

    Parameters
    ----------

        new_crs: a valid epsg number. check the link for more information between wkt to epsg string conversions

    Output
    ------

        rio_crs: a rasterio crs object
    
    Examples
    --------

    crs = change_to_mercator(4326)

    '''
    proj_crs = CRS.from_epsg(new_crs)
    if LooseVersion(rio.__gdal_version__) < LooseVersion("3.0.0"):
        rio_crs = proj_crs.to_wkt(WktVersion.WKT1_GDAL)
    else:
        rio_crs = proj_crs.to_wkt()
    return rio_crs

def perform_addition(data:np.array, dem_src:str) -> np.array:
    '''
    This function is similar to perform_restoration except it uses 0's as the column or row to be added.

    To correct the data values for a given array of map data by using the Digital Elevation Model (DEM) as a reference. 

    The function will read the DEM file and append or delete additional rows or columns depending on the size of your current data. 
    This function will add 0's as data values to the file.

    Parameters
    ----------

        data: a NumPy ndarray object.

        dem_src: a string containing the file and its file path to the DEM in the format: '../data/dem.tif'

    Output
    ------

        data: a NumPy ndarray object
    
    Examples
    --------

    data = perform_addition(np.array([1,2,3,4,5], [1,2,3,4,5]), '../data/dem/dem.tif')

    '''

    print("DEM SOURCE FILE: ", dem_src)
    with rio.open(dem_src) as src:
        row, cols = src.shape
        print("DEM DIMENSIONS: ", row, cols)

    data_row, data_cols = data.shape
    print("SOURCE FILE DIMENSIONS: ", data_row, data_cols)
    
    zero = np.array([0 for i in range(data.shape[1])])
    
    while data_row != row:

        if data_row > row:
            data = np.delete(data, -1, axis=0)
        else:
            zero = zero.reshape(1, -1)
            data = np.concatenate((data,zero), axis=0)

        data_row, data_cols = data.shape
        print("ROW MODIFICATION OCCURED FOR FILE. DIMENSIONS AFTER ROW MODIFICATIONS: ", data_row, data_cols)

    zero = np.array([0 for i in range(data.shape[0])])

    while data_cols != cols:
        if data_cols > cols:
            data = np.delete(data, -1, axis=1)

        else:
            zero = zero.reshape(-1,1)
            data = np.concatenate((data,zero), axis=1)
        data_row, data_cols = data.shape
        print("COLUMN MODIFICATION OCCURED FOR FILE. DIMENSIONS AFTER COLUMN MODIFICATIONS: ", data_row, data_cols)

    return data

def move_files(files:str, src:str, dest:str):
    '''
    A function to move files from one directory to another.

    Parameters
    ----------

        files: A list of files you want to move (Required)

        src: the source directory you want to move your files from (Required)

        dest: the name of the destination directory you want to move your files to (Required)

    Output
    ------

        None
    
    Examples
    --------

    move_files('../data/experiments/jakarta/exp5/', '../code/' 'landunits.tif', 'landunits_cubic.tif')

    '''

    if isinstance(files, list):
        for item in files:
            shutil.move(src+item, dest+item)
            shutil.move(src+item+'.aux.xml', dest+item+'.aux.xml')
    else:
        shutil.move(src+files, dest+files)
        
    print(f"successfully moved {files} from {src} to {dest}")
    return None

def perform_restoration(data:np.array, dem_src:str) -> np.array:
    '''
    This function is similar to perform_addition except it uses the last row or column in your current data as the row or column to be added.

    To correct the data values for a given array of map data by using the Digital Elevation Model (DEM) as a reference. 

    The function will read the DEM file and append or delete additional rows or columns depending on the size of your current data. 
    This function will add 0's as data values to the file.

    Parameters
    ----------

        data: a NumPy ndarray object. (Required)

        dem_src: a string containing the file and its file path to the DEM in the format: '../data/dem.tif' (Required)

    Output
    ------

        data: a NumPy ndarray object
    
    Examples
    --------

    data = perform_restoration(np.array([1,2,3,4,5], [1,2,3,4,5]), '../data/dem/dem.tif')

    '''

    dem_xds = rx.open_rasterio(
            dem_src,
            parse_coordinates=True,
            chunks=True,
            masked=False,
        )

    row, cols = dem_xds.rio.height, dem_xds.rio.width
    data_row, data_cols = data.shape
    
    last_row_val = data[-1]
    
    while data_row != row:
        # hack
        if data_row > row:
            data = np.delete(data, -1, 0)

        else:
            last_row_val = last_row_val.reshape(1, -1)
            data = np.concatenate((data,last_row_val), axis=0)
        data_row, data_cols = data.shape
        

        print("Number of Rows and Columns: ", data_row, data_cols)
        
    last_col_val = np.array([10e-5 for i in data])
        
    while data_cols != cols:
        if data_cols > cols:
            data = np.delete(data, -1, 1)

        else:
            last_col_val = last_col_val.reshape(-1,1)
            data = np.concatenate((data,last_col_val), axis=1)
        data_row, data_cols = data.shape
        
    print("Number of Rows and Columns: ", data_row, data_cols)
    return data


def create_folders(dest:str, folders:list = list()):
    '''
    To create the required folder names in the file directory: river, reprojection, fixed_raster, base_maps, and maps.

    Parameters
    ----------

        dest: a string containing the main working folder name ending with a slash /: '../data/' (Required)

    Output
    ------

        None
    
    Examples
    --------

    data = create_folders("../data/")

    '''    

    if len(folders) == 0:
        folders = ['river', 'reprojection', 'fixed_raster', 'base_maps', 'maps', 'building']

    created = list()
    for item in folders:
        if not os.path.isdir(dest+item):
            os.mkdir(dest+item)
            created.append(item)
            
    print(f"Folders {created} created")


def get_folder_files(incl:str, excl:str='', ext:str='.tif', dest:str='../data/landcover/epsg4326/generated_base/'):
    '''
    Get all files within the desired destination folder

    Parameters
    ----------
        
        incl: any characters in files you want to include in the output (Required)

        excl: any characters in files you want to exclude in the output (Optional)

        ext: the desired file extension you want to search for (Optional)

        dest: a string containing the main working folder name ending with a slash /: '../data/' (Optional)

    Output
    ------

        folder: the directory name you entered in as dest

        files: a list of file names matching incl, excl, and ext.
    
    Examples
    --------

    folder, files = get_folder_files("dem", "aux", ".map", "../data/")

    '''    

    folder = [i for i in os.walk(dest)][0]
    files = [j for j in folder[-1] if ext in j and incl in j]

    if len(excl) > 0:
        files = [j for j in files if excl not in j]


    print(folder[0], files)
    
    if ext == '.tif':
        for idx, item in enumerate(files):
            dataset = rio.open(folder[0] + "/" + item)
            print("Data Name {}, Data Rows: {}, Data Columns: {}".format(item, dataset.height, dataset.width))
    
    return folder[0], files


def open_and_plot_datasets(folder:str, files:list):
    '''
    Open a set of files and plot them

    Parameters
    ----------
        
        folder: a string containing the main working folder name ending with a slash /: '../data/' (Required)

        files: a list of file names you want to plot (Required)


    Output
    ------

        None. Plot is inline if you're running a notebook. 
    
    Examples
    --------

    open_and_plot_datasets("../data/experiments/jakarta/fixed_raster/", ['thetas1.tif', 
    'rr.tif', 'landunits.tif', 'chanmask.tif', 'ch.tif', 'n.tif', 'per.tif', 'lai.tif', 'ksat1.tif', 
    'psi1.tif', 'dem.tif', 'thetai1.tif'])

    '''    

    for idx, item in enumerate(files):
        try:
            print("Dataset plotted: ", item)
            dataset = rio.open(folder + "/" + item) # CHANGE THIS TO YOUR OWN FOLDER
        
            spatial_extent = rio.plot.plotting_extent(dataset)

            fig, ax = plt.subplots(figsize = (10,8))

            lidar_dem = dataset.read(1, masked=True)
            
            title = item[:item.find('.tif')]
            ep.plot_bands(lidar_dem, 
                 title=f"{title} Map for Indonesia", 
                 cmap='twilight_shifted',
                extent=spatial_extent,
                 ax=ax)

            ax.set_axis_off()

            plt.show()

        except ValueError:
            print("Cmax and Cmin value error for dataset: ", item)

def current_dir_files(dest:str) -> list:
    '''
    Get the files in the current directory

    Parameters
    ----------
        
        dest: a string containing the folder name /: '../data/' (Required)

    Output
    ------

        A list of files within the folder
    
    Examples
    --------

    folder, files = current_dir_files("../data/experiments/jakarta/exp11/fixed_raster/")

    ''' 

    dest = os.getcwd() or dest
    return [i for i in os.walk(dest)]

def set_current_folder() -> str:

    '''
    Returns the current folder. Used to actually set folder and have less imports in your environment

    Parameters
    ----------
        
        None


    Output
    ------

        Current working directory as given by os.getcwd()
    
    Examples
    --------

    curr_folder = set_current_folder()

    ''' 

    return os.getcwd()


def chg_to_folder(dest:str):
    '''
    changes your current working directory to the destination you provide

    Parameters
    ----------
        
        dest: a string containing the folder name /: '../data/' (Required)


    Output
    ------

       None
    
    Examples
    --------

    chg_to_folder('../data/experiments/jakarta/')

    ''' 

    os.chdir(dest)
    print(os.getcwd())

# function to rearrange the ordering of files so that DEM comes first and can be processed. 
def rearrange_dem(files:list):
    '''
    changes your current working directory to the destination you provide

    Parameters
    ----------
        
        files: a list of file names that ends in .tif 


    Output
    ------

       A list containing the file names with the DEM file on top
    
    Examples
    --------

    rearrange_dem(['thetas1.tif', 
    'rr.tif', 'landunits.tif', 'chanmask.tif', 'ch.tif', 'n.tif', 'per.tif', 'lai.tif', 'ksat1.tif', 
    'psi1.tif', 'dem.tif', 'thetai1.tif'])

    ''' 

    dem_idx = files.index('dem.tif')
    result = [files[dem_idx]] + files[:dem_idx] + files[dem_idx:]
    return result


def geojson_to_coordinates(file:str):
    '''
    Converts a geojson file into a set of coordinates

    Parameters
    ----------
        
        files: a geojson file you want to read (Required)


    Output
    ------

        file_folder: a list of files within the path that are of geojson format
    
    Examples
    --------

    geojson_to_coordinates('../data/experiments/indonesia/jakarta.geojson')

    '''
    jkt_gdf = gpd.read_file(file)
    all_coords = list(zip(*jkt_gdf.geometry[0].exterior.coords.xy))
    return all_coords

def extract_geojson_files_from_path(path, excl='jakarta'):
    '''
    Extract all files with the geojson extension from the path

    Parameters
    ----------
        
        path: a path you want to search. (Required)

        excl: a file type you want to exclude. By default the value is 'jakarta' (Optional)


    Output
    ------

        file_folder: a list of files within the path that are of geojson format
    
    Examples
    --------

    extract_geojson_files_from_path('../data/experiments/indonesia/', 'semarang')

    '''

    folders = [i for i in os.walk(path)][0]
    folder = [folders[0]+i for i in folders[1]]
    files = [list(os.walk(i)) for i in folder]
    file_folder = [[[j[0] + "/" + k for k in j[-1] if 'geojson' in k] for j in i if excl not in j[0]][0] for i in files ]
    return file_folder

def update_meta(meta:dict, p2:str, **kwargs):
    '''
    A function to update the metadata 

    Parameters
    ----------

        meta: a dictionary containing the metadata of a file (Required)

        p2: a string containing the coordinate reference system you want to update to. It must be a valid Rasterio CRS. (Required)

        **kwargs: additional items you want to update your metadata with. In the form of a kwarg dictionary (Optional)

    Output
    ------

        meta: a dictionary containing the metadata 
    
    Examples
    --------
    
    folder, item = "../data/experiments/base_maps/", "rr.tif"
    _, rio_crs, *unused = get_meta(folder, item)
    update_meta(meta = meta,  p2 = rio_crs, compress = "lzw", nodata = 0.0001, height = 348, width = 966)
    

    '''

    kwarg = kwargs
    print("EXTRA UPDATE ARGUMENTS", kwarg)
    if len(kwarg) > 0:
        for key, val in kwarg.items():
            if key == "dtype":
                meta.update(dtype = val)
            elif key == "driver":
                meta.update(driver = val)
            elif key == "nodata":
                meta.update(nodata = val)
            elif key == "width":
                meta.update(width = val)
            elif key == "height":
                meta.update(height = val)
            elif key == "crs":
                meta.update(crs = val)
            elif key == "transform":
                meta.update(transform = val)
            elif key == "compress":
                meta.update(compress = val)
    else:
        meta.update(dtype = rio.float32) # sample code for how you convert metadata
        meta.update(nodata = 10e-5)
    
    meta.update(crs = p2)
        
    return meta

def save_gpkg_to_layer(data_src:str, working_folder:str):
    '''
    This function is used to save a .gpkg data_source file to a layer. 
    More information here: https://stackoverflow.com/questions/56165069/can-geopandas-get-a-geopackages-or-other-vector-file-all-layers
    
    Parameters
    ----------

        data_src: a .gpkg data source string (Required)

        working_folder: a working folder directory in string format (Required)


    Output
    ------

        None
    
    '''
    layers, size = list(), list()
    for layername in fn.listlayers(data_src):
        print("LAYER NAME: ", layername)
        with fn.open(data_src, layer=layername) as src:
            print(layername, len(src))
            layers.append(layername)
            size.append(len(src))
            
    chosen_layer = layers[size.index(max(size))]
    print("MAX SIZE CHOSEN LAYER NAME", chosen_layer)
    
    folder = working_folder+'buildings/building.shp'
    
    gpkg_data = gpd.read_file(data_src, layer = chosen_layer)
    gpkg_data.to_file(folder, driver='ESRI Shapefile')
    
    print("GPKG DATA SUCCESSFULLY SAVED TO:", folder)

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

def generate_river_fill(data:np.array, zoom:tuple, meta:dict, grids:GridData, error_thres:float = 4e-1, algo:str = 'river'):
    '''
    This function is used to generate data points using a zoom factor, the metadata, and the data array
    
    PARAMETERS
    ----------
    
        data: a NumPy array with your data
        
        zoom: a tuple containing the x and y zoom factors in the form (x_zoom, y_zoom) (Required)

        meta: a dictionary containing metadata corresponding to the data array (Required)

        grids: a GridData class item (Required)
        
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
    generate_river_fill(data = building_data.data, zoom = (6, 6), meta = building_meta, error_thres = 2e-1)
    
    '''
    
    zoom_scale = zoom[0] if zoom[0] > zoom[1] else zoom[1]
    
    geos, points = list(), list()

    
    no_grid_cells = zoom[0] * zoom[1]

    # grids = GridData(data = hi_res_data.data, zoom = zoom)
    
    traverse_grid = generate_grids()
    
    # get the nodata value
    nd = meta['nodata']
    map_trans = meta['transform']
    
    print("---- PROCESSING ----")
    print("NODATA", nd)
    print("DATA TRANSFORM", map_trans)

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
    #                 print(tmp_col, tmp_row)

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




def generate_dataframe_from_points(points:list, geos:list):
    '''
    This is a function used to generate a dataframe from a set of data points and geometries.

    Parameters
    ----------

        points: a list of data points (Required)

        geos: a list of shapely.geometry.Point objects (Required)


    Output
    ------

        df: a dataframe containing all points and geometries

    Examples
    --------

    zoom = set_zoom_scale(hi_res_meta, building_meta)
    points, geos = generate_river_fill(building_data.data, zoom = zoom, meta = building_meta, error_thres = 2e-1, algo = "building")
    get_percentage_corrected(points, building_data)

    df = generate_dataframe_from_points(points = points, geos = geos)

    '''
    df = pd.DataFrame(columns=['geometry', 'data'])

    for idx, entry in tqdm(enumerate(points)):
        df_temp = pd.DataFrame({"geometry":geos[idx], "data":entry}, index = [i for i in range(len(entry))])
        frames = [df, df_temp]
        df = pd.concat(frames)

    df = df.astype({'data':'float32'})

    return df

def generate_rx_meta_data(src_file:str, rx_data, crs:str, write_data:np.array = np.array([]), nodata:float = 10e-5):
    '''
    This is a function to open a source file and update the metadata using rioxarray array data. 

    Parameters
    ----------

        src_file: the file you want to get the metadata from (Required)
        
        rx_data: a Rioxarray.xarray data (Required)

        crs: the crs to update the metadata to (Required)

        write_data: a NumPy array for your data (Optional)

        nodata: the no data replacement value. If no value, it is set to 10e-5 (Optional)


    Output
    ------

        None
    
    Examples
    --------

    folder, item = "../data/experiments/base_maps/", "rr.tif"

    xds = rx.open_rasterio(
                    folder+item,
                    parse_coordinates=True,
                    chunks=True,
                    masked=False,
                )
            
    _, rio_crs, nodata, *unused = get_meta(folder, item)        

    generate_rx_meta_data(src_file, rx_data = xds, write_data = xds[0].values, crs = rio_crs, nodata = nodata)


    '''
    map_data, meta = open_rio(src_file)

    meta.update(compress = 'lzw')
    meta.update(dtype = rio.float32)
    meta.update(crs = crs)

    if write_data.shape[0] == 0:
        meta.update(width = rx_data.rio.width)
        meta.update(height = rx_data.rio.height)
    else:
        meta.update(width = write_data.shape[0])
        meta.update(height = write_data.shape[1])

    meta.update(transform = rx_data.rio.transform())
    meta.update(dtype = 'float32')
    meta.update(nodata = nodata)

    return meta

def get_meta(folder:str, file:str):
    '''
    Checks whether a file exist in the path.

    Parameters
    ----------
        
        folder: a string containing the main working folder name ending with a slash /: '../data/' (Required)

        files: a file name you want to return the meta for (Required)


    Output
    ------

        xds: a Rioxarray xarray.core.dataarray.DataArray object

        xds.rio.crs: a rasterio.crs.CRS coordinate reference system object

        xds.rio.nodata: the no data value of the file
        
        xds.rio.width: the width (number of columns) of the file

        xds.rio.height: the height (number of rows) of the file

        xds.rio.transform(): an affine.Affine tuple that represents the geographical transformation in the format:
                                (x_increment, 0, x_starting_coordinate, 0, y_increment, y_starting_coordinate)
    
    Examples
    --------

    get_meta("../data/experiments/jakarta/reprojection", 'dem.tif')

    '''
    xds = rx.open_rasterio(
    folder + file,
    masked=False,
    parse_coordinates=True,
    chunks=True,
    )

    return xds, xds.rio.crs, xds.rio.nodata, xds.rio.bounds(), xds.rio.width, xds.rio.height, xds.rio.transform()


def open_rio(file:str):
    '''
    This function is used to open a raster file using rasterio and return the important information

    Parameters
    ----------

        file: a file name you want to open (with the backextension) (Required)

    Output
    ------

        map_data: mapdata in NumPy array format

        meta: a dictionary containing the metadata of the src_file
    
    Examples
    --------

    open_rio("../data/experiments/jakarta/reprojection/dem.tif")

    '''

    with rio.open(file) as src:
        # grab the meta data because we'll need it later
        meta = src.meta.copy()
        meta.update(compress='lzw')
        map_data = src.read(1, masked=True)

    print("---- OPENING SOURCE FILE ----")
    print("FILENAME, ", file)
    # View object dimensions
    print("DEM SHAPE:" , map_data.shape)

    # src = rio.open('output_alos.tif')
    print("METADATA:", src.meta)

    print("FILE CRS:",src.crs)
    print("FILE BOUNDS:",src.bounds)

    print("SIZE:", map_data.ravel().shape)
    
    return map_data, meta

def fix_nan_values(data:np.array, replace_val:float = 0):
    '''
    A function to convert all NaN values in a numpy array to a replacement value

    Parameters
    ----------

        data: a NumPy ndarray object. (Required)

        replace_val: the replacement value. If no value is provided, 0 is used. 


    Output
    ------

        clipped.values: a NumPy ndarray object.
    
    Examples
    --------

    fix_nan_values(np.array([1,2,4,np.nan,0]), 0)

    '''
    
    if len(data[np.isnan(data)]) > 0:
        clipped.values = np.where(np.isnan(data), replace_val, data)
        
    return clipped.values

def checkExist(file:str, path:str='default'):
    '''
    Checks whether a file exist in the path.

    Parameters
    ----------
        
        file: a file name (Required)

        path: a path you want to search. If not given, it just searches the current working directory using os.getcwd() (Optional)


    Output
    ------

       Boolean: True if the file exists, False if the file does not exist
    
    Examples
    --------

    checkExist('dem.tif', "../data/experiments/jakarta/reprojection")

    '''

    if path == "default": path = os.getcwd()
    else: path = path
    try:
        files = [i for i in os.walk(path)][0][-1]
        if file in files: return True
        else: return False      

    except IndexError as e:
        print("File does not exist")
        return True

# generate affine projection that creates long/lats
def affine_projection(file:str='../data/landcover/epsg4326/base_map/land_cover_jakarta.tif'):
    '''
    Create an affine projection for the dataset based on the source file's coordinate reference system.
    Taken from: https://gis.stackexchange.com/questions/129847/obtain-coordinates-and-corresponding-pixel-values-from-geotiff-using-python-gdal

    Parameters
    ----------
        
        file: a file name that ends in .tif (Required)


    Output
    ------

       longs: a NumPy ndarray of longitude coordinates

       lats: a NumPy ndarray of latitude coordinates

       p2: a string that represents the coordinate reference system. For example in the EPSG:900913 format:

       'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],
       AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]'
    
    Examples
    --------

    longs, lats, p2 = affine_projection('../data/landcover/epsg4326/lcm_Indonesia_Semarang.tif')

    '''

    # Read raster
    with rio.open(file) as r:
        T0 = r.transform  # upper-left pixel corner affine transform
        p1 = Proj(r.crs)
        A = r.read()  # pixel values

    # All rows and columns
    cols, rows = np.meshgrid(np.arange(A.shape[2]), np.arange(A.shape[1]))

    # Get affine transform for pixel centres
    T1 = T0 * Affine.translation(0.5, 0.5)
    # Function to convert pixel row/column index (from 0) to easting/northing at centre
    rc2en = lambda r, c: (c, r) * T1

    # All eastings and northings (there is probably a faster way to do this)
    eastings, northings = np.vectorize(rc2en, otypes=[np.float, np.float])(rows, cols)

    # Project all longitudes, latitudes
    p2 =  CRS.to_wkt(r.crs) #Proj(proj='latlong',datum='WGS84')
    longs, lats = transform(p1, p2, eastings, northings)
    
    return longs, lats, p2

def set_gdal_datatype_conv():
    '''
    This function is used to get the gdal datatype conversion data. Information from here:
    https://borealperspectives.org/2014/01/16/data-type-mapping-when-using-pythongdal-to-write-numpy-arrays-to-geotiff/

    Parameters
    ----------

        None

    Output
    ------

        The typemap and corresponding conversino number in GDAL
    
    Examples
    --------
    
    dtypes = set_gdal_datatype_conv()

    '''

    # typemap = {}
    # for name in dir(np):
    #     obj = getattr(np, name)
    #     if hasattr(obj, 'dtype'):
    #         try:
    #             npn = obj(0)
    #             nat = np.asscalar(npn)
    #             if gdal_array.NumericTypeCodeToGDALTypeCode(npn.dtype.type):
    #                 print(npn.dtype.name)
    #                 typemap[npn.dtype.name] = gdal_array.NumericTypeCodeToGDALTypeCode(npn.dtype.type)
    #         except:
    #             pass
    typemap = {
      "uint8": 1,
      "int8": 1,
      "uint16": 2,
      "int16": 3,
      "uint32": 4,
      "int32": 5,
      "float32": 6,
      "float64": 7,
      "complex64": 10,
      "complex128": 11,
    }
    return typemap

def write_array_to_raster_gdal(data:np.array, src_file:str, dest_file:str, title:str, meta:dict, target_crs:int = 4326, band:int = 1):
    '''
    This function is to raster a NumPy array into a destination file using GDAL. 

    Parameters
    ----------

        data: a geopandas GeoDataFrame containing target data (Required)

        src_file: a string containing the name of the source file: '../data/dem.tif' (Required)

        dst_file: a string containing the name of the destination file: '../data/dem.tif' (Required)
        
        meta: a dictionary containing the metadata for the source file (Required)

        target_crs: a crs to project to. If not provided, defaults to 4326 (Optional)

        band: the band to use for reading. If none provided defaults to 1 (Optional)    
    
    Output
    ------

        None
    
    Examples
    --------

    write_array_to_raster_gdal(data=np.array([1,2,3,4,5]), src_file="../data/experiments/jakarta/base_maps/rr.tif", 
    dest_file="../data/experiments/jakarta/fixed_raster/rr.tif")

    '''
    print("---- WRITING TO DESTINATION FILE ----")
    print("WRITING FILE FOR:", title)

    ds = gdal.Open(src_file)
    
    [cols, rows] = data.shape#meta['height'], meta['width']
    print("DATA DIMENSION:", cols, rows)
    print("DATA DIMENSION OF ARRAY TO WRITE:", data.shape)

    NP2GDAL_CONVERSION = set_gdal_datatype_conv()
    if data.dtype.name == "int64":
        gdaltype = 7
    else:
        gdaltype = NP2GDAL_CONVERSION[data.dtype.name]
    
    # setting up the destination file with options taken from: https://gis.stackexchange.com/questions/299059/osgeo-gdal-translate-how-to-set-compression-on-gdal-gtiff-driver
    # this uses LZW compression
    driver = gdal.GetDriverByName("GTiff")

    print("WRITING TO DESTINATION FILE:", dest_file)

    outdata = driver.Create(dest_file, rows, cols, 1, 7, options=['COMPRESS=LZW'])
    
    dest_transform = (meta['transform'][2], meta['transform'][0], meta['transform'][1], \
                     meta['transform'][5], meta['transform'][3], meta['transform'][4])
    

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(target_crs)

    outdata.SetProjection(srs.ExportToWkt())
    outdata.SetGeoTransform(dest_transform)
    outdata.GetRasterBand(1).WriteArray(data)

    if title == "dem": 
        nodata_repl = -10e38

    # TODO: Write wrapper functions to cover warnings like this. This is just for a quick fix
    if meta['nodata'] is not None:
        outdata.GetRasterBand(1).SetNoDataValue(meta['nodata'])
    else:
        logging.warning("SOURCE FILE NODATA VALUE IS NONE. This could potentially lead to future save rasters that result in problematic nodata conversions. If you think this is correct behavior, please ignore this warning.")

    outdata.FlushCache() ##saves to disk!!
    outdata = None
    band=None
    ds=None
    
    print(f"Wrote file {src_file} to {dest_file}")

    return None


def write_gdf_to_raster_gdal(data, selector:str, src_file:str, dest_file:str, meta:dict, target_crs:int = 4326, band:int = 1):
    '''
    This function is to raster a geopandas GeoDataFrame into a destination file using GDAL. 

    Parameters
    ----------

        data: a geopandas GeoDataFrame containing target data (Required)
        
        selector: a column name to select from the GeoDataFrame (Required)

        src_file: a string containing the name of the source file: '../data/dem.tif' (Required)

        dst_file: a string containing the name of the destination file: '../data/dem.tif' (Required)
        
        meta: a dictionary containing the metadata for the source file (Required)

        target_crs: a crs to project to. If not provided, defaults to 4326 (Optional)

        band: the band to use for reading. If none provided defaults to 1 (Optional)

    Output
    ------

        None
    
    Examples
    --------

    write_gdf_to_raster_gdal(data, 'dem', src_file = "../data/experiments/jakarta/base_maps/dem.tif", 
    dest_file = "../data/experiments/jakarta/fixed_raster/rr.tif", meta = meta)

    ''' 
    print("---- WRITING TO DESTINATION FILE ----")
    print("WRITING FILE FOR: ", selector)
    src_data = gdal.Open(src_file)
     
    # open the source file and get the dimensions
    # band = src_data.GetRasterBand(band)
    # arr = band.ReadAsArray()
    # [cols, rows] = arr.shape

    [cols, rows] = meta['height'], meta['width']
    print("DATA DIMENSION:", cols, rows)

    # reshape the geodataframe to match the correct resolution
    data_write_arr = np.array(data[selector]).reshape(cols, rows)


    # repurpose the dest_transform in the correct format as per https://gdal.org/tutorials/raster_api_tut.html
    dest_transform = (meta['transform'][2], meta['transform'][0], meta['transform'][1], \
                     meta['transform'][5], meta['transform'][3], meta['transform'][4])
    
    # https://gis.stackexchange.com/questions/158503/9999-no-data-value-becomes-0-when-writing-array-to-gdal-memory-file
    # an error arises for writing gdal files that truncates numbers
    # tmp = gdal.AutoCreateWarpedVRT(outFileRead, src_wkt, dst_wkt)

    NP2GDAL_CONVERSION = set_gdal_datatype_conv()

    if data_write_arr.dtype.name == "int64":
        gdaltype = 7
    else:
        gdaltype = NP2GDAL_CONVERSION[data_write_arr.dtype.name]

    driver = gdal.GetDriverByName("GTiff")

    # setting up the destination file with options taken from: https://gis.stackexchange.com/questions/299059/osgeo-gdal-translate-how-to-set-compression-on-gdal-gtiff-driver
    # this uses LZW compression
    outdata = driver.Create(dest_file, rows, cols, 1, gdaltype, options=['COMPRESS=LZW'])

    # generate the spatial reference using the crs you are using
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(target_crs)

    # outdata.SetProjection(meta['crs'])
    outdata.SetProjection(srs.ExportToWkt())
    outdata.SetGeoTransform(dest_transform)

    band = outdata.GetRasterBand(1)

    if meta['nodata'] is not None:
        band.SetNoDataValue(meta['nodata'])
        band.WriteArray(np.full((cols, rows), meta['nodata']))
    else:
        # TODO: Write wrapper functions to cover warnings like this. This is just for a quick fix
        logging.warning("SOURCE FILE NODATA VALUE IS NONE. This could potentially lead to future save rasters that result in problematic nodata conversions. If you think this is correct behavior, please ignore this warning.")
        band.SetNoDataValue(None)##if you want these values transparent
        band.WriteArray(np.full((cols, rows), None))

    band.WriteArray(data_write_arr)
    

    outdata.FlushCache() ##saves to disk!!
    outdata = None
    band = None
    ds = None


def write_gdf_to_raster_rio(data:gpd.GeoDataFrame, selector:str, dest_file:str, meta:dict, band:int = 1):
    '''
    This function is to raster a NumPy array into a destination file using rasterio. 

    Parameters
    ----------
        
        data: a geopandas GeoDataFrame containing target data (Required)
        
        select: a column name to select from the GeoDataFrame (Required)

        src_file: a string containing the name of the source file: '../data/dem.tif' (Required)

        dst_file: a string containing the name of the destination file: '../data/dem.tif' (Required)
        
        meta: a dictionary containing the metadata for the source file (Required)

        band: the band to use for reading. If none provided defaults to 1 (Optional)

    Output
    ------

        None
    
    Examples
    --------
    
    open_and_clip_datasets(data = gdf, selector = "rr", src_file = "../data/experiments/jakarta/base_maps/rr.tif", 
    dest_file = "../data/experiments/jakarta/fixed_raster/rr.tif", meta = meta)

    '''
    print("---- WRITING TO DESTINATION FILE ----")
    print("WRITING FILE FOR: ", selector)
    with rio.open(dest_file, 'w+', **meta) as out:
        out_arr = out.read(band)
    
        print(out_arr.shape)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(data.geometry, data[selector].astype(rio.float32)))

        burned = features.rasterize(shapes = shapes, fill = 0, out = out_arr, transform = out.transform, dtype = rio.float32)
        out.write_band(band, burned)
        out.close()
        
    print(f"Successfully wrote raster file for {selector} at {dest_file}")
    return None


def create_directory(path:str):
    '''
    Creates a new folder in the path

    Parameters
    ----------
        
        path: a directory name with the directory you want to create. (Required)


    Output
    ------

       None
    
    Examples
    --------

    create_directory('../data/experiments/jakarta/reprojection/.tif')

    '''

    if not os.path.isdir(path):
        path = path.lower()
        os.makedirs(path)
        print("created new directory at: {}".format(path))
    
    return

def show_sample(file:str):
    '''
    Show a plot of the first (1) band of the dataset

    Parameters
    ----------
        
        file: a file name with the full directory path (Required)


    Output
    ------

       None
    
    Examples
    --------

    show_sample("../data/experiments/jakarta/base_maps/mask.tif")

    '''

    raster = rio.open(file, 'r')

    red = raster.read(1)
    red = red.astype(float)

    show(red)

