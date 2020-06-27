import subprocess, os, shapely, fiona, logging, sys

import geopandas as gpd, pandas as pd, rasterio as rio, osmnx as ox

from shapely.geometry import Point
from osgeo import ogr, gdal
from disaster_preprocessing.code.common import *
from disaster_preprocessing.code.extract_data.griddata import GridData

def shapeify(inputVec:str, outputImg:str, \
             refImg:str, nodata:float = 10e-5):
    '''
    A script to rasterise a shapefile to the same projection & pixel resolution as a reference img.

    Inspiration from https://gis.stackexchange.com/questions/222394/how-to-convert-file-shp-to-tif-using-ogr-or-python-or-gdal


    Parameters
    ----------

        inputVec: a filename to a .shp file. (Required)

        outputImg: a filename containing the output file to save the new .tif to (Required)

        refImg: a reference file, in .tif format (Required)

        nodata: a nodata value to use. If none supplied, default is 10e-5 (Optional)
        
    Output
    ------

        None
    
    Examples
    --------

    shapeify(inputVec="../data/river/jkt_river_osm/n_jkt_rivers.shp", outputImg='../data/river/north_jkt_river_chan.tif',\
         refImg="../data/experiments/exp3/base_maps/land_cover_rr.tif")

    '''

    inputVector, outputimg, refimg = inputVec, outputImg, refImg
    print("REFERENCE IMAGE FILE:", refimg)
    gdalformat, datatype, burnVal = 'GTiff', gdal.GDT_Byte, 1
    
    img = gdal.Open(refimg, gdal.GA_ReadOnly)
    Imgband = img.GetRasterBand(1)
    
    # Open Shapefile
    shapefile = ogr.Open(inputVector)
    shapefile_layer = shapefile.GetLayer()

    # Rasterise
    print("Rasterising shapefile...")
    Output = gdal.GetDriverByName(gdalformat).Create(outputImg, img.RasterXSize, img.RasterYSize, 1, \
                                                     7, options=['COMPRESS=DEFLATE'])

    Output.SetProjection(img.GetProjectionRef())
    Output.SetGeoTransform(img.GetGeoTransform()) 

    # Write data to band 1
    Band = Output.GetRasterBand(1)
    nd = Imgband.GetNoDataValue()
    print("NO DATA VALUE", nd)

    if nd is not None:
        Band.SetNoDataValue(Imgband.GetNoDataValue())
    else:
        logging.warning("SOURCE FILE NODATA VALUE IS NONE")
        Band.SetNoDataValue(nodata)

    gdal.RasterizeLayer(Output, [1], shapefile_layer, burn_values=[burnVal])

    # Close datasets
    Band = None
    Output = None
    img = None
    Shapefile = None

    # Build img overviews
    subprocess.call("gdaladdo --config COMPRESS_OVERVIEW DEFLATE "+outputimg+" 2 4 8 16 32 64", shell=True)
    print("Done.")


def get_place_river(place:list = ['North Jakarta', 'Indonesia'], *args):
    '''
    this function is used to return the type of river features you want to extract from OSMnx. A list of available options
    can be found here: https://wiki.openstreetmap.org/wiki/Key:waterway


    Parameters
    ----------

        place: a list of names about the place you want to get information for (Required)

        *args: optional arguments to pass into OpenStreetMaps using OSMnx.  (Required)
        If you want to pass in arguments, pass in the arguments sequentially.
        
    Output
    ------

        A geopandas GeoDataFrame containing the information from OpenStreetMaps
    
    Examples
    --------

    shapeify(inputVec="../data/river/jkt_river_osm/n_jkt_rivers.shp", outputImg='../data/river/north_jkt_river_chan.tif',\
         refImg="../data/experiments/exp3/base_maps/land_cover_rr.tif")

    '''
    
    all_args = '|'.join([i.lower() for i in args])
    parsed_arg = f'way["waterway"~"{all_args}"]'
    
    place_name = ', '.join(place)
    G = ox.graph_from_place(place_name,
                       retain_all=True, truncate_by_edge=True, simplify=True,
                       network_type='none', infrastructure=parsed_arg)
    river_nodes, river_info = ox.graph_to_gdfs(G)
    return river_info
    
def preprocess_gdf(gdf:gpd.GeoDataFrame, d_type:str='numeric'):
    '''
    A function to select certain data types from the geodataframe.

    Parameters
    ----------

        gdf: a geopandas GeoDataFrame object (Required)

        d_type: a string that is either "numeric" or "object". If anything else is passed, it returns both "numeric" and "object" items.
        If no argument is passed, it defaults to "numeric". (Optional)
        
    Output
    ------

        A geopands GeoDataFrame object.
    
    Examples
    --------

    save_to_gj(gdf, '../data/expriments/jakarta/file.geojson')

    '''

    if d_type == "numeric":
        save_type = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    elif d_type == "object":
        save_type = ['object']
    else:
        save_type = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64', 'object']# iterable
    gdf_geom = gdf['geometry']
    gdf = gdf.select_dtypes(include=save_type)
    gdf['geometry'] = gdf_geom
    return gdf

def explode_df(gdf:gpd.GeoDataFrame):
    '''
    A function to call pandas transform on each element of a list-like to a row, replicating index values.
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.explode.html


    Parameters
    ----------

        gdf: a geopandas GeoDataFrame object (Required)
        
    Output
    ------

        A geopands GeoDataFrame object.
    
    Examples
    --------

    gdf = explode_df(gdf)

    '''

    df = pd.DataFrame(gdf)

    for item in df.columns.tolist():
        df = df.explode(item)

    gdf = gpd.GeoDataFrame(df)
    gdf.geometry = gdf['geometry']
    return gdf

def save_to_gj(gdf:gpd.GeoDataFrame, dest:str):
    '''
    A function to save a geodataframe to a geojson file.
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.explode.html


    Parameters
    ----------

        gdf: a geopandas GeoDataFrame object (Required)

        dest: the destination where you want to save your file. MUST end in .geojson or .json (Required)
        
    Output
    ------

        A geopands GeoDataFrame object.
    
    Examples
    --------

    save_to_gj(gdf, '../data/expriments/jakarta/file.geojson')

    '''

    if "geojson" not in dest or "json" not in dest:
        return None

    # building_info.drop(drop_cols)#.to_file("../data/river/n_jkt_rivers.geojson", driver="GeoJSON")
    with open(dest, 'w') as f:
        f.write(building_info.to_json())
    return None

def read_rivers(file:str):
    '''
    A function to read a vector based format such as geojson. To get more information on the available formats, read here:
    https://geopandas.org/io.html

    In general, geopandas.read_file() is pretty smart and should do what you want without extra arguments, but for more help, type:

    import fiona; help(fiona.open)

    Parameters
    ----------

        file: the file you want to read. (Required)

    Output
    ------

        A geopands GeoDataFrame object.
    
    Examples
    --------

    read_rivers("../data/experiments/jakarta/jakarta.geojson")

    '''
    return gpd.read_file(file)
    

def obtain_crs(file:str):
    '''
    A function to obtain the coordinate reference system of the file.

    To learn more about CRS, check here: https://pyproj4.github.io/pyproj/stable/crs_compatibility.html

    Hence, it requires GDAL to be installed and in the environemnt
    path variable.

    For a full list of options you can use to modify the command, go to https://gdal.org/programs/gdalwarp.html

    Parameters
    ----------

        file: the file you want to obtain the coordinate reference system for. It must be in .tif format . (Required)

    Output
    ------

        A string that contains the coordinate reference system information
    
    Examples
    --------

    set_up_cubic_spline('../data/experiments/jakarta/exp5/', '../code/' 'landunits.tif', 'landunits_cubic.tif')

    '''

    with rio.open(file) as src:
        meta = src.meta.copy()
        meta.update(compress='lzw')
    return meta['crs']

# %%time
def get_percentage_corrected(original_points:list, corrected_points:list):
    '''
    This is a function used to count the percentage of corrected cells.

    Parameters
    ----------

        original_points: a list of data points (Required)

        corrected_points: a list of corrected datapoints (Required)


    Output
    ------

        None

    Examples
    --------

    zoom = set_zoom_scale(hi_res_meta, ref_meta)
    grids = GridData(data = bldg_arr, zoom = zoom)

    points, geos = generate_river_fill(ref_data.data, zoom = zoom, meta = ref_meta, grids = grids, error_thres = 2e-1, algo = "building")
    get_percentage_corrected(points, ref_data)

    '''
    
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
    
    
def generate_res_data(folder:str, ref_file:str = 'dem_hi_res.tif', map_type:str = "river", folder_ext:str = ""):
    '''
    A function to run the river raster functions with reference files:

    Parameters
    ----------

        folder: the current working folder (Required)
        
        ref_file: a reference file to use. If none provided, uses dem_hi_res.tif (Optional)

        map_type: either "building" or "river". If none provided, uses "river". (Optional)

        folder_ext: a folder extension if any. (Optional)


    Output
    ------

        None
    
    Examples
    --------

    maps_for_generation = ["river", "buildings"]
    folder_ext = ""#'generated/'
    for map_item in maps_for_generation:
        # generate high resolution data.
        generate_res_data(folder = working_folder + generated_ext, \
                          ref_file = required_ext + 'datafiles/dem_hi_res.tif', map_type = map_item, folder_ext = folder_ext)
        
        # generate low resolution data.
        generate_res_data(folder = working_folder + generated_ext, \
                          ref_file = generated_ext + 'datafiles/dem_hi_res.tif', map_type = map_item, folder_ext = folder_ext)

    '''

    refFil = ref_file

    if map_type == "river":

        if "hi_res" in ref_file:
            dstFil, saveShpFil = folder + folder_ext + 'fixed_raster/chanmask_hi_res.tif', folder + folder_ext + 'river/river.shp'
        else:
            dstFil, saveShpFil = folder + folder_ext + 'fixed_raster/chanmask.tif', folder + folder_ext + 'river/river.shp'

        if not checkExist(folder + folder_ext + "river/", 'river.shp'):
            gdf = generate_osm_rast(['North Jakarta', 'Indonesia'], "river", 'river', 'stream', 'riverbank', 'tidal', 'channel')
            gdf.to_file(driver = 'ESRI Shapefile', filename = saveShpFil)

        else:
            gdf = gpd.read_file(saveShpFil)

    elif map_type == "building":
        if "hi_res" in ref_file:
            dstFil, saveShpFil = folder + folder_ext + 'fixed_raster/building_hi_res.tif', \
                                folder + folder_ext + 'building/building.shp'
        else:
            dstFil, saveShpFil = folder + folder_ext + 'fixed_raster/building.tif', \
                                folder + folder_ext + 'building/building.shp'

        if not checkExist(folder + folder_ext + "building/", 'building.shp'):
            gdf = generate_osm_rast(['North Jakarta', 'Indonesia'], "building")
            gdf.to_file(driver = 'ESRI Shapefile', filename = saveShpFil)

        else:
            gdf = gpd.read_file(saveShpFil)


    print("SHAPE FILE REFERENCE", saveShpFil)
    print("DESTINATION FILE", dstFil)
    
    
    shapeify(inputVec=saveShpFil, outputImg=dstFil,\
         refImg=refFil)
    
    print("SUCCESFULLY RASTERED", map_type, "TO .TIF FILE AT", dstFil, "USING", refFil)

def generate_maps(working_folder:str, folder:str, map_iter:list, ext:dict):
    '''
    A function to generate both the building and river maps.

    Parameters
    ----------

        working_folder: the overall working folder (Required)

        folder: the current working folder (Required)

        map_iter: a list of different maps to iterate over. (Required)

        ext: a dictionary of folder extensions. (Required)


    Output
    ------

        None
    
    Examples
    --------

    maps_for_generation = ["river", "building"]
    extensions = {"generated":generated_ext, "required":required_ext}

    generate_maps(working_folder = "../data/sample/", folder = "../data/sample/generated/raster_fix/", map_type = "river", folder_ext = extensions)
    '''
    
    dem_hi_res = 'dem_hi_res.tif'
    dem_folder = working_folder + ext['required'] + 'datafiles/'
    print("HIGH RES DEM FOLDER", dem_folder)
    
    if os.path.exists(dem_folder + 'modified_dem/' + dem_hi_res) == False:
        print("CLIPPED HIGH RES DEM DOES NOT EXIST. CLIPPING NOW")
        # clip the high resolution dem if it has not already been clipped
        open_and_clip_datasets(folder =  dem_folder, files = [dem_hi_res], \
                               working_folder = working_folder, dest_folder = dem_folder + 'modified_dem/', \
                               gjson = gj_ref_file)

    for map_item in map_iter:
        if map_item == "river": title = "chanmask"
        else: title = "building"

        print("********** GENERATING CORRECTION LAYER FOR", map_item, "**********")

        # generate high resolution data.
        print(" ---- GENERATING HIGH RES DATA ----")
        generate_res_data(folder = working_folder, \
                          ref_file = dem_folder + "modified_dem/" + dem_hi_res, map_type = map_item, folder_ext = ext['generated'])

        print(" ---- GENERATING LOW RES DATA ----")
        # generate low resolution data.
        generate_res_data(folder = working_folder, \
                          ref_file = working_folder + ext['generated'] + 'fixed_raster/dem.tif', map_type = map_item, folder_ext = ext['generated'])


        # open the low res and the high res file we just created
        print("OPENING REFERNCE LOW RES DEM FILE")
        reference_file = folder + "dem.tif"
        ref_data, ref_meta = open_rio(reference_file)

        print("OPENING GENERATED HIGH RES CHANNEL MASK FILE")
        src_file = folder + f'{title}_hi_res.tif'
        hi_res_data, hi_res_meta = open_rio(src_file)
        bldg_arr = hi_res_data.data

        print("SIZE OF NON ZERO CELLS (HIGH RES", "BUILDING", "DATA):", len(bldg_arr[bldg_arr != hi_res_meta['nodata']]))
        print("HIGH RESOLUTION GRID SIZE", hi_res_meta['width'] * hi_res_meta['height'])
        
        zoom = set_zoom_scale(hi_res_meta, ref_meta)
        grids = GridData(data = bldg_arr, zoom = zoom, nd = hi_res_meta['nodata'])

        points, geos = generate_river_fill(ref_data.data, zoom = zoom, meta = ref_meta, grids = grids, error_thres = 2e-1, algo = "building")
        get_percentage_corrected(points, ref_data)

        df = generate_dataframe_from_points(points = points, geos = geos)
        gdf = gpd.GeoDataFrame(df)
        
        print(gdf.head())
        src_file = folder + "dem.tif"
        
        write_gdf_to_raster_gdal(data = gdf, selector = "data", src_file = src_file, \
                                 dest_file = folder + f'{title}_2.tif', \
                                 meta = ref_meta)
        
        print("WROTE FILE FOR ITEM", map_item, "TO", folder + f'{map_item}_2.tif')


def generate_osm_rast(dest:str, map_type:str = 'river', *args):
    '''
    A function to read a vector based format such as geojson. To get more information on the available formats, read here:
    https://geopandas.org/io.html

    In general, geopandas.read_file() is pretty smart and should do what you want without extra arguments, but for more help, type:

    import fiona; help(fiona.open)

    Parameters
    ----------

        dest: a list containing the areas you want to extract from OpenStreetMaps using OSMnx. (Required)
        
        *args: optional arguments to pass into OpenStreetMaps using OSMnx.  (Required)
        If you want to pass in arguments, pass in the arguments sequentially.

    Output
    ------

        A geopands GeoDataFrame object.
    
    Examples
    --------

    river_rast(['North Jakarta', 'Indonesia'], 'river', 'stream', 'riverbank', 'tidal', 'channel')

    '''
    if map_type == "river":
        gdf = get_place_river(dest, *[i for i in args])

    elif map_type == "building":
        gdf = ox.footprints_from_place(dest)

    gdf = explode_df(gdf)
    gdf = preprocess_gdf(gdf, 'all')
    # gdf.crs = crs
    return gdf


if __name__ == '__main__':
    if len(os.args) < 1:
        refFil, dstFil, saveShpFil = '../data/experiments/semarang/base_maps/psi_fix.tif', '../data/experiments/semarang/base_maps/chanmask.tif', "../data/experiments/semarang/river/semarang_rivers.shp"
    else:
        refFil = os.argv[0]
        dstFil = os.argv[1]
        saveShpFil = os.argv[2]

        # crs, gj_file = ??

    # crs = obtain_crs(refFil)
    # gdf = get_place_river(['Indonesia', 'Semarang'], 'river', 'stream', 'riverbank', 'tidal', 'channel')

    gdf = river_rast(['Indonesia', 'Semarang'], 'river', 'stream', 'riverbank', 'tidal', 'channel')
    gdf.to_file(driver = 'ESRI Shapefile',\
                        filename= saveShpFil)

    shapeify(inputVec=saveShpFil, outputImg=dstFil,\
             refImg=refFil)
