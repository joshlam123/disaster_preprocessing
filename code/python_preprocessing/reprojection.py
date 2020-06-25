import rasterio as rio
from rasterio import features, Affine as A
from rasterio.plot import show
from rasterio.plot import plotting_extent
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from disaster_preprocessing.code import common

from shapely.geometry import Point

from osgeo import gdal, gdal_array, osr
from pyproj.crs import CRS
from distutils.version import LooseVersion
from pyproj.enums import WktVersion
from osgeo import osr

import pandas as pd, numpy as np, pandas as pd, rioxarray as rx, geopandas as gpd

import fiona, gdal, os
from tqdm import tqdm

def get_dest_affine(transform:list):
    '''
    This function is similar to get the destination affine transformation

    Parameters
    ----------

        transform: an affine transformation object containing the affine transform

    Output
    ------

        dest_trans: an affine transformation object
    
    Examples
    --------

    get_dest_affine(Affine(30, 0, 106, 0, -30, -6))

    '''
    dest_trans = A.translation(transform[2] + transform[2], transform[5] + transform[5]) * A.scale(transform[0], transform[4])


    return dest_trans

def return_meta(working_folder:str, dest_folder:str, p2:str, src_file:str, map_type:str = "mask", nodata:float = 10e-5):
    '''
    This takes in a CRS and sourcefile to get the calculated transform, map data, and the metadata.

    Parameters
    ----------
        
        working_folder: the working folder of your code (Required)

        dest_folder: a string containing the coordinate reference system of the file (Required)

        p2: a string containing the coordinate reference system of the file (Required)

        src_file: a string representing the source file .tif you want to return the meta data for (Required)

        map_type: a string which has the map type. If none supplied, it will default to mask (Optional)

        nodata: a float value you want for the no data fill value. Defaults to 10e-5 (Optional)

    Output
    ------
        
        meta: a dictionary with the source file's metadata

        map_src: the map data in NumPy array for
    
    Examples
    --------

    return_meta(, , CRS.from_wkt('PROJCS["Google_Maps_Global_Mercator",GEOGCS["WGS 84",DATUM["WGS_1984",
    SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]],PROJECTION["Mercator_1SP"],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],
    PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],
    AXIS["Northing",NORTH],EXTENSION["PROJ4","+proj = merc +a = 6378137 +b = 6378137 +lat_ts = 0 +lon_0 = 0 +x_0 = 0 +y_0 = 0 +k = 1 +units = m 
    +nadgrids = @null +wktext +no_defs"]]'), "../data/experiments/jakarta/dem.tif", "dem", 10e-5)

    '''

    if "dem" in map_type:
        with rio.open(src_file) as src:
            transform, _, _ = calculate_default_transform(
            src.crs, p2, src.width, src.height, *src.bounds)
            height, width = map_src.shape
            meta = src.meta.copy()
            meta.update(nodata = -10e38)
            map_src = src.read(1, masked = False)

    else:
        xds, crs, nodata, bounds, _, _, _ = get_meta('', src_file)

        *other_non_essential_data, width, height, transform = get_meta(working_folder+dest_folder, 'dem.tif')

        # grab the meta data because we'll need it later

        meta = dict()
        uniq_val = np.unique(xds[0].values)
        if len(uniq_val[uniq_val<=10e-5]) > 1:
            nodata = uniq_val[uniq_val<=10e-5][1]
        meta['nodata'] = nodata

        map_src = xds[0].values
    
    meta['compress'] = 'lzw'
    meta['dtype'] = 'float32'
    meta['transform'] = transform
    meta['width'] = width
    meta['height'] = height


    return meta, map_src




def create_points(working_folder:str, dest_folder:str, meta:dict, map_src:np.array, map_type:str = "mask", toggle:bool = False):
    '''
    This takes in a CRS and sourcefile to get the calculated transform, map data, and the metadata.

    Parameters
    ----------

        working_folder: the working folder of your code (Required)

        dest_folder: a string containing the coordinate reference system of the file (Required)
        
        meta: a dictionary with the metadata of the map (Required)

        map_src: a NumPy array with your map data (Required)

        map_type: a string which has the map type. If none supplied, it will default to mask (Optional)

        toggle: a boolean value to check whether you want to follow the dimensions of the source metadata or not. If empty defaults to False. (Optional)

    Output
    ------
        
        matching_geos: a NumPy array with geocoordinates of your map
        
        matching_points: a NumPy array with the points of your map
    
    Examples
    --------
    
    meta, dest_trans, map_src_data = return_meta(CRS.from_wkt('PROJCS["Google_Maps_Global_Mercator",GEOGCS["WGS 84",DATUM["WGS_1984",
    SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]],PROJECTION["Mercator_1SP"],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],
    PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],
    AXIS["Northing",NORTH],EXTENSION["PROJ4","+proj = merc +a = 6378137 +b = 6378137 +lat_ts = 0 +lon_0 = 0 +x_0 = 0 +y_0 = 0 +k = 1 +units = m 
    +nadgrids = @null +wktext +no_defs"]]'), "../data/experiments/jakarta/dem.tif")

    create_points("../data/experiments/jakarta/", "reprojection/", meta, map_src_data, "rr", "10e-5")

    '''

    trans = meta['transform']
    nodata = meta['nodata']
    print("NODATA value:", nodata)

    counter = 0 
    matching_geos = list()
    matching_points = list()

    print("---------------------")
    no_nan = len(map_src[map_src == nodata]) 
    print("NUMBER OF NO DATA ENTRIES FOR", map_type, ":", no_nan)
    print("TRANSFORM: ", trans)

    for idx, item in tqdm(enumerate(map_src)):
        
        geo, points = list(), list()

        top_left_x, top_left_y = trans[2], trans[5]+(counter * trans[4])


        for edx, entry in enumerate(item):
            

            if map_type == "mask":
                if entry !=  nodata:
                    entry = 1

            if np.isnan(entry) == True or entry == nodata:
                entry = nodata 
                    
            geo.append(Point([top_left_x, top_left_y]))
            points.append(entry)
            top_left_x += trans[0]
        
        if toggle:

            while len(geo) != meta['width']:
                if len(geo) < meta['width']:
                    top_left_x += trans[0]
                    geo.append(Point([top_left_x, top_left_y]))
                    points.append(nodata)

                elif len(geo) > meta['width']:
                    geo = geo[:meta['width']]
                    points = points[:meta['width']]

        matching_geos.append(geo)
        matching_points.append(points)

        counter += 1
            
    if toggle:
        while len(matching_geos) != meta['height']:

            if len(matching_geos) < meta['height']:

                top_left_x = trans[2]
                top_left_y = trans[5] + (counter * trans[4])

                matching_geos.append([Point(top_left_x + trans[4], top_left_y) for _ in range(meta['width'])])
                matching_points.append([nodata for _ in range(meta['width'])])
            
            elif len(matching_geos) > meta['height']:

                matching_geos = matching_geos[:meta['height']]
                matching_points = matching_geos[:meta['height']]

        print("AFTER COL MODIFICATIONS: ", len(points))
    # create the jakarta geopoly and return it
    
    return matching_geos, matching_points

def create_gdf(colname, map_crs, matching_geos, matching_points):
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

    df = pd.DataFrame(columns = ['geometry', colname])

    for idx, entry in tqdm(enumerate(matching_geos)):
        df_temp = pd.DataFrame({"geometry":entry, colname:matching_points[idx]})
        frames = [df, df_temp]
        df = pd.concat(frames)
    
    gdf = gpd.GeoDataFrame(df, geometry = df['geometry'])
    if colname !=  "dem":
        gdf[colname] = gdf[colname].fillna(10e-5)
    gdf.crs = map_crs

    print(df.shape)
    
    return gdf


    
def run_processing(folder:str, working_folder:str, dest_folder:str, target_crs:int, files:list, nodata:float = 10e-5):    
    '''
    This function is used to open up data from the source folder and clip them according to a geojson file with the polygon coordinates.

    Parameters
    ----------

        folder: the current directory to get the files from (Required)
        
        working_folder: the current working folder for your files (Required)

        dest_foler: the destination folder to save your files in (Required)

        target_crs: a crs to project to (Required)

        files: a list of files which you want to work with (Required)

        nodata: the no data replacement value. If no value, it is set to 10e-5 (Optional)

        gjson: the geojson file to clip from (Optional)

    Output
    ------

        None
    
    Examples
    --------

    open_and_clip_datasets("../data/jakarta/experiments/base_maps", ["rr.tif"], 
    "../data/jakarta/experiments/", "fixed_raster/", "../data/jkt_shape/epsg4326/north_jakarta.geojson")

    '''

    dem_src = folder + 'dem.tif'

    for idx, item in enumerate(files):
        
        # if 'mask' in item and 'chanmask' not in item:
        #     src_file = folder + 'dem.tif'
        # else:
            
        title = item[:item.find('.tif')]
        print("Map name:", title)

        if 'mask' in title and 'chanmask' not in title:
            src_file = folder + 'rr.tif'
        else:
            src_file = folder + item

        print("src_file:", src_file)
        dest_file = working_folder + dest_folder + item 
        print("destination file:", dest_file)
        

        if os.path.exists(dest_file) == False:
        
            rio_crs = change_to_mercator(target_crs)

            meta, map_src_data = return_meta(working_folder, dest_folder, rio_crs, src_file, title, nodata) # src_file

            if 'dem' in title:
                matching_geos, matching_points = create_points(working_folder, dest_folder, meta, map_src_data, title)
            else:
                matching_geos, matching_points = create_points(working_folder, dest_folder, meta, map_src_data, title, toggle = 1)

            gdf = create_gdf(title, rio_crs, matching_geos, matching_points)

            # write_gdf_to_raster_rio(gdf, title, src_file, destfile, meta)

            if 'dem' not in title:
                dem_file = working_folder + dest_folder + 'dem.tif'
                write_gdf_to_raster_gdal(data = gdf, selector = title, src_file = dem_file, dest_file = dest_file, meta = meta, target_crs = target_crs)

            else:
                write_gdf_to_raster_gdal(data = gdf, selector = title, src_file = src_file, dest_file = dest_file, meta = meta, target_crs = target_crs)
            
            
        else:
            print(f"File at {dest_file} exists!")
            
    return None