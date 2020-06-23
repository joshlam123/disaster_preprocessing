import pandas as pd, numpy as np
import os, json, math

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping
from preprocessing.common import *
from rasterio.mask import mask
from rasterio.plot import show, plotting_extent
from rasterio.merge import merge
from pyproj import Proj, transform
from osgeo import gdal,osr
from rasterio import features
from rasterio.crs import CRS

import geopandas as gpd, networkx as nx, earthpy.plot as ep, rasterio as rio, matplotlib.pyplot as plt


def convert_params(file:str = '../data/landcover/epsg4326/convert_files/params.csv'):
    '''
    A function to process to convert a csv file containing the parameters parameters into landcover.

    Parameters
    ----------

        file: a full path filename containing the .csv file with the parameters for landcover. If default, it will 
        look in the landcover/epsg4326/convert_files folder for params.csv (Optional)

    Output
    ------

        data_conv_dict: a dictionary mapping each number encoding to several other parameters associated with it.

        trans_dict: a dictionary used to map the string name for land use to a number encoding
    
    Examples
    --------
    
    data_conv_dict, trans_dict = convert_params(file = '../data/landcover/epsg4326/convert_files/params.csv')
    

    '''
    
    data_convert = pd.read_csv(file)
    # convert to uppercase for consistency
    data_convert['landunits'] = data_convert['landunits'].apply(lambda x: x.upper())
    
    data_conv_dict = data_convert.to_dict('id')
    trans_dict = {v['landunits']:k for k,v in data_conv_dict.items()}
    
    return data_conv_dict, trans_dict


# opening the mapping file and reading it
def read_category(category_file:str = "../data/landcover/epsg4326/convert_files/mapping.json"):
    '''
    This function is to raster a NumPy array into a destination file using GDAL. 

    Parameters
    ----------

        category_file: a full path filename containing the .json file with the full landcover parameter mapping converted from Google Earth. 
        If default, it will look in the landcover/epsg4326/convert_files folder for mapping.json (Optional)

        For more details: look here https://developers.google.com/earth-engine/datasets/catalog/ESA_GLOBCOVER_L4_200901_200912_V2_3
    
    Output
    ------

        category: a dictionary containing the mapping of category to land use. 
    
    Examples
    --------

    read_category(category_file = "../data/landcover/epsg4326/convert_files/mapping.json")

    '''

    with open(category_file, "r") as f:
        category = json.load(f)
    category = {int(k):v for k,v in category.items()}
    return category


def mapToCategory(item):
    '''
    This function is used to map a data entry within the landcover map to it's category. Meant to be used as a sort of pandas vectorisation
    inspired from: https://towardsdatascience.com/how-to-make-your-pandas-loop-71-803-times-faster-805030df4f06

    Parameters
    ----------

        item: any indexable object. (Required)
    
    Output
    ------

        category[item]: is a NumPy array containing the mapping of items.
    
    Examples (taken from generate_lcm.py)
    --------

    all_lidar_categories = list()
    all_info_cat = list()
    
    for idx, item in enumerate(lidar_dem.data):
        if idx % 1000 == 0: print(idx/len(lidar_dem)*100)
        all_lidar_categories.append(list(map(mapToCategory, item)))

    '''

    return category[item]

# mapping function that will be used to map a category to it's corresponding landcover values
def mapToLCM(item):
    '''
    This function is to raster a NumPy array into a destination file using GDAL. 

    Parameters
    ----------

        item: any indexable object. (Required)
    
    Output
    ------

        multiple return items: containing the mapping of either
    
    Examples
    --------

    all_lidar_categories = list()

    for idx, item in enumerate(all_lidar_categories):
        if idx % 1000 == 0: print(idx/len(lidar_dem)*100)
        all_info_cat.append(list(map(mapToLCM, item)))

    '''

    # If landcover is unclassified, it is treated as having 0 entry.
    if item == "UNCLASSIFIED":
        return 0, "UNCLASSIFIED", 0, 0, 0, 0, 0, 0, 0, 0
    else:
        translated_item = data_conv_dict[trans_dict[item]]
        return translated_item['id'], translated_item['landunits'], translated_item['n'], translated_item['rr'], translated_item['per'],\
                translated_item['ch'], translated_item['lai'], translated_item['ksat1'], translated_item['thetas1'],\
                 translated_item['thetai1'], translated_item['psi']

def generate_keys_dict(data_conv_dict:dict, all_info_cat:list):
    '''
    This function is to generate a dictionary of keys that map a category to it's corresponding values

    Parameters
    ----------

        data_conv_dict: a dictionary mapping each number encoding to several other parameters associated with it. (Required)

        all_info_cat: a list containing the mapping of each grid cell to a category (Required)
    
    Output
    ------

        all_keys: a list of the parameters in the dictinoary

        all_info_dict: an empty dictionary initialised for each parameter in all_keys, with each key having a list() as a value to
        hold the array map data.
    
    Examples
    --------

    write_array_to_raster_gdal(data=np.array([1,2,3,4,5]), src_file="../data/experiments/jakarta/base_maps/rr.tif", 
    dest_file="../data/experiments/jakarta/fixed_raster/rr.tif")

    '''

    all_keys = list(data_conv_dict[0].keys())
    all_info_dict = dict()
    for k in all_keys:
        all_info_dict[k] = list()
        
    for idx, k in enumerate(all_keys):
        for idx, item in enumerate(all_info_cat):
            all_info_dict[k].append(np.array([x[0] for x in item], dtype=np.float32))
            
    return all_keys, all_info_dict
    

# generate a geopandas dataframe by concatenating all items in the numpy array into a dataframe and turning geometries into Points
def generate_dataframe(all_info_cat:list, longs:np.array, lats:np.array):
    '''
    This function is to raster a NumPy array into a destination file using GDAL. 

    Parameters
    ----------

        all_info_cat: a list containing the mapping of each grid cell to a category (Required)

        longs: a NumPy array containing longitude data. (Required)

        lats: a NumPy array containing latitude data. (Required)  
    
    Output
    ------

        None
    
    Examples
    --------

    longs, lats, p2 = affine_projection(src_file)
    
    data = generate_dataframe(all_info_cat, longs, lats)

    '''

    df = pd.DataFrame(columns=["long", "lat", "id", "landunits", "n", "rr", "per", "ch", "lai", "ksat1", "thetas1", "thetai1", "psi1"])

    for idx, item in enumerate(all_info_cat):
        if idx % 500 == 0: print(idx/len(all_info_cat)*100)

        df1 = pd.DataFrame(item, columns=["id", "landunits", "n", "rr", "per", "ch", "lai", "ksat1", "thetas1", "thetai1", "psi1"])

        df1["Long"] = longs[idx]
        df1["Lat"] = lats[idx]

        frames = [df, df1]
        df = pd.concat(frames)
    
    geometry = [Point(xy) for xy in zip(df.Long, df.Lat)]
    df['geometry'] = geometry
    return df

def fix_df_na(data:pd.DataFrame, na_val:float = 0):
    '''
    This function is to fillna values with 0

    Parameters
    ----------

        data: a pandas DataFrame containing target data (Required)

        na_val: the NA fill value to replace the NaN with. If none supplied, defaults to 0. (Optional)
    
    Output
    ------

        None
    
    Examples
    --------

    longs, lats, p2 = affine_projection(src_file)
    
    df = generate_dataframe(all_info_cat, longs, lats)
    df = fix_df_na(df, 10e-5)

    '''
    for column in data.columns.tolist():
        if column == 'landunit':
            data[column] = data['id'].fillna(na_val = na_val)
        else:
            data[column] = data[column].fillna(naval = na_val)

    # data['n'] = data['n'].fillna(naval)
    # data['rr'] = data['rr'].fillna(naval)
    # data['per'] = data['per'].fillna(naval)
    # data['ch'] = data['ch'].fillna(naval)
    # data['lai'] = data['lai'].fillna(naval)
    # data['thetas1'] = data['thetas1'].fillna(naval)
    # data['thetai1'] = data['thetai1'].fillna(naval)
    # data['ksat1'] = data['ksat1'].fillna(naval)
    # data['psi1'] = data['psi1'].fillna(naval)

    return df

def generate_geo_df(data:pd.DataFrame, p2:str, drop_columns:list, drop_all:bool = False):
    '''
    This function is to raster a NumPy array into a destination file using GDAL. 

    Parameters
    ----------

        data: a pandas DataFrame containing target data (Required)

        p2: a string with the coordinate reference system of the geodataframe (Required)
        
        drop_columns: a list that contains the columns to be dropped from the dataframe. (Required)

        drop_all: an indicator flag whether to drop all columns from the dataframe. If empty, it defaults to False. (Optional)
    
    Output
    ------

        None
    
    Examples
    --------

    longs, lats, p2 = affine_projection(src_file)
    
    df = generate_dataframe(all_info_cat, longs, lats)
    df = fix_df_na(df, 10e-5)
    gdf = generate_geo_df(df = df, p2 = p2, drop_columns = ['landunits', 'id', 'Long', 'Lat'])

    '''

    if drop_all == True:
        drop_columns = ['landunits', 'id', 'long', 'lat', 'n', 'rr', 'per', 'ch', 'lai', 'ksat1', 'thetas1', 'thetai1', 'psi']

    df_copy = data.drop(columns = drop_columns, axis = 1)
    gdf = gpd.GeoDataFrame(df_copy, crs = p2, geometry = df['geometry'])
    return gdf

def get_geo_df_type_meta(gdf):
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

    d_type = dict()
    for item in gdf.columns:
        print(item)
        if item not in ["geometry"]:
            d_type[item] = rio.dtypes.get_minimum_dtype(gdf[item])
            print(f"Minimum Dtype for {item} : {d_type[item]}")
    return d_type




category = read_category()
data_conv_dict, trans_dict = convert_params()