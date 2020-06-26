from disaster_preprocessing.code.python_preprocessing.pp_raster import *
from disaster_preprocessing.code.common import *
from rasterio import features
import os

def generate_pre_dict(src_file:str, category:dict, data_conv_dict:dict, trans_dict:dict):
    '''
    A function to process to obtain a map a dictionary containing the category mapping to a numpy array

    Parameters
    ----------

        src_file: a string with the file name of your params file (Required)

        category: a dictionary mapping a land cover type to a numerical index (Required)

        data_conv_dict: a dictionary containing the mapping of numerical indices to parameters (Required)
        
        trans_dict: a dictionary containing the translation of parameters to extended parameters (Required)

    Output
    ------

        meta: a dictionary containing the metadata of your source file

        all_info_cat: a list containing the mapping of each grid cell to a category 
    
    Examples
    --------
    generate_pre_dict("../data/experiments/jakarta/land_cover.tif")

    '''

    rasterPrepare = PrepareRaster(category = category, data_conv_dict = data_conv_dict, trans_dict = trans_dict)

    lidar_dem, meta = open_rio(src_file)
    
    all_lidar_categories = list()
    all_info_cat = list()
    
    for idx, item in enumerate(lidar_dem.data):
        if idx % 1000 == 0: print(idx/len(lidar_dem)*100)
        all_lidar_categories.append(list(map(rasterPrepare.mapToCategory, item)))

    for idx, item in enumerate(all_lidar_categories):
        if idx % 1000 == 0: print(idx/len(lidar_dem)*100)
        all_info_cat.append(list(map(rasterPrepare.mapToLCM, item)))
    
    all_keys, all_info_dict = generate_keys_dict(data_conv_dict, all_info_cat)
    return meta, all_info_cat
    
def process_df(all_info_cat:list, src_file:str):
    '''
    A function to process to convert all_info_cat into a geodataframe

    Parameters
    ----------

        all_info_cat: a list containing the mapping of different categories of land cover (Required)

        src_file: a string with the file name of your land cover source file (Required)

    Output
    ------

        gdf: a geopandas GeoDataFrame object with the columns from all_info_cat

        d_types: a dictionary with the datatypes of each column in gdf

        p2: a string with the coordinate reference system of the geodataframe
    
    Examples
    --------
    
    meta, all_info_cat = do_preprocessing("../data/landcover/epsg4326/base_map/land_cover_jakarta.tif")
    gdf, d_types, crs = process_df(all_info_cat, "../data/landcover/epsg4326/base_map/land_cover_jakarta.tif")
    

    '''

    longs, lats, p2 = affine_projection(src_file)
    
    data = generate_dataframe(all_info_cat, longs, lats)
    data = fix_df_na(data = data, na_val = 10e-5)
    gdf = generate_geo_df(data = data, p2 = p2, drop_columns = ['id', 'long', 'lat'])
    
    d_types = get_geo_df_type_meta(gdf)
    d_types['landunits'] == 'float32'

    return gdf, d_types, p2
    
def do_preprocessing(srcfile:str, category:dict, data_conv_dict:dict, trans_dict:dict):
    '''
    A function to preprocess a sourcefile into a list with the mapped categories

    Parameters
    ----------

        srcfile: a string with the .tif map source file (Required)

        category: a dictionary mapping a land cover type to a numerical index (Required)

        data_conv_dict: a dictionary containing the mapping of numerical indices to parameters (Required)
        
        trans_dict: a dictionary containing the translation of parameters to extended parameters (Required)


    Output
    ------

        meta: a dictionary containing the metadata of the map source file 

        all_info_cat: a list containing the mapping of different categories of land cover
    
    Examples
    --------
    
    meta, all_info_cat = do_preprocessing("../data/landcover/epsg4326/base_map/land_cover_jakarta.tif")
    

    '''

    meta, all_info_cat = generate_pre_dict(srcfile, category, data_conv_dict, trans_dict)

    print("Successfully generated metadata and category mapping for all files!")
    return meta, all_info_cat

