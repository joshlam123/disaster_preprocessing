import rasterio as rio, numpy as np, geopandas as gpd, rioxarray as rx, os
from rasterio.plot import show, plotting_extent
from shapely.geometry import box, mapping
from osgeo import gdal, gdal_array
from disaster_preprocessing.code.common import *

def open_and_clip_datasets(folder:str, working_folder:str, dest_folder:str, files:list, 
                             nodata:float = 10e-5, gjson:str = "../data/jkt_shape/epsg4326/north_jakarta.geojson"):
    '''
    This function is used to open up data from the source folder and clip them according to a geojson file with the polygon coordinates.

    Parameters
    ----------

        folder: the current directory to get the files from (Required)
        
        working_folder: the current working folder for your files (Required)

        dest_foler: the destination folder to save your files in (Required)

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

    jkt_gdf = gpd.read_file(gjson)

    for idx, item in enumerate(files):

        dest_file = working_folder + dest_folder + item 

        if os.path.exists(dest_file) == False:
            print("----------- PROCESSING -----------")
            print("DEST FILE: ", dest_file)
            
            title = item[:item.find('.tif')]
            print("FILE TITLE: ", title)
            
            src_file = folder + item
            print("SOURCE FILE: ", src_file)
            
            dataset = rio.open(src_file) # CHANGE THIS TO YOUR OWN FOLDER
            xds = rx.open_rasterio(
                    dataset,
                    parse_coordinates=True,
                    chunks=True,
                    masked=False,
                )
            
            clipped = xds[0].rio.clip(jkt_gdf.geometry.apply(mapping), xds.rio.crs, all_touched = True, drop = True, invert = False)

            # d_typ = rio.dtypes.get_minimum_dtype(clipped.values)
            
            _, rio_crs, nodata, *unused = get_meta(folder, item)
            
            if 'dem' not in title:
                dem_src = working_folder + dest_folder + 'dem.tif'
                # data = perform_restoration(clipped, dem_src)
                data = perform_addition(clipped.values, dem_src)

                if 'landunits' in title: nodata = 0
                meta = generate_rx_meta_data(src_file, rx_data = xds, write_data = data, crs = rio_crs, nodata = nodata)

            else:

                meta = generate_rx_meta_data(src_file, rx_data = xds, crs = rio_crs, nodata = nodata)


            # rio_crs = change_to_mercator(target_crs)
            
            if 'dem' not in title:
                write_array_to_raster_gdal(data, src_file, dest_file, title, meta)
            else:
                write_array_to_raster_gdal(clipped.values, src_file, dest_file, title, meta)
            
    return None

if __name__ == '__main__':
    if len(os.argv) < 4:
        folder_dest = '../data/experiments/semarang/'
        incl = ''
        ctry, region = 'indonesia', 'Semarang'
    else:
        folder_dest = os.argv[0]
        incl = os.argv[1]
        ctry, region = os.argv[3], os.argv[4]

    b_map = folder_dest + 'base_maps/'
    d_map = folder_dest + 'fixed_raster/'
    gj_file = "../data/jkt_shape/epsg4326/" + f'{ctry}/{region}.geojson'

    folder, files = get_folder_files('', folder_dest)
    open_and_clip_datasets(folder, files, d_map, gj_file)