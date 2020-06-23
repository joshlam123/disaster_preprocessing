import rasterio as rio
from rasterio.plot import show
from rasterio.plot import plotting_extent
from shapely.geometry import box, mapping
from osgeo import gdal, gdal_array
import numpy as np
import rioxarray
import xarray
import earthpy.plot as ep
import matplotlib.pyplot as plt


def perform_restoration(data):
    dem_xds = rioxarray.open_rasterio(
            "../data/dem/dem_north_jakarta.tif",
        #     masked=True,
            parse_coordinates=True,
            chunks=True,
            masked=False,
        )
    row, cols = dem_xds.rio.height, dem_xds.rio.width
    data_cols, data_row = data.values.shape
    
    zero = np.array([0 for i in range(arr.shape[1])])
    
    while data_row != row:
        # hack
        zero = zero.reshape(1, -1)
        data.values = np.concatenate((data.values,zero), axis=0)
        data_cols, data_row = data.values.shape
    
    if data_cols != cols:
        zero = zero.reshape(-1,1)
        data.values = np.concatenate((data.values,zero), axis=1)
        data_cols, data_row = data.values.shape
        
    return data.values

def fix_nan_values(data):
    if len(data[np.isnan(data)]) > 0:
        clipped.values = np.where(np.isnan(data), 0, data)
        
    return clipped.values

def write_to_dataset(destination, dataset, clipped, write_type):
 
    file = dataset
    ds = gdal.Open(file)
    [cols, rows] = clipped.shape
    arr_min = clipped.min()
    arr_max = clipped.max()

    driver = gdal.GetDriverByName("GTiff")
    # outdata = driver.Create("../data/landcover/land_cover_coverc_3.tif", rows, cols, 1, gdal.GDT_UInt16)

    gdaltype = NP2GDAL_CONVERSION[clipped.dtype.name]

    # typecodeInteger = typemap[clipped.values.dtype.name]
    # typecode = gdal.GetDataTypeName(typecodeInteger)
    filename = f"{destination}land_cover_{dataset}_{write_type}.tif"
    
    outdata = driver.Create(filename, rows, cols, 1,\
                            gdaltype)
    outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    outdata.SetProjection(ds.GetProjection())##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(clipped)
    if write_type == "large":
        outdata.GetRasterBand(1).SetNoDataValue(10e5) #if you want these values transparent
    else:
        outdata.GetRasterBand(1).SetNoDataValue(-10e5) #if you want these values transparent
    outdata.FlushCache() ##saves to disk!!
    outdata = None
    band=None
    ds=None
    
    print(f"Wrote file {dataset} to {filename}")

    return None

def open_and_clip_raster_datasets(folder, files):
    wtypes = ['large', 'small']
    jkt_gdf = gpd.read_file("../data/jkt_shape/north_jakarta.geojson")
    destination = '../data/landcover/generated_maps/'
    
    NP2GDAL_CONVERSION = set_gdal_datatype_conv()
    
    for idx, item in enumerate(files):
        title = item[item.find('cover_')+6:item.find('tif')-5]
        
        dataset = rio.open(folder[0] + "/" + item) # CHANGE THIS TO YOUR OWN FOLDER
        xds = rioxarray.open_rasterio(
                dataset,
                parse_coordinates=True,
                chunks=True,
                masked=False,
            )
        
        clipped = xds[0].rio.clip(jkt_gdf.geometry.apply(mapping), xds.rio.crs, all_touched=True, drop=True, invert=False)
        d_typ = rio.dtypes.get_minimum_dtype(clipped.values)
        
        clipped.values = perform_hack_addition(clipped)
        
        clipped.values = fix_nan_values(clipped.values)
            # clip NAN VALUES
        
        for write_type in wtypes:
            write_to_dataset(destination, title, clipped.values, write_type)
            
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