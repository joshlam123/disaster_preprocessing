from ee import Image
from ee.batch import Task, Export
from disaster_preprocessing.code import common
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import fiona, json, ee
import geopandas as gpd, osmnx as ox


# Initialize the Earth Engine module.
ee.Initialize()

'''
Setup earth engine here: https://github.com/google/earthengine-api/tree/master/python
'''

def extract_countries(mapping:str, path:str, filename:str = None):
    '''
    This module is used to extract maps from the earth engine catalogue. For a list of the available datasets,
    go to: https://developers.google.com/earth-engine/datasets/catalog and put in the tag that corresponds to the map.
    For example, SRTM Digital Elevation Data 30m (https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003)
    corresponds to 'USGS/SRTMGL1_003'

    Parameters
    ----------

        mapping:
        
        path:
        
        filename:


    Output
    ------

        None
    
    Examples
    --------

    
    '''

    for k,v in mapping.items():
        for item in v:
            FILNAME = filename or item.replace(" ", "_")
            # if not isinstance(item, list):

            gdf_file = "".join(["{}".format(k),"_","{}".format(FILNAME),".geojson"])

            # else:

            #     gdf_file = "".join(["{}".format(k),"_","{}".format(),".geojson"])                

            gdf_path = path + k.lower()

            if checkExist(gdf_file, gdf_path):

                place = [f'{k}, {item}']
                print(place)


                create_folder(gdf_path)

                file = gdf_path+"/"+gdf_file
                print(file)

                place_gdf = ox.gdf_from_places(item)

                if isinstance(item, list):
                    
                    polygons = [i for i in place_gdf['geometry']]
                    place_gdf = gpd.GeoSeries(cascaded_union(polygons))
                
                else:
                    if len(place_gdf) > 1: 

                        polygons = [Polygon(list(zip(*i.exterior.coords.xy))) for i in jkt_gdf.geometry[0]]
                        place_gdf = pd.DataFrame({"geometry":polygons[0]}, index=[i for i in range(1)])


                place_gdf.to_file(file, driver="GeoJSON")

                print(f"file {gdf_file} saved at {gdf_path}!")
            else:
                print(f"file {gdf_file} exists at {gdf_path}!")
    return None

def create_ee_maps(file_folder, map_type:str="USGS/SRTMGL1_003", vec_scale:int=30):
    '''
    This module is used to extract maps from the earth engine catalogue. For a list of the available datasets,
    go to: https://developers.google.com/earth-engine/datasets/catalog and put in the tag that corresponds to the map.
    For example, SRTM Digital Elevation Data 30m (https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003)
    corresponds to 'USGS/SRTMGL1_003'

    Parameters
    ----------

        file_folder: a list that is a Shapely Polygon. (Required)

        map_type: a string that corresponds to map type. (e.g. 'USGS/SRTMGL1_003'). Default value = 'USGS/SRTMGL1_003' (Optional)

        vec_scale: an integer in meters that corresponds to the length to pixel ratio. Default value = 30 (Optional)

    Output
    ------

        export_task: a list containing all tasks run recently
    
    Examples
    --------

    
    '''

    for item in file_folder:
        print(item)
        for enum in item:
            print(enum)

            name = enum[enum.find("6/")+2:]
            ctry = name[:name.find("/")]

            file_name = name[name.find("/")+1:name.find(".geo")]

            j_gdf = gpd.read_file(enum)
            print(j_gdf.head())
            if type(j_gdf.geometry[0]) == Polygon:
                polygons = [list(zip(*j_gdf.geometry[0].exterior.coords.xy))]
            else:
                polygons = [list(zip(*i.exterior.coords.xy)) for i in j_gdf.geometry[0]]


            aoi = ee.Geometry.Polygon(polygons)

            img = extract_from_earth_engine(aoi, map_type=map_type, vec_scale=vec_scale)
        

            export_in_earth_engine(aoi, img)

def extract_from_earth_engine(aoi:list, map_type:str, vec_scale:float = 30):
    '''
    This module is used to extract maps from the earth engine catalogue. For a list of the available datasets,
    go to: https://developers.google.com/earth-engine/datasets/catalog and put in the tag that corresponds to the map.
    For example, SRTM Digital Elevation Data 30m (https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003)
    corresponds to 'USGS/SRTMGL1_003'

    Parameters
    ----------

        aoi : a list that is a Shapely Polygon. (Required)

        map_type: a string that corresponds to map type. (e.g. 'USGS/SRTMGL1_003'). Default value = 'USGS/SRTMGL1_003' (Optional)

        vec_scale: an integer in meters that corresponds to the length to pixel ratio. Default value = 30 (Optional)

    Output
    ------

        export_task: a list containing all tasks run recently
    
    Examples
    --------

    
    '''
    map_type = map_type
    VECTORIZATION_SCALE = vec_scale
    
    aoi = ee.Geometry.Polygon(aoi)
    # image = ee.Image('ESA/GLOBCOVER_L4_200901_200912_V2_3').clip(aoi)
    image = Image(map_type).clip(aoi)
    img = image.select('elevation')

    return img

def check_all_tasks():
    '''
    . No arguments required.

    Output
    ------

        active_tasks: 
    
    Examples
    --------
        
    '''

    tasks = ee.batch.Task.list()
    active_tasks = [ee.batch.Task.active(i) for i in tasks]
    print(active_tasks)

    return active_tasks

def cancel_tasks(active_tasks:list):
    '''
    

    Parameters
    ----------

        active_tasks: 

    Output
    ------

        export_task: a list containing all tasks run recently
    
    Examples
    --------
        
    '''
    for idx, item in enumerate(active_tasks):
            if item == True:
                ee.batch.Task.cancel(tasks[idx])
    return ee.batch.Task.list()

def get_ctry_list(file:str, path:str):
    '''
    

    Parameters
    ----------

        active_tasks: 

    Output
    ------

        export_task: a list containing all tasks run recently
    
    Examples
    --------
        
    '''

    # ctry_city_list = dict()
    # ctry_city_list = {"Cambodia":["Battambang"], "Indonesia": ["Sumbawa","Semarang"], "Laos": ["Louangphabang"], \
    #             "Malaysia": ["Terengganu"], "Myanmar": ["Amarapura"], "Philippines": ["Caranga"], "Thailand": ["Rayong","Pathumthani"], \
    #             "Vietnam": ["Hue"]}
    with open(path+file, 'r') as f:
        ctry_city_list = json.load(f)

    for k,v in ctry_city_list.items():

        if checkExist(k, path):
            with open(path+file, 'r') as f:
                ctry_map = json.load(f)
        else:
            with open(path+file, 'w', encoding='utf-8') as f:
                json.dump(ctry_city_list, f, ensure_ascii=False, indent=4)
    
    ctry_map = ctry_city_list
    return ctry_map

def export_in_earth_engine(aoi, img:Image, filename:str="dem_jakarta", dest_folder:str='/content/drive/My Drive/Capstone/', \
                            export_format:str='GeoTIFF'):
    '''
    This module is used to pass the vectorisation scale

    Parameters
    ----------

        img: a google EEImage instance. (Required)
        
        filename: the name of the file you want to save. (Required)
        
        dest_folder: a destination folder on your google drive. (Note earth engine must have been authorized to run on your 
        pc. (e.g. '/content/drive/My Drive'). Default value = '/content/drive/My Drive' (Optional)

        export_format: The string file format to which the image is exported. Currently only 'GeoTIFF' and 'TFRecord' are supported, 
        defaults to 'GeoTIFF'. (https://github.com/google/earthengine-api/blob/master/python/ee/batch.py) (Optional)

    Output
    ------

        export_task: a list containing all tasks run recently
    
    Examples
    --------
        


    
    '''

    export_task = Export.image.toDrive(**{
          'image': img,
          'description': filename,
          'scale': 30,
          'folder': dest_folder,
          'fileFormat': export_format,
          'region': aoi,
          'formatOptions': {
            'cloudOptimized': True
          }
        })
    Task.start(export_task)
    
    return export_task

def get_task_info():
    '''
    This function returns the current running tasks on google earth engine. No arguments required.

    Output
    ------

        export_task: a list containing all tasks run recently
    
    Examples
    --------
        
    '''

    return ee.batch.Task.list()
