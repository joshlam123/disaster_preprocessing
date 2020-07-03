from subprocess import Popen, PIPE
from .common import move_files, chg_to_folder
import os, re, pandas as pd

def gdal_set_up_cubic_spline(path:str, curr_folder:str=os.getcwd(), src:str='dem.tif', dest:str='dem_cubic.tif'):
    '''
    A function to use the built in gdal_warp and use the cubic spline algorithm to do interpolation for grid cells. This function
    invokes the command line using os.system and runs the command. It does this by changing to the destination directory, and 
    running the gdalwarp command there. 

    Hence, it requires GDAL to be installed and in the environemnt
    path variable.

    For a full list of options you can use to modify the command, go to https://gdal.org/programs/gdalwarp.html

    Parameters
    ----------

        path: the path that contains your .tif files to convert. (Required)

        curr_folder: your current working directory. If empty, it uses os.getcwd() (Optional)

        src: the name of the source file. If empty, it uses 'dem.tif' (Optional)
    
        dest: the name of the destination file. If empty, it uses 'dem_cubic.tif' (Optional)

    Output
    ------

        None
    
    Examples
    --------

    set_up_cubic_spline('../data/experiments/jakarta/exp5/', '../code/' 'landunits.tif', 'landunits_cubic.tif')

    '''

    if os.getcwd() != path:
        revert_to_folder(path)
    
    cmd = f"gdalwarp -ot Float32 -r cubicspline -order 3 -et 10e-3 -tr 30 30 {src} {dest}"
    print(cmd)
    os.system(cmd)
    
    print(f"successfully converted file {src} to cubic spline equivalent at {dest}")
    
    revert_to_folder(curr_folder)

def gdal_change_resolution(path, src, dest, res=30):
    '''
    A function to use the built in gdal_translate function to change from your desired resolution to your target resolution. This function
    invokes the command line using os.system and runs the command. It does this by changing to the destination directory, anbd 

    Hence, it requires GDAL to be installed and in the environemnt path variable.

    For a full list of options you can use to modify the command, go to https://gdal.org/programs/gdaltranslate.html

    NOTE: YOU HAVE TO MAKE SURE THE RESOLUTION YOU WANT IS IN THE CORRECT FORMAT AS YOUR SOURCE FILE (e.g. if you source file
    is in EPSG4326 and you want the target resolution to be 30m x 30m, you will need to do something like 30/111132 ~ 0.000026 degrees)

    Parameters
    ----------

        path: the path that contains your .tif files to convert. (Required)

        src: the name of the source file. If empty, it uses 'dem.tif' (Required)
    
        dest: the name of the destination file. If empty, it uses 'dem_cubic.tif' (Required)

        res: the target resolution you want to change your files to 

    Output
    ------

        None
    
    Examples
    --------

    set_up_cubic_spline('../data/experiments/jakarta/exp5/', '../code/' 'landunits.tif', 'landunits_cubic.tif')

    '''

    if os.getcwd() != path:
        chg_to_folder(path)
    tr = res/111132
    
    cmd = f"gdal_translate -ot Float32 -tr {tr} {tr} {src} {dest}"
    print(cmd)
    os.system(cmd)
    
    chg_to_folder(curr_folder)


def get_gdal_info(files:list, folder:str, curr_folder:str=os.getcwd()):
    '''
    A function to extract the text from metadata returned from get_gdal_info.
    
    Parameters
    ----------
    
        files: a list of files you want to get the gdal info for. They must all end with the .map extension (Required)

        folder: a folder directory that you want to look in for the files (Required)

        curr_folder: your current working directory. If not provided, it will default to os.getcwd() (Optional)
        
    Output
    ------
        
        meta_df: a pandas DataFrame consisting of the important columns: name, size, origin, pixel size, upper left, 
        lower left, upper right, lower right, center, block, min_max, and nodata values
    
    Example
    -------

    get_gdal_info(["mask.map", "rr.map", "chanmask.map"], "../data/experiments/jakarta/maps/", os.getcwd())

    '''
    return_msg = dict()
    if os.getcwd() != folder:
        chg_to_folder(folder)
    for item in files:
        p = Popen(["gdalinfo", item], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        return_msg[item[:item.find('.map')]] = output.decode("utf-8") 
    
    chg_to_folder(curr_folder)
    
    return return_msg

def extract_gdal_meta(gdal_return_msgs:dict):
    '''
    A function to extract the text from metadata returned from get_gdal_info().
    
    Parameters
    ----------
    
        gdal_return_msgs: a dictionary containing the text metadata in the format {filename:[strings of metadata]} (Required)
        
    Output
    ------
        
        meta_df: a pandas DataFrame consisting of the important columns: name, size, origin, pixel size, upper left, 
        lower left, upper right, lower right, center, block, min_max, and nodata values

    Examples
    --------

        return_msgs = get_gdal_info(["mask.map", "rr.map", "chanmask.map"], "../data/experiments/jakarta/maps/", os.getcwd())
        extract_gdal_meta(return_msgs)
    
    '''
    extracted_msgs = dict()
    for k, v in gdal_return_msgs.items():
        extracted_msgs[k] = [re.sub('^/(?!nan)[A-Za-z:,=/"()_\[\]]', ' ', i).strip() for i in v.split("\n")]

    filtered_msgs = {k:([j for j in v if len(j) > 2 or 'nan' in j]) for k,v in extracted_msgs.items()}

    col_headers = ['name','size', 'origin', 'pixel size', 'upper left', 'lower left', 'upper right', 'lower right', 'center',\
              'block', 'min_max', 'nodata']

    for k,val in filtered_msgs.items():
        # ////// IMPORTANT ///////
        # there are two types of messages.
        # lengths less than 12 usually are concise gdal generated metadata
        # lengths of 12 are usually self generated metadata 
        # 12 is more of a heuristic number from manual inspection, but it is consistent across those generated by my code
        # in the EPSG900913 format
        if len(val) > 12:
            if 'Size' in val[2]:
                new_val = [val[1]]+[val[2]]+val[-13:-11]+val[-8:]
            else:
                new_val = [val[1]]+[val[3]]+val[-14:-12]+val[-8:]

            filtered_msgs[k] = new_val


    meta_df = pd.DataFrame.from_records([v for k,v in filtered_msgs.items()], columns=col_headers)

    # to apply regex substitution on columns which contain the words with their header inside
    for v in col_headers:
        meta_df[v] = meta_df[v].apply(lambda x: re.sub('[()=,]', ' ', x).strip().lower().replace(v, ''))

    # custom substitution based on the values which I know from the column specifically
    meta_df['min_max'] = meta_df['min_max'].apply(lambda x: x.replace('min', '').replace('max', ''))
    meta_df['nodata'] = meta_df['nodata'].apply(lambda x: x.replace('value', ''))
    meta_df['size'] = meta_df['size'].apply(lambda x: x.replace('is', ''))
    meta_df['name'] = meta_df['name'].apply(lambda x: x.replace('files:', '').replace('.map', '').strip())

    return meta_df


def gdal_trans(folder:str, files:str, dest_folder:str, nodata_repl:float=10e-5, curr_folder:str=os.getcwd()):   
    '''
    A function to use the built in gdal_translate function to change a list of files into the .map format for PCRaster processin. 
    This function invokes the command line using os.system and runs the command. It does this by changing to the destination directory, anbd 
    
    At the end of translation, it calls the move_files() function from preprocessing.common and moves all the files to the destination

    Hence, it requires GDAL to be installed and in the environemnt path variable.

    For a full list of options you can use to modify the command, go to https://gdal.org/programs/gdaltranslate.html

    Parameters
    ----------
        
        folder: the directory name you entered in as dest (Required)

        files: a list of file names matching incl, excl, and ext. (Required)
        
        dest_folder: a destination folder where you want to save your files (Required)

        nodata_repl: the nodata value for the raster. If empty, uses 10e-5  (Optiona)

        curr_folder: your current working directory. If empty, it uses os.getcwd() (Optional)

    Output
    ------

        None
    
    Examples
    --------

    set_up_cubic_spline('../data/experiments/jakarta/exp5/', '../code/' 'landunits.tif', 'landunits_cubic.tif')

    '''

    new_file = list()
    if os.getcwd() != folder:
        chg_to_folder(folder)
    
    for idx, item in enumerate(files):

        title = item[:item.find('.tif')]
        print("DATASET: ", title)

        dest = title + '.map'
        new_file.append(dest)

        # if title == "dem":
        #     cmd = f"gdal_translate -ot Float32 -of PCRaster -mo PCRASTER_VALUESCALE=VS_SCALAR {item} {dest}"
        # else:
        #     
        cmd = f"gdal_translate -ot Float32 -of PCRaster -mo PCRASTER_VALUESCALE=VS_SCALAR {item} {dest}"
        print(cmd)
        os.system(cmd)
    
        print(f"successfully converted {title} to map file at {dest}")
    
    chg_to_folder(curr_folder)
    
    move_files(new_file, folder, dest_folder)

