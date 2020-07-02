from pcraster import * 

'''
The following setglobaloption functions are used to create a smoother Local Drainage Direction Map. 
Should more global options be necessary, please feel free to put them here.
'''

setglobaloption("lddin")
setglobaloption("lddfill")




def read_maps(src=''):
    
    '''
    This function is used to read and return the 14 .map files that will be used in the PCRaster preprocessing:
        DEM, n, lai, rr, ch, per, thetas1, thetai1, ksat1, landunit, mask, psi1, chanmask, buildings. 
    
    PARAMETERS
    ----------
    
        src: working directory for the folder that contains the maps
        
    OUTPUT
    ------
    
        The 14 maps that will be used for the PCRaster preprocessing
    
    EXAMPLES
    --------
    
    read_maps(src="C:/Users/John Doe/Documents/Maps")
    '''
    
    dem = readmap(src+"dem.map")
    n = readmap(src+"n.map")
    lai = readmap(src+"lai.map")
    rr = readmap(src+"rr.map")
    ch = readmap(src+"ch.map")
    per = readmap(src+"per.map")
    thetas1 = readmap(src+"thetas1.map")
    thetai1 = readmap(src+"thetai1.map")
    ksat1 = readmap(src+"ksat1.map")
    landunits = readmap(src+"landunits.map")
    mask = readmap(src+"mask.map")
    psi1 = readmap(src+"psi1.map")
    chanmask = readmap(src+"chanmask.map")
    buildings = readmap(src+"buildings.map")

    return dem, n, lai, rr, ch, per, thetas1, thetai1, ksat1, landunits, mask, psi1, chanmask, buildings


### === PART 1 Catchment Characteristics ====
def create_catchment(dem, landunits, mask, src = ""):
        
    '''
    This function is used to generate catchment characteristics maps (refer to documentation for more information).
    The output will be located in the user's desired output folder according to the src argument. If left void, the maps will be outputted
    in the same directory as this python script.
    
    PARAMETERS
    ----------
    
        Maps: Digital Elevation Map, Land Unit, Mask, and Fraction of cell covered by buildings
        
        Desired output location
        
    OUTPUT: saved maps
    ------
    
        Digital Elevation (masked), Land Unit (masked), and Fraction of cell covered by buildings (masked),
        Local Drainage Direction, Location of outlet and sub outlets, Slope Gradient
    
    EXAMPLES
    --------
    
    create_catchment(dem, landunits, mask, src="C:/Users/John Doe/Documents/Maps")
    '''
    # mask DEM
    dem = dem * mask
    # Generate Local Drainage Direction, Outlet, Gradient maps
    dem = lddcreatedem(dem,999999,999999,999999,999999)
    lddmap = lddcreate(dem, 1e20, 1e20, 1e20, 1e20)
    outlet = pit(lddmap)
    grad = max(sin(atan(slope(dem))), 0.001)
    
    # output maps
    report(dem, src+"dem.map")
    report(lddmap, src+"ldd.map")
    report(landunits*mask, src+"landunits.map")
    report(outlet, src+"outlet.map")
    report(grad, src+"grad.map")

### === PART 2 Vegetation and Soil Surface ====
def vegetation(n, lai, rr, ch, per, mask, buildings, src = ""):
    '''
    
    This function is used to generate vegetation and soil surface maps (refer to documentation for more information).
    The output will be located in the user's desired output folder according to the src argument. If left void, the maps will be outputted
    in the same directory as this python script.
    
    PARAMETERS
    ----------
    
        Maps: Manning's n, Leaf Area Index, Random Roughness, Crop Height, Fraction of Soil Covered by Vegetation,
        Mask, Fraction of Cell Covered by Buildings
        
        Desired output location
        
    OUTPUT: saved maps
    ------
    
        Manning's n (masked), Leaf Area Index (masked), Random Roughness (masked), Crop Height (masked), Fraction of Soil Covered by Vegetation (masked),
        Fraction of Cell Covered by Buildings (masked), Fraction of Cell with Hard Surface, Fraction of cell covered by stones, Fraction of ell covered by crust,
        Fraction of cell that is compacted, Fraction of cell with hard surface 

    
    EXAMPLES
    --------
    
    vegetation(n, lai, rr, ch, per, mask, buildings, src="C:/Users/John Doe/Documents/Maps")
    '''
    # create a map full of zeros 
    zeromask = 0 *mask
    
    # output maps
    report(n*mask, src+"n.map")
    report(lai*mask, src+"lai.map")
    report(rr*mask, src+"rr.map")
    report(ch*mask, src+"ch.map")
    report(per*mask, src+"per.map")
    report(buildings*mask, src+"buildings.map")
    report(zeromask, src+"drumstore.map")
    report(zeromask, "crust.map")
    report(zeromask, "stone.map")
    report(zeromask, "comp.map")
    report(zeromask, "bufferid.map")
    report(zeromask, "buffervol.map")

### === PART 3 Channel Data ===
def channels(mask, chanmask, dem, src = ""):
    
    '''
    This function is used to generate channels maps (refer to documentation for more information).
    The output will be located in the user's desired output folder according to the src argument. If left void, the maps will be outputted
    in the same directory as this python script.
    The channel parameters used here are default parameters as used by Dr. Victor Jetten in the GANSPOEL Project. Should users feel the need
    to change the parameters, please feel free to do so in the section below
    
    PARAMETERS
    ----------
    
        Maps: Mask, Channel Mask, Digital Elevation Map
        
        Desired Output Location
        
    OUTPUT: saved maps
    ------
    
        Channel network, Local drainage direction of channel network, Channel gradient, Manning’s n for channel network, Width of channel,
        Channel cross section shape, Saturated hydraulic conductivity for channel network

    EXAMPLES
    --------
    
    channels(mask, chanmask, dem, src="C:/Users/John Doe/Documents/Maps")
    '''
    
    # Channel parameters. Change as needed.
    chancoh = 10
    chanman = 0.2
    chanside = 0
    chanwidth = 2
    roadwidth = 6
    chanKsat = 20

    # generate necessary output maps
    lddchan = lddcreate(chanmask*dem, 1e20,1e20,1e20,1e20)
    changrad = max(0.001, sin(atan(slope(chanmask*dem))))
    chancoh = chanmask * scalar(chancoh)
    chanman = chanmask * scalar(chanman)
    chanside = chanmask * scalar(chanside)
    chanwidt = chanmask * scalar(chanwidth)
    hmxinit = mask * 0
    chanksat = chanmask *chanKsat

    # output maps
    report(lddchan, src+"lddchan.map")
    report(changrad, src+"changrad.map")
    report(chancoh, src+"chancoh.map")
    report(chanman, src+"chanman.map")
    report(chanside, src+"chanside.map")
    report(chanwidt, src+"chanwidt.map")
    report(hmxinit, src+"hmxinit.map")
    report(chanksat, src+"chanksat.map")


#### === PART 4 Green and Ampt (G&A) Layer ====
def GA(mask, ksat1, thetas1, thetai1, psi1, src = ""):
    
    '''
    This function is used to generate Green & Ampt maps (refer to documentation for more information).
    The output will be located in the user's desired output folder according to the src argument. If left void, the maps will be outputted
    in the same directory as this python script.
    The parameters used for the second G&A Layer here are default parameters as used by Dr. Victor Jetten in the GANSPOEL Project. Should users feel the need
    to change the parameters, please feel free to do so in the section below
    
    PARAMETERS
    ----------
    
        First G&A Layer Maps of: Saturated hydraulic conductivity, Saturated volume soil moisture content (porosity), Initial volumetric soil moisture content
        Soil water tension at wetting font, Soil depth
 
        
        Desired Output Location

        
        
    OUTPUT: saved maps
    ------
    
        First and Second G&A Layer Maps of: Saturated hydraulic conductivity, Saturated volume soil moisture content (porosity), Initial volumetric soil moisture content
        Soil water tension at wetting font, Soil depth
        
        

    EXAMPLES
    --------
    
    GA(mask, ksat1, thetas1, thetai1, psi1, src="C:/Users/John Doe/Documents/Maps")
    '''
    
    # second G&A Layer Parameters. Change as needed.
    soildep1 = scalar(200)
    soildep2 = scalar(1000)
    # the parameters defined here refer to the ratio between the second layer and first layer parameters
    # ksat_ratio of 0.9 means that the second layer of Ksat should have 0.9 the value of the first layer
    ksat_ratio = 0.9
    thetas_ratio = 1
    thetai_ratio = 1
    psi_ratio = 1
    
    # generte second G&A Layer
    ksat2 = ksat_ratio*ksat1*mask
    thetas2 = thetas_ratio*thetas1*mask
    thetai2 = thetai_ratio*thetai1*mask
    psi2 = psi_ratio*psi1*mask
    
    # output maps
    report(ksat1*mask, "ksat1.map")
    report(thetas1*mask, "thetas1.map")
    report(thetai1*mask, "thetai1.map")
    report(psi1*mask, "psi1.map")
    report(soildep1*mask, "soildep1.map")

    report(ksat2*mask, "ksat2.map")
    report(thetas2*mask, "thetas2.map")
    report(thetai2*mask, "thetai2.map")
    report(psi2*mask, "psi2.map")
    report(soildep2*mask, "soildep2.map")
    
"""
Main executable
"""

if __name__ == '__main__':
    # desired folder locations for each function
    base_maps_src = ""
    catchment_src = ""
    vegetation_src = ""
    channels_src = ""
    GA_src = ""
    
    dem, n, lai, rr, ch, per, thetas1, thetai1, ksat1, landunits, mask, psi1, chanmask, buildings = read_maps(base_maps_src)
    create_catchment(dem, landunits, mask, src = base_maps_src)
    vegetation(n, lai, rr, ch, per, mask, buildings, src = vegetation_src)
    channels(mask, chanmask, dem, src = channels_src)
    GA(mask, ksat1, thetas1, thetai1, psi1, src = GA_src)