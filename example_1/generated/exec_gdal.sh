for d in $(find ./reprojection -name "*.tif")
do 
    NAME=$(echo $d | cut -d'.' -f 2 | cut -d'/' -f 3)
    echo "Converting $NAME using GDAL"
    gdal_translate -ot Float32 -of PCRaster -mo PCRASTER_VALUESCALE=VS_SCALAR $d "${NAME}.map"
done


