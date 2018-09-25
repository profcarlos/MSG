import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

from sys import argv
from utils import *
import subprocess
import os
import numpy
from scipy.special import factorial

csv_file   = 'C:\\Users\\carlos.silveira\\Dropbox\\Gabriel\\pontos_coord_v2.csv'
input_dir  = 'T:\\MODIS\\'



exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
gdalIlwisDir = "C:\\Ilwis372\\Extensions\\Geonetcast-Toolbox\\GDAL\\bin\\" 
utilDir = 'C:\\ILWIS372\\Extensions\\GEONETCast-Toolbox\\util\\'
ilwDir = 'C:\\ilwis372\\'
osgeoDir = 'C:\\OSGeo4W\\bin\\'
tmpDir  = 'c:\\tmp\\'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'

class main():
  
    # Using A to NDVI band and B to QUAL band 
    fpar_seca  = '"(B&3==0)*(A>0)*(0.0040+0.9843*A*0.0001)"'
    fpar_chuva = '"(B&3==0)*(A>0)*(0.0564+0.8198*A*0.0001)"'
    for year in range(2000, 2016, 1):
        for doy in range(1, 366, 16):
            file_MOD13 = input_dir + '\\MOD13_250m\\' + str(year) + str('%03d'%doy) + '_MOD13.tif'
            file_FPAR  = input_dir + '\\FPAR_calc\\' + str(year) + str('%03d'%doy) + '_FPAR.tif'
            tmpfile = tmpDir + str(year) + str('%03d'%doy) + '_tmp.tif'
            if(os.path.isfile(file_MOD13) and not os.path.isfile(file_FPAR)):
                if (doy < 91 or doy > 301):
                    fpar_function = fpar_chuva
                else:
                    fpar_function = fpar_seca
                print('process file:' + file_MOD13)
                print('c:\\python34\\python.exe ' + osgeoDir + 'gdal_calc.py -A ' + file_MOD13 + ' --A_band=1 -B ' + file_MOD13 + ' --B_band=4  --type=Float32 --NoDataValue=0 --overwrite --outfile=' + tmpfile + ' --calc='+ fpar_function)
                subprocess.call('c:\\python34\\python.exe ' + osgeoDir + 'gdal_calc.py -A ' + file_MOD13 + ' --A_band=1 -B ' + file_MOD13 + ' --B_band=4 --type=Float32 --NoDataValue=0 --overwrite --outfile=' + tmpfile + ' --calc='+ fpar_function)
                print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile + ' ' + file_FPAR)
                subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile + ' ' + file_FPAR)
                delFile(tmpfile)
                #sys.exit(1)





