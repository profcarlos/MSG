import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

from sys import argv
from utils import *
import subprocess
import os
import numpy
from scipy.special import factorial

input_dir  = 'I:\\RASTER\\MODIS\\MOD17A2H\\'
output_dir = 'T:\\OUTROS\\STAT_GAB\\'
exeDir = 'C:\\osgeo4w\\bin\\'
osgeoDir = 'C:\\OSGeo4W\\bin\\'
modiscutfile = output_dir + 'mod_goias.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
tmpDir = 'c:\\tmp\\'

class main():
  

    for year in range(2014, 2015, 1):
        for doy in range(1, 366, 8):
            file_MOD17 = input_dir + 'MOD17A2H_006_Gpp_500m_' + str(year) + str('%03d'%doy) + '_brasil.tif'
            tmpfile1 = tmpDir + str(year) + str('%03d'%doy) + '_tmp1.tif'
            tmpfile2 = tmpDir + str(year) + str('%03d'%doy) + '_tmp2.tif'
            tmpfile3 = tmpDir + str(year) + str('%03d'%doy) + '_tmp3.tif'            
            outfile = output_dir + 'GPP_LAPIG\\'+ str(year) + str('%03d'%doy) + '_gpp.tif' 
            if(os.path.isfile(file_MOD17)):
                print('...cut input file')
                print(exeDir + 'gdalwarp -tr 0.001905476965202 0.001905476965202 -r "cubic" -cutline ' + modiscutfile + ' -crop_to_cutline ' + file_MOD17 + ' ' + tmpfile1) 
                subprocess.call(exeDir + 'gdalwarp -tr 0.001905476965202 0.001905476965202 -r "cubic" -cutline ' + modiscutfile + ' -crop_to_cutline ' + file_MOD17 + ' ' + tmpfile1) 

                print('...use modis cutfile')
                print(exeDir + 'gdal_translate.exe -projwin -53.2524035 -12.3932791 -45.9010733 -19.5007081 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
                subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2524035 -12.3932791 -45.9010733 -19.5007081 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)

                print('...process outfile')
                print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outfile)
                subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outfile)
 
                delFile(tmpfile1)
                delFile(tmpfile2)
                delFile(tmpfile3)
                #sys.exit()





