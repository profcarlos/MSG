import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import os
import glob
from utils import *
import subprocess
import sys

outputDir = 't:\\outros\\Indices_Solo\\'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
exeDir = 'C:\\OSGeo4W\\bin\\'

class main():


	tmpfile1 = outputDir + 'Indice_Medio_Solo.tif'
	output_filename = tmpfile1.replace('.tif', '_go.tif')
	tmpfile2 = outputDir + 'Indice_Medio_Solo' + '_tmpfile2.tif'
	tmpfile3 = outputDir + 'Indice_Medio_Solo' + '_tmpfile3.tif'
	print('\n...Cut file to Goias')
	print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "cubic" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
	subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "bilinear" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
	print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
	subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
	print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + output_filename)
	subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + output_filename)
	delFile(tmpfile2)
	delFile(tmpfile3)
	print('...Finish you win!!')