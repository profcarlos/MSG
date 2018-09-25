import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import os
import glob
from utils import *
import subprocess
import sys

#path = 'I:\\VETOR\\TRMM\\'
path = 'T:\\OUTROS\\TRMM_OTTO\\TRMM_DATA\\'
outputDir = 't:\\outros\\TRMM_OTTO\\'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
exeDir = 'C:\\OSGeo4W\\bin\\'

class main():

	for n_month in range(1,13,1):
		month = str('%02d'%n_month)
		output_filename = outputDir + 'pa_br_trmm_3000_M' + month + '.tif'
		filename = path + 'pa_br_trmm_3000_*_M' + month + '_lapig.tif'
		if(not os.path.isfile(output_filename)):
			files = glob.glob(filename)
			if(len(files) == 0):
				print('Error to read data: %s'%(filename))
				sys.exit()
			print('\n...Files to process month: %s [%d]' %(month, len(files)))
			start_data = True
			for file in files:
				print('Reading file: %s' %(file))
				if(start_data == True):
					DATA = readFileBand(file, 1)
					start_data = False
				else:
					DATA = DATA + readFileBand(file, 1)
			DATA = DATA/len(files)
			saveRASTERfile(output_filename, file, DATA)
		tmpfile1 = output_filename
		output_filename = output_filename.replace('br', 'go')
		if(not os.path.isfile(output_filename)):
			tmpfile2 = outputDir + 'pa_br_trmm_3000_M' + month + '_tmpfile2.tif'
			tmpfile3 = outputDir + 'pa_br_trmm_3000_M' + month + '_tmpfile3.tif'
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

	#c:\osgeo4w\bin\gdal_translate -b 1 -b 2 -b 3 -mask "none" "input_file.tif" "output_file.tif"
