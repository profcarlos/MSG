import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import os
import glob
from utils import *
import subprocess
import sys

outputDir = 't:\\outros\\TerraClass\\'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
exeDir = 'C:\\OSGeo4W\\bin\\'

def process_class(x, y, class_msg, CLASS_ORIG):

	'''
	3 km (MSG)
	Size is 211, 204
	Origin = (-53.271811612492662,-12.383945170790270)
	Pixel Size = (0.034840758326260,-0.034840758326260)
	Corner Coordinates:
	Upper Left  ( -53.2718116, -12.3839452) ( 53d16'18.52"W, 12d23' 2.20"S)                                                                 
																																																																																																																																  Lower Left  ( -53.2718116, -19.4914599) ( 53d16'18.52"W, 19d29'29.26"S)
	Upper Right ( -45.9204116, -12.3839452) ( 45d55'13.48"W, 12d23' 2.20"S)
	Lower Right ( -45.9204116, -19.4914599) ( 45d55'13.48"W, 19d29'29.26"S)
	'''

	'''
	Driver: GTiff/GeoTIFF
	Files: T:\OUTROS\TerraClass\TClass_maskVectorGoias.tif
		   T:\OUTROS\TerraClass\TClass_maskVectorGoias.tif.aux.xml
	Size is 27453, 26491
	Coordinate System is:
	GEOGCS["WGS 84",
		DATUM["WGS_1984",
			SPHEROID["WGS 84",6378137,298.257223563,
				AUTHORITY["EPSG","7030"]],
			AUTHORITY["EPSG","6326"]],
		PRIMEM["Greenwich",0],
		UNIT["degree",0.0174532925199433],
		AUTHORITY["EPSG","4326"]]
	Origin = (-53.270596718902418,-12.295126755674074)
	Pixel Size = (0.000276028422046,-0.000276028422046)
	Metadata:
	  AREA_OR_POINT=Area
	Image Structure Metadata:
	  INTERLEAVE=BAND
	Corner Coordinates:
	Upper Left  ( -53.2705967, -12.2951268) ( 53d16'14.15"W, 12d17'42.46"S)
	Lower Left  ( -53.2705967, -19.6073957) ( 53d16'14.15"W, 19d36'26.62"S)
	Upper Right ( -45.6927884, -12.2951268) ( 45d41'34.04"W, 12d17'42.46"S)
	Lower Right ( -45.6927884, -19.6073957) ( 45d41'34.04"W, 19d36'26.62"S)
	Center      ( -49.4816926, -15.9512612) ( 49d28'54.09"W, 15d57' 4.54"S)
	'''

	r3km_py  = -12.383945170790270
	r3km_px  = -53.271811612492662
	r3km_psy = -0.034840758326260
	r3km_psx = 0.034840758326260
	r3km_sy  = 204
	r3km_sx  = 211

	r250m_py  = -12.295126755674074
	r250m_px  = -53.270596718902418
	r250m_psy = -0.000276028422046
	r250m_psx = 0.000276028422046
	r250m_sy  = 26491
	r250m_sx  = 27453

	start_px = round(abs(((r3km_px + r3km_psx*x)-r250m_px)/r250m_psx))
	start_py = round(abs(((r3km_py + r3km_psy*y)-r250m_py)/r250m_psy))
	
	num_class   = 0
	num_outless = 0
	num_pixels = 125*125
	# abs(r3km_psx/r250m_psx - 1) = 125
	for y_class in range(start_py, start_py+125, 1):
		for x_class in range(start_px, start_px+125, 1):
			# Get equal class
			if(CLASS_ORIG[y_class][x_class] == class_msg):
				num_class = num_class + 1
			if(CLASS_ORIG[y_class][x_class] == 128):
				num_outless = num_outless + 1
	if((num_pixels - num_outless) > 0):
		perc_class = num_class/(num_pixels - num_outless)
	else:
		perc_class = 0
	return round(perc_class*100)


def histogram(filename):
	# Histogram  of file
	if(os.path.isfile(filename)):
		DATA = readFileBand(filename, 1)
		HIST = numpy.zeros(15)
		ds = gdal.Open(filename)
		y_band = ds.RasterYSize
		x_band = ds.RasterXSize
		
		for x in range(x_band):
			for y in range(y_band):
				if(DATA[y][x] == 128):
					continue
				else:
					HIST[DATA[y][x]] = HIST[DATA[y][x]] + 1

		print('Histogram of file: %s'%(filename))
		for i in range(15):
			print('%d\t%.0f' %(i, HIST[i]))
	else:
		print('Error to read input file')
	sys.exit()

class main():

	terraClass_orig = outputDir + 'TClass_maskVectorGoias.tif'
	terraClass_3km  = terraClass_orig.replace('.tif', '_go.tif')

	if(not os.path.isfile(terraClass_3km)):
		print('Create file TerraClass 3 km')
		tmpfile1 = terraClass_orig
		tmpfile2 = outputDir + 'TClass_maskVectorGoias' + '_tmpfile2.tif'
		tmpfile3 = outputDir + 'TClass_maskVectorGoias' + '_tmpfile3.tif'
		output_filename = terraClass_3km
		print('\n...Cut file to Goias')
		print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "mode" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
		subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "mode" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
		print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
		subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
		print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + output_filename)
		subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + output_filename)
		delFile(tmpfile2)
		delFile(tmpfile3)
	terraClass_perc = terraClass_3km.replace('.tif', '_go_perc.tif')
	if(not os.path.isfile(terraClass_perc)):
		ds = gdal.Open(terraClass_3km)
		y_MSG = ds.RasterYSize
		x_MSG = ds.RasterXSize
		ds = None
		CLASS_PERC = numpy.zeros((y_MSG, x_MSG))
		CLASS_MSG  = readFileBand(terraClass_3km,  1)
		CLASS_ORIG = readFileBand(terraClass_orig, 1)
		process = 0
		print('..Start process class percent in TerraClass Cerrado...')
		print('..0',end ="")
		for y in range(y_MSG):
			if(int(y/y_MSG*100) >= process + 5):
				process = int(y/y_MSG*100)
				print('..%d'%(process),end ="")
			for x in range(x_MSG):
				if(CLASS_MSG[y][x] != 15):
					CLASS_PERC[y][x] = process_class(x, y, CLASS_MSG[y][x], CLASS_ORIG)
					#if(CLASS_PERC[y][x] > 0.7):
					#	print('CLASS_MSG[%d][%d] = %d : %.2f'%(y, x, CLASS_MSG[y][x], CLASS_PERC[y][x]))
		print('..100 - done')	
		saveRASTERfile(terraClass_perc, terraClass_3km, [CLASS_MSG, CLASS_PERC])
	