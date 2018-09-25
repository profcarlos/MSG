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

	Size is 26627, 25525
	Coordinate System is:
	GEOGCS["GCS_SIRGAS_2000",
	    DATUM["Sistema_de_Referencia_Geocentrico_para_America_del_Sur_2000",
	        SPHEROID["GRS_1980",6378137,298.257222101]],
	    PRIMEM["Greenwich",0],
	    UNIT["degree",0.0174532925199433]]
	Origin = (-53.265785735546764,-12.310689655732755)
	Pixel Size = (0.000282251956231,-0.000282251956231)
	Metadata:
	  AREA_OR_POINT=Area
	Image Structure Metadata:
	  INTERLEAVE=BAND
	Corner Coordinates:
	Upper Left  ( -53.2657857, -12.3106897) ( 53d15'56.83"W, 12d18'38.48"S)
	Lower Left  ( -53.2657857, -19.5151708) ( 53d15'56.83"W, 19d30'54.62"S)
	Upper Right ( -45.7502629, -12.3106897) ( 45d45' 0.95"W, 12d18'38.48"S)
	Lower Right ( -45.7502629, -19.5151708) ( 45d45' 0.95"W, 19d30'54.62"S)
	Center      ( -49.5080243, -15.9129302) ( 49d30'28.89"W, 15d54'46.55"S)

	'''

	r3km_py  = -12.383945170790270
	r3km_px  = -53.271811612492662
	r3km_psy = -0.034840758326260
	r3km_psx = 0.034840758326260
	r3km_sy  = 204
	r3km_sx  = 211

	r250m_py  = -12.310689655732755
	r250m_px  = -53.265785735546764
	r250m_psy = -0.000282251956231
	r250m_psx = 0.000282251956231
	r250m_sy  = 25525
	r250m_sx  = 26627

	start_px = round(abs(((r3km_px + r3km_psx*x)-r250m_px)/r250m_psx))
	start_py = round(abs(((r3km_py + r3km_psy*y)-r250m_py)/r250m_psy))
	
	num_class_agro = 0
	num_class_past = 0
	num_class_natu = 0
	num_class_othe = 0
	perc_agro = 0
	perc_past = 0
	perc_natu = 0
	perc_othe = 0
	num_outless = 0
	num_pixels = 125*125
	# abs(r3km_psx/r250m_psx - 1) = 125
	for y_class in range(start_py, start_py+125, 1):
		for x_class in range(start_px, start_px+125, 1):
			# Get agriculture
			if(CLASS_ORIG[y_class][x_class] == 1 or CLASS_ORIG[y_class][x_class] == 2):
				num_class_agro = num_class_agro + 1
			elif(CLASS_ORIG[y_class][x_class] == 5):
				num_class_natu = num_class_natu + 1
			elif(CLASS_ORIG[y_class][x_class] == 11):
				num_class_past = num_class_past + 1
			elif(CLASS_ORIG[y_class][x_class] >= 1 and CLASS_ORIG[y_class][x_class] <= 14):
				num_class_othe = num_class_othe + 1
			else:
				num_outless = num_outless + 1
	if((num_pixels - num_outless) > 0):
		perc_agro = round(100*num_class_agro/(num_pixels - num_outless))
		perc_past = round(100*num_class_past/(num_pixels - num_outless))
		perc_natu = round(100*num_class_natu/(num_pixels - num_outless))
		perc_othe = round(100*num_class_othe/(num_pixels - num_outless))
	#sample = [class_msg, num_class_agro, num_class_natu, num_class_past, num_class_othe, perc_agro, perc_past, perc_natu, perc_othe, (num_pixels - num_outless)]
	#print(sample)
	return [perc_agro, perc_past, perc_natu, perc_othe]


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
	terraClass_perc = terraClass_3km.replace('.tif', '_go_perc_class.tif')
	if(not os.path.isfile(terraClass_perc)):
		ds = gdal.Open(terraClass_3km)
		y_MSG = ds.RasterYSize
		x_MSG = ds.RasterXSize
		ds = None
		CLASS_PERC = numpy.zeros((4, y_MSG, x_MSG))
		CLASS_MSG  = readFileBand(terraClass_3km,  1)
		CLASS_ORIG = readFileBand(terraClass_orig, 1)
		process = 0
		print('..Start process class percent in TerraClass Cerrado...')
		print('..0')
		for y in range(y_MSG):
			if(int(y/y_MSG*100) >= process + 5):
				process = int(y/y_MSG*100)
				print('..%d'%(process))
			for x in range(x_MSG):
				if(CLASS_MSG[y][x] != 15):
					[perc_agro, perc_past, perc_natu, perc_othe] = process_class(x, y, CLASS_MSG[y][x], CLASS_ORIG)
					CLASS_PERC[0][y][x] = perc_agro
					CLASS_PERC[1][y][x] = perc_past
					CLASS_PERC[2][y][x] = perc_natu
					CLASS_PERC[3][y][x] = perc_othe

					#if(CLASS_PERC[y][x] > 0.7):
					#	print('CLASS_MSG[%d][%d] = %d : %.2f'%(y, x, CLASS_MSG[y][x], CLASS_PERC[y][x]))
		print('..100 - done')	
		saveRASTERfile(terraClass_perc, terraClass_3km, CLASS_PERC)
	