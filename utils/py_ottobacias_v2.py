import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import os
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import datetime
import time
import glob
import numpy
from utils import *
import subprocess
import sys

outputDir = 't:\\outros\\Ottobacias\\'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
exeDir = 'C:\\OSGeo4W\\bin\\'

def saveRASTERfile(filename, fileMask, DATA):
	num_bands = len(DATA)
	# For DATA with one band only
	if(num_bands > 30):
		num_bands = 1
	for num_band in range (num_bands):  
		# Save georeference file
		if(num_band == 0):
			ds = gdal.Open(fileMask)
			band = ds.GetRasterBand(1)
			driver = gdal.GetDriverByName("GTiff")
			#print(band.DataType)
			if(num_bands > 1):
				dataType = type(DATA[0][0][0])
			else:
				dataType = type(DATA[0][0])
			# Test dataType if mask is different of DATA
			if((dataType is numpy.float32 or dataType is numpy.float64) and band.DataType <= 5 ):
				dataType = gdal.GDT_Float32
			else:
				dataType = band.DataType

			dsOut  = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, num_bands, dataType)
			CopyDatasetInfo(ds,dsOut)
			bandOut = dsOut.GetRasterBand(num_band+1)
			if(num_bands == 1):
				BandWriteArray(bandOut, DATA)
			else:    
				BandWriteArray(bandOut, DATA[num_band])
		else:
			bandOut = dsOut.GetRasterBand(num_band+1)
			try:
				BandWriteArray(bandOut, DATA[num_band])
			except:
				print('num_band: %d'%(num_band))
				return False
	return True



def readFileBand(file, band):
	if(os.path.isfile(file)):
		try:
			ds = gdal.Open(file)
		except:
			return []
		try:
			bd = ds.GetRasterBand(band)
		except:
			return []
		#print('Read file: %s Band: %d' %(file, band))
		return BandReadAsArray(bd)
	else:
		return []

def histogram(filename, interv, mult):
	# Histogram  of file
	if(os.path.isfile(filename)):
		DATA = readFileBand(filename, 1)
		HIST = numpy.zeros(interv)
		y_band = len(DATA)
		x_band = len(DATA[0])
#		print(x_band)
#		print(y_band)
#		print(range(x_band))
#		print(range(y_band))
		for x in range(x_band):
			for y in range(y_band):
				if(DATA[y][x] == 128 or DATA[y][x] == 0):
					continue
				else:
					#print(DATA[y][x])
					HIST[int(DATA[y][x])] = HIST[int(DATA[y][x])] + 1

		print('Histogram of file: %s'%(filename))
		for i in range(interv):
			print('%d\t%.0f' %(i, HIST[i]))
	else:
		print('Error to read input file')

def delFile(file):
	#print( '..delFile: %s' % str(glob.glob(file)))
	#print('Remove files:')
	for x in glob.glob(file):
		try:
			os.remove(x)
			#print(x)
		except: 
			print('Error in del file: %s' % str(glob.glob(file)))

def inc_date_doy(date, int_doy):
	for inc in range(int_doy):
		date = incDate(date)
	return date

class main():

	'''
	
	# Create ottobacias in 3km
	filename = 'T:\\OUTROS\\OTTOBACIAS\\ottobacia_nivel5.tif'

	print('Create file Ottobacias 3 km')
	tmpfile1 = filename
	tmpfile2 = outputDir + 'Ottobacias' + '_tmpfile2.tif'
	tmpfile3 = outputDir + 'Ottobacias' + '_tmpfile3.tif'
	tmpfile4 = outputDir + 'Ottobacias' + '_tmpfile4.tif'
	output_filename = filename.replace('.tif','_3km.tif')
	print('\n...Cut file to Goias')
	print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "mode" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
	subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "mode" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
	print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
	subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
	print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + output_filename)
	subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + tmpfile4)
	DATA = readFileBand(tmpfile4, 1)
	# get data only
	DATA[DATA>400] = 0
	# saveRASTERfile(filename, fileMask, DATA)
	saveRASTERfile(output_filename, tmpfile4, DATA)

	delFile(tmpfile2)
	delFile(tmpfile3)
	delFile(tmpfile4)

	'''

	# Transform data and do mean using otto basin
	ottofile = 't:\\OUTROS\\OTTOBACIAS\\ottobacia_nivel5_3km.tif'
	if(not os.path.isfile(ottofile)):
		print('...Error in path %s' %(ottofile))
		sys.exit()
	OTTO = readFileBand(ottofile,1)
	y_band = len(OTTO)
	x_band = len(OTTO[0])
	#OUT_DATA  = numpy.zeros((bands, y_band, x_band))
	OUT_DATA  = numpy.zeros((y_band, x_band))
	MEAN_DATA = numpy.zeros((y_band, x_band))
	STD_DATA = numpy.zeros((y_band, x_band))
	n_basin = numpy.amax(OTTO)
	print(n_basin)

	n_doy = 12
	dateStart = '20130101'
	dateEnd   = '20151231'
	date = dateStart
	while(date <= dateEnd):

		#datafile = 'T:\\OUTROS\\PREC_INT\\'+ str(date[0:4]) + '\\' + str(date[4:6]) + '\\'  + str(date) + '_PREC_INT.tif'
		datafile = 'T:\\EUMETSAT\\INTERP\\NDWI\\'+ str(date[0:4]) + '\\' + str(date[4:6]) + '\\'  + str(date) + '_interp_NDWI.tif'
		outputSubdir = 'T:\\EUMETSAT\\BASIN\\NDWI\\'+ str(date[0:4]) + '\\' + str(date[4:6]) + '\\'
		if(not os.path.isfile(datafile)):
			print('...file not exist: %s'%(datafile))
			sys.exit()
		else:
			#print('...file not exist: %s'%(datafile))
			if(not os.path.isdir(outputSubdir)):
				os.makedirs(outputSubdir)
			output_filename = outputSubdir +  str(date) + '_basin_NDWI.tif'
			DATA = readFileBand(datafile, 1)
			for basin in range(1,n_basin+1,1):
				sample = numpy.count_nonzero(DATA[OTTO == basin])
				#print('...process %d basin has %d'%(basin, sample))

				if(sample > 0):
					mean = numpy.mean(DATA[OTTO == basin])
					std = (DATA[OTTO == basin]).std()
					MEAN_DATA[OTTO == basin] = mean 
					STD_DATA[OTTO == basin]  = std
					
			saveRASTERfile(output_filename, ottofile, [MEAN_DATA, STD_DATA])
			print('...saving file %s'%(output_filename))
			OUT_DATA = MEAN_DATA 
			MEAN_DATA = numpy.zeros((y_band, x_band))
			STD_DATA = numpy.zeros((y_band, x_band))
			new_date = str(inc_date_doy(date, n_doy))
			if(int(str(new_date)[0:4]) > int(str(date)[0:4]) or '1227' in str(new_date)):
				date = (str(int(str(date)[0:4])+1) + '0101')
			else:
				date = str(new_date)
	'''
	# Generate NDDI = (NDVI - NDWI)/(NDVI + NDWI)
	ndvifile = 'T:\\EUMETSAT\OTTO\\ndvi_otto\\msg_ndvi_otto_m_all.tif' 
	ndwifile = 'T:\\EUMETSAT\OTTO\\ndwi_otto\\msg_ndwi_otto_m_all.tif' 
	nddifile = 'T:\\EUMETSAT\OTTO\\nddi_otto\\msg_nddi_otto_m_num.tif' 
	ottofile = 'T:\\OUTROS\\OTTOBACIAS\\ottobacia_nivel5_3km.tif'

	if(os.path.isfile(ottofile)):
		OTTO = readFileBand(ottofile, 1)
		y_band = len(OTTO)
		x_band = len(OTTO[0])
		n_basin = numpy.amax(OTTO)
	else:
		print('Error to read: %s'%(ottofile))
		sys.exit()
	NDVI = numpy.zeros((12, y_band,x_band))
	NDWI = numpy.zeros((12, y_band,x_band))
	if(os.path.isfile(ndvifile)):
		for band in range(1,12+1,1):
			NDVI[band-1] =  readFileBand(ndvifile, band)
	else:
		print('Error to read: %s'%(ndvifile))
		sys.exit()
	if(os.path.isfile(ndwifile)):
		for band in range(1,12+1,1):
			NDWI[band-1] =  readFileBand(ndwifile, band)
	else:
		print('Error to read: %s'%(ndwifile))
		sys.exit()

	NDDI = numpy.zeros((12,y_band, x_band))
	for band in range(12):
		NDDI[band] = (NDVI[band] - 2.5*NDWI[band])/(NDVI[band]+2.5*NDWI[band])
		saveRASTERfile(nddifile.replace('num', str(band+1)), ottofile, NDDI[band])
	saveRASTERfile(nddifile.replace('num', 'all'), ottofile, NDDI)
	'''

	'''
	# Get data table about otto basin
	bands = 12
	trmmfile = 'T:\\OUTROS\TRMM_OTTO\\trmm_otto_m_all.tif'
	ndvifile = 'T:\\EUMETSAT\OTTO\\ndvi_otto\\msg_ndvi_otto_m_all.tif' 
	ndwifile = 'T:\\EUMETSAT\OTTO\\ndwi_otto\\msg_ndwi_otto_m_all.tif' 
	ottofile = 'T:\\OUTROS\\OTTOBACIAS\\ottobacia_nivel5_3km.tif'
	outputDir = 'T:\\EUMETSAT\\ANALYSES\\OTTO_1010\\'
	percfile = 'T:\\OUTROS\\TerraClass\\TClass_maskVectorGoias_go_go_perc_class.tif'

	if(os.path.isfile(ottofile)):
		OTTO = readFileBand(ottofile, 1)
		y_band = len(OTTO)
		x_band = len(OTTO[0])
		n_basin = numpy.amax(OTTO)
	else:
		print('Error to read: %s'%(ottofile))
		sys.exit()
	TRMM = numpy.zeros((12, y_band,x_band))
	NDVI = numpy.zeros((12, y_band,x_band))
	NDWI = numpy.zeros((12, y_band,x_band))
	PERC_CLASS = numpy.zeros((4, y_band,x_band))
	print('Start get data in data bases')
	if(os.path.isfile(trmmfile)):
		for band in range(1,12+1,1):
			TRMM[band-1] = readFileBand(trmmfile, band)
	else:
		print('Error to read: %s'%(trmmfile))
		sys.exit()

	if(os.path.isfile(ndvifile)):
		for band in range(1,12+1,1):
			NDVI[band-1] =  readFileBand(ndvifile, band)
	else:
		print('Error to read: %s'%(ndvifile))
		sys.exit()

	if(os.path.isfile(ndwifile)):
		for band in range(1,12+1,1):
			NDWI[band-1] =  readFileBand(ndwifile, band)
	else:
		print('Error to read: %s'%(ndwifile))
		sys.exit()

	if(os.path.isfile(percfile)):
		for band in range(1,4+1,1):
			PERC_CLASS[band-1] =  readFileBand(percfile, band)
	else:
		print('Error to read: %s'%(percfile))
		sys.exit()
	DATA_SAMP = datafile('', outputDir, 'data_trmm_msg_ndvi_ndwi_perc.csv', 'BASIN SAMPLES CLASS1 CLASS2 CLASS3 CLASS4 MONTH TRMM NDVI NDWI NDVI_1 NDWI_1 NDVI_2 NDWI_2 NDVI_3 NDWI_3', 400*12, 16)

	for band in range(1,bands+1,1):
		print('Get data in band %d of hydrographic basin'%(band))
		for basin in range(1,n_basin+1,1):
			sample = numpy.sum(OTTO[OTTO == basin])
			if(sample != 0):
				trmm = numpy.mean(TRMM[band-1][OTTO == basin])
				ndvi = numpy.mean(NDVI[band-1][OTTO == basin]) 
				ndwi = numpy.mean(NDWI[band-1][OTTO == basin])
				class_1 = numpy.mean(PERC_CLASS[0][OTTO == basin])
				class_2 = numpy.mean(PERC_CLASS[1][OTTO == basin])
				class_3 = numpy.mean(PERC_CLASS[2][OTTO == basin])
				class_4 = numpy.mean(PERC_CLASS[3][OTTO == basin])
				ndvi_1 = numpy.mean(NDVI[(band+0)%12][OTTO == basin])
				ndwi_1 = numpy.mean(NDWI[(band+0)%12][OTTO == basin])
				ndvi_2 = numpy.mean(NDVI[(band+1)%12][OTTO == basin])
				ndwi_2 = numpy.mean(NDWI[(band+1)%12][OTTO == basin])
				ndvi_3 = numpy.mean(NDVI[(band+2)%12][OTTO == basin])
				ndwi_3 = numpy.mean(NDWI[(band+2)%12][OTTO == basin])

				DATA_SAMP.addAtrib(basin)
				DATA_SAMP.addAtrib(sample)
				DATA_SAMP.addAtrib(class_1)
				DATA_SAMP.addAtrib(class_2)
				DATA_SAMP.addAtrib(class_3)
				DATA_SAMP.addAtrib(class_4)
				DATA_SAMP.addAtrib(band)
				DATA_SAMP.addAtrib(trmm)
				DATA_SAMP.addAtrib(ndvi)
				DATA_SAMP.addAtrib(ndwi)
				DATA_SAMP.addAtrib(ndvi_1)
				DATA_SAMP.addAtrib(ndwi_1)
				DATA_SAMP.addAtrib(ndvi_2)
				DATA_SAMP.addAtrib(ndwi_2)
				DATA_SAMP.addAtrib(ndvi_3)
				DATA_SAMP.addAtrib(ndwi_3)
				print(DATA_SAMP.getLine())
			else:
				DATA_SAMP.addAtrib(basin)
				DATA_SAMP.addAtrib([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])	
			DATA_SAMP.addSample()
		DATA_SAMP.save()
	'''