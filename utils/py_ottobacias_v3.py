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

class datafile:

	def __init__(self, name, outputDir, filename, header, x, y):
		self.name = name
		self.outputDir = outputDir
		self.filename = filename.replace('.csv', name + '.csv')
		self.header  = header
		self.sample = 0
		self.pos = 0
		self.MAT = numpy.zeros((x, y))

	def addSample(self):
		self.sample = self.sample + 1
		self.pos = 0

	def addAtrib(self, data):
		try:
			n_atrib = len(data)
		except:
			#print('... Except ValueError')
			n_atrib = 1
			data = [data]
		for x in range(n_atrib):
			self.MAT[self.sample][self.pos] = data[x]
			#print('%s[%d] add: %.2f' %(self.name, self.pos, data[x]))
			self.pos = self.pos + 1

	def getAll(self):
		return self.MAT
	
	def getLine(self): 
		return self.MAT[self.sample]

	def save(self):
		# Clear data in MAT
		for id_x in range (len(self.MAT)):
			if(numpy.sum(self.MAT[id_x]) == 0):
				break
		'''
		for id_y in range (len(self.MAT[0])):
			if(self.MAT[0][id_y] < 0):
				break
		'''
		id_y = 	len(self.MAT[0])	
		# Make new array cutting zeros rows and cols
		MAT_OUT = numpy.zeros((id_x, id_y))
		for x in range(id_x):
			for y in range(id_y):
				MAT_OUT[x][y] = self.MAT[x][y]
		# Save data file
		file = self.outputDir + self.filename
		if(os.path.isfile(file)):
			delFile(file)
		#print('--- Save output file: %s'%(file))
		numpy.savetxt(file, MAT_OUT, fmt='%1.4f', comments='', header = self.header, delimiter = ' ', newline='\n')
		#self.savgol()
		
	def savgol(self):
		MAT_SAV = savgolFilter(self.MAT, 'MSG')
		filename = self.filename.replace('.csv', '_SAV.csv')
		print('--- Save output file: %s'%(self.outputDir + filename))
		numpy.savetxt(self.outputDir + filename, MAT_SAV, fmt='%1.4f', comments='', header = self.header, delimiter = ' ', newline='\n')

def delFile(file):
	#print( '..delFile: %s' % str(glob.glob(file)))
	#print('Remove files:')
	for x in glob.glob(file):
		try:
			os.remove(x)
			#print(x)
		except: 
			print('Error in del file: %s' % str(glob.glob(file)))

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
	'''
	# Transform data and do mean using otto basin
	#datafile = 'T:\\EUMETSAT\\PCA\\msg_pca_NDWI_20130101_20151231.tif'
	#datafile = 'T:\\OUTROS\TRMM\\pa_go_trmm_all_bands.tif'
	#datafile = 'T:\\OUTROS\TRMM_OTTO\\pa_go_trmm_all_bands.tif'
	datafile = 'T:\\OUTROS\\TerraClass\\TClass_maskVectorGoias_go_go_perc_class.tif'
	ottofile = 't:\\OUTROS\\OTTOBACIAS\\ottobacia_nivel5_3km.tif'
	#datafile =  't:\\EUMETSAT\\OTTO\\raster_idhm_2010.tif'
	#datafile = 'T:\\OUTROS\\Indices_Solo\\Indice_Medio_Solo_go.tif'
	
	#datafile = ottofile
	#output_filename = 'T:\\EUMETSAT\\OTTO\\msg_pca_NDWI_doy_num.tif'
	#output_filename = 'T:\\OUTROS\\TRMM_OTTO\\trmm_otto_m_num.tif'
	output_filename = 'T:\\OUTROS\\TerraClass\\princ_use_of_earth_basin.tif'
	#output_filename = 'T:\\OUTROS\\OTTOBACIAS\\otto_basin_len_all.tif'
	#output_filename = 't:\\EUMETSAT\\BASIN\\Indice_Medio_Solo_go_basin.tif'
	doy = 1
	#bands = 12
	bands = 4
	
	OTTO = readFileBand(ottofile,1)
	y_band = len(OTTO)
	x_band = len(OTTO[0])
	#OUT_DATA  = numpy.zeros((bands, y_band, x_band))
	OUT_DATA  = numpy.zeros((bands+1, y_band, x_band))
	n_basin = numpy.amax(OTTO)
	DATA = numpy.zeros((bands, y_band, x_band))

	print(n_basin)
	if(os.path.isfile(datafile)):
		for band in range(1,bands+1,1):
			DATA[band-1] =  readFileBand(datafile, band)
	else:
		print('Error to read: %s'%(datafile))
		sys.exit()
		
	for basin in range(1,n_basin+1,1):
		print('Get data in band %d of hydrographic basin'%(band))
		samples = []
		for band in range(0,bands,1):
			samples.append(numpy.mean(DATA[band][OTTO == basin]))

		if(numpy.sum(samples) <= 0):
			continue
		argmax = numpy.argmax(samples)
		if(samples[argmax] > 50):
			OUT_DATA[0][OTTO == basin] = argmax + 1
		else:
			OUT_DATA[0][OTTO == basin] = 4
		for band in numpy.arange(0,bands,1):
			OUT_DATA[band+1][OTTO == basin] = samples[band]
			
	saveRASTERfile(output_filename, ottofile, OUT_DATA)
	print('...saving file %s'%(output_filename))
	sys.exit()
	'''
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
	'''
	ottofile = 'T:\\OUTROS\\OTTOBACIAS\\ottobacia_nivel5_3km.tif'

	ndvifile = 'T:\\EUMETSAT\\STAND\\std_index_basin_NDVI_coef.tif'
	ndwifile = 'T:\\EUMETSAT\\STAND\\std_index_basin_NDWI_coef.tif'
	precfile = 'T:\\EUMETSAT\\STAND\\std_index_basin_PREC_coef.tif'
	solofile = 'T:\\EUMETSAT\\BASIN\\Indice_Medio_Solo_go_basin.tif'
	classfile = 'T:\\EUMETSAT\\OTTO\\princ_uso_terra_bacias.tif'
	percfile = 'T:\\OUTROS\\TerraClass\\TClass_maskVectorGoias_go_go_perc_class.tif'
	ndwisum = 'T:\\EUMETSAT\\OTTO\\ndwi_sum_year.tif'
	idhm = 'T:\\EUMETSAT\\OTTO\\raster_idhm_2010_go_basin.tif'
	trmmsum = 'T:\\EUMETSAT\\OTTO\\trmm_sum_year.tif'
	if(os.path.isfile(ottofile)):
		OTTO = readFileBand(ottofile, 1)
		y_band = len(OTTO)
		x_band = len(OTTO[0])
		PERC_CLASS = numpy.zeros((4, y_band,x_band))
		n_basin = numpy.amax(OTTO)
	else:
		print('Error to read: %s'%(ottofile))
		sys.exit()

	print('Start get data in data bases')

	if(os.path.isfile(ndvifile)):
		NDVI =  readFileBand(ndvifile, 1)
	else:
		print('Error to read: %s'%(ndvifile))
		sys.exit()
	if(os.path.isfile(ndwifile)):
		NDWI =  readFileBand(ndwifile, 1)
	else:
		print('Error to read: %s'%(ndwifile))
		sys.exit()
	if(os.path.isfile(precfile)):
		PREC = readFileBand(precfile, 1)
	else:
		print('Error to read: %s'%(precfile))
		sys.exit()
	if(os.path.isfile(solofile)):
		SOLO =  readFileBand(solofile, 1)
	else:
		print('Error to read: %s'%(solofile))
		sys.exit()
	if(os.path.isfile(classfile)):
		CLASS =  readFileBand(classfile, 1)
	else:
		print('Error to read: %s'%(classfile))
		sys.exit()
	if(os.path.isfile(percfile)):
		for band in range(1,4+1,1):
			PERC_CLASS[band-1] =  readFileBand(percfile, band)
	else:
		print('Error to read: %s'%(percfile))
		sys.exit()
	if(os.path.isfile(trmmsum)):
		TRMM_SUM =  readFileBand(trmmsum, 1)
	else:
		print('Error to read: %s'%(trmmsum))
		sys.exit()
	if(os.path.isfile(ndwisum)):
		NDWI_SUM =  readFileBand(ndwisum, 1)
	else:
		print('Error to read: %s'%(ndwisum))
		sys.exit()
	if(os.path.isfile(idhm)):
		IDHM =  readFileBand(idhm, 1)
	else:
		print('Error to read: %s'%(idhm))
		sys.exit()
	DATA_SAMP = datafile('', 't:\\EUMETSAT\\ANALYSES\\OTTO_2011\\', 'data_msg_ndvi_ndwi_perc_basin.csv', 'BASIN SAMPLES NDVI NDWI PREC SOLO TRMM_SUM NDWI_SUM IDHM CLASS_P CLASS1 CLASS2 CLASS3 CLASS4', 400, 15)

	for basin in range(1,n_basin+1,1):
		sample = numpy.sum(OTTO[OTTO == basin])
		if(sample != 0):
			ndvi = numpy.mean(NDVI[OTTO == basin]) 
			ndwi = numpy.mean(NDWI[OTTO == basin])
			prec = numpy.mean(PREC[OTTO == basin])
			class_p = numpy.mean(CLASS[OTTO == basin])
			class_1 = numpy.mean(PERC_CLASS[0][OTTO == basin])
			class_2 = numpy.mean(PERC_CLASS[1][OTTO == basin])
			class_3 = numpy.mean(PERC_CLASS[2][OTTO == basin])
			class_4 = numpy.mean(PERC_CLASS[3][OTTO == basin])
			solo    = numpy.mean(SOLO[OTTO == basin])
			trmm_sum   = numpy.mean(TRMM_SUM[OTTO == basin])
			ndwi_sum   = numpy.mean(NDWI_SUM[OTTO == basin])
			idhm  = numpy.mean(IDHM[OTTO == basin])
			DATA_SAMP.addAtrib(basin)
			DATA_SAMP.addAtrib(sample)
			DATA_SAMP.addAtrib(ndvi)
			DATA_SAMP.addAtrib(ndwi)
			DATA_SAMP.addAtrib(prec)
			DATA_SAMP.addAtrib(solo)
			DATA_SAMP.addAtrib(trmm_sum)
			DATA_SAMP.addAtrib(ndwi_sum)
			DATA_SAMP.addAtrib(idhm)
			DATA_SAMP.addAtrib(class_p)
			DATA_SAMP.addAtrib(class_1)
			DATA_SAMP.addAtrib(class_2)
			DATA_SAMP.addAtrib(class_3)
			DATA_SAMP.addAtrib(class_4)
		
			print(DATA_SAMP.getLine())
		else:
			DATA_SAMP.addAtrib(basin)
			DATA_SAMP.addAtrib([0,0,0,0,0,0,0,0,0,0,0,0,0])	
		DATA_SAMP.addSample()
	DATA_SAMP.save()
	'''


	munfile  = 'C:\\USERS\\CARLOS.SILVEIRA\\DROPBOX\\ANALYSES\\MUN_0112\\mun_goias_cod_ibge.tif'
	credfile = 'C:\\USERS\\CARLOS.SILVEIRA\\DROPBOX\\ANALYSES\\MUN_0112\\rec_pastagem_mun_2013_2015.tif'
	swifile  = 'T:\\EUMETSAT\\STAND\\std_index_NDWI_coef.tif'
	precfile = 'T:\\OUTROS\\TRMM_OTTO\\trmm_otto_m_all.tif'
	mesofile = 'C:\\USERS\\CARLOS.SILVEIRA\\DROPBOX\\ANALYSES\\MUN_0112\\mun_goias_cod_meso.tif'
	if(os.path.isfile(munfile)):
		MUN = readFileBand(munfile, 1)
		y_band = len(MUN)
		x_band = len(MUN[0])
		max_mun = 5222302
		min_mun = 5200050
		print('Start process...')
		print('num municipios: %s %s'%(str(max_mun), str(min_mun)))
	else:
		print('Error to read: %s'%(munfile))
		sys.exit()

	print('Start get data in data bases')

	if(os.path.isfile(credfile)):
		CRED =  readFileBand(credfile, 1)
	else:
		print('Error to read: %s'%(credfile))
		sys.exit()

	if(os.path.isfile(swifile)):
		SWI =  readFileBand(swifile, 1)
	else:
		print('Error to read: %s'%(swifile))
		sys.exit()
	if(os.path.isfile(precfile)):
		PREC =  readFileBand(precfile, 1)
	else:
		print('Error to read: %s'%(precfile))
		sys.exit()
	if(os.path.isfile(mesofile)):
		MESO =  readFileBand(mesofile, 1)
	else:
		print('Error to read: %s'%(mesofile))
		sys.exit()
	DATA_SAMP = datafile('', 'C:\\Users\\carlos.silveira\\Dropbox\\ANALYSES\\MUN_0112\\', 'data_mun_cred_n_pixels.csv', 'MUN MESO CRED PREC SWI_VERDE SWI_AZUL SWI_ALL', 250, 7)

	for mun in numpy.arange(min_mun,max_mun+1,1):
		sample = numpy.sum(MUN[MUN == mun])
		if(sample != 0):
			cred = numpy.mean(CRED[MUN == mun]) 
			pix_verde = numpy.sum(numpy.logical_and(SWI[MUN == mun]>-0.0008, SWI[MUN == mun] < 0))
			pix_azul  = numpy.sum(SWI[MUN == mun]>= 0)
			prec = numpy.mean(PREC[MUN == mun])
			meso = numpy.mean(MESO[MUN == mun])
			DATA_SAMP.addAtrib(mun)
			DATA_SAMP.addAtrib(meso)
			DATA_SAMP.addAtrib(cred)
			DATA_SAMP.addAtrib(prec)
			DATA_SAMP.addAtrib(pix_verde)
			DATA_SAMP.addAtrib(pix_azul)
			DATA_SAMP.addAtrib(pix_verde+pix_azul)
			print(DATA_SAMP.getLine())
			DATA_SAMP.addSample()
	DATA_SAMP.save()
