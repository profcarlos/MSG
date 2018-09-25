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

def saveRASTERfile(filename, fileMask, DATA):
	num_bands = len(DATA)
	# For DATA with one band only
	if(num_bands > 15):
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

	
	no_data = 15
	print('Start to generate file')
	outputDir = 'T:\\EUMETSAT\\ANALYSES\\PCA_2910\\'

	filename = 'T:\\EUMETSAT\\PCA\\maskGoias.tif'
	if(os.path.isfile(filename)):
		MASK = readFileBand(filename, 1)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()
	filename = 'T:\\OUTROS\\TerraClass\\TClass_maskVectorGoias_go_perc.tif'
	if(os.path.isfile(filename)):
		CLASS = readFileBand(filename, 1)
		PERC  = readFileBand(filename, 2)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()

	filename = 'T:\\EUMETSAT\\PCA\\msg_pca_3bands_20130101_20151231.tif'
	if(os.path.isfile(filename)):
		PC1_MSG = readFileBand(filename, 1)
		PC2_MSG = readFileBand(filename, 2)
		PC3_MSG = readFileBand(filename, 2)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()
	
	filename = 'T:\\EUMETSAT\\PCA\\mod_pca_3bands_v2010_20130101_20151231.tif'
	if(os.path.isfile(filename)):
		PC1_MOD = readFileBand(filename, 1)
		PC2_MOD = readFileBand(filename, 2)
		PC3_MOD = readFileBand(filename, 3)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()

	y_band = len(MASK)
	x_band = len(MASK[0])
	DATA_SAMP = datafile('', outputDir, 'class_pca_mod_msg_3axis.csv', 'CLASS PC1_MOD PC2_MOD PC3_MOD PC1_MSG PC2_MSG PC3_MSG', y_band*x_band, 7)
	for x in range(x_band):
		for y in range(y_band):
			if(CLASS[y][x] == no_data):
				continue
			else:
				DATA_SAMP.addAtrib(CLASS[y][x])
				DATA_SAMP.addAtrib(PC1_MOD[y][x])
				DATA_SAMP.addAtrib(PC2_MOD[y][x])
				DATA_SAMP.addAtrib(PC3_MOD[y][x])
				DATA_SAMP.addAtrib(PC1_MSG[y][x])
				DATA_SAMP.addAtrib(PC2_MSG[y][x])
				DATA_SAMP.addAtrib(PC3_MSG[y][x])
				print(DATA_SAMP.getLine())
				DATA_SAMP.addSample()


		DATA_SAMP.save()
	

	'''
	file = 'T:\\EUMETSAT\\PCA\\msg_pca_20130101_20151231.tif'
	num_bands = 30
	print('statistics in file: %s'%(file))
	print('BAND\tMEAN\tVAR')
	for band in range(num_bands):
		DATA = readFileBand(file, band+1)
		mean = numpy.mean(DATA[DATA>0])
		var  = DATA[DATA>0].var()
		print('%d\t%.4f\t%.4f'%(band+1, mean, var)) 

	file = 'T:\\EUMETSAT\\PCA\\mod_pca_v2010_20130101_20151231.tif'
	num_bands = 23
	print('statistics in file: %s'%(file))
	print('BAND\tMEAN\tVAR')
	for band in range(num_bands):
		DATA = readFileBand(file, band+1)
		mean = numpy.mean(DATA[DATA>0])
		var  = DATA[DATA>0].var()
		print('%d\t%.4f\t%.4f'%(band+1, mean, var)) 
	'''

	'''
	no_data = 15
	print('Start to generate file: class_mod_msg.csv')
	outputDir = 'T:\\EUMETSAT\\ANALYSES\\PCA_2010\\'
	filename = 'T:\\EUMETSAT\\PCA\\mod_kmean_4ids_pca_3bands_v2010_20130101_20151231_v3.tif'
	if(os.path.isfile(filename)):
		MOD = readFileBand(filename, 1)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()
	filename = 'T:\\OUTROS\\TerraClass\\TClass_maskVectorGoias_go_perc.tif'
	if(os.path.isfile(filename)):
		CLASS = readFileBand(filename, 1)
		PERC  = readFileBand(filename, 2)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()

	filename = 'T:\\EUMETSAT\\PCA\\msg_kmean_4ids_pca_3bands_20130101_20151231.tif'
	if(os.path.isfile(filename)):
		MSG = readFileBand(filename, 1)
	else:
		print('Error to read: %s'%(filename))
		sys.exit()
	

	y_band = len(MSG)
	x_band = len(MSG[0])
	DATA_SAMP = datafile('', outputDir, 'class_modis_msg_v2010_v2.csv', 'CLASS PERC LEV_PERC MOD MSG', y_band*x_band, 5)
	for x in range(x_band):
		for y in range(y_band):
			if(CLASS[y][x] == no_data):
				continue
			else:
				DATA_SAMP.addAtrib(CLASS[y][x])
				DATA_SAMP.addAtrib(PERC[y][x])
				DATA_SAMP.addAtrib(round(PERC[y][x]/10-5))
				DATA_SAMP.addAtrib(MOD[y][x])
				DATA_SAMP.addAtrib(MSG[y][x])
				print(DATA_SAMP.getLine())
				DATA_SAMP.addSample()


	DATA_SAMP.save()
	'''
	'''
	# View histogram of data file
	filename ='T:\\EUMETSAT\\ANALYSES\\PCA_2909\\mod13_kmean_4ids_pca_3bands_20130101_20151231.tif'
	histogram(filename, 16, 1)
	filename ='T:\\EUMETSAT\\ANALYSES\\PCA_2909\\msg_kmean_4ids_pca_3bands_20130101_20151231.tif'
	histogram(filename, 16, 1)	
	filename = 'T:\\OUTROS\\TerraClass\\TClass_maskVectorGoias_go_perc.tif'
	histogram(filename, 16, 1)
	sys.exit()
	'''

	'''
	#Create one band file (NDVI)
	#folder = 'T:\\EUMETSAT\\DAY_INT\\*.tif'
	folder = 'T:\\EUMETSAT\\ANALYSES\\PCA\\MSG_3KM_FILTER_123\\*DAY_INT.tif'
	files = glob.glob(folder)
	print(files)
	for file in files:
		print('new file: %s'%(file))
		BAND1 = readFileBand(file, 1)
		BAND1[BAND1 < 0] = 0 
		filename = file.replace('.tif', '_b1.tif')
		# saveRASTERfile(filename, fileMask, DATA)
		saveRASTERfile(filename, file, BAND1)
	'''
