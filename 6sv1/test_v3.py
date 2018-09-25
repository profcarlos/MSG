import sys
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')

import numpy
import time
import os
import itertools
from utils import * 
from sys import argv
import gc

otbDir = 'E:\\Install\\OTB-5.4.0-win64\\OTB-5.4.0-win64\\bin\\'
exeDir = 'C:\\OSGeo4W\\bin\\'
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoias123 = os.getcwd() + '\\shapes\\201401011300_123.tif'
fileGoiasHRV = os.getcwd() + '\\shapes\\201401011700_HRV.tif'
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMaskSHP = fileGoias
fileMaskHRV = fileGoiasHRV
fileMask123 = fileGoias123

usage = """\
Usage: %s [OPTIONS]
		-ty     type data (1day, 2days)
		-fi     folder input data
		-fo     folder output data
		-ds     date to start process (format: yyyyMMdd)
		-de     date to end process (format: yyyyMMdd)

""" % argv[0]

def importATMdata(inputDir, date):
	#print('ATMclass read')
	# Read ATM file
	# band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
	# band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
	# band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
	# band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)
	# ATM to 6S: Vector with atmospheric data [WV O AOT]

	ds = gdal.Open(fileMask123)
	y_123 = ds.RasterYSize
	x_123 = ds.RasterXSize
	# DATA[3][y_123][x_123] = [WV, OZONE, AOT550]
	DATA = numpy.zeros((3, y_123, x_123))

	inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
	filename = inputSubDir.replace('___', 'ATM') + str(date) + '_ATM.tif'
	if os.path.isfile(filename) == 0:
		#print('Error to read ATM file: %s' % (filename))
		return []
	else:
		#print('Read ATM file: %s' % (filename))
		'''
		# get Aerosol Data Band
		vet = readFileBand(filename, 4)
		if(not len(vet)):
			return [0, 0, 0]
		else:
			DATA[2] = vet[i][j]
		'''
		# get Water_Vapour  
		vet = readFileBand(filename, 10)
		if(not len(vet)):
			return []
		else:
			DATA[0] = vet

		# get Total Ozone Data Band
		vet = readFileBand(filename, 7)
		if(not len(vet)):
			return []
		else:
			DATA[1] = vet/1000 # Convert Dobson Unit to cm-atm
		# get Deep_Blue Data Band
		vet = readFileBand(filename, 1)
		if(not len(vet)):
			return []
		else:
			DATA[2] = vet
	return DATA

def importBRDFdata(inputDir, date):
	startHour = 1300
	endHour   = 1700    

	ds = gdal.Open(fileMask123)
	y_123 = ds.RasterYSize
	x_123 = ds.RasterXSize
	# DATA[17][8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
	DATA = numpy.zeros((17, 8, y_123, x_123))
	inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
	hour = startHour
	h = 0
	fileError = False
	while(hour <= endHour):
		# --------- Process ANGLE file
		filename = inputSubDir.replace('___', 'ANGLE') + str(date) + str(hour) + '_ANGLES.tif'
		if os.path.isfile(filename) == 0:
			print('..Error! Invalid data file: %s' %filename)
			fileError = True
		else:
			#print('.Read ANGLE file:%s'%(filename))
			# get VZA View Zenith Angle (msg_zenres.tif)
			vet = readFileBand(filename, 3)
			if(not len(vet)):
				fileError = True
			else:
				DATA[h][0] = vet
			# get VAZ View Azimuth Angle (msg_azres.tif)
			vet = readFileBand(filename, 1)
			if(not len(vet)):
				fileError = True
			else:
				DATA[h][1] = vet
			# get SZA Sun Zenith Angle (sun_zenres.tif)
			vet = readFileBand(filename, 4)
			if(not len(vet)):
				fileError = True
			else:
				DATA[h][2] = vet
			# get SAZ Sun Azimuth Angle (sol_azres.tif)
			vet = readFileBand(filename, 2)
			if(not len(vet)):
				fileError = True
			else:
				DATA[h][3] = vet
		# --------- Process 123 file
		if(fileError == False):
			#print('Process 123 file')
			filename = inputSubDir.replace('___', '123') + str(date) + str(hour) + '_123.tif'
			if os.path.isfile(filename) == 0:
				print('..Error! Invalid data file: %s' %filename)
				fileError = True
			else:
				#print('.Read 123 file:%s'%(filename))
				# get RED Data Band
				vet = readFileBand(filename, 1)
				if(not len(vet)):
					fileError = True
				else:
					DATA[h][4] = vet
				# get NIR Data Band
				vet = readFileBand(filename, 2)
				if(not len(vet)):
					fileError = True
				else:
					DATA[h][5] = vet
				# get WIR Data Band
				vet = readFileBand(filename, 3)
				if(not len(vet)):
					fileError = True
				else:
					DATA[h][6] = vet

		# --------- Process CLM file
		if(fileError == False):
			#print('Process CLM file')
			filename = inputSubDir.replace('___', 'CLM') + str(date) + str(hour) + '_CLM.tif'
			if os.path.isfile(filename) == 0:
				print('..Error! Invalid data file: %s' %filename)
				fileError = True
			else:
				#print('.Read CLM file:%s'%(filename))
				vet = readFileBand(filename, 1)
				if(not len(vet)):
					fileError = True
				else:
					DATA[h][7] = vet
		# Verify if error clear data, else inc position
		if(fileError == False):
			h = h + 1
		hour = incHour(hour)
	return DATA[:h]

class corr_atm():

	'''
		|Variable | Min value | Max value | Step   | Num Step| Mult |
		| VZA     | 55        | 65        | 2.5    | 4       | 10   |
		| VAZ     | 100       | 110       | 2.5    | 4       | 10   |
		| SZA     | 30        | 60        | 5.0    | 6       | 10   |
		| SAZ     | 120       | 230       | 20     | 6       | 10   |
		| WV      | 1.0       | 7.0       | 1.0    | 7       | 10   |
		| OZ      | 0.20      | 0.30      | 0.025  | 4       | 1000 |
		| AOT     | 0.00      | 0.4       | 0.1    | 4       | 1000 |
	'''
	def __init__(self):
		print('...__init__')
		self.VARS = ['VZA', 'VAZ', 'SZA', 'SAZ', 'WV', 'OZ', 'AOT', 'XA_RED', 'XB_RED', 'XC_RED', 'XA_NIR', 'XB_NIR', 'XC_NIR']
		self.START = [[550, 651, 25, 10], [1000, 1101, 25, 10], [300, 601, 50, 10], [1200, 2301, 200, 10], [10, 71, 10, 10],[200, 301, 25, 1000], [0, 401, 100, 1000]]
		# Order VZA VAZ SZA SAZ WV OZ AOT XA_RED XB_RED XC_RED XA_NIR XB_NIR XC_NIR
		# Order 550 1050 600 2200 30 250 200 516 10170 9516 710 6016 6533

		tac = time.clock()
		for vza in (self.drange(self.START[0][:3])):
			path = os.getcwd() + '\\LUT_6S\\LUT_6S_VZA_' + str(vza) + '.csv'
			print(vza)
			if(vza == 550):
				DATA = numpy.loadtxt(path, dtype='int32', delimiter = ' ', skiprows = 1)
			else:
				DATA = numpy.append(DATA, numpy.loadtxt(path, dtype='int32', delimiter = ' ', skiprows = 1), axis = 0)
			print('DATA [%d][%d]' %(len(DATA), len(DATA[0])))
		print('time to get data %.2f'%(time.clock()-tac))
		DATA = DATA.T
		for name, i in zip(self.VARS, range(len(self.VARS))):
			globals()[name] = numpy.vstack(DATA[i])
		DATA = []
		print('VZA [%d][%d]' %(len(VZA), len(VZA[0])))
		print('time to org data %.2f'%(time.clock()-tac))

	def adjust_parameters(self, PAR):
		print('...adjust_parameters')
		print('PAR input:')
		print(PAR) 
		# Put input parameters in range of data values of look-up tables
		for i in range(len(PAR)):
			par = PAR[i]*self.START[i][3]
			new_par = 0
			#print('i: %d par: %d'%(i, par))
			for lut in self.drange(self.START[i][:3]):
				#print('i: %d par: %.2f lut: %.2f'%(i, par, lut))
				#print('lut: %d (par - lut): %d dif: %d' %(lut, numpy.abs(par - lut), self.START[i][2]))
				if(numpy.abs(par - lut) <= self.START[i][2]):
					new_par = lut
					break
			if(i != 6 and new_par == 0):
				print(". Error! Verify data in look-up tables (variable %d): %(i)")
				print(PAR)
				sys.exit()
			PAR[i] = new_par
		print('PAR processed:')
		print(PAR)
		return(PAR)

	def calc_corr(self, PAR, red_rad, nir_rad):
		print('...calc_corr')
		tac = time.clock()
		VARS_sample = []
		for i in range(len(self.VARS[:7])):
			VARS_sample.append(self.VARS[i] + '_sample')
		print(PAR)
		for name, i in zip(VARS_sample, range(len(VARS_sample))):
			globals()[name] = PAR[i]
		paser = numpy.logical_and(VZA == VZA_sample, numpy.logical_and(VAZ == VAZ_sample, numpy.logical_and(SZA == SZA_sample, \
				 numpy.logical_and(SAZ == SAZ_sample, numpy.logical_and(WV == WV_sample, numpy.logical_and(OZ == OZ_sample, \
				 AOT  == AOT_sample))))))
		samples = len(paser[paser == True])
		print('samples: %d' %(samples))
		if(samples != 1):
			print('Error! not have sample to this data!')
		else:
			print('sample_data')
			print('RED_DATA: %d %d %d' %(XA_RED[paser], XB_RED[paser], XC_RED[paser]))
			print('NIR_DATA: %d %d %d' %(XA_NIR[paser], XB_NIR[paser], XC_NIR[paser]))
		print('time to search %.2f'%(time.clock()-tac))
		red_ref = self.calc_rad(red_rad, XA_RED[paser], XB_RED[paser], XC_RED[paser])
		nir_ref = self.calc_rad(nir_rad, XA_NIR[paser], XB_NIR[paser], XC_NIR[paser])
		print('[red, nir]: %.4f %.4f'%(red_ref, nir_ref))
		return [red_ref, nir_ref]
	
	def calc_rad(self, L, xa, xb, xc):
		L = L*1000
		xa = xa/100000
		xb = xb/100000
		xc = xc/100000
		print('L: %.4f xa: %.4f xb: %.4f xc: %.4f'%(L, xa, xb, xc))
		y = xa*L-xb
		acr = y/(1 + xc*y)
		return acr

	def drange(self, data):
		[start, stop, step] = data
		rng = start
		while rng < stop:
			yield rng
			rng += step
		return rng

class main():

	argDict = mapDict(argv, usage)
	gc.collect()
	
	if ("-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict):
		dateStart = int(argDict["-ds"])
		dateEnd = int(argDict["-de"])
		outputDir = argDict["-fo"]
		inputDir  = argDict["-fi"]
		typeData = argDict["-ty"]
	else:
		exit(usage)

	if(typeData != '1day' and typeData != '2days'):
		print('Error in typeData, verify!')
		exit(usage)

	ds = gdal.Open(fileMask123)
	y_123 = ds.RasterYSize
	x_123 = ds.RasterXSize

	# Process in range of data
	date = dateStart
	while(date < dateEnd):
		#Output file name 
		print('...Create 6S data: %s'%(date))  
		outputFolder = outputDir + '6S\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
		# Verify and create folders
		if(not os.path.isdir(outputFolder)):
			os.makedirs(outputFolder)

		outputFile = outputFolder + str(date) + '_6S.tif'    

		# Verify if files exist
		if(os.path.isfile(outputFile)):
			print('... Data file exist in: %s'%(outputFile))
			date = incDate(date)
			continue

		# DATA[17][8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
		DATA = importBRDFdata(inputDir, date)
		if(len(DATA) == 0):
			print('...Error to read BRDF data: %s'%(inputDir + str(date) + '___.tif'))
			date = incDate(date)
			continue

		# ATM = [WV OZ AOT][y_123][x_123]    
		ATM = importATMdata(inputDir, date)
		if(len(ATM) == 0):
			print('...Error to read ATM data: %s'%(inputDir + date + '_ATM.tif'))
			date = incDate(date)
			continue
		six = corr_atm()
		tac = time.clock()
		hours = len(DATA)
		#print('hours:%d [%d, %d]'%(hours, y_123, x_123))
		for x in range(x_123):
			for y in range(y_123):
				if(DATA[0][0][y][x] == 0):
					continue
				for h in range (hours):
					# PAR is parameters to process atmospheric correction
					# PAR = [VZA, VAZ, SZA, SAZ, WV, OZ, AOT]
					PAR = [DATA[h][0][y][x], DATA[h][1][y][x], DATA[h][2][y][x], DATA[h][3][y][x], ATM[0][y][x], ATM[1][y][x], ATM[2][y][x]]
					print('...PAR join')
					[red, nir] = six.calc_corr(six.adjust_parameters(PAR), DATA[h][4][y][x], DATA[h][5][y][x])
					sys.exit()