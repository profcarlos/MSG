import sys
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\MSG\\Dropbox\\newPython\\utils')
import numpy
import time
import os
import itertools
from utils import * 
from sys import argv
import gc
import random
import glob


otbDir = 'E:\\Install\\OTB-5.4.0-win64\\OTB-5.4.0-win64\\bin\\'
exeDir = 'C:\\OSGeo4W\\bin\\'
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoias123 = os.getcwd() + '\\shapes\\201301011300_123.tif'
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
	#print('ATM file: %s'%(filename))
	if os.path.isfile(filename) == 0:
		print('Error to read ATM file: %s' % (filename))
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

def importBRDFdata(inputDir, date, hour):
  
	ds = gdal.Open(fileMask123)
	y_123 = ds.RasterYSize
	x_123 = ds.RasterXSize
	# DATA[8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
	DATA = numpy.zeros((8, y_123, x_123))
	inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'

	fileError = False

	# --------- Process ANGLE file
	filename = inputSubDir.replace('___', 'ANGLE') + str(date) + str(hour) + '_ANGLES.tif'
	if os.path.isfile(filename) == 0:
		print('..Error! Invalid data file: %s' %filename)
		fileError = True
	else:
		print('.Read ANGL file: %s'%(filename))
		# get VZA View Zenith Angle (msg_zenres.tif)
		vet = readFileBand(filename, 3)
		if(not len(vet)):
			fileError = True
		else:
			DATA[0] = vet
		# get VAZ View Azimuth Angle (msg_azres.tif)
		vet = readFileBand(filename, 1)
		if(not len(vet)):
			fileError = True
		else:
			DATA[1] = vet
		# get SZA Sun Zenith Angle (sun_zenres.tif)
		vet = readFileBand(filename, 4)
		if(not len(vet)):
			fileError = True
		else:
			DATA[2] = vet
		# get SAZ Sun Azimuth Angle (sol_azres.tif)
		vet = readFileBand(filename, 2)
		if(not len(vet)):
			fileError = True
		else:
			DATA[3] = vet
	# --------- Process 123 file
	if(fileError == False):
		#print('Process 123 file')
		filename = inputSubDir.replace('___', '123_RAD') + str(date) + str(hour) + '_123.tif'
		if os.path.isfile(filename) == 0:
			print('..Error! Invalid data file: %s' %filename)
			fileError = True
		else:
			print('.Read 123 file: %s'%(filename))
			# get RED Data Band
			vet = readFileBand(filename, 1)
			if(not len(vet)):
				fileError = True
			else:
				DATA[4] = vet
			# get NIR Data Band
			vet = readFileBand(filename, 2)
			if(not len(vet)):
				fileError = True
			else:
				DATA[5] = vet
			# get WIR Data Band
			vet = readFileBand(filename, 3)
			if(not len(vet)):
				fileError = True
			else:
				DATA[6] = vet

	# --------- Process CLM file
	if(fileError == False):
		#print('Process CLM file')
		filename = inputSubDir.replace('___', 'CLM') + str(date) + str(hour) + '_CLM.tif'
		if os.path.isfile(filename) == 0:
			print('..Error! Invalid data file: %s' %filename)
			fileError = True
		else:
			print('.Read CLM file: %s'%(filename))
			vet = readFileBand(filename, 1)
			if(not len(vet)):
				fileError = True
			else:
				DATA[7] = vet
	if(fileError == True):
		return []
		print('...fileError')
	else:
		return DATA

class corr_atm():

	'''
		|Variable | Min value | Max value | Step   | Num Step| Mult |
		| VZA     | 55        | 65        | 2.5    | 4       | 10   |
		| VAZ     | 100       | 110       | 2.5    | 4       | 10   | 100 a 112.5 ?
		| SZA     | 0         | 60        | 10     | 6       | 10   |  0 a 700
		| SAZ     | 0         | 360       | 60     | 6       | 1    |  0 a 420
		| WV      | 0         | 7.5       | 1.5    | 5       | 10   |  0 a  90
		| OZ      | 0         | 0.30      | 0.06   | 5       | 1000 |  0 a 360 
		| AOT     | 0         | 0.9       | 0.1    | 9       | 1000 |  0 a 900
	'''

	def __init__(self):

		self.VARS = ['VZA', 'VAZ', 'SZA', 'SAZ', 'WV', 'OZ', 'AOT', 'XA_RED', 'XB_RED', 'XC_RED', 'XA_NIR', 'XB_NIR', 'XC_NIR', 'XA_WIR', 'XB_WIR', 'XC_WIR']
		self.START = [[550, 651, 25, 10], [1000, 1101, 25, 10], [0, 701, 100, 10], [0, 361, 60, 1], [0, 91, 15, 10],[0, 361, 60, 1000], [0, 901, 100, 1000]]		
		# Order VZA VAZ SZA SAZ WV OZ AOT XA_RED XB_RED XC_RED XA_NIR XB_NIR XC_NIR
		# Order 550 1050 600 2200 30 250 200 516 10170 9516 710 6016 6533
		tac = time.clock()
		start_table = True
		for vza in (self.drange(self.START[0][:3])):
			path = os.getcwd() + '\\LUT_6S_3bands\\LUT_6S_VZA_' + str(vza) + '*.csv'
			files = glob.glob(path)
			for file in files:
				print('.Reading data file: %s'%(file))
				if(start_table == True):
					DATA = numpy.loadtxt(file, dtype='int32', delimiter = ' ', skiprows = 1)
					start_table = False
				else:
					DATA = numpy.append(DATA, numpy.loadtxt(file, dtype='int32', delimiter = ' ', skiprows = 1), axis = 0)
				#print('DATA [%d][%d]' %(len(DATA), len(DATA[0])))
		print('..time to get data: %.2f'%(time.clock()-tac))
		DATA = DATA.T
		for name, i in zip(self.VARS, range(len(self.VARS))):
			globals()[name] = numpy.vstack(DATA[i])
		DATA = []
		#print('VZA [%d][%d]' %(len(VZA), len(VZA[0])))
		print('..time to org data: %.2f'%(time.clock()-tac))
		self.error = 0
		self.warn_aot = 0
		self.warn_aot_value = 0

	def adjust_parameters(self, PAR):

		# Put input parameters in range of data values of look-up tables
		for i in range(len(PAR)):
			par = PAR[i]*self.START[i][3]
			new_par = -9999
			#print('i: %d par: %d'%(i, par))
			for lut in self.drange(self.START[i][:3]):
				#print('i: %d par: %.2f lut: %.2f'%(i, par, lut))
				#print('lut: %d (par - lut): %d dif: %d' %(lut, numpy.abs(par - lut), self.START[i][2]))
				#print('i: %d par: %f lut: %f dif: %f comp: %f'%(i, par, lut, numpy.abs(par - lut), self.START[i][2]/2))
				if(numpy.abs(par - lut) <= self.START[i][2]/2):
					new_par = lut
					break
			if(i == 6 and par > 950):
				new_par = 900
				if(self.warn_aot_value == 0  or self.warn_aot_value != par):
					print('. Warning! AOT out of range [0 900] (%d)'%(par))
					self.warn_aot_value = par
				self.warn_aot = par
			if(new_par == -9999):
				print('. Error! Verify data in look-up tables (variable %d value: %f lut: %f): '%(i, PAR[i], lut))
				print(PAR)
				self.error = i
			PAR[i] = new_par
		return PAR

	def calc_corr(self, PAR, red_rad, nir_rad, wir_rad):
		#print('...calc_corr')
		VARS_sample = []
		for i in range(len(self.VARS[:7])):
			VARS_sample.append(self.VARS[i] + '_sample')
		for name, i in zip(VARS_sample, range(len(VARS_sample))):
			globals()[name] = PAR[i]
		paser = numpy.logical_and(VZA == VZA_sample, numpy.logical_and(VAZ == VAZ_sample, numpy.logical_and(SZA == SZA_sample, \
				 numpy.logical_and(SAZ == SAZ_sample, numpy.logical_and(WV == WV_sample, numpy.logical_and(OZ == OZ_sample, \
				 AOT  == AOT_sample))))))
		samples = len(paser[paser == True])
		#print('samples: %d' %(samples))
		if(samples != 1):
			print('Error! not have sample to this data!')
			return [-1, self.error, samples]

		'''
		else:
			print('sample_data')
			print('RED_DATA: %d %d %d' %(XA_RED[paser], XB_RED[paser], XC_RED[paser]))
			print('NIR_DATA: %d %d %d' %(XA_NIR[paser], XB_NIR[paser], XC_NIR[paser]))
		'''
		red_ref = self.calc_rad(red_rad, XA_RED[paser], XB_RED[paser], XC_RED[paser])
		nir_ref = self.calc_rad(nir_rad, XA_NIR[paser], XB_NIR[paser], XC_NIR[paser])
		wir_ref = self.calc_rad(wir_rad, XA_WIR[paser], XB_WIR[paser], XC_WIR[paser])
		#print('[red, nir]: %.4f %.4f'%(red_ref, nir_ref))
		return [red_ref, nir_ref, wir_ref]
	
	def calc_rad(self, L, xa, xb, xc):
		L = L
		xa = xa/100000
		xb = xb/100000
		xc = xc/100000
		#print('L: %.4f xa: %.4f xb: %.4f xc: %.4f'%(L, xa, xb, xc))
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
	
	if ("-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict):
		dateStart = int(argDict["-ds"])
		dateEnd = int(argDict["-de"])
		outputDir = argDict["-fo"]
		inputDir  = argDict["-fi"]
		typeData = argDict["-ty"]
	else:
		exit(usage)
	six = corr_atm()
	HOURS = [1300, 1315, 1330, 1345, 1400, 1415, 1430, 1445, 1500]
	ds = gdal.Open(fileMask123)
	y_123 = ds.RasterYSize
	x_123 = ds.RasterXSize
	ang_error = numpy.zeros((3000,8))
	warn_aot  = numpy.zeros((10000,1))
	# Process in range of data
	date = dateStart
	while(date < dateEnd):
		# ATM = [WV OZ AOT][y_123][x_123]    
		ATM = importATMdata(inputDir, date)
		six.warn_aot_value = 0
		if(len(ATM) == 0):
			print('...Error to read ATM data: %s'%(inputDir + str(date) + '_ATM.tif'))
			date = incDate(date)
			continue
		for hour in HOURS:
		#for hour in random.sample(HOURS, 2):
			#Output file name 
			outputFolder = outputDir + '6S_3bands\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
			outputFile1 = outputFolder + str(date) +  str(hour) + '_NDVI.tif'    
			outputFile2 = outputFolder + str(date) +  str(hour) + '_NDWI.tif'    
			# Verify if files exist
			if(os.path.isfile(outputFile1) and os.path.isfile(outputFile2)):
				print('...Data file exist in: %s'%(outputFile1))
				print('...Data file exist in: %s'%(outputFile2))
				continue
			else:
				print('...Create 6S data files: %s'%(str(date) + str(hour)))  

			# Verify and create folders
			if(not os.path.isdir(outputFolder)):
				os.makedirs(outputFolder)
			OUTDATA_1 = numpy.zeros((3, y_123, x_123))
			OUTDATA_2 = numpy.zeros((3, y_123, x_123))
			# DATA[8][y_123][x_123] = [VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
			DATA = importBRDFdata(inputDir, date, hour)
			#saveRASTERfile('c:\\tmp\\'+str(date)+'_test.tif', fileMask123, DATA)
			if(len(DATA) == 0):
				print('...Error to read 6S process data: %s'%(str(date)))
				continue
			#print('hours:%d [%d, %d]'%(hours, y_123, x_123))
			tac = time.clock()
			for x in range(x_123):
				for y in range(y_123):
					if(DATA[0][y][x] <= 0 or DATA[7][y][x] > 1.05):
						continue
					'''	
					saveRASTERfile('c:\\tmp\\teste_atm.tif', fileMask123, ATM)
					saveRASTERfile('c:\\tmp\\teste_angs.tif', fileMask123, DATA)
					print('... testa ai')
					sys.exit()
					'''
					# PAR is parameters to process atmospheric correction
					# PAR = [VZA, VAZ, SZA, SAZ, WV, OZ, AOT]
					#print('x: %d y: %d ATM: %f %f %f'%(x, y, ATM[0][y][x+1], ATM[1][y][x+1], ATM[2][y][x+1]))
					# Necessary correct angle position. It's necessary reprocess this data!
					PAR = [DATA[0][y][x], DATA[1][y][x], DATA[2][y][x], DATA[3][y][x], ATM[0][y][x], ATM[1][y][x], ATM[2][y][x]]
					#print(PAR)
					#print('PAR antes:')
					ADJ = six.adjust_parameters(PAR)
					[red, nir, wir] = six.calc_corr(ADJ, DATA[4][y][x], DATA[5][y][x], DATA[6][y][x])
					#print('RAD >> red: %f nir: %f' %(DATA[4][y][x+1], DATA[5][y][x+1]))
					#print('REF >> RED: %f NIR: %f'%(red, nir))
					#sys.exit()
					#print(ADJ, DATA[4][y][x], DATA[5][y][x])
					#print([red, nir])
					#sys.exit()
					if(red == -1):
						ang_error[int(DATA[nir][y][x]), int(nir)] = ang_error[int(DATA[nir][y][x]), int(nir)] + 1
						'''
						if(ang_error[int(DATA[nir][y][x]), int(nir)] > 10):
							numpy.savetxt("c:\\tmp\\6s_error.csv", ang_error, delimiter=",",fmt='%d')
							sys.exit()
						'''
						print('var: %d ang: %d samples: %d'%(nir, DATA[nir][y][x],wir))
						print('PAR antes:')
						print([DATA[0][y][x], DATA[1][y][x], DATA[2][y][x], DATA[3][y][x], ATM[0][y][x], ATM[1][y][x], ATM[2][y][x]])
						print('PAR depois:')
						print(ADJ)
					else:
						# Get NDVI data
						OUTDATA_1[0][y][x] = red
						OUTDATA_1[1][y][x] = nir
						OUTDATA_1[2][y][x] = calcNDVI(nir, red)
						# Get NDWI data
						# See in http://edo.jrc.ec.europa.eu/documents/factsheets/factsheet_ndwi.pdf
						OUTDATA_2[0][y][x] = nir
						OUTDATA_2[1][y][x] = wir
						OUTDATA_2[2][y][x] = calcNDVI(nir,wir)
					#sys.exit()
					if(six.warn_aot != 0):
						warn_aot[int(six.warn_aot)] = warn_aot[int(six.warn_aot)] + 1
			saveRASTERfile(outputFile1, fileMask123, OUTDATA_1)
			saveRASTERfile(outputFile2, fileMask123, OUTDATA_2)
			numpy.savetxt("c:\\tmp\\6s_warn_aot.csv", warn_aot, delimiter=",",fmt='%d')
			numpy.savetxt("c:\\tmp\\6s_error.csv", ang_error, delimiter=",",fmt='%d')
			print('..time to process: %d s'% int(time.clock()-tac))
		date = incDate(date)