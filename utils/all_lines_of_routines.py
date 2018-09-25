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
import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import subprocess
from subprocess import Popen, PIPE
import os
import glob
import datetime
import numpy
import gc
import shutil
import sys
import copy
import gdal
from osgeo import gdal, osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from time import  strftime, localtime, sleep
from utils import printlog, sleepTime, mapDict, sleepTime, delFile, infoFile, incHour, incDate
from osgeo import gdal, osr
import math
import copy

from sys import argv

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to save output files (optional, default is \\files\\)
        -ds     date to start process (format: yyyymmdd)
        -de     date to end process   (format: yyyymmdd) 
        -hs     hour to start process (format: hhMM)
        -he     hour to end process   (format: hhMM)
""" % argv[0]

# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMask = fileGoias
angleDir = os.getcwd() +  '\\angle\\'
#exeDir = 'C:\\Python34\\Lib\\site-packages\osgeo\\'
exeDir = 'C:\\OSGeo4W\\bin\\'
outputDir = os.getcwd() + '\\file\\'
tmpDir = 'C:\\tmp\\'



def timeAngle(day, hour, form):
    day = str(day)
    hour = str(hour)
    minuto = 10*int(hour[2])+int(hour[3]) 
    if(minuto == 00):
        hour = hour[:2]+'.00'
    else:
        if(minuto <= 15):
            hour = hour[:2]+'.15'
        else:
            if(minuto <= 30):    
                hour = hour[:2]+'.30'
            else:
                hour = hour[:2]+'.45'   
    
    if(form == 'SPACE'):
       string = day[:4]+' '+day[4]+day[5]+' '+day[6]+day[7]+' '+ hour
    if(form == 'UNDER'):
       string = '_' + day[:4]+'_'+day[4]+day[5]+'_'+day[6]+day[7]+'_'+ hour
    return string

    
def array2raster(outFile,array):

    #Xorigin = -66.96422963471274   # Oeste 
    #Xend    = +66.96364537142968   # Leste
    #Yend    = -62.27591109399421   # South 
    #Yorigin = +62.27979992238434   # North
    #Xpixel  = 0.034840758326259734
    #Ypixel  = 0.034840758326259734

    Xpixel = 0.45
    Ypixel = 0.45
 
    cols = 378
    rows = 378

    Xorigin = -85 #- Xpixel/2
    Yorigin = -85 #- Ypixel/2

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outFile, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((Xorigin, Xpixel, 0, Yorigin, 0, Ypixel))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    # this assumes the projection is Geographic lat/lon WGS 84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outRaster.SetProjection(srs.ExportToWkt())
    outband.FlushCache()
    outRaster = None
    infoFile(outFile)

def transformAngle(outputDir, day, hour, vetAngle):

    timeprocess = str(day) + str(hour)
    tmpfile1 = tmpDir + timeprocess + '_tmp1.tif'
    tmpfile2 = tmpDir + timeprocess + '_tmp2.tif'
    tmpfile3 = tmpDir + timeprocess + '_tmp3.tif'
    tmpfile4 = tmpDir + timeprocess + '_tmp4.tif'
    tmpfile5 = tmpDir + timeprocess + '_tmp5.tif'
    fileIn   = vetAngle[0]
    fileOut  = vetAngle[1]
    minValue = vetAngle[2]
    maxValue = vetAngle[3]
    resampling_method = vetAngle[4]
       
    print(True,'transform Angle: ' + fileIn)

    arrayTemp = numpy.fromfile(fileIn, dtype='>f')
    matrixTemp = numpy.reshape(arrayTemp,(-1,378))
    rows = 378
    cols = 378
    matrix = numpy.zeros((rows, cols), dtype = float)
    for i in range (0, rows, 1):
        for j in range (0, cols, 1):
            temp = float("{0:.5f}".format(matrixTemp[i][j]))
            if ('sunzen' in fileIn) and math.isnan(temp):
                #print(False,'...none: ' + str(temp))
                temp = 90.0001
            if temp > maxValue:
                printlog(False, '...maior: ' + str(temp))
                if ('sunzen' in fileIn): 
                    temp = 90.0000
                else:
                    temp = maxValue
            if temp < minValue:
                printlog(False, '...menor:'  + str(temp))
                temp = minValue
            matrix[i][j]=temp
    del arrayTemp
    del matrixTemp
    array2raster(tmpfile1, matrix)
    
    print(exeDir + 'gdalwarp -tr 0.034840758326259734 0.034840758326259734 -r  bilinear ' + ' ' + tmpfile1 + ' ' + tmpfile2)
    subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326259734 0.034840758326259734 -r  bilinear ' + ' ' + tmpfile1 + ' ' + tmpfile2)

    print(exeDir + 'gdalwarp -cutline ' + fileDataRetriever + ' -crop_to_cutline -overwrite '  + tmpfile2 + ' ' + tmpfile3)
    subprocess.call(exeDir + 'gdalwarp -cutline ' + fileDataRetriever + ' -crop_to_cutline -overwrite '  + tmpfile2 + ' ' + tmpfile3)
    
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)

    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + fileOut)
    subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + fileOut)

    #print(fileOut)
    #sys.exit()
    delFile(tmpDir + timeprocess + '_tmp' + '*')


def create1file(inputFile, outputFile):
    # Read input and output data
    printlog(False, '...create1file:')
    printlog(False, '..inputFile: %s' % inputFile)
    try:
        dsIn = gdal.Open(inputFile[0])
    except(RuntimeError, e):
        raise RuntimeError('...Error: file not found!')
    bandIn = dsIn.GetRasterBand(1)
    Nbands = len(inputFile)
    printlog(False, '..outputFile: %s' % outputFile)
    driver = gdal.GetDriverByName("GTiff")
    delFile(outputFile)
    dsOut  = driver.Create(outputFile, dsIn.RasterXSize, dsIn.RasterYSize, Nbands, bandIn.DataType)
    CopyDatasetInfo(dsIn,dsOut)
    # Write output file
    for x in range(0, len(inputFile), 1):
        try:
            dsIn = gdal.Open(inputFile[x])
        except(RuntimeError, e):
            raise RuntimeError('...Error: file %s not found!' % inputFile[x])
        bandIn = dsIn.GetRasterBand(1)
        arrayIn = BandReadAsArray(bandIn)
        bandOut = dsOut.GetRasterBand(x+1)
        BandWriteArray(bandOut, arrayIn)
        arrayIn = None
        bandIn = None
        dsIn = None
        delFile(inputFile[x])
        
class main():

    argDict = mapDict(argv, usage)

    if "-ds" in argDict and "-de" in argDict and "-hs" in argDict and "-he" in argDict:

        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        hourStart = int(argDict["-hs"])
        hourEnd = int(argDict["-he"])
        if "-fo" in argDict:
            outputDir = argDict["-fo"]
    else:
        exit(usage)
    gc.collect()
    gdal.UseExceptions()
    #printlog(True, 'os.getcwd() : ' + os.getcwd())
    inputFile = ['test', 'test', 'test', 'test']
    arrayAngle =  [['sataz','msg_azres',0, 360, 'average'], ['sunaz','sol_azres',0, 360, 'average'], ['satzen', 'msg_zenres', 0, 99, 'average'], ['sunzen', 'sol_zenres', 0, 98.9998, 'average']]
    day = dateStart
    
    while (day <= dateEnd):
        year  = str(day)[0:4]
        month = str(day)[4:6]
        subfolder = year + '\\' + month + '\\'
        if not(os.path.isdir(outputDir + subfolder)):
            os.makedirs(outputDir + subfolder)
        hour = hourStart
        while(hour <= hourEnd):
            outputFile = outputDir + subfolder + str(day) + str(hour) + '_angles.tif'
            if os.path.isfile(outputFile):
                printlog(True,'File exist: ' + str(day) + str(hour) + '_angles.tif')
                hour = incHour(hour)
                continue
            arrayAngle =  [['sataz','msg_azres',0, 360, 'average'], ['sunaz','sol_azres',0, 360, 'average'], ['satzen', 'msg_zenres', 0, 99, 'average'], ['sunzen', 'sol_zenres', 0, 98.9998, 'average']]
            matAngle = copy.copy(arrayAngle)
            #printlog(False, 'matAngle:')
            #printlog(False, matAngle)
            #printlog(False, 'arrayAngle:')
            #printlog(False, arrayAngle)

            for x in range(0, len(matAngle), 1):
                delFile(os.getcwd()  + '\\' + matAngle[x][0] + '_*')
                matAngle[x][0] += timeAngle(day, hour,'UNDER')
                matAngle[x][1] = str(day) + str(hour) + '_' + matAngle[x][1] + '.tif'
                inputFile[x] = os.getcwd() + '\\' + matAngle[x][1]
            
            #ans = check_output('java.exe -cp .;C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle\\operation.jar;C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle\\Jama-1.0.1.jar;C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle AngleMaps '+ timeAngle(day, hour,'SPACE'))   
            printlog(False, 'java -cp .;' + angleDir + 'operation.jar;' + angleDir  + 'Jama-1.0.1.jar;' + os.getcwd() +'\\angle AngleMaps ' +  timeAngle(day, hour,'SPACE'))
            try:
                subprocess.call('java -cp .;' + angleDir + 'operation.jar;' + angleDir  + 'Jama-1.0.1.jar;' + os.getcwd() +'\\angle AngleMaps ' +  timeAngle(day, hour,'SPACE'))
            except:
                printlog(True, '..Fail Generate Azimuth Files')
                sys.exit()
            for x in range(0, len(matAngle), 1):
                transformAngle(outputDir , day, hour, matAngle[x])
            printlog(False, 'inputFile: %s'  % inputFile)
            printlog(False, 'outputFile: %s' % outputFile)
            create1file(inputFile, outputFile)
            #infoFile(outputFile)
            for x in range(0, len(matAngle), 1):
                delFile(inputFile[x])
            hour = incHour(hour)
        day = incDate(day)
import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import subprocess
import os
import glob
import datetime
import gc
import subprocess
import shutil
import copy
import tarfile
from utils import *
from sys import argv
import numpy
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import threading
import queue
import time
import shutil

usage = """\
Usage: %s [OPTIONS]
        -ty     type of data to process (123, CLM, NDVI, MOD08, MYD08, MOD09, HRV, MOD13)
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)

""" % argv[0]


# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoiasGeo = os.getcwd() + '\\shapes\\201401011300_123.tif'
fileMaskRect   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMask = fileGoias

pathlog = os.getcwd() + '\\logs\\'
MRTDir =  os.getcwd() + '\\MRT\\'
dataRetrieve = os.getcwd() + '\\MSGDataRetriever\\'
pyModisDir = os.getcwd() + '\\\pyModis-2.0.5\\scripts'
exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
gdalIlwisDir = "C:\\Ilwis372\\Extensions\\Geonetcast-Toolbox\\GDAL\\bin\\" 
utilDir = 'C:\\ILWIS372\\Extensions\\GEONETCast-Toolbox\\util\\'
ilwDir = 'C:\\ilwis372\\'
osgeoDir = 'C:\\OSGeo4W\\bin\\'

tmpDir    = 'c:\\tmp\\'
logTarError = tmpDir + 'logTarError.txt'

threadList = ["Thread-1"]#, "Thread-2", "Thread-3"] #, "Thread-4", "Thread-5"]
queueLock = threading.Lock()
workQueue = queue.Queue(12)
threads = []
exitFlag = 0

class myThread (threading.Thread):
    def __init__(self, threadID, name, q, typeData, tmpDir):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        self.typeData = typeData
        self.tmpDir = tmpDir
    def run(self):
        print("Starting " + self.name)
        #print( "q: " + str(self.q) + " typeData: " + self.typeData + " donwloadDir: " + self.downloadDir)
        process_data(self.name, self.q, self.typeData, self.tmpDir)
        print("Exiting " + self.name)

def process_data(threadName, q, typeData, tmpDir):
    global queueLock
    global workQueue
    global exitFlag
    #print( "Enter in process_data. exitFlag: " + str(exitFlag))

    while not exitFlag:
        #print( "exitFlag is False")
        queueLock.acquire()
        #print( "\nqueueLock acquire")
        if not workQueue.empty():
            #print( "Queue is not empty")
            inputSubdir = q.get()
            outputSubdir = q.get()
            tiffile = q.get()
            queueLock.release()
            #print( "\nqueueLock release")
            print("\n%s processing %s" % (threadName, tiffile))
            # Generate file
            flagGenFile = False
            if(typeData == '123' or typeData == 'HRV'):
                flagGenFile = generateBands(tmpDir, inputSubdir, outputSubdir, tiffile, typeData)
            if(typeData == 'CLM'):        
                [flagGenFile, error] = generateCloudMask(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'NDVI'):
                flagGenFile = generateNDVI(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD08' or typeData == 'MYD08'):
                flagGenFile = generateMOD08(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD09' or typeData == 'MYD09'):
                flagGenFile = generateMOD09(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD13'):
                flagGenFile = generateMOD13(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD15'):
                flagGenFile = generateMOD15(tmpDir, inputSubdir, outputSubdir, tiffile)


            if(flagGenFile == True):    
                print( '... File generate type: ' + outputSubdir + tiffile)
            else:
                print( '... Error to generate file: ' + tiffile)
                print("FileExtractError: %s error: %d \n" % (tiffile, flagGenFile))
                logFileError = os.getcwd() + '\\logExtractError.txt'
                logfile = open(logFileError, "a+")
                logfile.write("FileExtractError: %s error: %d \n" % (tiffile, flagGenFile))
                logfile.close()
                sys.exit()
                if(typeData == 'CLM'):
                    print('Error ' + str(error) + ' in ' + typeData + 'file:'  + tiffile)        
                # Here is necessary start another process to get this data
                #sys.exit()
            
        else:
            #print( "Queue is empty")
            try:
                queueLock.release()
                #print( "\nqueueLock release")
            except:
                print( 'Fail to release queueLock !')
        time.sleep(1)
    #print( "exitFlag is True")

def generateBands(tmpDir, inputDir, outputDir, file, typeData):

    hour = file[8:12]
    if(hour[2:] == '45'): hour2 = hour[:2] + '57'
    if(hour[2:] == '30'): hour2 = hour[:2] + '42'
    if(hour[2:] == '15'): hour2 = hour[:2] + '27'
    if(hour[2:] == '00'): hour2 = hour[:2] + '12'
    filename = file[:file.rfind(hour)] + hour2

    #print('file: ' + file)
    #print('filename: ' + filename)
    # Find files
    files = glob.glob(inputDir + '*' + filename + '*.*')
    if(len(files) == 0):
        print( 'Not have this file in subfolder: ' + inputDir)
        return False
    filename = ''
    for x in range(len(files)):
        filename = files[x]
        # extract tar file
        if(tarfile.is_tarfile(filename)):
            print('. Trying extract file: ' + filename)
            if(extractTarfile(filename, tmpDir) == False):
                print('... File extract Error: ' + filename)
                print('save log in: %s' %logTarError)
                logfile = open(logTarError, "a+")
                logfile.write("FileExtractError: %s\n" % filename)
                logfile.close()
                if(x == (len(files)-1)):
                    return False
            else:
                break
        else: 
            print('... Not is tar file: ' + filename)
            logfile = open(logTarError, "a+")
            logfile.write("NotIsTarFile: %s\n" % filename)
            logfile.close()
            if(x == (len(files)-1)):
                return False
    # Process file
    tmpfile1   = tmpDir + 'temp1' + file + '.tif'
    tmpfile2   = tmpDir + 'temp2' + file + '.tif'
    delFile(tmpfile1)
    delFile(tmpfile2)    
    outputFile = outputDir + file + '_' + typeData + '.tif'
    print('. GenerateFile outputFile: ' + outputFile)
    if(typeData == '123'):

        # Radiance without cut
        # Use option L to radiance and T to reflectance 
        print(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',(1,2,3),Y,L,1,1) ' + tmpfile1)
        try:
            subprocess.call(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',(1,2,3),Y,L,1,1) ' + tmpfile1)
        except:
            return False

        print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        try:
            subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        except:
            return False

        print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        except:
            return False
   

    else:
        # Radiance without cut
        print(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',12,Y,L,1,1) ' + tmpfile1)
        try:
            subprocess.call(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',12,Y,L,1,1) ' + tmpfile1)
        except:
            return False

        print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        try:
            subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        except:
            return False

        print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        except:
            return False       
    delFile(tmpfile1)
    delFile(tmpfile2)
    delFile(tmpDir + '*' + file + '*')
    #sys.exit()
    return True

def generateCloudMask(tmpDir, inputDir, outputDir, file): # reviewed

    tmpfile1   = tmpDir + file + 'temp1_clm.tif'
    tmpfile2   = tmpDir + file + 'temp2_clm.tif'
    tmpfile3   = tmpDir + file + 'temp3_clm.tif'
    delFile(tmpfile1)
    delFile(tmpfile2)
    outputFile = outputDir + file + '_CLM.tif'
    files = glob.glob(inputDir + '*CLM*' + file + '*')
    print(inputDir + '*CLM*' + file + '*')
    if(len(files) > 1):
        print('Many files to time: ' + file)
    elif(len(files) == 0):
        print('Not have this file in subfolder: ' + inputDir)
        return False
    for x_files in range(len(files)):
        filename = files[0]      
        # Reproject image to WGS84
        print(exeDir + 'gdalwarp  -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" ' + filename + ' ' + tmpfile1)  
        try:
            subprocess.call(exeDir + 'gdalwarp -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" ' + filename + ' ' + tmpfile1)
        except:
            if(x_files < range(len(files))):
                continue
            else:
                return [False, -1]
        
        # Cut to dataRetrieve shape
        print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "bilinear" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
        try:
            subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "bilinear" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
        except:
            if(x_files < range(len(files))):
                delFile(filename)
                continue
            else:
                return [False, -2]
        # Translate data to rectangle Goias shape file                                      
        print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
        try:
            subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
        except:
            if(x_files < range(len(files))):
                continue
            else:
                return [False, -3]

        # Data cut using mask file
        print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outputFile)
        except:
            if(x_files < range(len(files))):
                continue
            else:
                return [False, -4]
    # Remove temp files
    delFile(tmpfile1)
    delFile(tmpfile2)
    delFile(tmpfile3)
    #sys.exit()
    return [True, 0]

def generateNDVI(tmpDir, inputDir, outputDir, file):
    files = glob.glob(inputDir + '*MSG3-SEVI-MSGNDVE*' + file + '*')
    if(len(files) > 1):
        print('Many files to time: ' + file)
    elif(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir + file)
        return False
    fileName = files[0]
    print('file: %s' % file)
    delFile(tmpDir + '*' + file + '*')
    print('fileName:' + fileName)
   
    print(gdalIlwisDir +  'gdal_translate.exe -of ILWIS HDF5:' + fileName + '://NDVImax '  +  tmpDir +  'tNDVImax'  +  file)
    try:
        subprocess.call(gdalIlwisDir +  'gdal_translate.exe -of ILWIS HDF5:' + fileName + '://NDVImax '  +  tmpDir +  'tNDVImax'  +  file)
    except:
        return False
    print(ilwDir + 'ilwis.exe -C ' + tmpDir + 'NDVImax_'   + file +'.mpr:=iff(' + tmpDir + 'tNDVImax'  + file +' le 100,' + tmpDir + 'tNDVImax'  + file + '/100,?);')
    try:
        subprocess.call(ilwDir + 'ilwis.exe -C ' + tmpDir + 'NDVImax_'   + file +'.mpr:=iff(' + tmpDir + 'tNDVImax'  + file +' le 100,' + tmpDir + 'tNDVImax'  + file + '/100,?);')
    except:
        return False
    # Convert study and define..
    print(ilwDir + 'ilwis.exe -C setgrf ' + tmpDir + 'NDVImax_'  + file + '.mpr ' + utilDir +'lritmsg;')
    try:
        subprocess.call(ilwDir + 'ilwis.exe -C setgrf ' + tmpDir + 'NDVImax_'  + file + '.mpr ' + utilDir +'lritmsg;')
    except:
        return False
    # Traduz imagem para coordenadas latitude e longitude
    print(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5568000 5568000 5568000 -5568000 ' + tmpDir + 'ndvimax_'  + file +'.mpr ' + tmpDir + '\\ndvimax_'+ file + '.tif')
    try:
        subprocess.call(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5568000 5568000 5568000 -5568000 ' + tmpDir + 'ndvimax_'  + file +'.mpr ' + tmpDir + file + '_tNDVImax.tif')
    except:
        return False

    # Reprojeta imagem para WGS84 near/bilinear/cubic
    print(exeDir  + 'gdalwarp  -tr 0.0348420636891 0.0348426244628 -r near -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask + ' -crop_to_cutline -dstalpha ' + tmpDir + file +'_tNDVImax.tif ' + outputDir +  file + '_NDVI.tif')
    try:
        subprocess.call(exeDir + 'gdalwarp  -tr 0.0348420636891 0.0348426244628 -r near -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask + ' -crop_to_cutline -dstalpha ' + tmpDir + file +'_tNDVImax.tif ' + outputDir +  file + '_NDVI.tif')
    except:
        return False


    # Remove temporal files
    delFile(tmpDir + '*' + file + '*')
    return True

def generateMOD08(tmpDir, inputDir, outputDir, file):

    tmpfile1 = tmpDir + file[:file.rfind('.')] + '_tmp1.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    outputFile = outputDir + file

    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')
    fileMODIS = file[0:4] + str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    
    files = glob.glob(inputDir + 'M?D08_D3.A' + fileMODIS + '*.hdf')

    if(len(files) > 1):
        print('Many files to time: ' + file)
    elif(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir)
        #print('findfile: '+ inputDir + '*' + file + '*')
        #print(files)
        #print(file)
        #sys.exit()
        return False
    filename = files[0]
    '''
    if('MOD08' in fileName):
        typeFile = '_MOD08.tif'
    else:
        typeFile = '_MYD08.tif'
    '''    
    #band1_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Surface Reflectance Band 1'
    #band2_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Surface Reflectance Band 2'
    #band3_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Ozone'

    band1_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57
    band2_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37
    band3_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Total_Ozone_Mean' # Band 829
    band4_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Atmospheric_Water_Vapor_Mean' # Band 853

    #print('-------')
    #print(band1_sub)
    #print('-------')
    #print(band2_sub)
    #print(fileName)
    #sys.exit(1)
        
    Xorigin = -180 
    Yorigin = 90
    Xpixel = 1
    Ypixel = -1
    cols = 180
    rows = 360
    
    # Read data
    try:
        ds1 = gdal.Open(filename)
    except:
        print('Error to open file %s' %(filename))
        return False
    subdata = ds1.GetSubDatasets()

    band1 = gdal.Open(band1_sub)
    band2 = gdal.Open(band2_sub)
    band3 = gdal.Open(band3_sub)
    band4 = gdal.Open(band4_sub)

    BAND1 = band1.ReadAsArray()
    BAND2 = band2.ReadAsArray()
    BAND3 = band3.ReadAsArray()
    BAND4 = band4.ReadAsArray()
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(tmpfile1, rows, cols, 4, gdal.GDT_Int16)
    outRaster.SetGeoTransform((Xorigin, Xpixel, 0, Yorigin, 0, Ypixel))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outRaster.SetProjection(srs.ExportToWkt())
    bandOut1 = outRaster.GetRasterBand(1)
    bandOut2 = outRaster.GetRasterBand(2)
    bandOut3 = outRaster.GetRasterBand(3)
    bandOut4 = outRaster.GetRasterBand(4)
    
    bandOut1.WriteArray(BAND1)
    bandOut2.WriteArray(BAND2[1])
    bandOut3.WriteArray(BAND3)
    bandOut4.WriteArray(BAND4)
    
    bandOut1.FlushCache()
    bandOut2.FlushCache()
    bandOut3.FlushCache()
    bandOut4.FlushCache()
    outRaster = None
    # Cut to dataRetrieve shape
    print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r near -srcnodata [-9999] -dstnodata [-9999] -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
    try:
        subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r near -srcnodata [-9999] -dstnodata [-9999] -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
    except:                                  # Dimension in 123 file: 0.034840758326260 0.034840758326260
        return False
    # Translate data to rectangle Goias shape file                                      
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
    except:
        return False

    # Data cut using mask file
    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outputFile)
    except:
        return False
    # Remove temp files
    delFile(tmpfile1)
    delFile(tmpfile2)
    delFile(tmpfile3)
    #sys.exit()
    return True

def generateMOD09(tmpDir, inputDir, outputDir, file):

    tmpfile1a = tmpDir + file[:file.rfind('.')] + '_tmp1a.tif'
    tmpfile1b = tmpDir + file[:file.rfind('.')] + '_tmp1b.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    tmpfile4 = tmpDir + file[:file.rfind('.')] + '_tmp4.tif'
    tmpfile5 = tmpDir + file[:file.rfind('.')] + '_tmp5.tif'
    outputFile = outputDir + file
    if (os.path.isfile(outputFile)):
        return True
    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')
    doy = str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    fileMODIS = file[0:4] + doy
    find_files = inputDir + file[0:4] + '\\' + doy + '\\M?D09GA.A' + fileMODIS + '*.hdf'
    files = glob.glob(find_files)
    print('--------------- MOD09')
    print(find_files)
    print(files)

    if(len(files) == 2):
        if(not 'h13v10' in files[0] and not 'h12v10' in files[1]):
            file_tmp = files[0]
            files[0] = files[1]
            files[1] = file_tmp
    if(len(files) == 0):
        print('...Error not have files to process')
        return(0)

    band = ['":MODIS_Grid_500m_2D:sur_refl_b01_1', '":MODIS_Grid_500m_2D:sur_refl_b02_1', '":MODIS_Grid_1km_2D:state_1km_1', '":MODIS_Grid_500m_2D:QC_500m_1']
    data_type = ['Int16', 'Int16', 'Int32', 'Int32']
    logModisError = os.getcwd() + '\\logModError.txt'
    
    delFile(tmpfile3)
    delFile(tmpfile4)
    delFile(tmpfile5)  
    for i in range(len(band)):
        delFile(tmpfile1a)
        delFile(tmpfile1b)
        delFile(tmpfile2)
        print('--------------- warp')
        print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] + ' HDF4_EOS:EOS_GRID:"' + files[0] + band[i] + ' ' + tmpfile1a)
        try:
            subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] + ' HDF4_EOS:EOS_GRID:"' + files[0] + band[i] + ' ' + tmpfile1a)
        except:
            return (-1)
        if(len(files) == 2):
            print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] +  ' HDF4_EOS:EOS_GRID:"' + files[1] + band[i] + ' ' + tmpfile1b)
            try:
                subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] +  ' HDF4_EOS:EOS_GRID:"' + files[1] + band[i] + ' ' + tmpfile1b)
            except:
                return (-1)    
        
        if(len(files) == 2): 
            print('--------------- merge')
            print('python ' + osgeoDir + 'gdal_merge.py  -n -28672  -ot ' + data_type[i] + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py  -n -28672 -ot ' + data_type[i] + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            except:
                return (-2)
        else:
            os.rename(tmpfile1a, tmpfile2)

        if('MODIS_Grid_1km_2D' in band[i]):
            print('------------- reproject')
            BAND_1km = readFileBand(tmpfile2, 1)
            ds = gdal.Open(tmpfile3)
            y_len = ds.RasterYSize
            x_len = ds.RasterXSize
            ds = None
            ds = gdal.Open(tmpfile2)
            y_len_1km = ds.RasterYSize
            x_len_1km = ds.RasterXSize
            ds = None
            print('tmptile2 [%d %d] tmpfile3 [%d %d]'%(y_len_1km, x_len_1km, y_len, x_len))
            BAND_1km_rep = numpy.zeros((y_len, x_len))
            for x in range(x_len):
                for y in range(y_len):
                    try:
                        BAND_1km_rep[y][x] = BAND_1km[int(y/2)][int(x/2)]
                    except:
                        return (-3)  
            delFile(tmpfile2)
            saveRASTERfile(tmpfile2, tmpfile3, BAND_1km_rep)
        if(not os.path.isfile(tmpfile3)):
            os.rename(tmpfile2, tmpfile3)
        else:
            print('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type[i] + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type[i] + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            except:
                return (-4)
            try:
                delFile(tmpfile3)
                os.rename(tmpfile4, tmpfile3)
            except:
                return (-5)
    # Translate data to rectangle Goias shape file                                      
    print('--------------- translate')
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    except:
        return (-6)
    print('--------------- warp')
    # Data cut using mask file
    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + tmpfile5)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + tmpfile5)
    except:
        return (-7)
    print('--------------- calc NDVI')
    RED = numpy.float32(readFileBand(tmpfile5, 1))*0.0001
    NIR = numpy.float32(readFileBand(tmpfile5, 2))*0.0001
    QUAL1 = numpy.int32(readFileBand(tmpfile5, 3))
    QUAL2 = numpy.int32(readFileBand(tmpfile5, 4))

    ds = gdal.Open(tmpfile5)
    y_len = ds.RasterYSize
    x_len = ds.RasterXSize
    ds = None
    NDVI = numpy.zeros((y_len, x_len))
    # Filter data, verify in Table 2: 500-meter QA Descriptions (32-bit)
    # https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod09ga
    for x in range(x_len):
        for y in range(y_len):
            if(RED[y][x] > 0 and NIR[y][x] > 0):
                NDVI[y][x] = (NIR[y][x] - RED[y][x])/(NIR[y][x] + RED[y][x])
                if(NDVI[y][x] > 1 or NDVI[y][x] < 0):
                    NDVI[y][x] = 0
            else:
                NDVI[y][x] = 0

    DATA  = [NDVI, RED, NIR, QUAL1, QUAL2] 
    saved = saveRASTERfile(outputFile, tmpfile5, DATA)
    if(not saved):
        return (-8)
    # Remove temp files
    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    delFile(tmpfile3)
    delFile(tmpfile4)
    delFile(tmpfile5)

    return True

def generateMOD15(tmpDir, inputDir, outputDir, file):

    tmpfile1a = tmpDir + file[:file.rfind('.')] + '_tmp1a.tif'
    tmpfile1b = tmpDir + file[:file.rfind('.')] + '_tmp1b.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    tmpfile4 = tmpDir + file[:file.rfind('.')] + '_tmp4.tif'

    doy = str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    fileMODIS = file[0:4] + doy
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_MOD15.tif'

    if (os.path.isfile(outputFile)):
        return True
    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')

    find_files = inputDir + 'MCD15A2H.A' + fileMODIS + '*.hdf'

    files = glob.glob(find_files)
    print('--------------- MCD15AH2')
    print(find_files)
    print(files)

    if(len(files) == 2):
        if(not 'h13v10' in files[0] and not 'h12v10' in files[1]):
            file_tmp = files[0]
            files[0] = files[1]
            files[1] = file_tmp
    if(len(files) == 0):
        print('...Error not have files to process')
        return(0)

    # Data name:HDF4_EOS:EOS_GRID:"MCD15A2H.A2014001.h12v10.006.2015273193555"

    band = '":MOD_Grid_MOD15A2H:Fpar_500m'
    data_type = 'Int16'
    logModisError = os.getcwd() + '\\logModError.txt'
    
    delFile(tmpfile3)
    delFile(tmpfile4)

    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    print('--------------- warp')
    print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + band + ' ' + tmpfile1a)
    try:
        subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + band + ' ' + tmpfile1a)
    except:
        return (-1)
    if(len(files) == 2):
        print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + band + ' ' + tmpfile1b)
        try:
            subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + band + ' ' + tmpfile1b)
        except:
            return (-1)    
    
    if(len(files) == 2): 
        print('--------------- merge')
        print('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
        try:
            subprocess.call('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
        except:
            return (-2)
    else:
        os.rename(tmpfile1a, tmpfile2)

    # Translate data to rectangle Goias shape file                                      

    print(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
    except:
        return False

    # Data cut using mask file
    print(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + tmpfile4)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + tmpfile4)
    except:
        return False
    print('--------------- scale factor')

    FPAR = numpy.float32(readFileBand(tmpfile4, 1))*0.01

    saved = saveRASTERfile(outputFile, tmpfile4, FPAR)
    if(not saved):
        return (-8)
    # Remove temp files
    #sys.exit()
    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    delFile(tmpfile3)
    delFile(tmpfile4)

    return True

def generateMOD15_lapig(tmpDir, inputDir, outputDir, file):

    tmpfile1 = tmpDir + file[:file.rfind('.')] + '_tmp1.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    band_names = ['fpar']
    fileMODIS = file[0:4] + str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    #outputFile = outputDir + file
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_MOD15.tif'


    num_band  = 1
    for pos_band in range(len(band_names)):
        subname    = inputDir + 'pa_br_'+ band_names[pos_band] +'_1000_lapig'
        filesearch = inputDir + 'pa_br_'+ band_names[pos_band] +'_1000_lapig' + '\\' + 'pa_br_'+ band_names[pos_band] +'_1000_'+ fileMODIS + '_lapig.tif'
        files = glob.glob(filesearch)
        print(filesearch)
        print(files)
        if(len(files) > 1):
            print('Many files to time: ' + file)
        elif(len(files) == 0):
            print('No have this file in subfolder: ' + filesearch)
            if(pos_band == 0):
                return False

        else:
            filename = files[0]


            print(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
            try:
                subprocess.call(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
            except:
                return False

            # Data cut using mask file
            print(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
            try:
                subprocess.call(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
            except:
                return False


            VET_DATA = readFileBand(tmpfile3, 1)

            if(pos_band == 0):
               # Save georeference file
                #print('------------------')
                #print(num_band)
                #print(outputDir + file)
                ds = gdal.Open(tmpfile3)
                band = ds.GetRasterBand(1)
                driver = gdal.GetDriverByName("GTiff")
                dsOut  = driver.Create(outputFile, ds.RasterXSize, ds.RasterYSize, len(band_names), band.DataType)
                CopyDatasetInfo(ds,dsOut)
                bandOut = dsOut.GetRasterBand(num_band)
                BandWriteArray(bandOut, VET_DATA)
            else:
                num_band = num_band + 1
                #print('------------------')
                #print(num_band)
                bandOut = dsOut.GetRasterBand(num_band)
                BandWriteArray(bandOut, VET_DATA) 
            #print('---- go out to test...')
            delFile(tmpDir + file[:file.rfind('.')] + '_tmp?.tif')
    return True

def generateMOD13(tmpDir, inputDir, outputDir, file):

    tmpfile1 = tmpDir + file[:file.rfind('.')] + '_tmp1.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    #band_names = ['NDVI', 'RED', 'NIR', 'MIR', 'composite_day', 'pixel_reliability']
    band_names = ['NDVI']
    fileMODIS = file[0:4] + str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    #outputFile = outputDir + file
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_MOD13.tif'


    num_band  = 1
    for pos_band in range(len(band_names)):
        subname    = inputDir + 'pa_br_'+ band_names[pos_band] +'_250_lapig'
        filesearch = inputDir + 'pa_br_'+ band_names[pos_band] +'_250_lapig' + '\\' + 'pa_br_'+ band_names[pos_band] +'_250_'+ fileMODIS + '_lapig.tif'
        files = glob.glob(filesearch)
        print(filesearch)
        print(files)
        if(len(files) > 1):
            print('Many files to time: ' + file)
        elif(len(files) == 0):
            print('No have this file in subfolder: ' + inputDir + subname)
            #print('findfile: '+ inputDir + '*' + file + '*')
            #print(files)
            #print(file)
            #sys.exit()
            return False
        filename = files[0]

        # don't transform spatial resolution 
        #tmpfile1 = filename

        '''
        # Cut to dataRetrieve shape and transform 3 x 3 km
        print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r cubic -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + filename + ' ' + tmpfile1) 
        try:
            subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r cubic -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + filename + ' ' + tmpfile1) 
        except:                                  # Dimension in 123 file: 0.034840758326260 0.034840758326260
            return False
        '''
        print(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
        except:
            return False

        # Data cut using mask file
        print(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
        except:
            return False

        VET_DATA = readFileBand(tmpfile3, 1)

        if(pos_band == 0):
           # Save georeference file
            #print('------------------')
            #print(num_band)
            #print(outputDir + file)
            ds = gdal.Open(tmpfile3)
            band = ds.GetRasterBand(1)
            driver = gdal.GetDriverByName("GTiff")
            dsOut  = driver.Create(outputFile, ds.RasterXSize, ds.RasterYSize, len(band_names), band.DataType)
            CopyDatasetInfo(ds,dsOut)
            bandOut = dsOut.GetRasterBand(num_band)
            BandWriteArray(bandOut, VET_DATA)
        else:
            num_band = num_band + 1
            #print('------------------')
            #print(num_band)
            bandOut = dsOut.GetRasterBand(num_band)
            BandWriteArray(bandOut, VET_DATA) 
        #print('---- go out to test...')
        delFile(tmpDir + file[:file.rfind('.')] + '_tmp?.tif')
    return True

class main():
    global exitFlag
 
    argDict = mapDict(argv, usage)

    if "-ty" in argDict and "-ds" in argDict and "-de" in argDict:
        typeData = str(argDict["-ty"])
        yearStart = int(argDict["-ds"][0:4])
        yearFinish = int(argDict["-de"][0:4])
        if('MOD' in typeData or 'MYD' in typeData):
            monthStart = int(argDict["-ds"][4:7])
            monthFinish = int(argDict["-de"][4:7])
        else:    
            monthStart = int(argDict["-ds"][4:6])
            monthFinish = int(argDict["-de"][4:6])
    else:
        exit(usage)
    if "-fo" in argDict and "-fi" in argDict:
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
    else:
        inputDir  = 'f:\\DADOS\\'
        outputDir = inputDir

    if (typeData not in ['CLM', '123', 'NDVI', 'MOD08', 'MYD08', 'MOD09', 'MYD09', 'HRV', 'MOD13', 'MOD15']):
        print('Error in typeData: %s'%(typeData))
        exit(usage)

    if(typeData == '123'):
        inputDir = inputDir + 'HRIT\\'
        outputDir = outputDir + '123\\'
    elif(typeData == 'HRV'):
        inputDir = inputDir + 'HRIT\\'
        outputDir = outputDir + 'HRV\\'        
    elif(typeData == 'CLM'):
        inputDir = inputDir + 'GRB\\'
        outputDir = outputDir + 'CLM\\'
    elif(typeData == 'NDVI'):
        inputDir = inputDir + 'HDF\\'
        outputDir = outputDir + 'NDVI\\'
    elif(typeData == 'MOD08'):
        inputDir = inputDir + 'MOD08\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MYD08'):
        inputDir = inputDir + 'MYD08\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MOD09'):
        inputDir = inputDir + 'MOD09\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MYD09'):
        inputDir = inputDir + 'MYD09\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MOD13'):
        inputDir = inputDir       
        outputDir = outputDir + 'MOD13\\' 
    elif(typeData == 'MOD15'):
        inputDir = inputDir + 'MCD15\\'      
        outputDir = outputDir + 'MOD15\\' 
    gc.collect()
    gdal.UseExceptions()
    lst = []
    while(1):
        print('....................... Generate data in folder: ' + inputDir)
        for year in range(yearStart, yearFinish+1):
            print('...........Generate data in year: ' + str(year))
            folder = outputDir + '\\' + str(year)
            if not(os.path.isdir(folder)):
                os.makedirs(folder)
            if('MOD' in typeData or 'MYD' in typeData):
                if(yearStart == yearFinish):
                    monthInic = monthStart
                    monthEnd  = monthFinish                    
                elif(year == yearStart):
                    monthInic = monthStart
                    monthEnd  = 365
                elif(year == yearFinish):
                    monthInic = 1
                    monthEnd  = monthFinish
                else:
                    monthInic = 1
                    monthEnd  = 365

            else:
                if(yearStart == yearFinish):
                    monthInic = monthStart
                    monthEnd  = monthFinish  
                elif(year == yearStart):
                    monthInic = monthStart
                    monthEnd  = 12
                elif(year == yearFinish):
                    monthInic = 1
                    monthEnd  = monthFinish
                else:
                    monthInic = 1
                    monthEnd  = 12                   

            for month in range(monthInic, monthEnd+1):
                str_month = ''
                if('MOD' in typeData or 'MYD' in typeData):
                    str_month = str('%03d'%month) 
                    outputSubdir = outputDir + str(year) + '\\' + str('%02d'%calcMonth(year, month)) + '\\'
                    if(typeData == 'MOD08' or typeData == 'MYD08'):
                        inputSubdir  = inputDir + str(year) + '\\' + str_month + '\\'
                    else:
                        inputSubdir  = inputDir
                else:
                    str_month = str('%02d'%month)
                    outputSubdir = outputDir + str(year) + '\\'  + str_month + '\\'
                    inputSubdir  = inputDir + str(year) + '\\'  + str_month + '\\'
                #print('str_month: ' + str_month)
                #print('typeData: ' + typeData)
                #print('month: ' + str(month))
                #sys.exit()
                print('....... Generate data in month/doy: ' + str_month)

                print('... Generate files in folder: ' + outputSubdir)
                findfile = inputDir+ str(year) + '\\' + str_month + '\\' + '*-' + str(year) + str_month + '*.*'
                if(typeData == 'HRV'):
                    findfile = findfile.replace('*.*', '??19*.*')
                if(typeData == 'CLM'):
                    findfile = findfile.replace('*.*', '*.grb')
                if(typeData == 'NDVI'):
                    findfile = findfile.replace('*.*', '*.h5')
                if(typeData == 'MOD08' or typeData == 'MYD08'):
                    findfile = inputDir + str(year) + '\\' + str_month + '\\' + '*A' + str(year) + str_month + '*.hdf'
                if(typeData == 'MOD09' or typeData == 'MYD09'):
                    findfile = inputDir + str(year) + '\\' + str_month + '\\' + '*A' + str(year) + str_month + '*h13v10*.hdf'
                if(typeData == 'MOD13'):
                    findfile = inputDir +'pa_br_ndvi_250_lapig\\pa_br_ndvi_250_' + str(year) + str_month + '_lapig.tif'
                if(typeData == 'MOD15'):
                    #findfile = inputDir +'pa_br_fpar_1000_lapig\\pa_br_fpar_1000_' + str(year) + str_month + '_lapig.tif'
                    findfile = inputDir + 'MCD15A2H.A' + str(year) + str_month + '*h13v10*.hdf'
                   

                print('... Finding file in folder: ' + findfile)
                files = glob.glob(findfile)
                print(findfile)
                print(files)
                #print('findfile: ' + findfile)
                rows = len(files)
                if(rows == 0):
                    print('... No files to folder: ' + outputSubdir)
                    continue
                if not(os.path.isdir(outputSubdir)):
                    os.makedirs(outputSubdir)
                id = len(inputDir)
                for i in range (0, rows, 1):
                    print('... file files['+str(i)+']: '+ files[i])
                    # Verify if tar file has files to this year and month
                    if not ('MOD' in typeData or 'MYD' in typeData):
                        id = files[i].index('-201')
                        day  = files[i][id+7:id+9]
                        hour = files[i][id+9:id+13]
                        if(typeData == '123' or typeData == 'HRV'):
                            if(hour[2:] == '57'): hour = hour[:2] + '45'
                            if(hour[2:] == '42'): hour = hour[:2] + '30'
                            if(hour[2:] == '27'): hour = hour[:2] + '15'
                            if(hour[2:] == '12'): hour = hour[:2] + '00'

                    print('. Verify file: ' + files[i])
                    
                    if('MOD' in typeData or 'MYD' in typeData):
                        tiffile = str(year) + str('%02d'%calcMonth(year, month)) + str('%02d'%calcDay(year, month))
                        extension = '_' + typeData + '.tif'
                        tiffile = tiffile + extension
                        testfile = outputSubdir + tiffile 
                    else:
                        tiffile = str(year) + str_month + str(day) + str(hour)
                        extension = '_' + typeData + '.tif'
                        testfile = outputSubdir + tiffile + extension
                    if(os.path.isfile(testfile)):
                        print('... File exist in: ' + testfile)
                        if(typeData == 'MOD09'):
                            if(numpy.sum(readFileBand(testfile, 5)) == 0 or numpy.sum(readFileBand(testfile, 4)) == 0):
                                print('... Error in file band, I''m sorry, I delete file!')
                                delFile(testfile)
                                lst.append(inputSubdir)
                                lst.append(outputSubdir)
                                lst.append(tiffile)
                        continue
                    else:
                        print('... File not exist in: ' + testfile)
                        lst.append(inputSubdir)
                        lst.append(outputSubdir)
                        lst.append(tiffile)

        if (lst == []):
            print('Waiting new files')
            sleepTime(1)
            continue
        print( "----- Create new threads")
        threadID = 0
        for tName in threadList:
            thread = myThread(threadID, tName, workQueue, typeData, tmpDir)
            thread.start()
            threads.append(thread)
            threadID += 1

        print( "----- Fill the queue and wait case many files")
        
        x = 0
        print('----- lst')
        for i in range(len(lst)):
            while(x < len(lst)):
                queueLock.acquire()
                if not workQueue.full():
                    word =  lst[x]
                    #print( 'x: ' + str(x) + '\tfilename: ' + word)
                    workQueue.put(word)
                    x = x + 1
                    queueLock.release()
                    
                else:
                    queueLock.release()
                    time.sleep(5)


        print( "----- Wait for queue to empty")
        while not workQueue.empty():
            #print('Files in Queue: ' + str(workQueue.qsize()))
            #time.sleep(5)
            pass
        print( '----- Queue is empty')

        # Notify threads it's time to exit
        exitFlag = 1

    # Wait for all threads to complete
        for t in threads:
            t.join()
        print( "Exiting Main Thread")
import sys
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\SMAC')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\SMAC')

import numpy
import matplotlib.pyplot as plt
import itertools
import kernels
import os
import gc
import glob
import sys
from utils import * 
from sys import argv
from kernels import Kernels # Import ther kernels
import psutil
from Py6S import *
from SMAC_CESBIO import *

numpy.set_printoptions(threshold=np.inf)
usage = """\
Usage: %s [OPTIONS]
        -ps     process (BRDF or 6S)
        -ty     type data (1day, 2days)
        -fi     folder input data
        -fo     folder output data
        -py     pixel y to process
        -px     pixel x to process
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
        -class  data class (culture or pasture)

Points to verify data:
        1. [143, 48] [-51.557, -17.376] [Rio Verde] Agricultura Anual  (Terraclass 1)
        2. [178, 14] [-52.751, -18.624] [Chapadão do Céu]
        3. [155, 32] [-52.105, -17.820] [Jatai]
        4. [144, 38] [-52.084, -17.668] [Jatai]
        5. [146, 33] [-52.068, -17.510] [Perolandia]
        6. [146, 58] [-51.202, -17.504] [Montividiu]
        7. [139, 57] [-51,256, -17.265] [Montividiu]
        8. [144, 67] [-50.916, -17.439] [Rio Verde]
        9. [134, 79] [-50.501, -17.089] [Parauna]
        10. [130, 63] [-51.028, -16.949] [Parauna]
        [84, 122] [-48.990, -15.330] Agricultura Perene (Terraclass 2)
        [16, 172] Natural Campestre  (Terraclass 5)
        [32, 159] Natural Florestal  (Terraclass 6)
        11. [35,  81] [-50.427, -13.634] [Novo Mundo] Pastagem           (Terraclass 11)
        12. [8 ,  83] [-50.361, -12.707] [São Miguel do Araguaia] 
        13. [32, 118] [-49.122, -13.531] [Porangatu]
        14. [34, 123] [-48.969, -13.607] [Santa Tereza de Goias]
        15. [50, 83]  [-50.352, -14.161] [Nova Crixas]
        16. [123, 81] [-50.426, -16.701] [Aurilandia]
        17. [130, 88] [-50.163, -16.923] [Jandaia]
        18. [156, 108][-49.460, -17.830] [Joviania]
        19. [15, 80]  [-50.436, -12.941] [Sao Miguel do Araguaia]
        20. [44, 188] [-46,698, -13,956] [Iaciara]
        [43, 177] Natural Savânica   (Terraclass 12) *** EXCLUIR DAS AMOSTRAS DE PASTAGEM ***
""" % argv[0]

CULTURE = [[143,48],[178,14],[155,32],[144,38],[146,33],[146,58],[139,57],[144,67],[134,79],[130,63]]
PASTURE = [[35,81],[8,83],[32,118],[34,123],[50,83],[123,81],[130,88],[156,108], [15, 80], [44,188]]


def process6Sdata(RED, NIR, ATM, GEOM):

    if(ATM[0] != -1 and ATM[1] != -1 and ATM[2] != -1):
        # process 6S data
        print('WV: %.4f OZONE: %.4f AOT: %.2f'% (ATM[0], ATM[1], ATM[2]))
        print('SZA: %.2f SAZ %.2f VZA: %.2f VAZ: %.2f'%(GEOM[0], GEOM[1], GEOM[2], GEOM[3]))
        OUT_6S = sixsv1(RED, NIR, ATM, GEOM)
        if(OUT_6S[0] != -1 and OUT_6S[1] != -1):
            return OUT_6S
        else:
            print('...Error in OUT_6S data')
            sys.exit(1)
    else:
        print('...Error in ATM read data')
        sys.exit(1)

def sixsv1(RED, NIR, ATM, GEOM):
    
    # Config 6SV1
    [WV, OZONE, AOT] = ATM
    [SZA, SAZ, VZA, VAZ] = GEOM
    try:
        s = SixS(os.getcwd() + "\\sixsV1_1_lab.exe")
    except:
        printlog(True, ' Except in sixv1 path.')
        s.produce_debug_report()

    s.altitudes.set_sensor_satellite_level()
    s.altitudes.set_target_sea_level()
    s.geometry = Geometry.User()
    s.geometry.solar_z = SZA
    s.geometry.solar_a = SAZ
    s.geometry.view_z = VZA
    s.geometry.view_a = VAZ
    s.geometry.day = 1
    s.geometry.month = 1
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(WV, OZONE)
    s.aot550 = AOT


    # Calc to RED band
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(RED) 
    s.wavelength = Wavelength (0.56, 0.71)  

    try:
        s.run()
    except:
        printlog(True, ' Error in Py6S (RED)')
        OUT_6S =  [-1, 0, 0]
        return OUT_6S
    
    RED = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out RED: %.2f'%(RED))

    # Calc to NEAR INFRARED band    
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(NIR) 
    s.wavelength = Wavelength (0.74, 0.88)
    try:
        s.run()
    except:
        printlog(True, ' Error in Py6S (NIR)')
        OUT_6S =  [0, -1, 0]
        return OUT_6S
    #s.produce_debug_report()
    #time.sleep(random.randint(1,10))
    NIR = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out NIR: %.2f'%(NIR))
    NDVI = calcNDVI(NIR, RED)
    OUT_6S = [RED, NIR, NDVI]
    #printlog(True, OUT_6S)
    return OUT_6S

class ATMclass:

    def __init__(self):
        #print('ATMclass init')
        self.date = 0
        self.i = 0
        self.j = 0
        self.inputDir = ''
        self.DATA = [0, 0, 0]

    def read(self, inputDir, date, i, j):
        #print('ATMclass read')
        # Read ATM file
        # band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
        # band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
        # band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
        # band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)
        # ATM to 6S: Vector with atmospheric data [WV O AOT]
        if(self.date == date and self.i == i and self.j == j):
            print('----------- return self.DATA')
            return self.DATA
        
        self.inputDir = inputDir
        self.date = date
        self.i = i
        self.j = j
        inputSubDir = self.inputDir + '___\\' + str(self.date)[0:4] + '\\' + str(self.date)[4:6] + '\\'
        strFile = inputSubDir.replace('___', 'ATM') + str(self.date) + '_ATM.tif'
        DATA = numpy.zeros((3))
        if os.path.isfile(strFile) == 0:
            print('Error to read ATM file: %s' % (strFile))
            DATA = [-1, -1, -1]
        else:
            #print('Read ATM file: %s' % (strFile))
            '''
            # get Aerosol Data Band
            vet = readFileBand(strFile, 4)
            if(not len(vet)):
                return [0, 0, 0]
            else:
                DATA[2] = vet[i][j]
            '''
            # get Water_Vapour  
            vet = readFileBand(strFile, 10)
            if(not len(vet)):
                DATA[0] = -1
            else:
                DATA[0] = vet[i][j] 

            # get Total Ozone Data Band
            vet = readFileBand(strFile, 7)
            if(not len(vet)):
                DATA[1] = -1
            else:
                DATA[1] = vet[i][j]
                DATA[1] = DATA[1]/1000 # Convert Dobson Unit to cm-atm
            # get Deep_Blue Data Band
            vet = readFileBand(strFile, 1)
            if(not len(vet)):
                DATA[2] = -1
            else:
                DATA[2] = vet[i][j]
        self.DATA = DATA
        return self.DATA


def importBRDFdata(inputDir, outputDir, dateStart, dateEnd, i, j):
    startHour = 1300
    endHour   = 1700    
    date = dateStart
    num_days = calcDays(dateStart, dateEnd)
    if(num_days == -1):
        print(' Error in dates, please verify.')
        return []
    data = numpy.zeros((int(num_days*(endHour - startHour)/20) , 11))
    h = 0

    while(date <= dateEnd):

        print(' Reading data files: %s' %(date))
        inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        hour = startHour
        while(hour <= endHour):
            # Flag to read error file
            fileError = False
            #print('Reading hour: %d' %hour)
            data[h][0] = date
            data[h][1] = hour
            
            # --------- Process ANGLE file
            #print('Process ANGLE file')
            strFile = inputSubDir.replace('___', 'ANGLE') + str(date) + str(hour) + '_ANGLES.tif'
            if os.path.isfile(strFile) == 0:
                #print('..Error! Invalid data file: %s' %strFile)
                fileError = True
            else:
                # get VAZ View Azimuth Angle (msg_azres.tif)
                vet = readFileBand(strFile, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][4] = vet[i][j]
                # get SAZ Sun Azimuth Angle (sol_azres.tif)
                vet = readFileBand(strFile, 2)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][6] = vet[i][j]
                # get VZA View Zenith Angle (msg_zenres.tif)
                vet = readFileBand(strFile, 3)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][3] = vet[i][j]
                # get SZA Sun Zenith Angle (sun_zenres.tif)
                vet = readFileBand(strFile, 4)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][5] = vet[i][j]

            # --------- Process 123 file
            if fileError == False:
                #print('Process 123 file')
                strFile = inputSubDir.replace('___', '123') + str(date) + str(hour) + '_123.tif'
                if os.path.isfile(strFile) == 0:
                    #print('..Error! Invalid data file: %s' %strFile)
                    fileError = True
                else:
                    # get RED Data Band
                    vet = readFileBand(strFile, 1)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][7] = vet[i][j]
                    # get NIR Data Band
                    vet = readFileBand(strFile, 2)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][8] = vet[i][j]
                    # get WIR Data Band
                    vet = readFileBand(strFile, 3)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][9] = vet[i][j]

            # --------- Process CLM file
            if fileError == False:
                #print('Process CLM file')
                strFile = inputSubDir.replace('___', 'CLM') + str(date) + str(hour) + '_CLM.tif'
                if os.path.isfile(strFile) == 0:
                    #print('..Error! Invalid data file: %s' %strFile)
                    fileError = True
                else:
                    vet = readFileBand(strFile, 1)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][2] = vet[i][j] 
            # Verify if error clear data, else inc position
            if(fileError == False):
                h = h + 1
            hour = incHour(hour)
        date = incDate(date)
    # data[1:16] = [VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM]
    filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(i) +'_'+ str(j) +  '_BRDF.csv'
    print(' Save output file: %s'%(outputDir + filename))
    numpy.savetxt(outputDir + filename, data, fmt='%1.4f', comments = '', header = 'DATE HOUR CLM VZA VAZ SZA SAZ RED NIR WIR', delimiter = ' ', newline='\n')
    return data

def importMOD13data(inputDir, outputDir, dateStart, dateEnd, i, j):
    
    date = dateStart
    num_days = calcDays(dateStart, dateEnd)
    if(num_days == -1):
        print(' Error in dates, please verify.')
        return []
    data = numpy.zeros((int(num_days/16+1) , 4))
    h = 0
    # MOD13 band_names = ['NDVI', 'RED', 'NIR', 'MIR', 'day_of_the_year', 'pixel_reliability']
    while(date <= dateEnd):
        inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        strFile = inputSubDir.replace('___', 'MOD13') + str(date) + '_MOD13.tif'
        if os.path.isfile(strFile):
            print(' Reading data files: %s' %(date))
            # get DOY data band
            vet = readFileBand(strFile, 5)
            if(not len(vet)):
                return []
            else:
                data[h][0] = vet[i][j]
            # get RED data band
            vet = readFileBand(strFile, 2)
            if(not len(vet)):
                return []
            else:
                data[h][2] = vet[i][j]
            # get NIR data band  
            vet = readFileBand(strFile, 3)
            if(not len(vet)):
                return []
            else:
                data[h][3] = vet[i][j]
            # get NDVI data band
            vet = readFileBand(strFile, 1)
            if(not len(vet)):
                return []
            else:
                if(vet[i][j] != 0):
                    data[h][1] = vet[i][j]
                # if NDVI = 0
                else:
                    if(data[h][3] + data[h][2] != 0):
                        data[h][1] = (data[h][3] - data[h][2])/(data[h][3] + data[h][2])
                    else:
                       data[h][1] = 0 
            h = h + 1
            #print(data)
        #else:
            #print('...Error to read: %s'%(strFile))
        date = incDate(date)
    # data = [DOY, NDVI, RED, NIR]
    filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(i) +'_'+ str(j) +  '_MOD13.csv'
    print(' Save output file: %s'%(outputDir + filename))
    numpy.savetxt(outputDir + filename, data, fmt='%1.4f', comments = '', header = 'DOY NDVI RED NIR', delimiter = ' ', newline='\n')
    return data


class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    
    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict and (("-px" in argDict and "-py" in argDict) or "-class" in argDict):
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData = argDict["-ty"]
        if("-class" in argDict):
            classData = argDict["-class"]
        else:            
            MSG_pix_x = int(argDict["-px"])
            MSG_pix_y  = int(argDict["-py"])
            classData = None
    else:
        exit(usage)
    if(typeData != '1day' and typeData != '2days'):
        print('Error in typeData, verify!')
        exit(usage)

    if(classData == 'culture'):
        CLASS = CULTURE
    elif(classData == 'pasture'):
        CLASS = PASTURE
    elif(classData == None):
        CLASS = [[MSG_pix_x, MSG_pix_y]]
    else:
        exit(usage)

    # Go out if folder not exist
    if(os.path.isdir(inputDir) == False or os.path.isdir(outputDir) == False):
        print('... Error in folder adress!')
        sys.exit()
    # Create class ATM to get atmospheric data
    ATMdata = ATMclass()
    num_days = calcDays(dateStart, dateEnd)
    filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) +  '_class_' + str(classData) + '_typ_' + str(typeData) + '_NDVI_ALL_.csv'
    header = 'DATE DOY HOUR_1330 MAX_DAY MIN_DAY MEAN_DAY MAX_SPARSE_THICK MIN_SPARSE_THICK MEAN_SPARSE_THICK KERN_FAIL_SPARCE_THICK 6S_SPARCE_THICK MAX_GAO MIN_GAO MEAN_GAO KERN_FAIL_GAO 6S_GAO'
    data_ALL = datafile('BRDF', outputDir, filename, header, num_days*10, 16)

    for pos_class in range(len(CLASS)):
        MSG_pix_x = CLASS[pos_class][0]
        MSG_pix_y = CLASS[pos_class][1]
        
        # Create BRDF data base
        filenameBRDF = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) +  '_typ_' + str(typeData) + '_NDVI_.csv'
        data_BRDF = datafile('BRDF', outputDir, filenameBRDF, header, num_days, 16)
 
        # Verify and/or process BRDF file
        filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_BRDF.csv'
        if(os.path.isfile(outputDir + filename)):
            print('... Load data file: %s'%(filename))
            DATA = numpy.loadtxt(outputDir + filename, delimiter = ' ', skiprows = 1)
        else:
            print('... Create data file: %s'%(filename))
            DATA = importBRDFdata(inputDir, outputDir, dateStart, dateEnd, MSG_pix_x, MSG_pix_y)
        
        # Read or create MOD13 data file
        filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_MOD13.csv'
        if(os.path.isfile(outputDir + filename)):
            print('... Load data file: %s'%(filename))
            MOD13 = numpy.loadtxt(outputDir + filename, delimiter = ' ', skiprows = 1)
        else:
            print('... Create data file: %s'%(filename))
            MOD13 = importMOD13data(inputDir, outputDir, dateStart, dateEnd, MSG_pix_x, MSG_pix_y)
            # To calc only MOD13 data
        #continue

        # Verify if files exist
        if(os.path.isfile(outputDir + filename.replace('.csv', 'BRDF.csv'))):
            print('... Data files NDVI exist! Verify!')
            continue

        DATE = DATA[:,0]
        HOUR = DATA[:,1]
        CLM  = DATA[:,2]
        VZA  = DATA[:,3]
        SZA  = DATA[:,5]
        RAA  = DATA[:,4] - DATA[:,6]
        VAZ  = DATA[:,4]
        SAZ  = DATA[:,6]
        BAND = DATA[:,7:10]
        RED  = DATA[:,7]
        NIR  = DATA[:,8]

        date = dateStart
        while(date <= dateEnd):
            doy_today = calcDoy(int(str(date)[0:4]), int(str(date)[4:6]), int(str(date)[6:8]))
            # Select number of days to create NDVI
            date_up   = date    
            if(typeData == '2days'):
                date_up   = incDate(date)
            passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))
            passer2= numpy.logical_and (CLM < 1.05, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(HOUR == 1330, numpy.logical_and(RED > 0, NIR > 0))))
            samples = len(passer[passer == True])
            if(samples < 5):
                date = incDate(date)
                continue
            VZA_filt = VZA[passer]
            SZA_filt = SZA[passer]
            RAA_filt = RAA[passer]
            VAZ_filt = VAZ[passer]
            SAZ_filt = SAZ[passer]
            RED_filt = RED[passer]
            NIR_filt = NIR[passer]
            BAND_filt = BAND[passer]
            print('............ doy: %d' %(doy_today))
            # Mount RESULT matrix
            data_BRDF.addAtrib([date, doy_today])
            data_ALL.addAtrib([date, doy_today])
            ATM = ATMdata.read(inputDir, date, MSG_pix_x, MSG_pix_y);
            # NDVI Sample 1330
            sample_1330 = len(passer2[passer2 == True])
            if(sample_1330 == 1):
                NDVI_1330 = calcNDVI(NIR[passer2], RED[passer2])
                data_BRDF.addAtrib(NDVI_1330)
                data_ALL.addAtrib(NDVI_1330)
                '''
                print('NDVI_1330')
                print(NDVI_1330)
                print(NIR[passer2])
                print(RED[passer2])
                input()
                '''
            else:
                data_BRDF.addAtrib(0)
                data_ALL.addAtrib(0)           
            # NDVI Mean day
            NDVI = numpy.zeros((3, samples))    #[RED NIR NDVI]
            NDVI[0] = RED_filt
            NDVI[1] = NIR_filt
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])

            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])])
            '''
            print('NDVI: %.4f %.4f %.4f'%(numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])))  
            print('... Array NDVI')
            print(NDVI)
            input()
            '''
            print('... Ross-Li Sparce-Thick Normal') 
            geo_kernel = 'Sparse'
            vol_kernel = 'Thick'
            # Generate the kernels, only bother for obs where the QA is OK
            K_obs =  Kernels( VZA_filt, SZA_filt, RAA_filt, \
                LiType=geo_kernel, doIntegrals=False, \
                normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                RossType= vol_kernel )
            kern = numpy.ones (( numpy.sum(passer==True), 3 )) # Store the kernels in an array
            kern[ :, 1 ] = K_obs.Ross
            kern[ :, 2 ] = K_obs.Li
            '''
            print ("%s\t%s\t%s\t%s\t%s" % ( "NUM", "PASSER", "CLM", "DATE", "KERN"))
            for x in range(len(passer[passer == True])):
                if(passer[x] == True):
                    print ("%2d\t%s\t%d\t%d\t%6f\t%6f\t%6f" % (x, passer[x], CLM[x], DATE[x], kern[x][0], kern[x][1], kern[x][2]))
            '''    
            # Calc to RED and NIR
            NDVI = numpy.zeros((3, samples))
            COV  = [0,0]
            kern_fails = 0
            for band in range(2):
                #print('passer: %d  kern: %d'%(len(passer), len(kern)))
                tmp_obs = BAND_filt
                obs_col = numpy.zeros((len(tmp_obs),1))
                obs_lin = numpy.zeros((len(tmp_obs)))
                for x in range(len(tmp_obs)):
                    obs_col[x] = tmp_obs[x][band]
                    obs_lin[x] = tmp_obs[x][band]

                if(len(kern) >= len(passer)):
                    K = kern[passer, :]
                else:
                    K = kern
                (f, rmse, rank, svals ) = numpy.linalg.lstsq( K, obs_col )
                print('----- Results Band[%d]' %(band))
                print ("%-20s %20s" % ( "Kernel", "Value"))
                for i, k in enumerate( ["Isotropic", "Ross-" + vol_kernel, "Li-" + geo_kernel] ):
                    print ("%-20s %20f" % ( k, f[i] ))
                if(f[0] < 0):
                    kern_fails = (band + 1) + kern_fails
                fwd = K.dot(f)
                fwd_lin = numpy.zeros(len(fwd))
                for x in range(len(fwd)):
                    fwd_lin[x] = fwd[x]
                COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
                NDVI[band] = fwd_lin
            print('kern_fails: %d'%(kern_fails))
            # Calc NDVI = (NIR - RED)/(NIR + RED)
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_BRDF.addAtrib(kern_fails)
            data_ALL.addAtrib(kern_fails)
            error = rmse
            '''
            print('--- NDVI BRDF')
            print(NDVI)
            '''
            print(' error: \t%f'%(error))
            # Process 6S correction in BRDF data
            id = numpy.argmax(NDVI[2]) 
            RED_BRDF  = NDVI[0][id]
            NIR_BRDF  = NDVI[1][id]
            NDVI_BRDF = NDVI[2][id]

            GEOM = [SZA_filt[id], SAZ_filt[id], VZA_filt[id], VAZ_filt[id]]
            [RED_6S, NIR_6S, NDVI_6S] = [0, 0, 0] #process6Sdata(RED_BRDF, NIR_BRDF, ATM, GEOM)
            data_BRDF.addAtrib(NDVI_6S)
            data_ALL.addAtrib(NDVI_6S)

            print('... Ross-Li Sparce-Thick Filter Data') 
            print('all samples: %d'%samples)
            
            # passer to +/- 20% NDVI variation
            NDVI = numpy.zeros((3, samples))    #[RED NIR NDVI]
            NDVI[0] = RED_filt
            NDVI[1] = NIR_filt
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
            NDVI = NDVI[2]
            NDVI_temp = np.partition(-NDVI, 4)
            NDVI_max = -NDVI_temp[:4]
            NDVI_mean = np.mean(NDVI_max)
            print('NDVI:')
            print(NDVI)
            print('NDVI_max:')
            print(NDVI_max)
            print('NDVI_mean: %.4f' %NDVI_mean) 
            
            NDVI = calcNDVI(NIR, RED)
            #print('NDVI_all')
            #print(NDVI)
            passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(NDVI > 0.2*NDVI_mean, numpy.logical_and(NDVI < 1.8*NDVI_mean, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))))

            SAZ_mean = numpy.mean(SAZ_filt)
            #passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(SAZ > 0.4*SAZ_mean, numpy.logical_and(SAZ < 1.6*SAZ_mean, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))))
            passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(NDVI > 0.8*NDVI_mean, numpy.logical_and(NDVI < 1.2*NDVI_mean, numpy.logical_and(SAZ > 0.2*SAZ_mean, numpy.logical_and(SAZ < 1.8*SAZ_mean,numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))))))

            samples = len(passer[passer == True])
            if(samples < 5):
                print('samples: %d'%(samples))
                date = incDate(date)
                data_BRDF.addAtrib([0,0,0,-1,0]) 
                data_ALL.addAtrib([0,0,0,-1,0]) 
                data_BRDF.addSample() 
                data_ALL.addSample()
                #input()
                continue
            VZA_filt = VZA[passer]
            SZA_filt = SZA[passer]
            RAA_filt = RAA[passer]
            VAZ_filt = VAZ[passer]
            SAZ_filt = SAZ[passer]
            RED_filt = RED[passer]
            NIR_filt = NIR[passer]
            BAND_filt = BAND[passer]

            '''
            geo_kernel = 'Sparse'
            vol_kernel = 'Thick'
            # Generate the kernels, only bother for obs where the QA is OK
            K_obs =  Kernels( VZA_filt, SZA_filt, RAA_filt, \
                LiType=geo_kernel, doIntegrals=False, \
                normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                RossType= vol_kernel )
            kern = numpy.ones (( numpy.sum(passer==True), 3 )) # Store the kernels in an array
            kern[ :, 1 ] = K_obs.Ross
            kern[ :, 2 ] = K_obs.Li

            print ("%s\t%s\t%s\t%s\t%s" % ( "NUM", "PASSER", "CLM", "DATE", "KERN"))
            for x in range(len(passer[passer == True])):
                if(passer[x] == True):
                    print ("%2d\t%s\t%d\t%d\t%6f\t%6f\t%6f" % (x, passer[x], CLM[x], DATE[x], kern[x][0], kern[x][1], kern[x][2]))
   
            # Calc to RED and NIR
            NDVI = numpy.zeros((3, samples))
            COV  = [0,0]
            kern_fails = 0
            for band in range(2):
                #print('passer: %d  kern: %d'%(len(passer), len(kern)))
                tmp_obs = BAND_filt
                obs_col = numpy.zeros((len(tmp_obs),1))
                obs_lin = numpy.zeros((len(tmp_obs)))
                for x in range(len(tmp_obs)):
                    obs_col[x] = tmp_obs[x][band]
                    obs_lin[x] = tmp_obs[x][band]

                if(len(kern) >= len(passer)):
                    K = kern[passer, :]
                else:
                    K = kern
                (f, rmse, rank, svals ) = numpy.linalg.lstsq( K, obs_col )
                print('----- Results Band[%d]' %(band))
                print ("%-20s %20s" % ( "Kernel", "Value"))
                for i, k in enumerate( ["Isotropic", "Ross-" + vol_kernel, "Li-" + geo_kernel] ):
                    print ("%-20s %20f" % ( k, f[i] ))
                if(f[0] < 0):
                    kern_fails = (band + 1) + kern_fails
                fwd = K.dot(f)
                fwd_lin = numpy.zeros(len(fwd))
                for x in range(len(fwd)):
                    fwd_lin[x] = fwd[x]
                COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
                NDVI[band] = fwd_lin
            print('kern_fails: %d'%(kern_fails))
            # Calc NDVI = (NIR - RED)/(NIR + RED)
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_BRDF.addAtrib(kern_fails)
            data_ALL.addAtrib(kern_fails)
            error = rmse

            print('--- NDVI BRDF')
            print(NDVI)

            print(' error: \t%f'%(error))
            # Process 6S correction in BRDF data
            id = numpy.argmax(NDVI[2]) 
            RED_BRDF  = NDVI[0][id]
            NIR_BRDF  = NDVI[1][id]
            NDVI_BRDF = NDVI[2][id]

            '''
            # Calc NDVI using Gao et al. (2002)
            #print('---------------- Calc NDVI using Gao et al. (2002)')
            print('... Gao (2002)')
            NDVI_calc = numpy.zeros((3, samples))
            NDVI_calc[0] = RED[passer]
            NDVI_calc[1] = NIR[passer]
            NDVI_calc[2] = (NDVI_calc[1] - NDVI_calc[0])/(NDVI_calc[1] + NDVI_calc[0])
            NDVI_mean = numpy.mean(NDVI_calc[2])
            W_ini   = (NDVI_calc[2]/NDVI_mean)**2 
            geo_kernel = 'Sparce'
            vol_kernel = 'Thick'
            error = 1
            loop = 0
            while(error > 0.0001 and loop <= 5):

                K_obs =  Kernels( VZA_filt, SZA_filt, RAA_filt, \
                    LiType=geo_kernel, doIntegrals=False, \
                    normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                    RossType= vol_kernel )
                # Store the kernels in an array
                kern = numpy.ones (( numpy.sum(passer==True), 3 )) 
                kern[ :, 1 ] = K_obs.Ross
                kern[ :, 2 ] = K_obs.Li
                # Calc to RED and NIR
                NDVI = numpy.zeros((3, samples))
                COV  = [0,0]
                kern_fails = 0
                for band in range(2):
                    #print('passer: %d  kern: %d'%(len(passer), len(kern)))
                    tmp_obs = BAND_filt
                    obs_col = numpy.zeros((len(tmp_obs),1))
                    obs_lin = numpy.zeros((len(tmp_obs)))
                    for x in range(len(tmp_obs)):
                        obs_col[x] = tmp_obs[x][band]
                        obs_lin[x] = tmp_obs[x][band]

                    if(len(kern) >= len(passer)):
                        K = kern[passer, :]
                    else:
                        K = kern
                    # Using Weight Inversion
                    K = K *numpy.sqrt(W_ini)[:,None]
                    #print('--- obs_col')
                    #print(obs_col)
                    #print('--- W_ini')
                    #print(W_ini)
                    #print('--- K')
                    #print(K)
                    obs_col = obs_col*numpy.sqrt(W_ini)[:,None]
                    #print('--- obs_col')
                    #print(obs_col)
                    (f, rmse, rank, svals ) = numpy.linalg.lstsq( K, obs_col )
                    print('----- Results Band[%d]' %(band))
                    print ("%-20s %20s" % ( "Kernel", "Value"))
                    for i, k in enumerate( ["Isotropic", "Ross-" + vol_kernel, "Li-" + geo_kernel] ):
                        print ("%-20s %20f" % ( k, f[i] ))
                    fwd = K.dot(f)
                    fwd_lin = numpy.zeros(len(fwd))
                    for x in range(len(fwd)):
                        fwd_lin[x] = fwd[x]
                    COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
                    NDVI[band] = fwd_lin
                    if(f[0] < 0):
                        kern_fails = (1 + band) + kern_fails
                # Calc NDVI = (NIR - RED)/(NIR + RED)
                NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
                W_obs = (NDVI[2]/NDVI_calc[2])**2
                error = 0
                #print('---- W_ini')
                #print(W_ini)
                #print('---- W_obs')
                #print(W_obs)
                for i in range(len(W_ini)):
                    error = error + (W_obs[i] - W_ini[i])**2
                W_ini = W_obs
                print(' error: %f in loop %d'%(error, loop))
                loop = loop + 1
            print('kern_fails: %d'%(kern_fails))
            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])])
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])])
            data_BRDF.addAtrib(kern_fails)
            data_ALL.addAtrib(kern_fails)
            # Process 6S correction in BRDF data
            id = numpy.argmax(NDVI[2]) 
            RED_BRDF  = NDVI[0][id]
            NIR_BRDF  = NDVI[1][id]
            NDVI_BRDF = NDVI[2][id]
            GEOM = [SZA_filt[id], SAZ_filt[id], VZA_filt[id], VAZ_filt[id]]
            [RED_6S, NIR_6S, NDVI_6S] = [0, 0, 0] #process6Sdata(RED_BRDF, NIR_BRDF, ATM, GEOM)
            data_BRDF.addAtrib(NDVI_6S)
            data_ALL.addAtrib(NDVI_6S)

            data_BRDF.addSample() 
            data_ALL.addSample()
            #input()   
            date = incDate(date)
            if(type == '2days'):
                date = incDate(date)
        data_BRDF.save()
    data_ALL.save()
import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 't:\\carlos\\newpython\\utils')

import numpy
import math
import os
from time import  strftime, localtime, sleep
import sys
import gc
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from time import strftime, localtime
from sys import argv
from utils import mapDict, incHour, incDate, infoFile, printlog,  readFileBand, readMaskGoias
import time
import threading
import queue
from Py6S import *
import psutil
import random


#def __debug__

usage = """\
Usage: %s [OPTIONS]
        -ty     type of process BRDF (to create...)
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
""" % argv[0]

# Constantes
# To print lib version
# lib.__version__
fileGoiasGeo = os.getcwd() + '\\shape\\201408221300.tif'

# Variables of threads
NUM_THREADS = 4
queueLock = threading.Lock()
workQueue = queue.Queue(4*NUM_THREADS)
threads = []
exitFlag = 0

# Process time data
startHour = 1300
endHour   = 1700

# Constants of BRDF Process
br = 1       # b/r used in K_geo, to calc SZA_line and VZA
hb = 1.5       # h/b used in K_geo, to calc cos(t)       


def sixsv1(i,j,VET_BRDF, VET_ATM):
    # VET: Matrix with data [REF VZA VAZ SZA SAZ] in 3 bands data
    # ATM: Vector with atmospheric data [WV O AOT]
    
    # Config 6SV1
    try:
        s = SixS(os.getcwd() + "\\sixsV1_1_dell.exe")
        #s.test()
    except:
        s = SixS(os.getcwd() + "\\sixsV1_1_lab.exe")
        printlog(True, '----- Except in sixv1 path')
        #s.test()
    #s.produce_debug_report()

    # Set the atmosphere profile to be base
    # Set the wavelength to be that of
    # Data file Py6S_artigo_ex3.py
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(VET_ATM[0], VET_ATM[1])
    s.aot550 = VET_ATM[2]

    # Calc to RED band
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(VET_BRDF[0,0]) 
    s.wavelength = Wavelength (0.56, 0.71)  
    s.geometry = Geometry.User()
    s.geometry.view_z  = VET_BRDF[0,1]
    s.geometry.view_a  = VET_BRDF[0,2]
    s.geometry.solar_z = VET_BRDF[0,3]
    s.geometry.solar_a = VET_BRDF[0,4]
    try:
        s.run()
    except:
        printlog(True, '----- Error in Py6S (RED)')
        printlog(True, 'Using Memory: %.2f CPU: %.2f' %(psutil.virtual_memory().percent, psutil.cpu_percent()))
        DATA_BRDF =  [ 0, -1, 0, 0]
        return DATA_BRDF
    
    RED = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out RED: %.2f'%(RED))
    # Calc to NEAR INFRARED band    
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(VET_BRDF[1,0]) 
    s.wavelength = Wavelength (0.74, 0.88)
    s.geometry = Geometry.User()
    s.geometry.view_z  = VET_BRDF[1,1]
    s.geometry.view_a  = VET_BRDF[1,2]
    s.geometry.solar_z = VET_BRDF[1,3]
    s.geometry.solar_a = VET_BRDF[1,4]
    try:
        s.run()
    except:
        printlog(True, '----- Error in Py6S (NIR)')
        printlog(True, 'Using Memory: %.2f CPU: %.2f' %(psutil.virtual_memory().percent, psutil.cpu_percent()))
        DATA_BRDF =  [ 0, 0, -1, 0]
        return DATA_BRDF
    #s.produce_debug_report()
    #time.sleep(random.randint(1,10))
    NIR = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out NIR: %.2f'%(NIR))
    # Calc to SHORT INFRARED band    
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(VET_BRDF[2,0]) 
    s.wavelength = Wavelength (1.50, 1.78)
    s.geometry = Geometry.User()
    s.geometry.view_z  = VET_BRDF[2,1]
    s.geometry.view_a  = VET_BRDF[2,2]
    s.geometry.solar_z = VET_BRDF[2,3]
    s.geometry.solar_a = VET_BRDF[2,4]
    try:
        s.run()
    except:
        printlog(True, '----- Error in Py6S (WIR)')
        printlog(True, 'Using Memory: %.2f CPU: %.2f' %(psutil.virtual_memory().percent, psutil.cpu_percent()))
        DATA_BRDF =  [ 0, 0, 0, -1]
        return DATA_BRDF
    #s.produce_debug_report()
    #time.sleep(random.randint(1,10))
    WIR = float(s.outputs.pixel_reflectance)
    if((NIR + RED) < 0.0001 and (NIR + RED) > -0.0001):
        printlog(True,'... Zero division warning')
        NDVI = 9999
        printlog(True, [i, j, NDVI, RED, NIR, WIR])
    else:
        NDVI = float((NIR - RED)/(NIR + RED))

    #printlog(True, '----- Out WIR: %.2f'%(WIR))
    '''
    printlog(False, 'RED: %.4f' %(RED))
    printlog(False, 'NIR: %.4f' %(NIR))
    printlog(False, 'NIR: %.4f' %(WIR))
    printlog(False, 'NDVI: %.4f' %(NDVI))
    '''
    DATA_BRDF = [NDVI, RED, NIR, WIR]
    #printlog(True, DATA_BRDF)
    return DATA_BRDF

class myThread (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        printlog (False, "Starting " + self.name)
        process_data(self.name, self.q)
        printlog (False, "Exiting " + self.name)

def process_data(threadName, q):
    global queueLock
    global workQueue
    global exitFlag
    global DATA_6S
    #printlog(True, "Enter in process_data. exitFlag: " + str(exitFlag))
    if __debug__: print('#')
    
    while not exitFlag:
        #printlog(True, "exitFlag is False")
        queueLock.acquire()
        #printlog(True, "\nqueueLock acquire")
        if not workQueue.empty():
            #printlog(True, "Queue is not empty")
            i = q.get()
            j = q.get()
            VET_BRDF = q.get()
            VET_ATM  = q.get()
            queueLock.release()
            DATA_6S[0,i,j] = -1
            if __debug__: print('l')
            while(DATA_6S[0,i,j] == -1):
                # Return NDVI, RED, NIR, WIR
                [DATA_6S[0,i,j], DATA_6S[1,i,j], DATA_6S[2,i,j], DATA_6S[3,i,j]]  = sixsv1(i, j, VET_BRDF, VET_ATM) 
                if __debug__: print('!')
        else:
            queueLock.release()
        time.sleep(1)
    
def calc_K_vol(VZA, SZA, RAA):
    # Pre-process variable
    EPS = math.acos(math.cos(SZA)*math.cos(VZA) + math.sin(SZA)*math.sin(VZA)*math.cos(RAA))
    # K_vol calc
    #printlog(False,'calc_k_vol:')
    #printlog(False,'SZA: %.2f\tVZA: %.2f\tRAA: %.2f\tEPS: %.2f' %(SZA, VZA, RAA, EPS))
    K_tmp = (math.pi/2 - EPS)*math.cos(EPS) + math.sin(EPS)
    K_vol  = K_tmp/(math.cos(SZA) + math.cos(VZA)) - math.pi/4
    #printlog(False,'K_tmp: %.2f\tK_vol: %.2f' %(K_tmp, K_vol))    
    #print '... BRDF >> K_vol = ' + str(K_vol)
    return K_vol
                         
def calc_K_geo(VZA, SZA, RAA):
    # Pre-process variable
    SZA = math.atan(br*math.tan(SZA))
    VZA = math.atan(br*math.tan(VZA))
    EPS = math.acos(math.cos(SZA)*math.cos(VZA) + math.sin(SZA)*math.sin(VZA)*math.cos(RAA))
    #printlog(False,'calc_k_geo:')
    #printlog(False,'SZA: %.2f\tVZA: %.2f\tRAA: %.2f\tEPS: %.2f' %(SZA, VZA, RAA, EPS))
    D = math.sqrt(math.pow(math.tan(SZA), 2) + math.pow(math.tan(VZA), 2) - 2*math.tan(SZA)*math.tan(VZA)*math.cos(RAA))
    t_tmp = hb*math.sqrt(math.pow(D, 2) + math.pow((math.tan(SZA)*math.tan(VZA)*math.sin(RAA)), 2))
    cost = t_tmp/(1/math.cos(SZA) + 1/math.cos(VZA))
    if(cost > 1 or cost < -1):
        #printlog(False, '---- Warning! cost not in range [-1, 1]. SZA: %.2f VZA: %.2f D: %.2f t_tmp: %.2f cost: %.2f' %(SZA, VZA, D, t_tmp, cost))
        if(cost > 1):
            t = math.acos(1)
            cost = 1
        else:
            t = math.acos(-1)
            cost = -1
    else:    
        t = math.acos(cost)
    O = (1/math.pi)*(t - math.sin(t)*cost)*(1/math.cos(SZA) + 1/math.cos(VZA))
    #printlog(False,'D:  %.2f\tO:  %.2f\tt_tmp: %.2f\tt: %.2f' %(D, O, t_tmp, t))

    # K_geo calc                     
    K_geo = O - 1/math.cos(SZA) - 1/math.cos(VZA) + 1/2*(1 + math.cos(EPS))*(1/math.cos(SZA)*1/math.cos(VZA)) 
    #print '... BRDF >> K_geo = ' + str(K_geo)
    return K_geo
                                    
def calc_BRDF(vet):
    # vet = [VZA VAZ SZA SAZ RED NIR WIR NDVI]
    row = len(vet)
    if (row < 5):
        printlog(True,'Error! Little samples to calc BRDF')
        return 0
    K = numpy.zeros((row,3))
    R = numpy.zeros((4,row))
    BRDF = numpy.zeros((3,row,3))
    reflect = numpy.zeros((3,row))
    COEFS = numpy.zeros((3,3))
    VZA = numpy.zeros((row))
    SZA = numpy.zeros((row))
    RAA = numpy.zeros((row))
    for i in range (0, row, 1):
        #printlog(False,'------- LOOP %d' %(i))
        VZA[i] = math.radians(vet[i][0])
        SZA[i] = math.radians(vet[i][2])
        RAA[i] = math.radians(vet[i][1] - vet[i][3])
        #printlog(True,'VZA: %.2fº (%.2f)\tSZA: %.2fº (%.2f)\tRAA: %.2fº (%.2f)' %(vet[i][0], VZA[i], vet[i][2], SZA[i], vet[i][1] - vet[i][3], RAA[i]))

        R[0][i] = vet[i][4]
        R[1][i] = vet[i][5]
        R[2][i] = vet[i][6]
        R[3][i] = vet[i][7]
        K[i][0] = 1
        K[i][1] = calc_K_vol(VZA[i], SZA[i], RAA[i])
        K[i][2] = calc_K_geo(VZA[i], SZA[i], RAA[i])
    #printlog(True, 'VZA_var: %.2f\t\tSZA_var: %.2f\t\tRAA_var: %.2f' %(numpy.var(VZA), numpy.var(SZA), numpy.var(RAA)))

    W_inic = numpy.zeros((row))
    neg_x = 0
    W_SZA = 1
    W_RAA = 1

    for i in range(0, row, 1):
        W_inic[i] = 1.6 - W_SZA*math.sin(abs(SZA[i])) - W_RAA*math.sin(abs(RAA[i]))
    #printlog(True, 'W_inic:')
    #printlog(True, W_inic)

    for loop in range (0, 4, 1):
        for i in range(0, row, 1):
            K[i][0] = W_inic[i]
        '''
        printlog(False,'----- K [K_iso K_vol K_geo]:')
        printlog(False, K)
        printlog(False,'----- R [R_vis R_nir R_wir]:')
        printlog(False, R)
        printlog(False,'----- lstsq:')
        printlog(False,'x[n]:')
        '''
        for i in range(0, 3, 1):
            x,resid,rank,s = numpy.linalg.lstsq(K, R[i])
            COEFS[i] = x
            '''    
            printlog(False,'----- calc to K and R[%d]:'%(i))   
            printlog(False, x)
            printlog(False,'resid:')
            printlog(False, resid)
            printlog(False,'rank:')
            printlog(False, rank)
            printlog(False,'s:')
            printlog(False, s)
            '''
            #printlog(False,'----- BRDF in R[%d]:' %(i))
            
            if(x[0] < 0):
                neg_x = neg_x + 1
            for j in range(0, row, 1):
                BRDF[i,j] = K[j]*x
                reflect[i,j] = numpy.sum(BRDF[i,j])
                #printlog(False,'BRDF: [%.4f %.4f %.4f]\treflect: [%.4f]' %(BRDF[i,j,0], BRDF[i,j,1], BRDF[i,j,2], reflect[i,j]))
            #printlog(False,'BRDF MEAN: [%.4f %.4f %.4f]\treflect: [%.4f]' %(numpy.mean(BRDF[i,:,0]), numpy.mean(BRDF[i,:,1]), numpy.mean(BRDF[i,:,2]), numpy.mean(reflect[i,:])))
        if(neg_x != 0):
            if(numpy.var(SZA) > numpy.var(RAA)):
                W_SZA = W_SZA - 0.2
                W_RAA = W_RAA + 0.2
                
            else:
                W_SZA = W_SZA + 0.2
                W_RAA = W_RAA - 0.2

            if(W_SZA < 0): W_SZA = 0
            if(W_RAA < 0): W_RAA = 0

            for i in range(0, row, 1):
                W_inic[i] = 1.6 - W_SZA*math.sin(abs(SZA[i])) - W_RAA*math.sin(abs(RAA[i]))
            '''
            printlog(False, 'W_SZA: %.2f\t\t W_RAA: %.2f neg_x: %d' %(W_SZA, W_RAA, neg_x)) 
            printlog(False, 'W_inic:')
            printlog(False, W_inic)
            '''
            neg_x = 0
        else:
            break
    # Calc NDVI to return result
    '''
    # Calc NDVI using 'num1' values with 10 percent variation of NDVI 
    NDVI_calc = numpy.zeros((row))
    NDVI_mean  = 0
    NDVI_mean2 = 0
    for i in range (0, row, 1):
        NDVI_calc[i] = (reflect[1,i] - reflect[0,i])/(reflect[1,i] + reflect[0,i])

    printlog(False, 'NDVI_calc:')
    printlog(False, NDVI_calc)
    printlog(False, 'NDVI_mean : %.4f' %(NDVI_mean))
    '''
    '''
    num1 = 0
    for i in range (0, row, 1):
        if NDVI_calc[i] < 1.1*NDVI_mean and NDVI_calc[i] > 0.9*NDVI_mean:
            NDVI_mean2 = NDVI_mean2 + NDVI_calc[i]
            num1 = num1 + 1
    NDVI_mean2 = NDVI_mean2/num1
    printlog(False, 'NDVI_mean2: %.4f (samples: %d)' %(NDVI_mean2, num1))
    '''
    # mat = [VZA VAZ SZA SAZ RED NIR WIR NDVI]
    reflect_mean = numpy.zeros(3)
    num = numpy.zeros(3)
    # vet = [VZA VAZ SZA SAZ RED NIR WIR NDVI]
    # Matrix to return data to 6S algorithm [REF VZA VAZ SZA SAZ] in 3 bands data
    # printlog(False, '----- MAT_BRDF:')
    MAT_BRDF = numpy.zeros((3,5))
    for i in range (0, 3, 1):
        reflect_mean[i] = numpy.mean(reflect[i,:])
        for j in range (0, row, 1):
            if(reflect[i,j] < 1.15*reflect_mean[i] and reflect[i,j] > 0.85*reflect_mean[i]):
                MAT_BRDF[i,0] = MAT_BRDF[i,0] + reflect[i,j]
                MAT_BRDF[i,1] = MAT_BRDF[i,1] + vet[i][0] # VZA
                MAT_BRDF[i,2] = MAT_BRDF[i,2] + vet[i][1] # VAZ
                MAT_BRDF[i,3] = MAT_BRDF[i,3] + vet[i][2] # SZA
                MAT_BRDF[i,4] = MAT_BRDF[i,4] + vet[i][3] # SAZ
                num[i] = num[i] + 1

        # Try against case data problem
        if(num[i] == 0):
            '''
            printlog(True, '. Problem (%d, %d):'%(i,j))
            printlog(True, 'reflect[i]:')
            printlog(True,  reflect[i])
            printlog(True, 'reflect_mean[i]:')
            printlog(True,  reflect_mean[i])
            printlog(True, 'vet[i]:')
            printlog(True,  vet[i])
            '''
            MAT_BRDF[i,:] = [0, 0, -1, -1, -1]
        else:
            MAT_BRDF[i,:] = MAT_BRDF[i,:]/num[i]
    
    RED  = MAT_BRDF[0,0]
    NIR  = MAT_BRDF[1,0]
    WIR  = MAT_BRDF[2,0]

    if(NIR + RED < 0.0001 and NIR + RED > -0.0001):
        NDVI = 9999 
        printlog(True, '...RuntimeWarning. RED: %.2f NIR: %.2f '%(RED, NIR))
    else:
        NDVI = (NIR - RED)/(NIR + RED)

    '''
    printlog(False, 'reflect_mean:')
    printlog(False, reflect_mean)
    printlog(False, 'num:')
    printlog(False, num)
    printlog(False, 'NDVI_mean_return: %.4f' %(NDVI_return))
    printlog(False, 'MAT_BRDF/num: ')
    printlog(False, MAT_BRDF)
    '''    
    
    '''
    printlog(True, 'COEFS:')
    printlog(True, COEFS)
    # wait to test routine
    wait()
    # return NDVI median of all values
    '''
    return NDVI, RED, NIR, WIR, COEFS, MAT_BRDF

def createMaskGoias(filename, band):
    # command to create maskGoias file: filename and band to read
    #createMaskGoias(inputDir + '123\\' + str(day) + str(hour) + '.tif', 4)
    if os.path.isfile(filename) == 0:
        printlog(True,'..Erro! Arquivo Invalido: %s' %strFile)
        return
    infoFile(filename)
    vetMask = readFileBand(filename, 4)
    vetMask = vetMask/255
    numpy.savetxt(os.getcwd() + '\\shape\\maskGoias.txt', vetMask, fmt='%1.0d')
    
def geraMascara(inputFile, outputFile):
    printlog(True,'...geraMascara [x][y]')
    ds1 = gdal.Open(inputFile, GA_ReadOnly )
    printlog(True,'inputFile: ' + inputFile)
    band1 = ds1.GetRasterBand(1)
    array = BandReadAsArray(band1)
    rows = len(array)
    cols = len(array[0])
    #printlog(False,'rows: ' + str(rows))
    #printlog(False,'cols: ' + str(cols))
    MASK1 = numpy.zeros((rows, cols))
    MASK2 = numpy.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            MASK1[i,j] = i
            MASK2[i,j] = j
    # Armazena em arquivo o resultado de selecao dos melhores pixels do dia
    driver = gdal.GetDriverByName("GTiff")
    dsOutMASK = driver.Create(outputFile, ds1.RasterXSize, ds1.RasterYSize, 2, band1.DataType)
    CopyDatasetInfo(ds1,dsOutMASK)
    bandOut1=dsOutMASK.GetRasterBand(1)
    bandOut2=dsOutMASK.GetRasterBand(2)
    BandWriteArray(bandOut1, MASK1)
    BandWriteArray(bandOut2, MASK2)
    #printlog(False,'outputFile: ' + outputFile)
      
def saveCSVfile(matrix, filename):
    timeFile = strftime("%Y%m%d%H%M%S", localtime()) 
    saveFile = os.getcwd()+ '\\logs\\' + filename + '_' + timeFile + '.csv'
    numpy.savetxt(saveFile, matrix, fmt='%.4f',delimiter=";")
    f = open(saveFile,'r')
    filedata = f.read()
    newdata = filedata.replace(".",",")
    f.close()
    f = open(saveFile,'w')
    f.write(newdata)
    f.close()
    
def generateBRDF(inputDir, outputDir, day):
    hour = startHour
    row = 204
    col = 211
    # data = [HOUR][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][X, Y]
    data = numpy.zeros((20, 8, 204, 211), dtype=numpy.float)
    DATA_BRDF = numpy.zeros((4, row, col), dtype=numpy.float)
    global DATA_6S 
    DATA_6S = numpy.zeros((4, row, col), dtype=numpy.float)
    COEFS = numpy.zeros((row, col, 3, 3), dtype=numpy)
    global workQueue
    global exitFlag

    lst = []
    time_brdf = 0
    cont_time_brdf = 0
    time_6S = 0
    cont_time_6S = 0
    start_time = 0
    #printlog(True,data.shape)
    # Read maskGoias
    vetMask = readMaskGoias()
    #printlog(True, 'generate BRDF: ' + str(day))
    #printlog(True,'vetMask: [%d] [%d] ' %(len(vetMask), len(vetMask[1])))
    h = 0 # mark number file in time 15 min
    errors = 0 # Get number errors to read data files
    printlog(True,'----- Reading data files')
    while(hour <= endHour):
        # Flag to read error file
        fileError = False
        printlog(False,'Reading hour: %d' %hour)
        
        # --------- Process ANGLE file
        printlog(False,'Process ANGLE file')
        strFile = inputDir.replace('___', 'ANGLE') + str(day) + str(hour) + '_ANGLES.tif'
        if os.path.isfile(strFile) == 0:
            printlog(False,'..Error! Invalid data file: %s' %strFile)
            fileError = True
        else:
            # get VAZ View Azimuth Angle (msg_azres.tif)
            vet = readFileBand(strFile, 1)
            if(not len(vet)):
                fileError = True
            else:
                data[h][1][:][:] = vet[:][:]
            # get SAZ Sun Azimuth Angle (sol_azres.tif)
            vet = readFileBand(strFile, 2)
            if(not len(vet)):
                fileError = True
            else:
                data[h][3][:][:] = vet[:][:]
            # get VZA View Zenith Angle (msg_zenres.tif)
            vet = readFileBand(strFile, 3)
            if(not len(vet)):
                fileError = True
            else:
                data[h][0][:][:] = vet[:][:]
            # get SZA Sun Zenith Angle (sun_zenres.tif)
            vet = readFileBand(strFile, 4)
            if(not len(vet)):
                fileError = True
            else:
                data[h][2][:][:] = vet[:][:]

        # --------- Process 123 file
        if fileError == False:
            printlog(False,'Process 123 file')
            strFile = inputDir.replace('___', '123') + str(day) + str(hour) + '_123.tif'
            if os.path.isfile(strFile) == 0:
                printlog(False,'..Error! Invalid data file: %s' %strFile)
                fileError = True
            else:
                # get RED Data Band
                vet = readFileBand(strFile, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][4][:][:] = vet[:][:]
                # get NIR Data Band
                vet = readFileBand(strFile, 2)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][5][:][:] = vet[:][:]
                # get WIR Data Band
                vet = readFileBand(strFile, 3)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][6][:][:] = vet[:][:]

        # --------- Process CLM file
        if fileError == False:
            printlog(False,'Process CLM file')
            strFile = inputDir.replace('___', 'CLM') + str(day) + str(hour) + '_CLM.tif'
            if os.path.isfile(strFile) == 0:
                printlog(False,'..Error! Invalid data file: %s' %strFile)
                fileError = True
            else:
                vet = readFileBand(strFile, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][7][:][:] = vet[:][:] 
        # For next loop    
        h = h + 1
        hour = incHour(hour)
        if(fileError == True):
            errors = errors + 1
    if(errors > 10):
        printlog(False,'..Error! Invalid data file: %s' %strFile)
        return False
    else:
        printlog(True, '----- Now generate BRDF %s'%(day))
    # Read ATM file
    # band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
    # band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
    # band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
    # band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)
    printlog(False,'Process ATM file')
    VET_ATM = numpy.zeros((4, row, col), dtype=numpy)
    strFile = inputDir.replace('___', 'ATM') + str(day) + '_ATM.tif'
    if os.path.isfile(strFile) == 0:
        printlog(False,'..Error! Invalid data file: %s' %strFile)
        return False
    else:
        # get Deep_Blue Data Band
        vet = readFileBand(strFile, 1)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[0][:][:] = vet[:][:]
        # get Aerosol Data Band
        vet = readFileBand(strFile, 2)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[1][:][:] = vet[:][:]
        # get Total Ozone Data Band
        vet = readFileBand(strFile, 3)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[2][:][:] = vet[:][:]
        # get Water_Vapour  
        vet = readFileBand(strFile, 3)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[3][:][:] = vet[:][:]  
    # data = [HOUR][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][X, Y]
    # mat = [VZA VAZ SZA SAZ RED NIR WIR]
    percent = 0
    start_time = time.time()
    for i in range(0, row, 1):
        
        if (float(100*i)/float(row) - 5) > percent :
            percent = float(100*i)/float(row)
            printlog(True, "... Process %.2f percent in %.2f minutes"%(percent, (time.time() - start_time)/60) )
            '''
            try:
                printlog(True, 'Time BRDF: %.4f \tTime 6S: %.4f' %(time_brdf/cont_time_brdf, time_6S/cont_time_6S))
            except:
                print('')
            '''
        for j in range(0, col, 1):
            # Verify if pixel is inside Goias area
            if(vetMask[i][j] == 1 and data[0,0,i,j]*data[0,4,i,j] !=0):
                # Verify cloud mask in pixel data
                cont_cloud_free = 0
                cont_images = 0
                #printlog(True, '------- data:')
                for x in range(0, 20, 1):
                    # Verify if pixel had cloud
                    # Cloud Mask: 0 - Clear sky over water, 1 - Clear sky over land, 2 - Cloud, 3 - No data
                    #printlog(True, data[x,0:8,i,j])
                    if(data[x,7,i,j] <= 1 and abs(data[x,1,i,j] - data[x,3,i,j]) < 100 and data[x,0,i,j] != 0):
                        cont_cloud_free = cont_cloud_free + 1
                    else:
                        cont_images = cont_images + 1
                if(cont_cloud_free > 5):
                    #printlog(True,'[%d, %d] Pixel generate BRDF (%d/%d)'%(i, j, cont_cloud_free, cont_cloud_free+cont_images))
                    mat = numpy.zeros((cont_cloud_free, 8))
                    ccf = 0 # count cloud free to loop
                    for x in range(0, 20, 1):
                        # Verify if pixel had cloud
                        if(data[x,7,i,j] <= 1 and abs(data[x,1,i,j] - data[x,3,i,j]) < 100  and data[x,0,i,j] != 0):
                            # To put NDVI in data vector for analyse
                            try:
                                data[x,7,i,j] = (data[x,5,i,j]-data[x,4,i,j])/(data[x,5,i,j]+data[x,4,i,j])
                            except ZeroDivisionError:
                                data[x,7,i,j] = 0
                            mat[ccf] = data[x,0:8,i,j]
                            ccf = ccf + 1
                    '''
                    printlog(False,'-------- mat:')
                    printlog(False, '[VZA, VAZ, SZA, SAZ, RED, NIR, WIR, NDVI]')
                    printlog(False, mat)
                    '''
                    #start_time = time.time()
                    DATA_BRDF[0,i,j], DATA_BRDF[1,i,j], DATA_BRDF[2,i,j],DATA_BRDF[3,i,j],COEFS[i,j,:], MAT_BRDF = calc_BRDF(mat)
                    #print(".time do BRDF: %.4f s" %(time.time() - start_time))

                    if (DATA_BRDF[0,i,j] + DATA_BRDF[1,i,j] + DATA_BRDF[2,i,j] + DATA_BRDF[3,i,j])!= 0:
                        #time_brdf = time_brdf + (time.time() - start_time)
                        #cont_time_brdf = cont_time_brdf + 1
                        #start_time = time.time()
                        #DATA_6S[i,j]= sixsv1(VET, 'angles')
                        #print(".time do 6SV1: %.4f s" %(time.time() - start_time))
                        #time_6S = time_6S + (time.time() - start_time)
                        #cont_time_6S = cont_time_6S + 1
                        if(workQueue.full()):
                            if __debug__: print('f')
                            while(workQueue.full()):
                                pass
                        #Send VET_ATM = [Water_Vapour, Total Ozone, Deep_Blue, Aerosol]
                        ATM = [VET_ATM[3,i,j], VET_ATM[2,i,j], VET_ATM[0,i,j], VET_ATM[1,i,j]]
                        queueLock.acquire()
                        workQueue.put(i)
                        workQueue.put(j)
                        workQueue.put(MAT_BRDF)
                        workQueue.put(ATM)
                        queueLock.release()
                        if __debug__: print('.')
                '''
                else:
                    printlog(True,'[%d, %d] Pixel insuficient observations'%(i,j))
            else:
                if(vetMask[i][j] == 0):
                    printlog(True,'[%d, %d] Pixel out area'%(i,j))
                else:
                    printlog(True,'[%d, %d] Pixel border area'%(i,j))
            '''
    printlog(True, "----- Wait for queue to empty")
    while not workQueue.empty():
        #print('Files in Queue: ' + str(workQueue.qsize()))
        #time.sleep(5)
        pass
    printlog(True, '----- Queue is empty')

    # Notify threads it's time to exit
    exitFlag = 1

    # Wait for all threads to complete
    for t in threads:
        if (t.isAlive()):
            t.join(5.0)
                
    printlog(True, "----- Exiting Main Thread")
    exitFlag = 0

    # Save DATA BRDF in georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly )
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")    
    filename = outputDir + str(day) + '_BRDF.tif'
    dsOut = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 4, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    for i in range(0, 4, 1):
        bandOut = dsOut.GetRasterBand(i+1)
        BandWriteArray(bandOut, DATA_BRDF[i])
   
    # Save DATA 6S in georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly )
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")    
    filename = outputDir + str(day) + '_6S.tif'
    dsOut = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 4, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    for i in range(0, 4, 1):
        bandOut = dsOut.GetRasterBand(i+1)
        BandWriteArray(bandOut, DATA_6S[i])

    # Save COEFs in georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly )
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")    
    filename = outputDir + str(day) + '_COEFS.tif'
    dsOut = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 9, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    for i in range(0, 3, 1):
        for j in range(0, 3, 1):
            bandOut = dsOut.GetRasterBand(3*i+j+1)
            BandWriteArray(bandOut, COEFS[:,:,i,j])
    return True

class main():

    argDict = mapDict(argv, usage)
    lst = []
    global exitFlag
    global threads
    global workQueue

    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict:
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
    else:
        exit(usage)

    date = dateStart
    gc.collect()
    gdal.UseExceptions()
    start_time = time.time()
    exitFlag = 0

    numpy.set_printoptions(formatter={'float': '{: 0.4f}'.format})  
    date = dateStart
    while (date <= dateEnd):
        datestr = str(date)
        printlog(True, '----- Verify BRDF: %s' %(date))
        outputSubdir = outputDir + 'BRDF\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        inputSubdir  = inputDir + '___\\'   + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        if(os.path.isfile(outputSubdir + str(date) + '_BRDF.tif')):
            printlog(True, '... File exist in: %s' %(outputSubdir))
            date = incDate(date)
            continue
        else:
            printlog(True, '... Generate file: %s' %(date))
        if not(os.path.isdir(outputSubdir)):
            os.makedirs(outputSubdir)
        printlog(True, "----- Create new threads")
        threadID = 0
        for num_thread in range(0,NUM_THREADS,1):
            thread = myThread(threadID, 'T-' + str(num_thread), workQueue)
            thread.start()
            threads.append(thread)
            threadID += 1
        if(generateBRDF(inputSubdir, outputSubdir, date) == False):
            printlog(True, '----- Error to generate BRDF %s'%(date))
        date = incDate(date)
    printlog(True, '...Finish process')
    # Notify threads it's time to exit
    exitFlag = 1
import sys
sys.path.insert(0, 'C:\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
import numpy
import matplotlib.pyplot as plt
import itertools
import kernels
import os
import gc
import glob
import sys
from utils import * 
import pandas as pd

usage = """\
Usage: %s [OPTIONS]
        -ps     process (BRDF or 6S)
        -ty     type data (1day, 2days)
        -fi     folder input data
        -fo     folder output data
        -py     pixel y to process
        -px     pixel x to process
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
        -class  data class (culture or pasture)

Points to verify data:
        1. [143, 48] [-51.557, -17.376] [Rio Verde] Agricultura Anual  (Terraclass 1)
        2. [178, 14] [-52.751, -18.624] [Chapadão do Céu]
        3. [155, 32] [-52.105, -17.820] [Jatai]
        4. [144, 38] [-52.084, -17.668] [Jatai]
        5. [146, 33] [-52.068, -17.510] [Perolandia]
        6. [146, 58] [-51.202, -17.504] [Montividiu]
        7. [139, 57] [-51,256, -17.265] [Montividiu]
        8. [144, 67] [-50.916, -17.439] [Rio Verde]
        9. [134, 79] [-50.501, -17.089] [Parauna]
        10. [130, 63] [-51.028, -16.949] [Parauna]
        [84, 122] [-48.990, -15.330] Agricultura Perene (Terraclass 2)
        [16, 172] Natural Campestre  (Terraclass 5)
        [32, 159] Natural Florestal  (Terraclass 6)
        11. [35,  81] [-50.427, -13.634] [Novo Mundo] Pastagem           (Terraclass 11)
        12. [8 ,  83] [-50.361, -12.707] [São Miguel do Araguaia] 
        13. [32, 118] [-49.122, -13.531] [Porangatu]
        14. [34, 123] [-48.969, -13.607] [Santa Tereza de Goias]
        15. [50, 83]  [-50.352, -14.161] [Nova Crixas]
        16. [123, 81] [-50.426, -16.701] [Aurilandia]
        17. [130, 88] [-50.163, -16.923] [Jandaia]
        18. [156, 108][-49.460, -17.830] [Joviania]
        19. [15, 80]  [-50.436, -12.941] [Sao Miguel do Araguaia]
        20. [44, 188] [-46,698, -13,956] [Iaciara]
        [43, 177] Natural Savânica   (Terraclass 12) 
""" % argv[0]

# E O HISTOGRAMA DO HORÁRIO DA MELHOR AMOSTRA DO DIA

CULTURE = [[143,48],[178,14],[155,32],[144,38],[146,33],[146,58],[139,57],[144,67],[134,79],[130,63]]
PASTURE = [[35,81],[8,83],[32,118],[34,123],[50,83],[123,81],[130,88],[156,108], [15, 80], [44,188]]
TYPE_FILE = ['BRDF2', '6S', 'SMAC']
use_cols = ['DOY', 'HOUR_1330', 'MAX_DAY', 'MEAN_DAY','MAX_SPARSE_THIN', 'MAX_SPARSE_THICK']
#use_cols = ['DOY', 'MAX_SPARSE_THIN', 'MAX_SPARSE_THICK', 'MAX_DENSE_THIN',  'MAX_DENSE_THICK', 'MAX_ROUJEAN_THIN', 'MAX_ROUJEAN_THICK', 'MAX_GAO', 'MAX_PROUD', 'MAX_SILVEIRA']

def insert_in_cols (type, cols):
    # get end of cols
    end = len(cols)
    if(type == 'BRDF2'): 
        type = 'BRDF'
    # start in third col
    for x in range(end):
        if(cols[x] != 'DOY'):
            cols[x] = str(type) + '_NDVI_' + str(cols[x])
    return cols


class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    
    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict and (("-px" in argDict and "-py" in argDict) or "-class" in argDict):
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData = argDict["-ty"]
        if("-class" in argDict):
            classData = argDict["-class"]
        else:            
            MSG_pix_x = int(argDict["-px"])
            MSG_pix_y  = int(argDict["-py"])
            classData = None
    else:
        exit(usage)
    if(typeData != '1day' and typeData != '2days'):
        print('Error in typeData, verify!')
        exit(usage)

    if(classData == 'culture'):
        CLASS = CULTURE
    elif(classData == 'pasture'):
        CLASS = PASTURE
    elif(classData == None):
        CLASS = [[MSG_pix_x, MSG_pix_y]]
    else:
        exit(usage)

    for typeFile in TYPE_FILE:
        ALL = []
        for pos_class in range(len(CLASS)):
            MSG_pix_x = CLASS[pos_class][0]
            MSG_pix_y = CLASS[pos_class][1]

            # Create name files
            mod_filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_MOD13.csv'
            msg_filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_typ_'+ str(typeData) + '_NDVI_' + str(typeFile) + '.csv'
            out_filename = 'info_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_typ_'+ str(typeData) + '_NDVI_' + str(typeFile) + '.csv'
            cmp_filename = out_filename.replace('.csv', '_CMP.csv')
            all_filename =  'info_' + str(dateStart) +'_'+ str(dateEnd) + '_class_' + str(classData) + '_typ_'+ str(typeData) + '_NDVI_' + str(typeFile) + '_ALL.csv'
            
            # Read data files
            # MOD13  = [DOY NDVI RED NIR]
            if(os.path.isfile(inputDir + mod_filename)):
                MOD13 = pd.read_csv(inputDir + mod_filename, sep = " ", index_col = False)
            else:
                print("Error! Verify data to read: %s"%(inputDir + mod_filename))
                sys.exit() 

            if(os.path.isfile(inputDir + msg_filename)):
                type_cols = insert_in_cols (typeFile, use_cols)
                MSG_IN = pd.read_csv(inputDir + msg_filename, sep = " ", index_col = False, usecols = type_cols)
            else:
                print("Error! Verify data to read: %s"%(inputDir + msg_filename))
                sys.exit()
        
            # Process data
            print("Process data: %s"%(msg_filename))
            MSG_OUT = MSG_IN[MSG_IN.DOY > 0]
            MSG_OUT = MSG_IN[MSG_IN.DOY < 366]
            MSG_OUT.to_csv(outputDir + out_filename, sep=',', index = False, encoding='utf-8')
            MSG_COMP = MSG_OUT[MSG_OUT.DOY == 0]
            mod_doy_ant = 0
            for mod_doy in MOD13["DOY"]:
                # Verify most near value in msg_doy
                dif = 30
                if(mod_doy_ant > mod_doy):
                    break
                for msg_doy in MSG_OUT["DOY"]:
                    if(mod_doy < msg_doy):
                        break
                    if(mod_doy - msg_doy < dif):
                        dif = mod_doy - msg_doy
                        if (dif == 0):
                            break
                mod_doy_ant = mod_doy
                #print('mod_doy: %d msg_doy: %d'%(mod_doy, msg_doy))
                #print(MSG_SAV[MSG_SAV.DOY == msg_doy])
                MSG_COMP = MSG_COMP.append(MSG_OUT[MSG_OUT.DOY == msg_doy], ignore_index=True)
            #print(MSG_COMP)
            MSG_COMP = MSG_COMP.rename(columns={'DOY': 'DOY_MSG'})
            MOD13 = MOD13.rename(columns={'DOY': 'DOY_MOD'})
            #print('-----------COMP')
            #print(MOD13)
            #print(MSG_COMP)
            MOD13 = pd.concat([MSG_COMP, MOD13], axis = 1)
            MOD13 = MOD13.sort_index(axis=1)
            #print(MOD)
            MOD13.to_csv(outputDir + cmp_filename, sep=',', index = False, encoding='utf-8')
            if(len(ALL) == 0):
                ALL = MOD13
            else:
                ALL = ALL.append(MOD13, ignore_index=True)
            id = ALL.shape[0]
            ALL = ALL[:id-1]
        ALL.to_csv(outputDir + all_filename, sep=',', index = False, encoding='utf-8')
        del ALL
import sys
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\SMAC')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\SMAC')

import numpy
import matplotlib.pyplot as plt
import itertools
import kernels
import os
import gc
import glob
import sys
from utils import * 
from sys import argv
from kernels import Kernels # Import ther kernels
from Py6S import *
from SMAC_CESBIO import *
import time

otbDir = 'E:\\Install\\OTB-5.4.0-win64\\OTB-5.4.0-win64\\bin\\'
exeDir = 'C:\\OSGeo4W\\bin\\'
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoias123 = os.getcwd() + '\\shapes\\201401011300_123.tif'
fileGoiasHRV = os.getcwd() + '\\shapes\\201401011700_HRV??.tif' # Colocar arquivo atualizado
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMaskSHP = fileGoias
fileMaskHRV = fileGoiasHRV
fileMask123 = fileGoias123

usage = """\
Usage: %s [OPTIONS]
        -ty     type data (1day, 2days)
        -ps     type input data (123, 6S)
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


def importBRDFdata(inputDir, date, type):
    startHour = 1300
    endHour   = 1500    

    ds = gdal.Open(fileMask123)
    y_123 = ds.RasterYSize
    x_123 = ds.RasterXSize
    # DATA[17][8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
    DATA = numpy.zeros((17, 8, y_123, x_123))
    inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
    hour = startHour
    h = 0

    while(hour <= endHour):
        fileError = False
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
            filename = inputSubDir.replace('___', type) + str(date) + str(hour) + '_' + type + '.tif'
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

def calc_lstsq(BAND, kern):
    # Calc to RED and NIR
    samples = len(kern)
    len_data= len(BAND)
    NDVI = numpy.zeros((3, samples))
    COV  = [0,0]
    # Receive n_sample to process NDVI
    s_sample = 99
    for band in range(2):
        #print('passer: %d  kern: %d'%(len(passer), len(kern)))
        tmp_obs = BAND
        obs_col = numpy.zeros((len(tmp_obs),1))
        obs_lin = numpy.zeros((len(tmp_obs)))
        for x_obs in range(len(tmp_obs)):
            obs_col[x_obs] = tmp_obs[x_obs][band]
            obs_lin[x_obs] = tmp_obs[x_obs][band]

        if(samples > len_data):
            K = kern[len_data, :]
        else:
            K = kern
        r_sample = [samples]
        for x in range(samples-2, 4, -2):
            if(x >= 5):
                r_sample.append(x)
            else:
                break     
        
        #print('r_sample:')
        #print(r_sample)
        for n_sample in r_sample:
            (f, rmse, rank, svals ) = numpy.linalg.lstsq(K[:n_sample], obs_col[:n_sample])
            '''
            print('n_samples: %d'%(n_sample))
            print('K:')
            print(K[:n_sample])
            print('obs_col:')
            print(obs_col[:n_sample])

            print('----- Results Band[%d]' %(band))
            print ("%-20s %20s" % ( "Kernel", "Value"))
            for i, k in enumerate( ["Isotropic", "Ross-Thick", "Li-Sparce"] ):
                print ("%-20s %20f" % ( k, f[i] ))
            '''
            if(f[0] > 0):
                break

        # Error if k_iso <= 0
        if(f[0] <= 0):
            #print('...return [0] in else')       
            return [0]
        else:
            if(s_sample > n_sample):
                s_sample = n_sample

        fwd = K.dot(f)
        fwd_lin = numpy.zeros(len(fwd))
        for x_fwd in range(len(fwd)):
            fwd_lin[x_fwd] = fwd[x_fwd]
        COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
        NDVI[band] = fwd_lin
    # Calc NDVI = (NIR - RED)/(NIR + RED)
    NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
    #print('...return NDVI (s_sample: %d of %d):'%(s_sample, len(NDVI[0])))
    NDVI_filt = numpy.zeros((3,s_sample))
    for x in range(3):
        NDVI_filt[x] = NDVI[x][:s_sample]
    #print(NDVI_filt)
    return NDVI_filt

class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    
    if ("-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict and "-ps" in argDict):
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData = argDict["-ty"]
        processData = argDict["-ps"]
    else:
        exit(usage)
    if(typeData != '1day' and typeData != '2days'):
        print('Error in typeData, verify!')
        exit(usage)

    if(processData != '123' and processData != '6S'):
        print('Error in processData, verify!')
        exit(usage)


    # Go out if folder not exist
    if(os.path.isdir(inputDir) == False or os.path.isdir(outputDir) == False):
        print('... Error in folder adress!')
        sys.exit()
    # Create output data matrix
    ds = gdal.Open(fileMask123)
    y_123 = ds.RasterYSize
    x_123 = ds.RasterXSize
    KERNEL = numpy.zeros((y_123, x_123, 5, 5))
    ds = None   
    dateSaveKernel = dateStart
    # Process in range of data
    date = dateStart
    while(date < dateEnd):
        # Output file name   
        outputFile = outputDir + 'BRDF_' + processData + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        # Verify and create folders
        if(not os.path.isdir(outputFile)):
            os.makedirs(outputFile)
        kernelFile = outputFile + str(date) + '_KERN.npy'
        if(date == dateStart):
            if(os.path.isfile(kernelFile)):
                KERNEL = numpy.load(kernelFile)
                print('Reading Kernel file data: %s'%(kernelFile))
        outputFileBRDF = outputFile + str(date) + '_BRDF.tif'    

        # Verify if files exist
        if(os.path.isfile(outputFileBRDF)):
            print('...File exist in: %s'%(outputFileBRDF))
            date = incDate(date)
            continue

        NDVI_BRDF = numpy.zeros((5, y_123, x_123))

        # DATA[17][8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
        DATA = importBRDFdata(inputDir, date, processData)
        if(len(DATA) == 0):
            print('..Error! Don\'t is possible create %s_BRDF.tif' %(date))
            date = incDate(date)
            continue
        print('...Create BRDF data: %s'%(date))
        tac = time.clock()
        hours = len(DATA)
        #print('hours:%d [%d, %d]'%(hours, y_123, x_123))
        for x in range(x_123):
            for y in range(y_123):
                if(DATA[0][0][y][x] == 0):
                    continue

                # Copy data to vector process
                VZA = numpy.zeros((hours,1))
                VAZ = numpy.zeros((hours,1))
                SZA = numpy.zeros((hours,1))
                SAZ = numpy.zeros((hours,1))
                RAA = numpy.zeros((hours,1))
                RED = numpy.zeros((hours,1))
                NIR = numpy.zeros((hours,1))
                CLM = 2*numpy.ones((hours,1))
                BAND = numpy.zeros((hours,2))
                HOUR = numpy.zeros((hours,1)) 
                NDVI = numpy.zeros((hours,1))
                hour = 1300

                for h in range (hours):
                    if(DATA[h][0][y][x] != 0):
                        VZA[h]  = DATA[h][0][y][x]
                        VAZ[h]  = DATA[h][1][y][x]
                        SZA[h]  = DATA[h][2][y][x]
                        SAZ[h]  = DATA[h][3][y][x]
                        RED[h]  = DATA[h][4][y][x]
                        NIR[h]  = DATA[h][5][y][x]
                        CLM[h] = DATA[h][7][y][x]
                        HOUR[h] = hour
                        BAND[h,0] = RED[h]
                        BAND[h,1] = NIR[h]
                        NDVI[h] = calcNDVI(NIR[h], RED[h])
                        # WIR = DATA[h][0][y][x]
                        hour = incHour(hour)
                RAA  = SAZ - VZA
                # Verify number of samples to process BRDF
                if(processData == '123'):
                    passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(RED > 0, NIR > 0))
                else:
                    passer = numpy.logical_and(RED > 0, NIR > 0)
                samples = len(passer[passer == True])
                if(samples == 0):
                    continue
                '''    
                NDVI_mean = numpy.mean(NDVI)
                SAZ_mean  = numpy.mean(SAZ[passer1])
                print('NDVI_mean: %.4f SAZ_mean: %.4f'%(NDVI_mean, SAZ_mean)) 
                if(processData == '123'): 
                    passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(RED > 0, numpy.logical_and(NIR> 0, \
                         numpy.logical_and (NDVI > 0.8*NDVI_mean, numpy.logical_and(NDVI < 1.2*NDVI_mean, \
                         numpy.logical_and (SAZ > 0.2*SAZ_mean, SAZ < 1.6*SAZ_mean))))))
                else:
                    passer = numpy.logical_and(RED > 0, numpy.logical_and(NIR> 0, \
                         numpy.logical_and (NDVI > 0.8*NDVI_mean, numpy.logical_and(NDVI < 1.2*NDVI_mean, \
                         numpy.logical_and (SAZ > 0.2*SAZ_mean, SAZ < 1.6*SAZ_mean)))))
                samples = len(passer[passer == True])
                if(samples == 0):
                    print('...continue in samples == 0 second test')
                    date = incDate(date)
                    continue
                '''
                VZA = VZA[passer]
                SZA = SZA[passer]
                RAA = RAA[passer]
                VAZ = VAZ[passer]
                SAZ = SAZ[passer]
                RED = RED[passer]
                NIR = NIR[passer]
                CLM = CLM[passer]
                NDVI = NDVI[passer]
                HOUR = HOUR[passer]
                if (samples >= 5):                  
                    BAND = numpy.zeros((len(RED),2))
                    BAND[:,0] = RED
                    BAND[:,1] = NIR
                else:
                    BAND = numpy.zeros((5,2))
                    for i in range(len(RED)):
                        BAND[i,0] = RED[i]
                        BAND[i,1] = NIR[i]                              

                id_brdf   = 0 
                red_brdf  = 0
                nir_brdf  = 0
                ndvi_brdf = 0
                hour_brdf = 0

                geo_kernel = 'Sparse'
                vol_kernel = 'Thick'
                flag_kernel = False
                if(samples >= 5):
                    #print('--------- Calc NDVI using Ambrals (2000)')
                        
                    # Generate the semiempirical kernels
                    K_obs =  Kernels( VZA, SZA, RAA, LiType=geo_kernel, doIntegrals=False, \
                        normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, RossType= vol_kernel )
                    kern = numpy.ones (( numpy.sum(passer==True), 3 )) # Store the kernels in an array
                    #print(K_obs)
                    kern[ :, 1 ] = K_obs.Ross
                    kern[ :, 2 ] = K_obs.Li
                    NDVI = calc_lstsq(BAND, kern)
                    # If kernel is negative not use in the calc
                    if(len(NDVI) > 1):
                        # Did calc BRDF in this pixel
                        flag_kernel = True
                        samples = len(NDVI[0])
                        # Calc best NDVI
                        id_brdf   = numpy.argmax(NDVI[2]) 
                        red_brdf  = NDVI[0, id_brdf]
                        nir_brdf  = NDVI[1, id_brdf]
                        ndvi_brdf = NDVI[2, id_brdf]
                        hour_brdf = HOUR[id_brdf]
                        qual_brdf = 0
                        # save BRDF KERNEL data to process future day
                        id_ord = numpy.argsort(NDVI[2])
                        id_ord = id_ord[::-1]
                        n = 0
                        for i in id_ord[:5]:
                            KERNEL[y,x,n] = [kern[i, 0], kern[i, 1], kern[i,2], RED[i], NIR[i]]
                            n = n + 1
                    else:
                        continue
                if(flag_kernel == False):
                    if(samples > 5):
                        samples = 3
                    if(numpy.mean(KERNEL[y,x,0]) == 0 or numpy.mean(KERNEL[y,x,1]) == 0):
                        continue
                    # Generate the semiempirical kernels
                    K_obs =  Kernels( VZA, SZA, RAA, LiType=geo_kernel, doIntegrals=False, \
                        normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, RossType= vol_kernel )
                    kern = numpy.ones ((5, 3)) # Store the kernels in an array
                    #print(K_obs)
                    for i in range(len(RED[:samples])):
                        kern[i,1] = K_obs.Ross[i]
                        kern[i,2] = K_obs.Li[i]
                    for i in range(5):
                        if(samples <= i):
                            kern[i] = KERNEL[y,x,i,:3]
                            BAND[i] = KERNEL[y,x,i,3:]
                    NDVI = calc_lstsq(BAND, kern)
                     # If kernel is negative not use in the calc
                    if(len(NDVI) > 1):
                        #samples = len(NDVI[0])
                        # Calc best NDVI in this day
                        id_brdf   = numpy.argmax(NDVI[2, :samples])
              
                        red_brdf  = NDVI[0, id_brdf]
                        nir_brdf  = NDVI[1, id_brdf]
                        ndvi_brdf = NDVI[2, id_brdf]
                        hour_brdf = HOUR[id_brdf]
                        qual_brdf = 1
                    else:
                        continue
                '''        
                print('samples: %d id_brdf: %d'%(samples, id_brdf))
                print('y: %d x: %d date %s'%(y, x, date))
                print('BAND')
                print(BAND[:samples])
                print('kern:')
                print(kern[:samples]) 
                print('samples:')
                print(samples)
                print('NDVI')
                print(NDVI)
                print('HOUR:')
                print(HOUR[:samples])
                '''
                NDVI_BRDF[0][y][x] = ndvi_brdf
                NDVI_BRDF[1][y][x] = red_brdf
                NDVI_BRDF[2][y][x] = nir_brdf
                NDVI_BRDF[3][y][x] = hour_brdf
                NDVI_BRDF[4][y][x] = qual_brdf
        #print('.Calc data time: %.2f s' %(time.clock() - tac))
        saveRASTERfile(outputFileBRDF, fileMask123, NDVI_BRDF)
        # Save kernel file data in 5 to 5 days
        if(calcDays(dateSaveKernel, date) > 5):
            numpy.save(kernelFile, KERNEL)
            dateSaveKernel = date
        #print('.Save data time: %.2f s' %(time.clock() - tac))
        #sys.exit()
        date = incDate(date)
# ftp-ex.py
# http://www.blog.pythonlibrary.org/2012/07/19/python-101-downloading-a-file-with-ftplib/
# http://www.pythonforbeginners.com/code-snippets-source-code/how-to-use-ftp-in-python
# https://docs.python.org/2/library/ftplib.html
# https://docs.python.org/2/tutorial/errors.html

import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import os
from ftplib import FTP
import os.path
import tarfile
import gc
from utils import printlog, sleepTime, mapDict
from time import strftime, localtime
from sys import argv
import threading
import queue
import time

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to output data
        -ty     type input data (in FTP server: GBR, HRIT, HDF)
""" % argv[0]

threadList = ["Thread-1", "Thread-2", "Thread-3"] #, "Thread-4", "Thread-5"]
queueLock = threading.Lock()
workQueue = queue.Queue(10)
threads = []
exitFlag = 0

def ftpCopy(downloadDir, typeData):
    try:
        ftp = FTP("ftp.lapig.iesa.ufg.br", "lapig-ftp", "lapig123ftp")
    except:
        printlog(True, 'Error! Don''t have Connection to Internet.')
        return(-1)
    try:
        ftp.login()
        printlog(True, 'Logged in ftp server')
        printlog(ftp.getwelcome())
    except:
        if(ftp.getwelcome()!= None):
            printlog(True, 'Logged in ftp server')
            printlog(True, ftp.getwelcome())
        else:
            printlog(True, 'Error! Verify user name and password')
            return(-2)
        printlog(False, 'Finally in login')

    # To cloud mask filtes
    if(typeData == 'GRB'):
        filter1 = 'MSG3-SEVI-MSGCLMK-0100-0100'
        filter2 = '.grb'
    elif(typeData == 'HRIT'):
        filter1  = 'MSG3-SEVI-MSG15-0100-NA-201'
        filter2  = '.tar'
    elif(typeData == 'HDF'):
        filter1 = 'MSG3-SEVI-MSGNDVE-0100-0100-201'
        filter2 = '.h5'
    print(downloadDir)
    folder = typeData
    print('--------------')
    #ftp.retrlines("LIST")
    ftp.cwd("eumetsat")
    ftp.cwd(folder)
    # ftp.cwd("subFolder")
    # or ftp.cwd("folderOne\\subFolder")
    print('--------------')
    #ftp.retrlines("LIST")
    listing = []
    lstFiles = []
    ftp.retrlines("LIST", listing.append)
    poslist = -1
    filename = []
    print('--------------')
    print('Number of files: ' + str(len(listing)))
    #print(listing)
    while(poslist < len(listing)):

        try:
            while True:
                poslist = poslist + 1
                print('.... File: ' + str(poslist))
                words = listing[poslist].split(None, 8)
                filename = words[-1].lstrip()
                # Filter to year data
                #print(filename)
                if(filter1 in str(filename) and  filter2 in str(filename)):
                    #print('!! pass in filters !!')
                    break
        except:
            if (poslist > len(listing)):
                printlog(True,'Finish copy files')
                break
        # download the file
        fileType = True
        #print(filename)
        try:
            id = filename.index('-201')
            #print('id in file')
        except:
            fileType = False
            printlog(True, 'File type different!')
        if(fileType == True):
            year  = filename[id+1:id+5]
            month = filename[id+5:id+7]
            subfolder = year + '\\' + month + '\\'
            #print(subfolder)
            if not(os.path.isdir(downloadDir + subfolder)):
                os.makedirs(downloadDir + subfolder)
            if(os.path.isfile(downloadDir + subfolder + filename) == True):
                printlog(True, 'File yet copy: ' + filename)
            else:
                printlog(True, 'File to copy: ' +  filename)
                lstFiles.append(filename)

    ftp.quit()
    return (lstFiles)

class myThread (threading.Thread):
    def __init__(self, threadID, name, q, typeData, downloadDir):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        self.downloadDir = downloadDir
        self.typeData = typeData
    def run(self):
        printlog (True, "Starting " + self.name)
        #printlog(True, "q: " + str(self.q) + " typeData: " + self.typeData + " donwloadDir: " + self.downloadDir)
        process_data(self.name, self.q, self.typeData, self.downloadDir)
        printlog (True, "Exiting " + self.name)

def process_data(threadName, q, typeData, downloadDir):
    global queueLock
    global workQueue
    global exitFlag
    #printlog(True, "Enter in process_data. exitFlag: " + str(exitFlag))

    while not exitFlag:
        #printlog(True, "exitFlag is False")
        queueLock.acquire()
        #printlog(True, "\nqueueLock acquire")
        if not workQueue.empty():
            #printlog(True, "Queue is not empty")
            filename = q.get()
            queueLock.release()
            #printlog(True, "\nqueueLock release")
            printlog (True, "\n%s processing %s" % (threadName, filename))
            copyFTPfile(typeData, filename, downloadDir)
        else:
            #printlog(True, "Queue is empty")
            try:
                queueLock.release()
                #printlog(True, "\nqueueLock release")
            except:
                printlog(True, 'Fail to release queueLock !')
        time.sleep(1)
    printlog(True, "exitFlag is True")

def copyFTPfile(typeData, inputFile, outputDir):
    try:
        ftp = FTP("ftp.lapig.iesa.ufg.br", "lapig-ftp", "lapig123ftp")
    except:
        printlog(True, 'Error! Don''t have Internet Connection.')
        return(-1)
    try:
        ftp.login()
        printlog(False, 'Logged in ftp server')
        printlog(False, ftp.getwelcome())
    except:
        if(ftp.getwelcome()!= None):
            printlog(False, 'Logged in ftp server')
            printlog(False, ftp.getwelcome())
        else:
            printlog(False, 'Error! Verify user name and password')
            return(-2)
        printlog(False, 'Finally in login')
    folder = typeData
    ftp.cwd("eumetsat")
    ftp.cwd(folder)
    id = inputFile.index('-201')
    year  = inputFile[id+1:id+5]
    month = inputFile[id+5:id+7]
    subfolder = year + '\\' + month + '\\'
    outputFile = outputDir + subfolder + inputFile
    #printlog(True, "Create FTP data file")
    lf = open(outputFile, "wb")
    ftp.retrbinary("RETR " + inputFile, lf.write, 8*1024)
    lf.close()

class main:

    global exitFlag
    global threads
    global workQueue
    threadID = 1

    argDict = mapDict(argv, usage)

    if "-ty" in argDict and "-fo" in argDict:
        typeData = argDict["-ty"]
        downloadDir = argDict["-fo"]
    else:
        exit(usage)     
    if not(typeData == 'GRB' or typeData == 'HRIT' or typeData == 'HDF'):
        exit(usage)
    #print('-------- teste')
    #print(downloadDir)

    while(1):
        gc.collect()
        lst = ftpCopy(downloadDir, typeData)
        print('List files: ' + str(len(lst)))
        #print(lst)

        try:
            if(lst == -1):
                printlog(True, 'Waiting internet connection')
                sleepTime(0.1)
        except:
            break
        try:
            if(lst == []):    
                printlog(True, 'Waiting new files')
                sleepTime(1)
        except:
            break
        printlog(True, "----- Create new threads")
        for tName in threadList:
            thread = myThread(threadID, tName, workQueue, typeData, downloadDir)
            thread.start()
            threads.append(thread)
            threadID += 1

        printlog(True, "----- Fill the queue and wait case many files")
        
        x = 0
        while(x < len(lst)):
            queueLock.acquire()
            if not workQueue.full():
                filename =  lst[x]
                #printlog(True, 'x: ' + str(x) + '\tfilename: ' + filename)
                workQueue.put(filename)
                x = x + 1
                
            else:
                time.sleep(30)
            queueLock.release()

        printlog(True, "----- Wait for queue to empty")
        while not workQueue.empty():
            #print('Files in Queue: ' + str(workQueue.qsize()))
            #time.sleep(5)
            pass
        printlog(True, '----- Queue is empty')

        # Notify threads it's time to exit
        exitFlag = 1

    # Wait for all threads to complete
        for t in threads:
            t.join()
        printlog(True, "Exiting Main Thread")
        printlog(True, 'Waiting new files')
        sleepTime(1)
    
import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import sys
import subprocess
import os
import glob
from utils import *

inputDir = 'z:\\OUTROS\\INMET\\'
subfolder = 'BCK_ANO_year\\'
filename  = inputDir + subfolder + 'BCK_ANO_year_EST_station.txt' 
outfile   = inputDir + 'data_all_year_all_stations.csv' 
stat_csv  = inputDir + 'stations_goias.csv'
outputDir = 'z:\\OUTROS\\PREC_DAILY\\'
grid_stations ='C:\\Users\\carlos.silveira\\Dropbox\\newPython\\kriging\\grid_stations.csv'

exeDir = 'C:\\osgeo4w\\bin\\'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'

def calcDoy(Y,M,D):
	""" given year, month, day return day of year
		Astronomical Algorithms, Jean Meeus, 2d ed, 1998, chap 7 """
	if is_leap_year(Y):
		K = 1
	else:
		K = 2
	N = int((275 * M) / 9.0) - K * int((M + 9) / 12.0) + D - 30
	return N

def is_leap_year(year):
	""" if year is a leap year return True
		else return False """
	if year % 100 == 0:
		return year % 400 == 0
	return year % 4 == 0

def calcMonth(Y,N):
	""" given year = Y and day of year = N, return year, month, day
		Astronomical Algorithms, Jean Meeus, 2d ed, 1998, chap 7 """    
	if is_leap_year(Y):
		K = 1
	else:
		K = 2
	M = int((9 * (K + N)) / 275.0 + 0.98)
	if N < 32:
		M = 1
	D = N - int((275 * M) / 9.0) + K * int((M + 9) / 12.0) + 30
	return M
	#return Y, M, D

def calcDay(Y, N):
	""" given year = Y and day of year = N, return year, month, day
		Astronomical Algorithms, Jean Meeus, 2d ed, 1998, chap 7 """    
	if is_leap_year(Y):
		K = 1
	else:
		K = 2
	M = int((9 * (K + N)) / 275.0 + 0.98)
	if N < 32:
		M = 1
	D = N - int((275 * M) / 9.0) + K * int((M + 9) / 12.0) + 30
	return D
	#return Y, M, D 

class main():
	#STATIONS = [NAME, LONG, LAT]
	STATIONS = numpy.loadtxt(stat_csv, delimiter = ',', skiprows = 1)
	n_stations = len(STATIONS)
	# HEADER has all variables 
	HEADER = ['STATION','YEAR','MON','DAY','OBS.','BAT','TEMP_CPU','AIR_INST.','AIR_MAX','AIR_MIN','UMID_INST','UMID_MAX','UMID_MIN', 'DP_INST','DP_MAX','DP_MIN','PRES_INST', 'PRES_MAX','PRES_MIN', 'WIND_Speed','WIND_Dir','WIND_Gust','RAD','PREC',	'CLOUD_TOT','CLOUD_CODE','CLOUD_BASE','CLOUD_VISIB']
	# VARS has variables to get
	VARS = ['YEAR','MON','DAY','OBS.','PREC']
	data_header = ['STAT','YEAR', 'MON','DAY','PREC']
	usecols = []
	for i, var in zip(range(len(HEADER)), HEADER):
		if var in VARS:
			usecols.append(i)
	print('n_stations: %d usecols of data: %s' %(n_stations, usecols))
	if(not os.path.isfile(outfile)):
		DATA = numpy.zeros((366*3*n_stations,len(VARS)))
		n_data = 0
		for year in range(2013, 2016, 1):
			for x in range(n_stations):
				file_station = filename
				file_station = file_station.replace('year', str(year))
				file_station = file_station.replace('station', 'A' + str('%03d'%int(STATIONS[x][0])))
				print('. reading file: %s'%(file_station))
				TAB = numpy.genfromtxt(file_station, delimiter = ' ', usecols = usecols, invalid_raise=False, unpack=True, missing_values=['//////', '/////', '//', '/', '=', ''])
				#TAB = numpy.transpose(TAB)
				#print(TAB)
				for i in range(len(VARS)):
					globals()[VARS[i]] = TAB[i]
				
				for month in range(1,13,1):
					for day in range(1,32,1):
						paser = numpy.logical_and(MON==month,DAY==day)
						if (numpy.sum(paser) > 0):
							DATA[n_data] = [STATIONS[x][0], year, month, day,numpy.sum(PREC[paser])]
							#print(DATA[n_data])
							n_data = n_data + 1
		header = 'STAT'
		for var in VARS:
			if(var != 'OBS.'):
				header = header + ' ' + str(var) 
		print(header)
		numpy.savetxt(outfile, DATA[:n_data], delimiter=" ", header = header, comments = '', fmt='%.1f')
	else:	
		DATA = numpy.loadtxt(outfile, delimiter = ' ', skiprows = 1, unpack = True)

	for i in range(len(VARS)):
		globals()[VARS[i]] = None
	for i in range(len(data_header)):
		globals()[data_header[i]] = DATA[i]

	for year in range(2013, 2016,1):
		for month in range(1,13,1):
			for day in range(1,32,1):
				# Verify if file exist
				outputSubfolder = outputDir + str('%04d'%int(year)) + '\\' + str('%02d'%int(month)) + '\\' 
				if(not os.path.isdir(outputSubfolder)):
					os.makedirs(outputSubfolder)
				outfile = outputSubfolder + str('%04d'%int(year)) + str('%02d'%int(month)) + str('%02d'%int(day)) + '_prec.tif'
				tmpfile1 = outfile.replace('.tif', '_tmpfile1.tif')
				tmpfile2 = outfile.replace('.tif', '_tmpfile2.tif')
				tmpfile3 = outfile.replace('.tif', '_tmpfile3.tif')
				if(os.path.isfile(outfile)):
					#print('file exist: %s'%(outfile))
					TEST = readFileBand(outfile, 1)
					has_nan = numpy.any(numpy.isnan(TEST))
					if(has_nan):
						os.remove(outfile)
						print('...REMOVING [has nan] file %s'%(outfile))
					else:
						continue
				print('process file: %s'%(outfile))
				print('...create grid stations file')
				if(os.path.isfile(grid_stations)):
					os.remove(grid_stations)
					#print('deleting grid_stations file %s'%(grid_stations))
				try:
					file = open(grid_stations,'w')
				except:
					print('... Error to open grid_stations file')
					sys.exit()
				file.write('LONG,LAT,PREC\n')
				if(os.path.isfile(outfile)):
					os.remove(outfile)
				n_error = 0
				for i, stat in zip(range(n_stations), numpy.transpose(STATIONS)[0]):
					#print('search in %d of %d %d %d'%(stat, year, month, day))
					#print(STAT[0])
					#print(YEAR[0])
					#print(MON[0])
					paser = numpy.logical_and(STAT==stat, numpy.logical_and(numpy.logical_and(YEAR==year, MON==month), DAY==day))
					if(numpy.sum(paser) == 0):
						n_error = n_error + 1
						#print('error...')
					else:
						prec_station = numpy.sum(PREC[paser])
						if(not numpy.isnan(prec_station) or prec_station < 0):
							prec_day = str('%.6f,%.6f,%.2f'%(STATIONS[i][1], STATIONS[i][2], prec_station))
							print(prec_day)
							if(i != n_stations-1):
								prec_day = prec_day + '\n'
							file.write(prec_day)
						else:
							n_error = n_error + 1
				file.close()
				if(n_error > 0.1*n_stations):
					print('...error to generate file !! [n_error: %d]'%(n_error))
					continue
				print('...process create image')
				
				print(exeDir + 'gdal_grid -a invdist:power=2.0:smoothing=1.0 -zfield PREC -l grid_stations  grid_stations.vrt ' + ' ' + tmpfile1)
				subprocess.call(exeDir + 'gdal_grid -a invdist:power=2.0:smoothing=1.0 -zfield PREC -l grid_stations grid_stations.vrt '+ ' ' + tmpfile1)

				#print(exeDir + 'gdalwarp.exe -t_srs "+proj=latlong +datum=WGS84" ' + tmpfile1 + ' ' + tmpfile2)
				#subprocess.call(exeDir + 'gdalwarp.exe -t_srs "+proj=latlong +datum=WGS84" ' + tmpfile1 + ' ' + tmpfile2)
				
				print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "cubic" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
				subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "bilinear" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 

				print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
				subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)

				print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outfile)
				subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outfile)
				#sys.exit()
				delFile(tmpfile1)
				delFile(tmpfile2)
				delFile(tmpfile3)
				#sys.exit()
'''
As described corresponding section of the manual - if one wants to perform pansharpening using GUI (Monteverdi), one have to do it in two steps:
1. Use Superimpose module to make multi-spectral imagery to have the same extent and resolution as the panchormatic one.
2. Use Simple RCS Pansharpening module on panchromatic and resampled multi-spectral imagery.

ECHO Bundle command

E:\Install\OTB-5.4.0-win64\OTB-5.4.0-win64\bin\otbcli_BundleToPerfectSensor.bat -inp E:\DADOS\HRV\2013\01\201301011700_HRV.tif -inxs E:\DADOS\123\2013\01\201301011300_123.tif -elev.dem E:\DADOS\SRTM\ -mode default -out E:\DADOS\PAN\2013\01\201301011300_PAN1.tif 

ECHO Pansharp command
REM E:\Install\OTB-5.4.0-win64\OTB-5.4.0-win64\bin\otbcli_Pansharpening.bat -inp E:\DADOS\HRV\2013\01\201301011700_HRV.tif   -inxs E:\DADOS\PAN\2013\01\201301011300_PAN1.tif -out E:\DADOS\PAN\2013\01\201301011300_PAN2.tif -method rcs
E:\Install\OTB-5.4.0-win64\OTB-5.4.0-win64\bin\otbcli_Pansharpening.bat -inp E:\DADOS\HRV\2013\01\201301011700_HRV.tif   -inxs E:\DADOS\PAN\2013\01\201301011300_PAN1.tif -out E:\DADOS\PAN\2013\01\201301011300_PAN3.tif -method lmvm 

'''

import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

from sys import argv                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
from utils import *
import subprocess
import os
import numpy
from scipy.special import factorial
import random
from scipy.signal import savgol_filter

usage = """\
Usage: %s [OPTIONS]
		-ty     type of data to process (DAY, TAB, PAN, PRO, NVDI, VIEW, MEAN, HIST)
				DAY: Generate day data using range hour 1300 to 1330
				DAY_INT: Do NDVI composite of n days
				PREC_INT: Do PREC sum of n days
				TAB_ALL  : Generate table all days getting DAY_123, DAY_6S, BRDF_123, BRDF_6S, MOD13
				TAL_MOD13: Generate table with MOD13 days getting DAY_123, DAY_6S, BRDF_123, BRDF_6S, MOD13
				PAN: Generate PAN data using HRV file >> It's necessary update!!
				PRO: Create PAN data using HRV image day (1700 to 1900) 
				NDVI: Fusion data BRDF to PAN
				VIEW: Generate table data (month, day, n_brdf, n_pro, n_ndvi)
				MEAN: Calculate month's NDVI data using BRDF data
				HIST: Calculate year's histogram data
		-ps     select process type (123 or 6S)
		-fi     folder input data
		-fo     folder output data
		-ds     date to start process (format: yyyyMMdd)
		-de     date to end process (format: yyyyMMdd)

""" % argv[0]

otbDir = 'C:\\OTB-5.4.0-win64\\bin\\'                                                                                                                                                                                                                                                                                                                     
exeDir = 'C:\\OSGeo4W\\bin\\'
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoias123 = os.getcwd() + '\\shapes\\201401011300_123.tif'
fileGoiasHRV = os.getcwd() + '\\shapes\\201301011700_HRV.tif'
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileMaskPAN   = os.getcwd() + '\\shapes\\201301011700_PAN.tif'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
fileMaskMOD13 = os.getcwd() + '\\shapes\\maskGoiasMOD13.tif'

fileMaskSHP = fileGoias
fileMaskHRV = fileGoiasHRV
fileMask123 = fileGoias123

def linreg(X, Y):
    """
    return a,b in solution to y = ax + b such that root mean square distance between 
    trend line and original points is minimized
    """
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in zip(X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    return (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det

def conv_3km_to(res_data, r3km_y, r3km_x):

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

	500 m (MOD09)
	Size is 1929, 1865
	Origin = (-53.252403461635794,-12.393279067399060)
	Pixel Size = (0.003810953930405,-0.003810953930405)
	Corner Coordinates:
	Upper Left  ( -53.2524035, -12.3932791) ( 53d15' 8.65"W, 12d23'35.80"S)
	Lower Left  ( -53.2524035, -19.5007081) ( 53d15' 8.65"W, 19d30' 2.55"S)
	Upper Right ( -45.9010733, -12.3932791) ( 45d54' 3.86"W, 12d23'35.80"S)
	Lower Right ( -45.9010733, -19.5007081) ( 45d54' 3.86"W, 19d30' 2.55"S)
	Center      ( -49.5767384, -15.9469936) ( 49d34'36.26"W, 15d56'49.18"S)
	
	250 m (MOD13)
	Size is 3515, 3399
	Origin = (-53.252606259537075,-12.393583326851415)
	Pixel Size = (0.002091294669953,-0.002091294669953)
	Corner Coordinates:
	Upper Left  ( -53.2526063, -12.3935833) ( 53d15' 9.38"W, 12d23'36.90"S)
	Lower Left  ( -53.2526063, -19.5018939) ( 53d15' 9.38"W, 19d30' 6.82"S)
	Upper Right ( -45.9017055, -12.3935833) ( 45d54' 6.14"W, 12d23'36.90"S)
	Lower Right ( -45.9017055, -19.5018939) ( 45d54' 6.14"W, 19d30' 6.82"S)
	'''

	r3km_py  = -12.383945170790270
	r3km_px  = -53.271811612492662
	r3km_psy = -0.034840758326260
	r3km_psx = 0.034840758326260
	r3km_sy  = 204
	r3km_sx  = 211

	r500m_py  = -12.393279067399060
	r500m_px  = -53.252403461635794
	r500m_psy = -0.003810953930405
	r500m_psx = 0.003810953930405
	r500m_sy  = 1865
	r500m_sx  = 1929

	r250m_py  = -12.393583326851415
	r250m_px  = -53.252606259537075
	r250m_psy = -0.002091294669953
	r250m_psx = 0.002091294669953
	r250m_sy  = 3399
	r250m_sx  = 3515

	if(res_data == '250m'):
		ret_x = round(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r250m_px - r250m_psx/2)/abs(r250m_psx))
		#print(r3km_px + r3km_x*r3km_psx + r3km_psx/2 )
		#print(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r250m_px - r250m_psx/2))
		#print(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r250m_px - r250m_psx/2)/abs(r250m_psx))
		#print(round(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r250m_px - r250m_psx/2)/abs(r250m_psx)))

		ret_y = round(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r250m_py - r250m_psy/2)/abs(r250m_psy))
		#print(r3km_py + r3km_y*r3km_psy + r3km_psy/2 )
		#print(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r250m_py - r250m_psy/2))
		#print(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r250m_py - r250m_psy/2)/abs(r250m_psy))
		#print(round(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r250m_py - r250m_psy/2)/abs(r250m_psy)))

		if(ret_y > r250m_sy or ret_y > r250m_sy):
			print('Error exceding array 250m [%d, %d]'%(r250m_y, r250m_y))

	if(res_data == '500m'):
		ret_x = round(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r500m_px - r500m_psx/2)/abs(r500m_psx))
		#print(r3km_px + r3km_x*r3km_psx + r3km_psx/2 )
		#print(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r500m_px - r500m_psx/2))
		#print(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r500m_px - r500m_psx/2)/abs(r500m_psx))
		#print(round(abs(r3km_px + r3km_x*r3km_psx + r3km_psx/2 - r500m_px - r500m_psx/2)/abs(r500m_psx)))

		ret_y = round(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r500m_py - r500m_psy/2)/abs(r500m_psy))
		#print(r3km_py + r3km_y*r3km_psy + r3km_psy/2 )
		#print(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r500m_py - r500m_psy/2))
		#print(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r500m_py - r500m_psy/2)/abs(r500m_psy))
		#print(round(abs(r3km_py + r3km_y*r3km_psy + r3km_psy/2 - r500m_py - r500m_psy/2)/abs(r500m_psy)))

		if(ret_y > r500m_sy or ret_y > r500m_sy):
			print('Error exceding array 500m [%d, %d]'%(r500m_y, r500m_y))

	if(res_data == '3km'):
		ret_x = r3km_px + r3km_x*r3km_psx + r3km_psx/2
		ret_y = r3km_py + r3km_y*r3km_psy + r3km_psy/2

	return [ret_y, ret_x]

def process_mod_clm(x, y, MOD09_CLM, MOD09_LIM):

	r3km_py  = -12.383945170790270
	r3km_px  = -53.271811612492662
	r3km_psy = -0.034840758326260
	r3km_psx =  0.034840758326260

	r500m_py  = -12.393279067399060
	r500m_px  = -53.252403461635794
	r500m_psy = -0.003810953930405
	r500m_psx =  0.003810953930405
	r500m_sy  = 1865
	r500m_sx  = 1929

	start_px = round(abs(((r3km_px + r3km_psx*x + r500m_psx/2)-r500m_px)/r500m_psx))
	start_py = round(abs(((r3km_py + r3km_psy*y + r500m_psy/2)-r500m_py)/r500m_psy))
	
	num_mod_clr = 0
	num_mod_clm = 0
	num_mod = 0
	for y_mod in range(start_py, start_py+8, 1):
		for x_mod in range(start_px, start_px+8, 1):
			# Get data if inside Goias
			if(MOD09_LIM[y_mod][x_mod] == 1):
				# If CLM equal zero count clear pixels
				num_mod = num_mod + 1
				if(MOD09_CLM [y_mod][x_mod] == 0 or MOD09_CLM [y_mod][x_mod] == 2 or  MOD09_CLM [y_mod][x_mod] == 4):
					num_mod_clr = num_mod_clr + 1
				else:
					num_mod_clm = num_mod_clm + 1
	return [num_mod, num_mod_clm]

class main():
  
	argDict = mapDict(argv, usage)

	if("-ty" in argDict and "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict):
		typeData    = str(argDict["-ty"])
		yearStart   = int(argDict["-ds"][0:4])
		yearFinish  = int(argDict["-de"][0:4])
		monthStart  = int(argDict["-ds"][4:6])
		monthFinish = int(argDict["-de"][4:6])
		outputDir = argDict["-fo"]
		inputDir  = argDict["-fi"]
		dateStart = int(argDict["-ds"])
		dateEnd   = int(argDict["-de"])
		if(typeData == 'DAY'):
			if("-ps" in argDict):
				processData = argDict["-ps"]
				if("-class" in argDict):
					classData = argDict["-class"]
				else:
					print('Error in data class, verify!')
					exit(usage)                    
				if(processData not in ['123','6S','NDVI','NDWI']):
					print('Error in data process, verify!')
					exit(usage)                    
			else:
				print('Error in data process, verify!')
				exit(usage)
		if(typeData == 'MOD_SAMP'):
			if("-ps" in argDict):
				processData = argDict["-ps"]
				if(processData not in ['MOD09', 'MYD09']):
					print('Verify process data options (MOD09 or MYD09)')
					sys.exit(usage)          
			else:
				print('Necessary process data options (MOD09 or MYD09)')
				sys.exit(usage) 
		if(typeData == 'TAB_CLM'):
			if("-ps" in argDict):
				perc = int(argDict["-ps"])
			else:
				print('Necessary percent data options [-ps perc]')
				sys.exit(usage)
		if(typeData in ['PCA', 'DAY_INT', 'PREC_INT', 'STD_INDEX']):
			if("-ps" in argDict):
				processData = argDict["-ps"]
			else:
				print('Necessary type data options [-ps processData]')         
	else:
		sys.exit(usage)

	if(typeData not in ['PCA', 'STD_INDEX', 'DAY', 'DAY_INT', 'PREC_INT' ,'MOD09_CLM', 'MOD09_CLM_MSG', 'TAB_CLM', 'MOD_SAMP', 'TAB_ALL', 'TAB_MOD13', 'DAY_AN', 'PAN', 'PRO', 'NDVI', 'VIEW', 'MEAN', 'HIST']):
		print('Error in data type, verify! You used: %s'%(typeData))
		exit(usage)

	# --------------------------------------
	# DAY generate daily data using range 1300 to 1330 
	if(typeData == 'DAY'):
		outputSubfolder = outputDir + processData + '_6S' + '\\'
		if not(os.path.isdir(outputSubfolder)):
			os.makedirs(outputSubfolder)

		ALL_HOURS = [1300, 1315, 1330, 1345, 1400, 1415, 1430, 1445, 1500, 1515, 1530, 1545, 1600, 1615, 1630, 1645, 1700]
		try:
			startClass = int(classData[0:4])
			finishClass= int(classData[6:10])
			id_start = ALL_HOURS.index(startClass+300)
			id_finish= ALL_HOURS.index(finishClass+300)+1
			HOURS = ALL_HOURS[id_start:id_finish]
			print(HOURS)
		except:
			print('Error in data class, verify!')
			exit(usage)

		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		OUT_MONTH = numpy.zeros((y_123, x_123))
		# Process in range of data
		date = dateStart
		while(date < dateEnd):
			print('. Create data to: %s'%(date))
			#for hour in HOURS:
			CLM  = numpy.zeros((len(HOURS), y_123, x_123))
			# OUT_NDVI = [IND B1 B2 HOUR]
			OUT_IND = numpy.zeros((4, y_123, x_123))
			# IND = [IND B1 B2]
			IND = numpy.zeros((len(HOURS), 3, y_123, x_123))
			cont_files = 0
			for i in range(len(HOURS)):
				hour = HOURS[i]
				print('hour: %s'%(hour))
				# Read bands RED and NIR
				inputSubDir = inputDir + '6S' +'\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
				filename = inputSubDir.replace('___', processData) + str(date) + str(hour) + '_' + processData + '.tif'
		  
				if os.path.isfile(filename) == 0:
					print('..Warning! Invalid data file: %s' %filename)

				else:
					B1 = readFileBand(filename, 1)
					B2 = readFileBand(filename, 2)
					B3 = readFileBand(filename, 3)

					if(len(B1) != 0 or len(B2) != 0):
					#if(len(NIR) != 0 or len(WIR) != 0):
						cont_files = cont_files + 1
						for x in range(x_123):
							for y in range(y_123):
								if(B1[y,x] > 0 and B2[y,x] > 0):
									IND[i][0][y,x] = B3[y,x]
									IND[i][1][y,x] = B1[y,x]
									IND[i][2][y,x] = B2[y,x]
					else:
						print('...Error process data file: %s'%(filename))
					# Read CLM data
					filename = inputSubDir.replace('6S', 'CLM') + str(date) + str(hour) + '_CLM.tif'
					if os.path.isfile(filename) == 0:
						print('..Error! Invalid data file: %s' %filename)
						sys.exit()
					else:
						CLM[i] = readFileBand(filename, 1)
			if(cont_files >= 1):
				for x in range(x_123):
					for y in range(y_123):
						for i in range(len(HOURS)):
							if(IND[i][1][y,x] > 0 and IND[i][2][y,x] > 0 and CLM[i][y,x] <= 1.0 and CLM[i][y,x] > 0):
								OUT_IND[0][y,x] = IND[i][0][y,x]
								OUT_IND[1][y,x] = IND[i][1][y,x]
								OUT_IND[2][y,x] = IND[i][2][y,x]
								OUT_IND[3][y,x] = HOURS[i]  
								break
			outputSubfolder2 = outputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
			if not(os.path.isdir(outputSubfolder2)):
				os.makedirs(outputSubfolder2)
			output_filename = outputSubfolder2 + str(date) + '_' + processData +'.tif'
			saveRASTERfile(output_filename, fileMask123, OUT_IND)
			new_date = incDate(date)
			# If add month save data
			if(str(date)[4:6] != str(new_date)[4:6]):
				output_filename = outputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'
				saveRASTERfile(output_filename, fileMask123, OUT_MONTH)         
				OUT_MONTH = numpy.zeros((y_123, x_123))
			else:
				OUT_MONTH = OUT_MONTH + 1*(OUT_IND[3]!=0)
			date = new_date
		output_filename = outputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'
		if(not os.path.isfile(output_filename)):
			saveRASTERfile(output_filename, fileMask123, OUT_MONTH) 

	# --------------------------------------
	# DAY generate combination image in data interval (12 days)
	if(typeData == 'DAY_INT'):
		fileMask123 = fileMaskPAN
		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		ds = None
		DATA = numpy.zeros((3, y_123, x_123))
		OUT_DATA = numpy.zeros((5, y_123, x_123))
		# Process in range of data

		interv = 12
		cont_interv = 0

		date = dateStart
		for doy in range(1, 366, interv):
			print('doy to verify: %d'% doy)
			if((calcDoy_date(date) >= doy) and (calcDoy_date(date) - doy < interv)):
				break

		while(date < dateEnd):
			cont_interv = cont_interv + 1
			if(cont_interv == 1):
				if(doy == 349):
					interv = 16
				if(doy > 360):
				   date = int(str(int(str(date)[0:4])+1) + '0101')
				   doy = 1 
				print('. Create data to: doy %d starting in %.0f'%(doy, date))
				outputSubfolder = outputDir + processData +'_INT' + '\\' + str(date)[0:4] + '\\' + str(date)[4:6]
				if not(os.path.isdir(outputSubfolder)):
					os.makedirs(outputSubfolder)
				output_filename = outputSubfolder + '\\' + str(date) + '_' + processData + '_INT' + '.tif'
			
			filename = inputDir + processData + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date) + '_'+ 'NDVI.tif' #processData +'.tif'
			if(os.path.isfile(filename)):
				print('Read file: %s'%(filename))
				DATA[0] = readFileBand(filename, 1)
				DATA[1] = readFileBand(filename, 1)#2)
				DATA[2] = readFileBand(filename, 1)#3)

				for  x in range(x_123):
					for y in range(y_123):
						if(DATA[0][y][x] == 0):
							continue
						else:
							if(OUT_DATA[0][y][x] == 0):
								OUT_DATA[0][y][x] = DATA[0][y][x]
								OUT_DATA[1][y][x] = DATA[1][y][x]
								OUT_DATA[2][y][x] = DATA[2][y][x]
								OUT_DATA[3][y][x] = doy
								OUT_DATA[4][y][x] = int(str(date)[4:8])
							else:
								if(DATA[0][y][x] > OUT_DATA[0][y][x]):
									OUT_DATA[0][y][x] = DATA[0][y][x]
									OUT_DATA[1][y][x] = DATA[1][y][x]
									OUT_DATA[2][y][x] = DATA[2][y][x]
									OUT_DATA[3][y][x] = doy
									OUT_DATA[4][y][x] = int(str(date)[4:8])
			else:
				print('Error to read: %s'%(filename))

			doy = doy + 1
			date = incDate(date)
			# If get data in interval dates
			if(cont_interv == interv):
				print('Saving file: %s'%(output_filename))
				saveRASTERfile(output_filename, fileMask123, OUT_DATA)
				OUT_DATA = numpy.zeros((5, y_123, x_123))
				cont_interv = 0 
				interv = 12
				#sys.exit()
	# --------------------------------------
	# DAY generate combination image in data interval (12 days)
	if(typeData == 'PREC_INT'):
		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		ds = None
		DATA = numpy.zeros((y_123, x_123))
		# Process in range of data

		interv = 12
		cont_interv = 0
		n_error = 0
		date = dateStart
		for doy in range(1, 366, interv):
			print('doy to verify: %d'% doy)
			if((calcDoy_date(date) >= doy) and (calcDoy_date(date) - doy < interv)):
				break

		while(date < dateEnd):
			cont_interv = cont_interv + 1
			if(cont_interv == 1):
				if(doy == 349):
					interv = 16
				if(doy > 360):
				   date = int(str(int(str(date)[0:4])+1) + '0101')
				   doy = 1 
				print('. Create data to: doy %d starting in %.0f'%(doy, date))
				outputSubfolder = outputDir + processData[0:4] +'_INT' + '\\'  + str(date)[0:4] + '\\' + str(date)[4:6]
				if not(os.path.isdir(outputSubfolder)):
					os.makedirs(outputSubfolder)
				output_filename = outputSubfolder + '\\' + str(date) + '_' + processData[0:4] + '_INT' + '.tif'
			
			filename = inputDir + processData + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date) + '_'+ 'PREC.tif' #processData +'.tif'
			if(os.path.isfile(filename)):
				print('Read file: %s'%(filename))
				DATA_DAY = readFileBand(filename, 1)
				if(numpy.sum(DATA_DAY)>0 and not numpy.any(numpy.isnan(DATA_DAY))):
					DATA = DATA + DATA_DAY
				#DATA[1] = readFileBand(filename, 1)#2)
				#DATA[2] = readFileBand(filename, 1)#3)
			else:
				print('Error to read: %s'%(filename))
				n_error = n_error + 1

			doy = doy + 1
			date = incDate(date)
			# If get data in interval dates
			if(cont_interv == interv):
				if(n_error > 0.8*interv):
					print('Error to generate: %s'%(output_filename))
				else:
					print('Saving file: %s'%(output_filename))
					saveRASTERfile(output_filename, fileMask123, DATA)
				DATA = numpy.zeros((y_123, x_123))
				cont_interv = 0
				n_error = 0
				interv = 12
				#sys.exit()

	# --------------------------------------
	# Process Standard Index  
	if(typeData == 'STD_INDEX'):
		date = dateStart

		int_doy = 12
		year_doys  = 30
		num_file = round(calcDays(date, dateEnd)/int_doy)
		'''
		if(processData in ['NDVI', 'NDWI']):
			inputfile = inputDir + 'INTERP\\' + processData + '\\year\\month\\yearmonthday_interp_' + processData + '.tif'
		else:
			inputfile = 't:\\outros\\' + processData + '\\year\\month\\yearmonthday_'+ processData + '.tif'
		'''
		inputfile = inputDir + 'BASIN\\' + processData + '\\year\\month\\yearmonthday_basin_' + processData + '.tif'
		outputfile = outputDir + 'STAND\\' 'std_index_basin_'+ processData + '.tif'
		# Create variable to process   
		MASK = readFileBand(fileMask123,1)
		ds = gdal.Open(fileMask123)
		y_band = ds.RasterYSize
		x_band = ds.RasterXSize
		ds = None
		DATA = numpy.zeros((num_file, y_band, x_band))
		DATA_OUT = numpy.zeros((num_file, y_band, x_band))
		n_file = 0
		# Read data num_file 
		print('... Read data files')  
		while(date < dateEnd):
			filename = inputfile.replace('year', str(date)[0:4])
			filename = filename.replace('month',str(date)[4:6])
			filename = filename.replace('day',  str(date)[6:8])
			if(os.path.isfile(filename) == False):
				print('error to read: %s'%(filename))
				#sys.exit()
			else:
				print('reading: %s'%(filename))

				DATA[n_file] = readFileBand(filename,1)
				n_file = n_file + 1				
			new_date = inc_date_doy(date, int_doy)

			#print(new_date)
			#print(date)
			if(int(str(new_date)[0:4]) > int(str(date)[0:4]) or '1227' in str(new_date)):
				date = int(str(int(str(date)[0:4])+1) + '0101')
			else:
				date = new_date
		# Process Standard Index
		print('... Process Standard Index')
		if(year_doys > n_file):
			doys_to_out = n_file
		else:
			doys_to_out = year_doys
		for x in range(x_band):
			for y in range(y_band):
				if(MASK[y][x] != 0):
					# Calc MEAN and STD
					MEAN = numpy.zeros((doys_to_out))
					STD  = numpy.zeros((doys_to_out))
					for n_samp in range(doys_to_out):
						#print('.. n_sample: %d'%(n_samp))
						samples = []
						for n_doy in range(n_samp, n_file, year_doys):
							samples.append(DATA[n_doy][y][x])
						mean_samples = numpy.mean(samples)
						std_samples  = numpy.std(samples)
						samples.append(mean_samples)
						samples.append(std_samples)
						#print(samples)
						for n_doy in range(n_samp, n_file, year_doys):
							if(abs(std_samples) < 10**-5):
								DATA_OUT[n_doy][y][x] = (DATA[n_doy][y][x] - mean_samples)
								samples.append(DATA_OUT[n_doy][y][x])	
								#print(samples)
							else:
								DATA_OUT[n_doy][y][x] = (DATA[n_doy][y][x] - mean_samples)/std_samples
		# Save standard index file
		print('...Save standard file')
		yearStart = int(str(dateStart)[0:4])
		yearEnd   = int(str(dateEnd)[0:4]) + 1
		for n_year, year in zip(range(yearEnd - yearStart), numpy.arange(yearStart, yearEnd+1,1)):
			for n_doy in range(year_doys):
				if(n_year*year_doys + n_doy + 1 > n_file):
					continue
				outputSubfolder = outputDir + 'STAND\\' + processData + '\\' + str(year) + '\\'
				if(not os.path.isdir(outputSubfolder)): 
					os.makedirs(outputSubfolder)
				outputfile = outputSubfolder + str('%d%03d'%(year, int_doy*n_doy+1)) + '_std_index_basin_'+ processData + '.tif'
				print('saving: %s'%(outputfile))
				saveRASTERfile(outputfile, fileMask123, DATA_OUT[n_year*year_doys + n_doy])
		# Process corrcoef (angular coef) of data
		outputfile = outputDir + 'STAND\\' 'std_index_basin_'+ processData + '.tif'
		COEF = numpy.zeros((2, y_band, x_band))
		print('... Calc angular coef of data')
		for x in range(x_band):
			for y in range(y_band):
				temp_series = numpy.zeros((n_file))
				if(MASK[y][x] != 0):
					for n_doys in range(n_file):
						temp_series[n_doy] = DATA_OUT[n_doy][y][x]
					if(numpy.count_nonzero(temp_series)>0):
						#linreg(range(len(x)),x)
						COEF[0][y][x], COEF[1][y][x] = numpy.polyfit(range(len(temp_series)), temp_series,1) #linreg(range(len(temp_series)), temp_series)
						#print('y [%d] x [%d]'%(y, x))
						#print(temp_series)
						#print(range(n_file))
						#print(COEF[y][x])
						#os.system("pause")
						

		saveRASTERfile(outputfile.replace('.tif','_coef.tif'), fileMask123, COEF)
		# Save output file
		print('... Save output data file')
		num_band = 1
		for pos_band in range(n_file):
			if(pos_band == 0):
				# Save georeference file
				#print('------------------')
				#print(num_band)
				#print(outputDir + file)
				ds = gdal.Open(fileMask123)
				band = ds.GetRasterBand(1)
				driver = gdal.GetDriverByName("GTiff")
				dsOut  = driver.Create(outputfile, ds.RasterXSize, ds.RasterYSize, n_file, band.DataType)
				CopyDatasetInfo(ds,dsOut)
				bandOut = dsOut.GetRasterBand(num_band)
				BandWriteArray(bandOut, DATA_OUT[num_band-1])
			else:
				num_band = num_band + 1
				#print('------------------')
				#print(num_band)
				bandOut = dsOut.GetRasterBand(num_band)
				BandWriteArray(bandOut, DATA_OUT[num_band-1]) 
				#print('---- go out to test...')


	# --------------------------------------
	# Process PCA in MSG and MOD13 data  
	if(typeData == 'PCA'):
		date = dateStart
				
		if(processData == 'MOD13'):
			inputSubfolder = inputDir + 'MOD13_3KM' + '\\'
			int_doy = 16
			year_doys  = 23
			num_file = round(calcDays(date, dateEnd)/int_doy)
			outputFile = inputDir + 'PCA\\mod_pca_'+ processData + '_' + str(date) + '_' + str(dateEnd) + '.tif'
			
		else:
			inputSubfolder = inputDir + processData + '_6S\\'
			int_doy = 12
			year_doys  = 30
			num_file = round(calcDays(date, dateEnd)/int_doy)
			year_start = int(str(dateStart)[0:4])
			year_end   = int(str(dateEnd)[0:4])
			outputFile = inputDir + 'PCA\\msg_pca_'+ processData + '_' + str(date) + '_' + str(dateEnd) + '.tif'
		# Create variable to process   
		MASK = readFileBand(fileMask123,1)
		ds = gdal.Open(fileMask123)
		y_band = ds.RasterYSize
		x_band = ds.RasterXSize
		ds = None
		DATA = numpy.zeros((num_file, y_band, x_band))
		DATA_OUT = numpy.zeros((year_doys, y_band, x_band))
		#ONES = 7*numpy.ones((y_band, x_band))
		
		name_file = []
		n_file = 0
		# Read data num_file 
		print('... Read data files')  
		while(date < dateEnd):
			filename = inputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date)[0:4] + str(date)[4:8] + '_' + processData + '.tif'
			if(os.path.isfile(filename) == False):
				print('error to read: %s'%(filename))
				sys.exit()
			else:
				print('reading: %s'%(filename))

			SAMP = readFileBand(filename,1)
			if(processData == 'MOD13'):
				CLM = readFileBand(filename,6)
			name_file.append(str(date)[0:4] + str(date)[4:8] + '_' + processData + '.tif')
			#print('. filter nan data')
			#print(n_file)
			DATA_FILT = SAMP
			if(processData == 'NDVI'):
				DATA_FILT[DATA_FILT <= 0] = numpy.nan
			else:
				DATA_FILT[DATA_FILT < -1] = numpy.nan
				DATA_FILT[DATA_FILT >  1] = numpy.nan
				DATA_FILT[DATA_FILT == 0] = numpy.nan
				if(processData == 'MOD13'):
					DATA_FILT[CLM != 0] = numpy.nan
			DATA[n_file] = DATA_FILT
			n_file = n_file + 1
			print('n_file: %d'%(n_file))
			new_date = inc_date_doy(date, int_doy)

			print(new_date)
			print(date)
			if(int(str(new_date)[0:4]) > int(str(date)[0:4]) or '1227' in str(new_date)):
				date = int(str(int(str(date)[0:4])+1) + '0101')
			else:
				date = new_date
		input("Press Enter to continue...")
		# Filter nan data and process sav golae     
		print('... Apply interpolate and savgol_filter') 
		for y in range(y_band):
			for x in range(x_band):
				if(MASK[y][x] <= 0):
					continue
				time = numpy.zeros(n_file)
				for n_samp in range(n_file):
					time[n_samp] = DATA[n_samp][y][x]
				#print('data[%d][%d] in n_file %d:'%(y, x, n_file))
				#print(time)
				nans, lamb = numpy.isnan(time), lambda z: z.nonzero()[0]

				if(numpy.sum(~nans) == 0):
					time_new = numpy.nan_to_num(time)
				else:
					if(numpy.sum(nans) > 0):
						time[nans] = numpy.interp(lamb(nans), lamb(~nans), time[~nans])
					time_new   = savgol_filter(time, 5, 2, mode ='nearest')
				#print(time_new)
				for n_samp in range(n_file):
					DATA[n_samp][y][x] = time_new[n_samp]
				#sys.exit()
		DATA = numpy.nan_to_num(DATA)
		# Save interpolated data files
		'''
		if(processData != 'MOD13'):
			n_doy = 1
			band  = 1
			date = dateStart
			print('Save interpolated files')
			while(date < dateEnd):
				interp_outputsubfolder = outputDir + 'INTERP\\' + processData + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
				if(not os.path.isdir(interp_outputsubfolder)):
					os.makedirs(interp_outputsubfolder)
				interp_outputfile = interp_outputsubfolder + str(date)[0:4] + str(date)[4:8] + '_interp_' + processData + '.tif'
				print('Saving interpolated file: %s'%(interp_outputfile))
				saveRASTERfile(interp_outputfile, filename, DATA[band])
				band = band + 1
				new_date = inc_date_doy(date, int_doy)
				print('date: %s new_date: %s' %(date,new_date))
				if(int(str(new_date)[0:4]) > int(str(date)[0:4]) or '1227' in str(new_date)):
					date = int(str(int(str(date)[0:4])+1) + '0101')
				else:
					date = new_date
		'''
		# Calc mean over data generate in same days
		print('... Calc means to years')
		if(year_doys > n_file):
			doys_to_out = n_file
		else:
			doys_to_out = year_doys
		print('doys_to_out: %d year_doys: %d n_file %d'%(doys_to_out, year_doys, n_file))
		for n_samp in range(doys_to_out):
			print('.. n_samp: %d'%(n_samp))
			n_div = 0
			for n_doy in range(n_samp, n_file, year_doys):
				#print(n_doy)
				print('name_file[%d]: %s'%(n_doy, name_file[n_doy]))
				DATA_OUT[n_samp] = DATA_OUT[n_samp] + DATA[n_doy]
				n_div = n_div + 1
			if(n_div > 0):
				DATA_OUT[n_samp] = DATA_OUT[n_samp]/n_div
			else:
				print('ERRROOOOOOORRR')
				print(n_samp)
				print(n_file)
				print(n_doy)
				print(n_div)
				sys.exit()


		# Save output file
		print('... Save output data file')
		num_band = 1
		for pos_band in range(doys_to_out):
			if(pos_band == 0):
				# Save georeference file
				#print('------------------')
				#print(num_band)
				#print(outputDir + file)
				ds = gdal.Open(fileMask123)
				band = ds.GetRasterBand(1)
				driver = gdal.GetDriverByName("GTiff")
				dsOut  = driver.Create(outputFile, ds.RasterXSize, ds.RasterYSize, year_doys, band.DataType)
				CopyDatasetInfo(ds,dsOut)
				bandOut = dsOut.GetRasterBand(num_band)
				BandWriteArray(bandOut, DATA_OUT[num_band-1])
			else:
				num_band = num_band + 1
				#print('------------------')
				#print(num_band)
				bandOut = dsOut.GetRasterBand(num_band)
				BandWriteArray(bandOut, DATA_OUT[num_band-1]) 
				#print('---- go out to test...')

	# --------------------------------------
	# DAY generate daily data using range 1300 to 1330 
	if(typeData == 'MOD_SAMP'):
		first_file = True
		# Process in range of data
		date = dateStart
		inputSubfolder = inputDir + 'MODIS' + '\\'
		n_month = 0
		while(date < dateEnd):
			filename = inputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date)[0:4] + str(date)[4:8] + '_' + processData + '.tif'
			if(os.path.isfile(filename) == False):
				print('...error to read: %s'%(filename))
			else:
				print('...reading: %s'%(filename))
				if(first_file == True):
					first_file = False
					ds = gdal.Open(filename)
					y_123 = ds.RasterYSize
					x_123 = ds.RasterXSize
					ds = None
					OUT_MONTH = numpy.zeros((y_123, x_123))
					OUT_YEAR = numpy.zeros((y_123, x_123))
					ONES = 7*numpy.ones((y_123, x_123))
					fileMask123 = filename
				DATA = readFileBand(filename, 1)
				CLM  = readFileBand(filename, 4)
				FILT = numpy.logical_not(numpy.bitwise_and(CLM.astype(numpy.int32), ONES.astype(numpy.int32)))
				SAMP = numpy.logical_and(DATA > 0, FILT > 0)
			new_date = incDate(date)

			# If add month save data
			if(str(date)[4:6] != str(new_date)[4:6]):
				outputSubfolder = outputDir + '\\' + str(date)[0:4] + '\\'
				output_filename = outputSubfolder + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'
				if(not os.path.isfile(output_filename)):
					if(not os.path.isdir(outputSubfolder)):
						os.makedirs(outputSubfolder)
					saveRASTERfile(output_filename, fileMask123, OUT_MONTH)    
					print('save file : %s' %(outputSubfolder + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'))
				else:
					print('save exist : %s' %(outputSubfolder + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'))
				OUT_MONTH = numpy.zeros((y_123, x_123))
				OUT_YEAR = OUT_YEAR + OUT_MONTH
				n_month = n_month + 1

				'''
				if(str(date)[:4] != str(new_date)[:4]):
					output_filename = outputSubfolder + str(date)[0:4] + '_MEAN.tif'
					saveRASTERfile(output_filename, fileMask123, OUT_YEAR)
					OUT_YEAR = numpy.zeros((y_123, x_123))
				'''                
			else:
				OUT_MONTH = OUT_MONTH + 1*SAMP
				#output_filename = outputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[0:4] + str(date)[4:8] + '_SAMP.tif'
				#saveRASTERfile(output_filename, fileMask123, SAMP) 
			date = new_date
		output_filename = outputSubfolder + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'
		if(not os.path.isfile(output_filename)):
			saveRASTERfile(output_filename, fileMask123, OUT_MONTH)
			OUT_YEAR = OUT_YEAR + OUT_MONTH
		output_filename = outputSubfolder + str(date)[0:4] + '_MEAN.tif'
		if(not os.path.isfile(output_filename)):
			saveRASTERfile(output_filename, fileMask123, OUT_YEAR)
		print('n_month: %d'%(n_month))
		output_filename = outputSubfolder + str(date)[0:4] + '_MEAN_div.tif'
		if(not os.path.isfile(output_filename)):
			saveRASTERfile(output_filename, fileMask123, OUT_YEAR/n_month)

	# -------------------------------------- DAY ANALYSES DATA
	if(typeData == 'DAY_AN'):
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()
		NAMES = ['MEAN', 'STD', 'MEDIAN', 'MAX', 'MIN', 'NZERO']
		SUB_NAMES = [['DAY_123_1000as1100', 'SAMP'],['DAY_123_1000as1200', 'SAMP'],['DAY_123_1000as1400', 'SAMP'],['DAY_123_1030as1130', 'SAMP'],['DAY_123_1030as1230', 'SAMP'],['DAY_123_1100as1200','SAMP'],['DAY_123_1100as1300','SAMP'],['DAY_123_1130as1230','SAMP'],['DAY_123_1130as1330','SAMP'],['DAY_123_1200as1300','SAMP'],['DAY_123_1200as1400','SAMP'],['DAY_123_1230as1330','SAMP'],['DAY_123_1300as1400','SAMP']]
		n_val = len(NAMES)*len(SUB_NAMES) + 1
		header = 'MONTH'
		for x in range(len(SUB_NAMES)):
			for name in NAMES:
				header = header + ' ' + name + SUB_NAMES[x][0][7:]
		print(header)
		print(num_days)
		print(n_val)
		DATA_SAMP = datafile('', outputDir, 'samp_day_analyses.csv', header, num_days, n_val)
		
		numpy.set_printoptions(precision=4)
		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		ds = None
	 
		
		for year in range(yearStart, yearFinish+1):
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				DATA_SAMP.addAtrib(month)
				for k in range(len(SUB_NAMES)):
					file = inputDir + SUB_NAMES[k][0] + '\\' + str(year) + '\\' + str(year) + str_month + '_' + SUB_NAMES[k][1] + '.tif' 
					print('Reading: ' + file)
					# All files have month's samples number 
					DATA = readFileBand(file, 1)
					if (DATA != []):
						LIST = []
						for x in range(x_123):
							for y in range(y_123):
								if(DATA[y,x] > 0):
									LIST.append(DATA[y,x])
						n_zero = 24072 - len(LIST)
						for x in range(n_zero):
							LIST.append(0)
						# MEAN STD MEDIAN MAX MIN
						DATA_SAMP.addAtrib(numpy.mean(LIST))
						DATA_SAMP.addAtrib(numpy.std(LIST))
						DATA_SAMP.addAtrib(numpy.median(LIST))
						DATA_SAMP.addAtrib(numpy.max(LIST))
						DATA_SAMP.addAtrib(numpy.min(LIST))
						DATA_SAMP.addAtrib(n_zero)
				DATA_SAMP.addSample()

		DATA_SAMP.save()

	# -------------------------------------- CLM Analyses           
	if(typeData == 'MOD09_CLM'):

		if(yearStart > yearFinish):
			inc = 1
		else:
			inc = 1
		
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()
		for year in range(yearStart, yearFinish+inc):
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				for day in range (1, 32):
					doy = calcDoy(year,month,day)
					str_day = str('%02d'%day)
					folder_out = inputDir + 'MOD09_CLM' + '\\' + str(year) + '\\' + str_month + '\\'
					if not (os.path.isdir(folder_out)):
						os.makedirs(folder_out)
					file_out = folder_out + str(year) + str_month + str_day + '_MOD09_CLM.tif' 
					if(os.path.isfile(file_out)):
						continue
					print(' Process day: %s%s%s'%(str(year),str_month, str_day))
					'''
					file = inputDir + CLM + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '_CLM.tif'                  
					if(os.path.isfile(file)):
						MSG = readFileBand(file, 1)
						print('... Reading  file: %s'%(file))
					else:
						print('... Error to read: %s'%(file))
						continue
					'''
					file = inputDir + 'MODIS' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '_MOD09.tif' 
					if(os.path.isfile(file)):
						MOD09_STATE_1KM = readFileBand(file, 4) 
						MOD09_LIM = readFileBand(file, 2) 
						print('... Reading  file: %s'%(file))
					else:
						print('... Error to read: %s'%(file))
						continue
					ds = gdal.Open(file)
					y_123 = ds.RasterYSize
					x_123 = ds.RasterXSize
					ONES = 7*numpy.ones((y_123, x_123))
					ds = None
					MOD09_CLM = numpy.bitwise_and(MOD09_STATE_1KM.astype(numpy.int16), ONES.astype(numpy.int16))
					MOD09_LIM = numpy.logical_not(numpy.logical_not(MOD09_LIM>-1))
					saveRASTERfile(file_out, file, [MOD09_CLM, MOD09_LIM])
					print(' Save file: %s'%(file_out))
					#print('test...')
					#sys.exit()

	# -------------------------------------- CLM Analyses           
	if(typeData == 'MOD09_CLM_MSG'):

		if(yearStart > yearFinish):
			inc = 1
		else:
			inc = 1
		
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()
		header = 'YEAR MONTH DOY N_FT N_FT_25 N_FT_50 N_FT_75 N_FT_100 N_TT N_TT_25 N_TT_50 N_TT_75 N_TT_100'
		DATA = datafile('', outputDir, typeData + '_1315_to_1330_year2015' + '_.csv', header, num_days, 14)
		for year in range(yearStart, yearFinish+inc):
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				MSG_T_MOD_T = [0, 0, 0, 0]
				MSG_F_MOD_T = [0, 0, 0, 0]
				for day in range (1, 32):
					doy = calcDoy(year,month,day)
					str_day = str('%02d'%day)
					print(' Process day: %s%s%s'%(str(year),str_month, str_day))

					file = inputDir + 'CLM' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '1315_CLM.tif'                  
					if(os.path.isfile(file)):
						MSG_CLM1 = readFileBand(file, 1)
						ds = gdal.Open(file)
						y_MSG = ds.RasterYSize
						x_MSG = ds.RasterXSize
						print(' reading  file: %s'%(file))
					else:
						print(' error to read: %s'%(file))
						continue

					file = inputDir + 'CLM' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '1330_CLM.tif'                  
					if(os.path.isfile(file)):
						MSG_CLM2 = readFileBand(file, 1)
						print(' reading  file: %s'%(file))
					else:
						print(' error to read: %s'%(file))
						continue
					
					file = inputDir + 'MOD09_CLM' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '_MOD09_CLM.tif' 
					if(os.path.isfile(file)):
						MOD09_CLM = readFileBand(file, 1) 
						MOD09_LIM = readFileBand(file, 2) 
						print(' reading  file: %s'%(file))
					else:
						print(' error to read: %s'%(file))
						continue
					for y in range(y_MSG):
						for x in range(x_MSG):
							if((MSG_CLM1[y][x] < 1 or MSG_CLM1[y][x] > 2) and (MSG_CLM2[y][x] < 1 or MSG_CLM2[y][x] > 2)):
								continue

							[num_mod, num_mod_clm] = process_mod_clm(x, y, MOD09_CLM, MOD09_LIM)
							#wait()
							# If MSG CLM > 1.05 then pixel cloud
							if(num_mod < 50):
								continue
							printa = 0
							clm_perc = int(round(100*(num_mod_clm/num_mod))/25)
							if(clm_perc == 4):
								clm_perc = 3

							if(MSG_CLM1[y][x] < 1.05 and MSG_CLM2[y][x] < 1.05):
								MSG_F_MOD_T[clm_perc] = MSG_F_MOD_T[clm_perc] + 1
								#print('MSG_F_MOD_T')
								#print(MSG_F_MOD_T) 
								#printa = 1

							if(MSG_CLM1[y][x] >= 1.05 and MSG_CLM2[y][x] >= 1.05):
								MSG_T_MOD_T[clm_perc] = MSG_T_MOD_T[clm_perc] + 1
								#print('MSG_T_MOD_T')
								#print(MSG_T_MOD_T) 
								#printa = 1
							if(printa == 1):
								print('-------------------------------------------')
								print('MSG[%d][%d] in %.4f %.4f has 1315: %.2f 1330: %.2f MOD[%d][%d] in %.4f %.4f'%(y,x, r3km_py+r3km_psy*y, r3km_px + r3km_psx*x, MSG_CLM1[y][x], MSG_CLM2[y][x], start_py, start_px, r3km_py + r3km_psy*y + r500m_psy/2, r3km_px + r3km_psx*x + r500m_psx/2))
								print('MOD[%d][%d] num_mod_clm: %d num_mod_clr: %d num_mod: %d' %(y_mod,x_mod, num_mod_clm, num_mod_clr, num_mod))
								wait()

							
				DATA.addAtrib([year, month, doy])
				DATA.addAtrib(numpy.sum(MSG_F_MOD_T))
				DATA.addAtrib(MSG_F_MOD_T/numpy.sum(MSG_F_MOD_T))
				DATA.addAtrib(numpy.sum(MSG_T_MOD_T))
				DATA.addAtrib(MSG_T_MOD_T/numpy.sum(MSG_T_MOD_T))
				line = DATA.getLine()
				print(' %d %2d %3d MSG_F_MOD_T: %d %.2f %.2f %.2f %.2f MSG_T_MOD_T: %d %.2f %.2f %.2f %.2f'%(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12]))
				DATA.addSample()
				DATA.save()
				#wait()
				#print('test...')
				#sys.exit()
	# -------------------------------------- Generate Table to process CLM           
	if(typeData == 'TAB_CLM'):

		if(yearStart > yearFinish):
			inc = 1
		else:
			inc = 1
		
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()

		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		ds = None
		# MSG BANDS (1,2,3,4,5,9)
		BANDS = numpy.zeros((6, y_123, x_123))
		header    = 'DATE LAT LON BAND_1 BAND_2 BAND_3 BAND_4 BAND_5, BAND_9, CLM MSG_CLM'
		DATA = datafile('', outputDir, typeData + '_perc_1_to_' + str(perc) +'_date_' + str(dateStart) + '_' + str(dateEnd) + '.csv', header, x_123*y_123*num_days/10, 11)

		for year in range(yearStart, yearFinish+inc):
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				for day in range (1, 32):
					# Select data using random between 1/10
					if(random.randint(1, perc) != 1):
						continue
					doy = calcDoy(year,month,day)
					str_day = str('%02d'%day)
					date = int(str(year) + str_month + str_day)

					print(' Process day: %s%s%s'%(str(year),str_month, str_day))

					file = inputDir + 'BANDS' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '1330_123.tif'                  
					if(os.path.isfile(file)):
						for band in range(1,7):
							BANDS[band-1] = readFileBand(file, band)
						print(' reading  file: %s'%(file))
					else:
						print(' error to read: %s'%(file))
						continue
					
					file = inputDir + 'MOD09_CLM' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '_MOD09_CLM.tif' 
					if(os.path.isfile(file)):
						MOD09_CLM = readFileBand(file, 1) 
						MOD09_LIM = readFileBand(file, 2) 
						print(' reading  file: %s'%(file))
					else:
						print(' error to read: %s'%(file))
						continue
					file = inputDir + 'CLM' + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '1330_CLM.tif' 
					if(os.path.isfile(file)):
						MSG_CLM = readFileBand(file, 1) 
						print(' reading  file: %s'%(file))
					else:
						print(' error to read: %s'%(file))
						continue

					for y in range(y_123):
						for x in range(x_123):
							if(BANDS[0][y][x] == 0):
								continue
							[num_mod, num_mod_clm] = process_mod_clm(x, y, MOD09_CLM, MOD09_LIM)
							if(num_mod < 60):
								continue
							
							clm_perc = int(round(100*(num_mod_clm/num_mod))/25)
							if(clm_perc == 4):
								clm_perc = 3
							[lat, lon] = conv_3km_to('3km', y, x)
							DATA.addAtrib([date, lat, lon])
							for band in range(0,6):
								DATA.addAtrib(BANDS[band][y][x])
							DATA.addAtrib(clm_perc)
							DATA.addAtrib(MSG_CLM[y][x])
							DATA.addSample()
					DATA.save()
					#print('............ test now!!')   
					#sys.exit()


	# -------------------------------------- TAB ANALYSES DATA
	if(typeData == 'TAB_ALL' or typeData == 'TAB_MOD13'):
		'''
		Points to verify data:
		1.  [49, 144] [815,  2402] Agricultura Anual  (Terraclass 1)
		2.  [15, 179] [249,  2986]
		3.  [33, 156] [548,  2602]
		4.  [50, 156] [831,  2602]
		5.  [12, 144] [198,  2402]
		6.  [59, 146] [982,  2436]
		7.  [58, 140] [965,  2336]
		8.  [69, 144] [1148, 2402]
		9.  [79, 134] [1315, 2236]
		10. [64, 130] [1065, 2168]

		11. [81,  35] [1348, 586 ] Pastagem (Terraclass 11)
		12. [83 ,  8] [1380, 137 ] 
		13. [112, 45] [1864, 753 ]
		14. [154, 39] [2564, 653 ]
		15. [84, 50 ] [1398, 836 ] 
		16. [81, 123] [1348, 2052]
		17. [89, 130] [1482, 2170]
		18. [109,156] [1815, 2602]
		19. [81, 15 ] [1348, 254 ] 
		20. [186, 39] [3097, 654 ]

		21. [ 8, 162] [131,  2703] Natural Campestre  (Terraclass 5)
		22. [173, 96] [2881, 1603]
		23. [164, 33] [2731, 553 ]
		24. [87, 35 ] [1447, 587 ] Natural Florestal  (Terraclass 6) 
		25. [66, 54 ] [1098, 903 ]
		26. [60, 86 ] [998,  1436]
		27. [182, 33] [3030,  553] Natural Savânica   (Terraclass 12)
		28. [1, 148 ] [15,   2469]
		29. [130, 56] [2164,  937] Natural (Terraclass 13)
		30. [137, 72] [2281, 1203] Natural Savânica   (Terraclass 12)
		'''

		CULTURE = [[49,144],[15,179],[33,156],[50,156],[12,144],[59,146],[58,140],[69,144],[79,134],[64,130]]
		PASTURE = [[81,35],[83,8],[112,45],[154,39],[84,50],[81,123],[89,130],[109,156], [81,15], [186,39]]
		NATURE =  [[8,162],[173,96],[164,33],[87,35],[66,54],[60,86],[182,33],[1,148],[130,56],[137,72]]
		
		#CULTURE_MOD = [[815, 2402],[249,2986],[548,2602],[831,2602],[198,2402],[982,2436],[965,2336],[1148,2402],[1315,2236],[1065,2168]]
		#PASTURE_MOD = [[1348,586],[1380,137],[1864,753],[2564,653],[1398,836],[1348,2052],[1348,2052],[1482,2170],[1815,2602],[1348,254],[3097,654]]
		#NATURE_MOD =  [[131,  2703],[2881, 1603],[2731, 553 ],[1447, 587 ],[1098, 903 ],[998,  1436],[3030,  553], [15,   2469], [2164,  937],[2281, 1203]]
		
		#r3km_y = CULTURE[k][0]
		#r3km_x = CULTURE[k][1]
		#[r250m_y, r250m_x] = conv_3km_to('250m', r3km_y, r3km_x)
		#[r500m_y, r500m_x] = conv_3km_to('500m', r3km_y, r3km_x)
	

		if(yearStart > yearFinish):
			inc = 1
		else:
			inc = 1
		
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()
		
		numpy.set_printoptions(precision=4)
		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		ds = None
		# Data to compare BRDF types
		#n_val  = 21
		#header    = 'TYPE SAMP YEAR MONTH DAY DAY_123_NDVI DAY_123_RED DAY_123_NIR DAY_6S_NDVI DAY_6S_RED DAY_6S_NIR BRDF_123_NDVI BRDF_123_RED BRDF_123_NIR BRDF_6S_NDVI BRDF_6S_RED BRDF_6S_NIR MOD13_3KM_NDVI MOD13_3KM_RED MOD13_3KM_NIR MOD13_250M_NDVI'
		#SUB_NAMES = [['DAY_123', 'DAY'],['DAY_6S', 'DAY'],['BRDF_123', 'BRDF'],['BRDF_6S', 'BRDF'],['MOD13_3KM','MOD13'], ['MOD13_250M','MOD13']]

		#Data to compare DAYS

		#SUB_NAMES = [['DAY_123_1000as1100', 'DAY'],['DAY_123_1000as1200', 'DAY'],['DAY_123_1000as1400', 'DAY'],['DAY_123_1030as1130', 'DAY'],['DAY_123_1030as1230', 'DAY'],['DAY_123_1100as1200','DAY'],['DAY_123_1100as1300','DAY'],['DAY_123_1130as1230','DAY'],['DAY_123_1130as1330','DAY'],['DAY_123_1200as1300','DAY'],['DAY_123_1200as1400','DAY'],['DAY_123_1230as1330','DAY'],['DAY_123_1300as1400','DAY'],['MODIS', 'MOD09'],['MOD13_3KM','MOD13'], ['MOD13_250M','MOD13']]
		#SUB_NAMES = [['DAY_123_1000as1200', 'DAY'], ['DAY_6S_1000as1200', 'DAY'], ['BRDF_123', 'BRDF'], ['BRDF_6S', 'BRDF'], ['MODIS', 'MOD09'],['MOD13_3KM','MOD13'], ['MOD13_250M','MOD13']]
		SUB_NAMES = [['DAY_123_1000as1200', 'DAY'],['DAY_6S_1000as1200', 'DAY'], ['BRDF_123', 'BRDF'],['BRDF_6S', 'BRDF'],['MODIS', 'MOD09'], ['MOD13_3KM','MOD13'], ['MOD13_250M','MOD13']]
		NAMES = ['NDVI', 'RED', 'NIR']
		n_val = len(NAMES)*len(SUB_NAMES) + 5
		header    = 'TYPE SAMP YEAR MONTH DAY'
		
		for x in range(len(SUB_NAMES)):
			'''
			for name in NAMES:
				if(SUB_NAMES[x][0] != 'MOD13_250M' and SUB_NAMES[x][0] != 'MOD13_3KM' and SUB_NAMES[x][0] != 'MODIS'):
					header = header + ' ' + name + SUB_NAMES[x][0][7:]
				else:
					header = header + ' ' + name + '_'+ SUB_NAMES[x][0]
					if(SUB_NAMES[x][0] == 'MOD13_250M'):
						break
			'''
			for name in NAMES:
				header = header + ' ' + name + '_'+ SUB_NAMES[x][0]
		print(header)

		DATA_CULT = datafile('CULTURE', outputDir, typeData + '_.csv', header, len(CULTURE)*num_days*2, n_val)
		DATA_PAST = datafile('PASTURE', outputDir, typeData + '_.csv', header, len(PASTURE)*num_days*2, n_val)
		DATA_NATU = datafile('NATURE', outputDir,  typeData + '_.csv',  header, len(NATURE)*num_days*2, n_val)
		# [FOLDER _NAME]
 
		if(typeData == 'TAB_ALL'):
			row_past = numpy.zeros((n_val))
			row_cult = numpy.zeros((n_val))
			row_natu = numpy.zeros((n_val))
			ant_past = numpy.zeros((n_val))
			ant_cult = numpy.zeros((n_val))
			ant_natu = numpy.zeros((n_val))
		else:
			row_past = numpy.zeros((len(CULTURE), n_val))
			row_cult = numpy.zeros((len(PASTURE), n_val))
			row_natu = numpy.zeros((len(NATURE), n_val))
		for year in range(yearStart, yearFinish+inc):
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				for day in range (1, 32):
					doy = calcDoy(year,month,day)
					MOD13_250M = []
					MOD09_NDVI = []
					str_day = str('%02d'%day)
					print(' Process day: %s%s%s'%(str(year),str_month, str_day))
					# DATA_TMP = [DAY_123 DAY_6S BRDF_123 BRDF_6S MOD13]
					DATA_TMP = numpy.zeros((3*len(SUB_NAMES), y_123, x_123))
					for k in range(len(SUB_NAMES)):
						file = inputDir + SUB_NAMES[k][0] + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '_' + SUB_NAMES[k][1] + '.tif' 
						print('data: %s\tfile: %s'%(SUB_NAMES[k][0], file))

						# All files have NDVI in fist band
						#print('k: %d file: %s'%(k, file))
						tmp_file = readFileBand(file, 1)
						# Register fail if don't possible read data file, exception in MOD13
						if (tmp_file != []):
							# Case MOD09
							if(SUB_NAMES[k][0] == 'MODIS'):
								MOD09_NDVI = readFileBand(file, 1)
								MOD09_QUAL = readFileBand(file, 4)
							else:
								# Case find MOD13_3KM then too read MOD13_250M
								if(SUB_NAMES[k][0] == 'MOD13_250M'):
									file = inputDir + SUB_NAMES[k][0] + '\\' + str(year) + '\\' + str_month + '\\' + str(year) + str_month + str_day + '_' + SUB_NAMES[k][1] + '.tif' 
									print('data: %s\tfile: %s'%(SUB_NAMES[k][0], file))
									MOD13_250M = readFileBand(file, 1)
									if(MOD13_250M == []):
										print('Fail to read: %s'%(file))
								else:
									DATA_TMP[3*k] = tmp_file
									DATA_TMP[3*k+1] = readFileBand(file, 2)
									DATA_TMP[3*k+2] = readFileBand(file, 3)

						else:
							if(SUB_NAMES[k][0] != 'MOD13_3KM'):
								print('Fail to read: %s'%(file))
								break
					if(typeData == 'TAB_ALL'):
						#print('---- culture')
						for sample in range(len(CULTURE)):
							#print('sample:%d'%(sample))
							x = CULTURE[sample][0]
							y = CULTURE[sample][1]
							
							[r250m_x, r250m_y] = conv_3km_to('250m', y, x)
							[r500m_x, r500m_y] = conv_3km_to('500m', y, x)
							for data_in in range(len(SUB_NAMES)):
								if(SUB_NAMES[data_in][0] == 'MODIS'):
									if(MOD09_NDVI != []):
										# If quality data MOD09 is good
										if((int(MOD09_QUAL[r500m_y][r500m_x]) & 7 == 0) and (row_cult[5+3*data_in+0] == 0 or (MOD09_NDVI[r500m_y][r500m_x] < 1.2*row_cult[5+3*data_in+0] and MOD09_NDVI[r500m_y][r500m_x] > 0.8 *row_cult[5+3*data_in+0]))):
											row_cult[5+3*data_in+0] = MOD09_NDVI[r500m_y][r500m_x]
								else:
									if(ant_cult[sample] == 0 and DATA_TMP[3*data_in][y][x] > 0.01):
										ant_cult[sample] = DATA_TMP[3*data_in][y][x]
									if(DATA_TMP[3*data_in][y][x] < 1.20*ant_cult[sample] and DATA_TMP[3*data_in][y][x] > 0.80*ant_cult[sample]):
										ant_cult[sample] = DATA_TMP[3*data_in][y][x]
										row_cult[5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
										row_cult[5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
										row_cult[5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
							data_in = data_in + 1
							if(MOD13_250M != []):
								row_cult[5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
								row_cult[5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
								row_cult[5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
								row_cult[5+3*data_in+3] = MOD13_250M[r250m_y][r250m_x]
							#print(DATA_CULT.getLine())
							row_cult[0:5] = [0, sample, month, day, doy]
							DATA_CULT.addAtrib(row_cult)
							row_cult = numpy.zeros((n_val))
							DATA_CULT.addSample()
						#print('---- pasture')
						for sample in range(len(PASTURE)):
							x = PASTURE[sample][0]
							y = PASTURE[sample][1]
							
							[r250m_x, r250m_y] = conv_3km_to('250m', y, x)
							[r500m_x, r500m_y] = conv_3km_to('500m', y, x)
							for data_in in range(len(SUB_NAMES)):
								if(SUB_NAMES[data_in][0] == 'MODIS'):
									if(MOD09_NDVI != []):
										# If quality data MOD09 is good
										if((int(MOD09_QUAL[r500m_y][r500m_x]) & 7 == 0) and (row_cult[5+3*data_in+0] == 0 or (MOD09_NDVI[r500m_y][r500m_x] < 1.2*row_past[5+3*data_in+0] and MOD09_NDVI[r500m_y][r500m_x] > 0.8 *row_past[5+3*data_in+0]))):
											row_past[5+3*data_in+0] = MOD09_NDVI[r500m_y][r500m_x]
								else:
									if(ant_past[sample] == 0 and DATA_TMP[3*data_in][y][x] > 0.01):
										ant_past[sample] = DATA_TMP[3*data_in][y][x]
									if(DATA_TMP[3*data_in][y][x] < 1.20*ant_past[sample] and DATA_TMP[3*data_in][y][x] > 0.80*ant_past[sample]):
										ant_past[sample] = DATA_TMP[3*data_in][y][x]
										row_past[5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
										row_past[5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
										row_past[5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
							data_in = data_in + 1
							if(MOD13_250M != []):
								row_past[5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
								row_past[5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
								row_past[5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
								row_past[5+3*data_in+3] = MOD13_250M[r250m_y][r250m_x]
							#print(DATA_PAST.getLine())
							row_past[0:5] = [1, sample, month, day, doy]
							DATA_PAST.addAtrib(row_past)
							row_past = numpy.zeros((n_val))
							DATA_PAST.addSample()
						#print('---- nature')
						for sample in range(len(NATURE)):
							x = NATURE[sample][0]
							y = NATURE[sample][1]
							
							[r250m_x, r250m_y] = conv_3km_to('250m', y, x)
							[r500m_x, r500m_y] = conv_3km_to('500m', y, x)
							for data_in in range(len(SUB_NAMES)):
								if(SUB_NAMES[data_in][0] == 'MODIS'):
									if(MOD09_NDVI != []):
										# If quality data MOD09 is good
										if((int(MOD09_QUAL[r500m_y][r500m_x]) & 7 == 0) and (row_cult[5+3*data_in+0] == 0 or (MOD09_NDVI[r500m_y][r500m_x] < 1.2*row_natu[5+3*data_in+0] and MOD09_NDVI[r500m_y][r500m_x] > 0.8 *row_natu[5+3*data_in+0]))):
											row_natu[5+3*data_in+0] = MOD09_NDVI[r500m_y][r500m_x]
								else:
									if(ant_natu[sample] == 0 and DATA_TMP[3*data_in][y][x] > 0.01):
										ant_natu[sample] = DATA_TMP[3*data_in][y][x]
									if(DATA_TMP[3*data_in][y][x] < 1.20*ant_natu[sample] and DATA_TMP[3*data_in][y][x] > 0.80*ant_natu[sample]):
										ant_natu[sample] = DATA_TMP[3*data_in][y][x]
										row_natu[5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
										row_natu[5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
										row_natu[5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
							data_in = data_in + 1
							if(MOD13_250M != []):
								row_natu[5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
								row_natu[5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
								row_natu[5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
								row_natu[5+3*data_in+3] = MOD13_250M[r250m_y][r250m_x]
							#print(DATA_PAST.getLine())
							row_natu[0:5] = [2, sample, month, day, doy]
							DATA_NATU.addAtrib(row_natu)
							row_natu = numpy.zeros((n_val))
							DATA_NATU.addSample()
					else:
						#print('---- culture')
						for sample in range(len(CULTURE)):

							#print('sample:%d'%(sample))
							x = CULTURE[sample][0]
							y = CULTURE[sample][1]
							[r250m_x, r250m_y] = conv_3km_to('250m', CULTURE[sample][1], CULTURE[sample][0])
							[r500m_x, r500m_y] = conv_3km_to('500m', CULTURE[sample][1], CULTURE[sample][0])
							print('[%d %d] [%d %d] [%d %d]'%(y,x, r250m_y, r250m_x, r500m_y, r500m_x))
							for data_in in range(len(SUB_NAMES)-1):
								if(SUB_NAMES[data_in][0] == 'MODIS'):
									if(MOD09_NDVI != []):
										# If quality data MOD09 is good
										if((int(MOD09_QUAL[r500m_y][r500m_x]) & 7) == 0 and MOD09_NDVI[r500m_y][r500m_x] > row_cult[sample][5+3*data_in+0]):
											row_cult[sample][5+3*data_in+0] = MOD09_NDVI[r500m_y][r500m_x]
								elif(DATA_TMP[3*data_in][y][x] > row_cult[sample][5+3*data_in]):
									row_cult[sample][5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
									row_cult[sample][5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
									row_cult[sample][5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
							if(MOD13_250M != []):
								data_in = data_in + 1
								row_cult[sample][5+3*data_in] = MOD13_250M[r250m_y][r250m_x]

								row_cult[sample][0:5] = [0, sample, year, month, day]
								DATA_CULT.addAtrib(row_cult[sample])
								row_cult[sample] = numpy.zeros((n_val))
								DATA_CULT.addSample()
						
						#print('---- pasture')
						for sample in range(len(PASTURE)):
							x = PASTURE[sample][0]
							y = PASTURE[sample][1]
							[r250m_y, r250m_x] = conv_3km_to('250m', PASTURE[sample][1], PASTURE[sample][0])
							[r500m_y, r500m_x] = conv_3km_to('500m', PASTURE[sample][1], PASTURE[sample][0])
							for data_in in range(len(SUB_NAMES)-1):
								if(SUB_NAMES[data_in][0] == 'MODIS'):
									if(MOD09_NDVI != []):
									# If quality data MOD09 is good
										if((int(MOD09_QUAL[r500m_y][r500m_x]) & 7) == 0 and MOD09_NDVI[r500m_y][r500m_x] > row_past[sample][5+3*data_in+0]):
											row_past[sample][5+3*data_in+0] = MOD09_NDVI[r500m_y][r500m_x]                                
								elif(DATA_TMP[3*data_in][y][x] > row_past[sample][5+3*data_in]):
									row_past[sample][5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
									row_past[sample][5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
									row_past[sample][5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
							if(MOD13_250M != []):
								data_in = data_in + 1
								row_past[sample][5+3*data_in+0] = MOD13_250M[r250m_y][r250m_x]

								row_past[sample][0:5] = [1, sample, year, month, day]
								DATA_PAST.addAtrib(row_past[sample])
								row_past[sample] = numpy.zeros((n_val))
								DATA_PAST.addSample()
						#print('---- nature')
						for sample in range(len(NATURE)):
							x = NATURE[sample][0]
							y = NATURE[sample][1]
							[r250m_y, r250m_x] = conv_3km_to('250m', NATURE[sample][1], NATURE[sample][0])
							[r500m_y, r500m_x] = conv_3km_to('500m', NATURE[sample][1], NATURE[sample][0])
							for data_in in range(len(SUB_NAMES)-1):
								if(SUB_NAMES[data_in][0] == 'MODIS'):
									if(MOD09_NDVI != []):
										# If quality data MOD09 is good
										if((int(MOD09_QUAL[r500m_y][r500m_x]) & 7) == 0 and MOD09_NDVI[r500m_y][r500m_x] > row_natu[sample][5+3*data_in+0]):
											row_natu[sample][5+3*data_in+0] = MOD09_NDVI[r500m_y][r500m_x]                               
								elif(DATA_TMP[3*data_in][y][x] > row_natu[sample][5+3*data_in]):
									row_natu[sample][5+3*data_in+0] = DATA_TMP[3*data_in+0][y][x]
									row_natu[sample][5+3*data_in+1] = DATA_TMP[3*data_in+1][y][x]
									row_natu[sample][5+3*data_in+2] = DATA_TMP[3*data_in+2][y][x]
							if(MOD13_250M != []):
								data_in = data_in + 1
								row_natu[sample][5+3*data_in+0] = MOD13_250M[r250m_y][r250m_x]
								row_natu[sample][0:5] = [2, sample, year, month, day]
								DATA_NATU.addAtrib(row_natu[sample])
								row_natu[sample] = numpy.zeros((n_val))
								DATA_NATU.addSample()

				DATA_CULT.save()
				DATA_PAST.save()
				DATA_NATU.save()
	# Generate fusion file 
	# ------------------------------------------------------------------------------------------------------------
	if(typeData == 'PAN'):
		if(yearStart == yearFinish):
			inc = 1
		else:
			inc = 0
		for year in range(yearStart, yearFinish+inc):
			print('...........Create PAN data in year: ' + str(year))
			outputSubfolder = outputDir + 'PAN' + '\\' + str(year)
			if not(os.path.isdir(outputSubfolder)):
				os.makedirs(outputSubfolder)
				print('.Create ' + outputSubfolder)
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				print('......Create HRV data in month: ' + str_month)
				outputSubfolder = outputDir + 'PAN' + '\\' + str(year) + '\\'  + str_month
				print('. Create PAN files in folder: ' + outputSubfolder)
				if not(os.path.isdir(outputSubfolder)):
					os.makedirs(outputSubfolder)
					print('.Create ' + outputSubfolder)
				inputSubfolder = inputDir + 'HRV' '\\' + str(year) + '\\'  + str_month
				print(inputSubfolder + '\\*' + str(year) + str_month + '*_HRV.tif')
				files = glob.glob(inputSubfolder + '\\*' + str(year) + str_month + '*_HRV.tif')
				#print(files)
				for hrvfile in files:
					bandfile = hrvfile.replace('HRV', '123')
					tmpfile  = 'c:\\tmp\\tmp_' + str(year) + str_month +'_hrv.tif'
					tmpfile1 = 'c:\\tmp\\tmp1_' + str(year) + str_month +'_hrv.tif'
					tmpfile2 = 'c:\\tmp\\tmp2_' + str(year) + str_month +'_hrv.tif'
					panfile = hrvfile.replace('HRV', 'PAN')
					if(os.path.isfile(panfile)):
						print('. PAN file found: %s'%(panfile))
						continue
					if not(os.path.isfile(bandfile)):
						print ('. Band file not found: %s' %(bandfile))
						continue
					print('. PAN band to generate: %s' %(panfile))
					
					ds = gdal.Open(fileMask123)
					y_123 = ds.RasterYSize
					x_123 = ds.RasterXSize
					DATA = numpy.zeros((2, y_123, x_123))
					DATA[0] = readFileBand(bandfile, 1)
					DATA[1] = readFileBand(bandfile, 2)
					saveRASTERfile(tmpfile, fileMask123, DATA)
					#print(otbDir + 'otbcli_BundleToPerfectSensor.bat -inp ' + hrvfile + ' -inxs ' + bandfile + ' -mode default -out ' + tmpfile1)
					subprocess.call(otbDir + 'otbcli_BundleToPerfectSensor.bat -inp ' + hrvfile + ' -inxs ' + tmpfile + ' -mode default -out ' + tmpfile1)
					#subprocess.call(exeDir + 'gdalwarp -tr 0.011613584620619 0.011613584620619 -r cubic '+ tmpfile + ' ' + tmpfile1 + ' -cutline ' + fileMaskSHP + ' -crop_to_cutline ')
					
					#print(otbDir + 'otbcli_Pansharpening.bat -inp ' + hrvfile + ' -inxs ' + tmpfile1 + ' -out ' + tmpfile2 + ' -method lmvm') 
					subprocess.call(otbDir + 'otbcli_Pansharpening.bat -inp ' + hrvfile + ' -inxs ' + tmpfile1 + ' -out ' + tmpfile2 + ' -method lmvm -method.lmvm.radiusx 2 -method.lmvm.radiusy 2') 

					#print(exeDir + 'gdalwarp '+ tmpfile2 + ' ' + panfile + ' -cutline ' + fileMask + ' -crop_to_cutline ')
					subprocess.call(exeDir + 'gdalwarp '+ tmpfile2 + ' ' + panfile + ' -cutline ' + fileMaskSHP)
					#sys.exit()    
					delFile(tmpfile)
					delFile(tmpfile1)
					delFile(tmpfile2)

	# Calc NDVI using pancromatic band (12)
	# --------------------------------------------------------------------------------------------------------------------------
	if(typeData == 'PRO'):
		if(yearStart == yearFinish):
			inc = 1
		else:
			inc = 0
		for year in range(yearStart, yearFinish+inc):
			print('...........Create PRO data in year: ' + str(year))
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				print('...........Create PRO data in month: %s' %(str_month))
				for day in range (1, 32):
					str_day = str('%02d'%day)
					ds = gdal.Open(fileMask123)
					y_123 = ds.RasterYSize
					x_123 = ds.RasterXSize
					ds = gdal.Open(fileMaskPAN)
					y_hrv = ds.RasterYSize
					x_hrv = ds.RasterXSize
					# PRO = [NDVI RED NIR HOUR]
					pro_bands = 4 
					hours = [1700, 1715, 1730, 1745, 1800, 1815, 1830, 1845, 1900]
					files = [1700, 1715, 1730, 1745, 1800, 1815, 1830, 1845, 1900]
					CLM = 2*numpy.ones((len(hours), y_hrv, x_hrv))
					PAN = numpy.zeros((3*len(hours), y_hrv, x_hrv))
					PRO = numpy.zeros((pro_bands,  y_hrv, x_hrv))
					pos = 0     
					for hour in hours:
						tmp_clmfile = 'c:\\tmp\\' + str(year) + str_month + str_day + str(hour) + '_CLM.tif'
						panfile = inputDir + 'PAN' '\\' + str(year) + '\\'  + str_month + '\\'  + str(year) + str_month + str_day + str(hour) + '_PAN.tif'
						profile = panfile.replace(str(hour) + '_PAN.tif', '_PRO.tif')
						profile = profile.replace('PAN', 'PRO')
						clmfile = panfile.replace('PAN', 'CLM')
						outputSubfolder = outputDir + 'PRO' '\\' + str(year) + '\\'  + str_month + '\\'
						if not(os.path.isdir(outputSubfolder)):
							os.makedirs(outputSubfolder)
						if not(os.path.isfile(panfile)):
							print('. Error in PAN file: %s' %(panfile))
							files[pos] = []
						else:
							print('. Readind PAN file: %s' %(panfile))
							# PAN = [RED NIR]
							PAN[2*pos]   = readFileBand(panfile, 1)
							PAN[2*pos+1] = readFileBand(panfile, 2)
							if not(os.path.isfile(clmfile)):
								print('. Error in CLM file: %s' %(clmfile))
							else:
								tmpfile1  = 'c:\\tmp\\tmp1_' + str(year) + str_month +'_pro.tif'
								tmpfile2 = 'c:\\tmp\\tmp2_' + str(year) + str_month +'_pro.tif'
								tmpfile3 = 'c:\\tmp\\tmp3_' + str(year) + str_month +'_pro.tif'
								saveRASTERfile(tmpfile1, fileMaskPAN, numpy.ones((y_hrv, x_hrv)))                                
								subprocess.call(otbDir + 'otbcli_BundleToPerfectSensor.bat -inp ' + tmpfile1 + ' -inxs ' + clmfile + ' -mode default -out ' + tmpfile2)
								subprocess.call(exeDir + 'gdalwarp '+ tmpfile2 + ' ' + tmpfile3 + ' -cutline ' + fileMaskSHP)
								CLM[pos] = readFileBand(tmpfile3, 1)
								#print('test......................')
								#sys.exit()
								#saveRASTERfile('c:\\tmp\\test_' + str(year) + str_month + str(day) + '_pro.tif', fileMaskPAN, CLM)
								delFile(tmpfile1)
								delFile(tmpfile2)
								delFile(tmpfile3)
								
						pos = pos + 1
					for  x in range(x_hrv):
						for y in range(y_hrv):
							for sample in range(5):
								if ((CLM[sample, y, x] < 1.05 and CLM[sample, y, x] > 0) and (hours[sample] in files)):
									if(not nearzero(PAN[2*sample, y, x]) and not nearzero(PAN[2*sample+1, y, x])):
										
										# NDVI = calcNDVI(NIR, RED)
										PRO[0, y, x] = calcNDVI(PAN[2*sample+1, y, x], PAN[2*sample, y, x])
										PRO[1, y, x] = PAN[2*sample,   y, x]
										PRO[2, y, x] = PAN[2*sample+1, y, x]
										PRO[3, y, x] = hours[sample]
										# Print results
										#print('[y,x]: %d, %d CLM: %.2f PRO: %.4f %.4f %.4f %.4f %.0f BND: %.4f %.4f '%(y, x, CLM[sample, y, x], PRO[0,3*y+py,3*x+px], PRO[1,3*y+py,3*x+px], PRO[2,3*y+py,3*x+px], PRO[3,3*y+py,3*x+px], PRO[4,3*y+py,3*x+px], PRO[5,3*y+py,3*x+px], PRO[6,3*y+py,3*x+px], BND[3*sample, y, x], BND[3*sample+1, y, x]))
					print('. Save output in: ' + profile)
					saveRASTERfile(profile, fileMaskPAN, PRO)
					#print('test......................')
					#sys.exit()       

	# Calculate NDVI using PRO data
	# -----------------------------------------------------------------------------------------------------------
	if(typeData == 'NDVI'):
		if(yearStart == yearFinish):
			inc = 1
		else:
			inc = 0
		for year in range(yearStart, yearFinish+inc):
			print('......Create NDVI data in year: ' + str(year))
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				print('...Create NDVI data in month: %s' %(str_month))
				for day in range (1, 32):
					str_day = str('%02d'%day)
					ds = gdal.Open(fileMask123)
					y_123 = ds.RasterYSize
					x_123 = ds.RasterXSize
					ds = gdal.Open(fileMaskHRV)
					y_hrv = ds.RasterYSize
					x_hrv = ds.RasterXSize
					OUT = numpy.zeros((y_hrv, x_hrv))
					profile = inputDir + 'PRO' '\\' + str(year) + '\\'  + str_month + '\\'  + str(year) + str_month + str_day + '_PRO.tif'
					dayfile = profile.replace('_PRO', '_DAY')
					dayfile = dayfile.replace('PRO', 'DAY_123_1000as1200')
					outfile = profile.replace('PRO', 'NDVI')
					profile = inputDir + 'PRO_1700_to_1900' '\\' + str(year) + '\\'  + str_month + '\\'  + str(year) + str_month + str_day + '_PRO.tif'
					outputSubfolder = outputDir + 'NDVI' '\\' + str(year) + '\\'  + str_month + '\\'
					if not(os.path.isdir(outputSubfolder)):
						os.makedirs(outputSubfolder)
					#PRO = [NDVI RED NIR HOUR]
					if not(os.path.isfile(profile)):
						print('. Error in PRO file: %s' %(profile))
						continue
					else:
						PRO = readFileBand(profile, 1)

					if not(os.path.isfile(dayfile)):
						print('. Error in DAY file: %s' %(dayfile))
						continue
					else:
						DAY = readFileBand(dayfile, 1)

					# Process Fusion Band using HPF (High Pass Filter): 1 pixel 3 km is the area to process 9 pixels 1 km      
					for  x in range(x_123):
						for y in range(y_123):
							if(DAY[y,x] == 0):
								continue
							pro_vet = numpy.zeros((9))
							for px in range(3):
								for py in range(3):
									if (3*x+px  >= x_hrv or 3*y+py >= y_hrv):
										continue
									pro_vet[3*py+px] = PRO[3*y + py, 3*x + px]
								pro_mean = pro_vet.mean()
								pro_std = pro_vet.std()
							if(pro_mean == 0):
								continue
							for px in range(3):
								for py in range(3):
									if (3*x+px  >= x_hrv or 3*y+py >= y_hrv):
										continue
									OUT[3*y + py, 3*x + px] = DAY[y,x] + (PRO[3*y + py, 3*x + px] - pro_mean)
					saveRASTERfile(outfile, fileMaskHRV, OUT)
					#sys.exit()


	# Calculate month's NDVI data processed in BRDF folder
	# -----------------------------------------------------------------------------------------------------------
	if(typeData == 'MEAN'):
		if(yearStart == yearFinish):
			inc = 1
		else:
			inc = 0
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()
		print('num_days:  %d'%(num_days))
		#SUB_NAMES = [['DAY_123_1000as1200', 'DAY']]
		SUB_NAMES = [['MODIS', 'MOD09']]
		n_val = 6
		header  = 'YEAR MONTH DOY N_DOY DAY'#1000as1200 1000as1400'
		n_step = 20
		n_start = 1
		# Number pixels not zero in image [DAY BRDF] to DAY_6S is 24055 and to BRDF_6S is 24036
		# Number pixels DAY | BRDF
		#n_pixels = [24054, 24042]
		n_pixels = [100, 100, 100, 100]
		num_row = 0
		for n_doy in range(n_start,n_step+1,1):
			num_row = num_row + int(num_days/n_doy)
		
		FILE = datafile('', outputDir, 'data_n_doys_mod09_1_to_20_rain.csv', header, num_row, n_val)
		first_file = False
		ds = None
		for n_doy in range(n_start,n_step+1,1):
			print('...generate to doy interval: %d'%(n_doy))
			date  = dateStart
			cont_doy = n_doy
			while(date < dateEnd):
				doy = calcDoy_date(date)
				# out date if in rain wheater
				#if(doy < 91 or doy > 301):
				# out date if in wet wheater
				if(doy >= 91 and doy <= 301):
					date = incDate(date)
					continue
				else:
					print('.get doy %d'%(doy))
				for n_file in range(len(SUB_NAMES)):

					filename =  inputDir + SUB_NAMES[n_file][0] + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date) + '_' + SUB_NAMES[n_file][1] + '.tif' 
					print('Read: ' + filename)
					TMP = readFileBand(filename, 1)
					if('09' in SUB_NAMES[n_file][1]):
						CLM = readFileBand(filename, 4)
					if(TMP != []):
						if(first_file == False):
							first_file = True
							ds = gdal.Open(filename)
							y_123 = ds.RasterYSize
							x_123 = ds.RasterXSize
							ONES = 7*numpy.ones((y_123, x_123))
							DATA = numpy.zeros((len(SUB_NAMES), y_123,x_123))
							ds = None
						#print('readind file: %s'%(filename))
						if('09' in SUB_NAMES[n_file][1]):
							FILT = numpy.logical_not(numpy.bitwise_and(CLM.astype(numpy.int32), ONES.astype(numpy.int32))) 
							TMP = numpy.logical_and(TMP > 0, FILT > 0)
						DATA[n_file] = numpy.logical_or(numpy.logical_or(DATA[n_file],0), numpy.logical_or(TMP,0))
					else:
						print('fail to read: %s'%(filename))
				cont_doy = cont_doy - 1
				#saveRASTERfile('c:\\tmp\\'+str(date) +'_MEAN_'+SUB_NAMES[n_file][0]+'.tif', filename, DATA[n_file])

				# Verify if give number of doys to generate data
				if(cont_doy == 0):
					cont_doy = n_doy
					n_NDVI = []
					data_atrib = [int(str(date)[0:4]), int(str(date)[4:6]), doy, n_doy]
					for n_file in range(len(SUB_NAMES)):
						data_atrib.append(int(numpy.sum(DATA[n_file])))
						#data_atrib.append(int(100*numpy.sum(DATA[n_file])/n_pixels[n_file]))

					print(data_atrib)
					FILE.addAtrib(data_atrib)
					FILE.addSample()

					#output_filename = outputDir + 'DOY\\n_doy_' + str(n_doy) + '_date_' + str(calcDoy_date(date)) + '.tif'
					#saveRASTERfile(output_filename, fileMask123, DATA) 
					DATA = numpy.zeros((len(SUB_NAMES), y_123,x_123))

				date = incDate(date)
				#print('go out.........')
				#sys.exit()
			FILE.save()
		print('...save file')


	# Calculate year's histogram data
	# ------------------------------------------------------------------------------------------------------------ 
	if(typeData == 'HIST'):

		if(yearStart == yearFinish):
			inc = 1
		else:
			inc = 0
		
		num_days = calcDays(dateStart, dateEnd)
		if(num_days == -1):
			print(' Error in dates, please verify.')
			sys.exit()
		ds = gdal.Open(fileMask123)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		ANG = numpy.zeros((360,5))
		ATM = numpy.zeros((2000,4))
		for i in range(360):
			ANG[i,0] = i
		for i in range(2000):
			ATM[i,0] = i
		for year in range(yearStart, yearFinish+inc):
			print('......Create CMP data')
			if(year != yearStart):
				monthStart = 1
			if(year == yearFinish):
				monthEnd = monthFinish
			else:
				monthEnd = 12
			for month in range(monthStart, monthEnd+1):
				str_month = str('%02d'%month)
				for day in range (1, 32):
					str_day = str('%02d'%day)
					for hour in (1300, 1315, 1330, 1345, 1400, 1415, 1430, 1445, 1500, 1515, 1530, 1545, 1600, 1615, 1630, 1645, 1700):
						str_hour = str(hour)
						anglefile   = inputDir + 'ANGLE\\' + str(year) + '\\'  + str_month + '\\'  + str(year) + str_month + str_day + str_hour + '_ANGLES.tif'
						if(not os.path.isfile(anglefile)):
							print('... Error in data file: %s' %(anglefile))
							continue
						else: 
							print('... Read data file: %s' %(anglefile))                                     
							# get VZA View Zenith Angle (msg_zenres.tif)
							VZA = readFileBand(anglefile, 3)
							# get VAZ View Azimuth Angle (msg_azres.tif)
							VAZ = readFileBand(anglefile, 1)
							# get SZA Sun Zenith Angle (sun_zenres.tif)
							SZA = readFileBand(anglefile, 4)
							# get SAZ Sun Azimuth Angle (sol_azres.tif)
							SAZ = readFileBand(anglefile, 2)
							for  x in range(x_123):
								for y in range(y_123):
									if(VZA[y,x] == 0):
										continue
									else:
										ANG[int(VZA[y,x]), 1] = ANG[int(VZA[y,x]), 1] + 1
										ANG[int(VAZ[y,x]), 2] = ANG[int(VAZ[y,x]), 2] + 1
										ANG[int(SZA[y,x]), 3] = ANG[int(SZA[y,x]), 3] + 1
										ANG[int(SAZ[y,x]), 4] = ANG[int(SAZ[y,x]), 4] + 1

					atmfile = inputDir + 'ATM\\' + str(year) + '\\'  + str_month + '\\'  + str(year) + str_month + str_day + '_ATM.tif'
					if(not os.path.isfile(atmfile)):
						print('... Error in data file: %s' %(atmfile))
						continue
					else:
						print('... Read data file: %s' %(atmfile)) 
						# get Water_Vapour  
						WV = readFileBand(atmfile, 10)*100
						# get Total Ozone Data Band, convert Dobson Unit to cm-atm
						OZ = readFileBand(atmfile, 7)
						# get Deep_Blue Data Band
						DB = readFileBand(atmfile, 1)*1000
						for  x in range(x_123):
							for y in range(y_123):
								if(VZA[y,x] == 0):
									continue
								else:
									ATM[int(WV[y,x]), 1] = ATM[int(WV[y,x]), 1] + 1
									ATM[int(OZ[y,x]), 2] = ATM[int(OZ[y,x]), 2] + 1
									ATM[int(DB[y,x]), 3] = ATM[int(DB[y,x]), 3] + 1

		#Save data file
		print('--- Save output file: %s'%(outputDir + 'atmfile_hist.csv'))
		numpy.savetxt(outputDir + 'atmfile_hist.csv', ATM, fmt='%1.4f', comments='', header = "N WV OZ DB", delimiter = ' ', newline='\n')
		print('--- Save output file: %s'%(outputDir + 'angfile_hist.csv'))
		numpy.savetxt(outputDir + 'angfile_hist.csv', ANG, fmt='%1.4f', comments='', header = "N VZA VAZ SZA SAZ", delimiter = ' ', newline='\n')

import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import math
import os
from time import  strftime, localtime, sleep
import sys
import gc
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from time import strftime, localtime
from sys import argv
from utils import mapDict, incHour, incDate, decDate, infoFile, printlog,  readFileBand, readPixel, readMaskGoias, wait
import time
import threading
import queue
import psutil
import random
from sklearn.preprocessing import Imputer
import subprocess

usage = """\
Usage: %s [OPTIONS]
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
""" % argv[0]

# Constantes
# To print lib version
# lib.__version__
fileGoiasGeo = os.getcwd() + '\\shapes\\201301011300_123.tif'
exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
# Variables of threads
NUM_THREADS = 12
queueLock = threading.Lock()
workQueue = queue.Queue(4*NUM_THREADS)
threads = []
exitFlag = 0

class myThread (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        printlog (True, "Starting " + self.name)
        process_data(self.name, self.q)
        printlog (True, "Exiting " + self.name)

def process_data(threadName, q):
    global queueLock
    global workQueue
    global exitFlag
    #print("Enter in process_data. exitFlag: " + str(exitFlag))

    while not exitFlag:
        #print("exitFlag is False")
        queueLock.acquire()
        #print("\nqueueLock acquire")
        if not workQueue.empty():
            #print("Queue is not empty")
            # put data in input parameters
            queueLock.release()
            #process name
        else:
            queueLock.release()
        time.sleep(1)


def modis_atm(inputDir, outputDir, date):
    
    row = 204
    col = 211
    neighbor  = 3
    num_band  = 4
    NO_DATA   = -9999
    LIM_BAND  = [800, 800, 4000, 8000] 
    MAT_BACK  = numpy.zeros((num_band, neighbor, row, col), dtype=numpy.float) 
    MAT_FORW  = numpy.zeros((num_band, neighbor, row, col), dtype=numpy.float)

    VALUE      = numpy.zeros((num_band, row, col), dtype=numpy.float)
    CONT_FORW  = numpy.zeros((num_band, row, col), dtype=numpy.float)
    CONT_BACK  = numpy.zeros((num_band, row, col), dtype=numpy.float)
    VALUE_RE   = numpy.zeros((num_band, row, col), dtype=numpy.float)

    mask = readFileBand(fileGoiasGeo, 1)
    n_zero = 0
    # Calc n_zero in mask
    for i in range(0, row, 1):
        for j in range(0, col, 1):
            if(mask[i,j] <= 0):
                n_zero = n_zero + 1

    if(len(mask) == 0):
        print('Error to read Mask file')
        return 0
    print('----- Reading Files')
    mat_cont = 0
    search_date = date
    filename = ''
    print('----- Read back files')
    for x in range (neighbor):
        if(x == 0):
            filename = inputDir + str(search_date)[0:4] + '\\' + str(search_date)[4:6] + '\\'+ str(search_date) + '_MOD08.tif'
        else:
            if(filename.find('_MOD08.tif') != -1):
                search_date = decDate(search_date)
                filename =  inputDir + str(search_date)[0:4] + '\\' + str(search_date)[4:6] + '\\'+ str(search_date) + '_MYD08.tif' 
            else:
                filename =  inputDir + str(search_date)[0:4] + '\\' + str(search_date)[4:6] + '\\'+ str(search_date) + '_MOD08.tif'
        if(len(readFileBand(filename, 1))):
            for band in range(num_band):
                MAT_BACK[band, mat_cont, :, :] = readFileBand(filename, band + 1)
            mat_cont = mat_cont + 1
            #print('Reading file: ' + filename)
        else:
            print('Error to read file: ' + filename)
    # Test if get many files to process
    if(mat_cont < (0.7*neighbor)):
        print('Error in back neighbor (%d) to generate: %d' %(mat_cont, date))
        return 0
    print('----- Read forward files')        
    mat_cont = 0
    search_date = date
    for x in range (neighbor):
        if(x == 0):
            filename = inputDir + str(search_date)[0:4] + '\\' + str(search_date)[4:6] + '\\'+ str(search_date) + '_MYD08.tif'
        else:
            if(filename.find('_MYD08.tif') != -1):
                search_date = incDate(search_date)
                filename =  inputDir + str(search_date)[0:4] + '\\' + str(search_date)[4:6] + '\\'+ str(search_date) + '_MOD08.tif' 
            else:
                filename =  inputDir + str(search_date)[0:4] + '\\' + str(search_date)[4:6] + '\\'+ str(search_date) + '_MYD08.tif'
        if(len(readFileBand(filename, 1))):
            for band in range(num_band):
                MAT_FORW[band, mat_cont, :, :] = readFileBand(filename, band + 1)
            mat_cont = mat_cont + 1
            #print('Reading file: ' + filename)
        else:
            print('Error to read file: ' + filename)
    # Test if get many files to process
    if(mat_cont < (0.7*neighbor)):
        print('Error in forward neighbor (%d) to generate: %d' %(mat_cont, date))
        return 0
    print('----- Process Band')
    #print('mat_cont_forw: ' + str(mat_cont))
    for band in range(num_band):
        print('----- Process Band: %i' %(band))
        for i in range(0, row, 1):
            for j in range(0, col, 1):
                # Verify pixels in maskGoias
                if(mask[i,j] <= 0):
                    continue
                value_back = 0
                value_forw = 0
                cont_back = 0
                cont_forw = 0
                for x in range(neighbor):
                    #print(x)
                    value = MAT_BACK[band, x, i, j]
                    if(value != NO_DATA and value != 0):
                        cont_back = x
                        value_back = value
                        break
                #print(value_back)
                #print('---')
                for x in range(neighbor):
                    #print(x)
                    value = MAT_FORW[band, x, i, j]
                    if(value != NO_DATA and value != 0):
                        cont_forw = x
                        value_forw = value
                        break
                #print(value_forw)
                if(value_back <= 0 and value_forw <= 0):
                    VALUE[band, i,j] = numpy.nan
                elif(value_back <= 0 or value_back > 10*value_forw):
                    VALUE[band, i,j] = value_forw
                elif(value_forw <= 0 or value_forw > 10*value_back):
                    VALUE[band, i,j] = value_back
                else:
                    VALUE[band, i,j] = (value_back + value_forw)/2
                if(VALUE[band, i,j] > LIM_BAND[band]):
                    print('problem in band %d lim %d data is %d'%(band, LIM_BAND[band],VALUE[band, i,j]))
                    VALUE[band, i,j] = numpy.nan

                CONT_BACK[band, i,j] = cont_back
                CONT_FORW[band, i,j] = cont_forw

    # Complete data in image using near data
    imp = Imputer(strategy="mean")
    for band in range(num_band):
        VALUE[band] = imp.fit_transform(VALUE[band])

    print('----- Process reasample data')
    for i in range(0, row, 1):
        for j in range(0, col, 1):
            # Verify pixels in maskGoias
            if(mask[i,j] > 0):
                for band in range(num_band):
                    sum_data = 0
                    cont_data = 0
                    for i_k in range(-5, 5, 1):
                        for j_k in range(-5, 5, 1):
                            # Verify limits
                            if((i + i_k) > 0 and (i + i_k) < row and (j + j_k) > 0 and (j + j_k) < col):
                                # Verify data mask
                                if(mask[i+i_k,j+j_k] > 0 and VALUE[band, i+i_k,j+j_k] > 0):
                                    sum_data = sum_data + VALUE[band, i+i_k,j+j_k]
                                    cont_data = cont_data + 1
                    if(cont_data == 0):
                        if(VALUE[band, i,j] > 0):
                            VALUE_RE[band, i,j] = VALUE[band, i,j]
                            CONT_FORW[band, i,j] = 1
                        else:
                            VALUE_RE[band, i,j] = numpy.sum(VALUE[band])/(col*row-n_zero)
                            CONT_FORW[band, i,j] = 2
                            CONT_BACK[band, i,j] = (col*row-n_zero)
                    else:    
                        if(VALUE[band, i,j] <= 0):
                            VALUE_RE[band, i,j] = sum_data/cont_data
                            CONT_FORW[band, i,j] = 3
                        else:
                            VALUE_RE[band, i,j] = sum_data/cont_data
                            CONT_FORW[band, i,j] = cont_data
    VALUE = VALUE_RE
    # Save georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly)
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")
    suboutputDir = outputDir + str(date)[0:4] + '\\' + str(date)[4:6] 
    if not(os.path.isdir(suboutputDir)):
        os.makedirs(suboutputDir)
    filename = suboutputDir + '\\' + str(date) + '_ATM.tif'
    dsOut  = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 3*num_band, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    # Correct band value
    # band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
    # band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
    # band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
    # band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)

    print('----- Saving file')
    for band in range(num_band):
        print('----- Saving band: %d' %(band))
        bandOut=dsOut.GetRasterBand(3*band + 1)
        if(band != 2):
            BandWriteArray(bandOut, VALUE[band, :, :]*0.001)
        else:
            BandWriteArray(bandOut, VALUE[band, :, :]*0.1)
        bandOut=dsOut.GetRasterBand(3*band + 2)
        BandWriteArray(bandOut, CONT_FORW[band, :, :])
        bandOut=dsOut.GetRasterBand(3*band + 3)
        BandWriteArray(bandOut, CONT_BACK[band, :, :])


class main():

    argDict = mapDict(argv, usage)
    global exitFlag
    gc.collect()
    gdal.UseExceptions()    
    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict:
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
    else:
        exit(usage)
    if(inputDir.find('MODIS') == -1):
        inputDir = inputDir + 'MODIS\\'
    if(outputDir.find('ATM') == -1):
        outputDir = outputDir + 'ATM\\'
    date = dateStart
    while(date < dateEnd):
        if(modis_atm(inputDir, outputDir, date) == 0):
            print("----- Error to generate MODIS: %d"%(date))
        else:
            print("----- Generate MODIS: %d"%(date))
        date = incDate(date)

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
		print('..0')
		for y in range(y_MSG):
			if(int(y/y_MSG*100) >= process + 5):
				process = int(y/y_MSG*100)
				print('..%d'%(process))
			for x in range(x_MSG):
				if(CLASS_MSG[y][x] != 15):
					CLASS_PERC[y][x] = process_class(x, y, CLASS_MSG[y][x], CLASS_ORIG)
					#if(CLASS_PERC[y][x] > 0.7):
					#	print('CLASS_MSG[%d][%d] = %d : %.2f'%(y, x, CLASS_MSG[y][x], CLASS_PERC[y][x]))
		print('..100 - done')	
		saveRASTERfile(terraClass_perc, terraClass_3km, [CLASS_MSG, CLASS_PERC])
	import numpy


class main:
	folder = 'C:\\Users\\carlos\\Dropbox\\ANALYSES\\OTTO_1010\\'
	tab_data = folder + 'data_trmm_msg_ndvi_ndwi_v2.csv'
	NAME = ['CLASS', 'MONTH', 'TRMM', 'NDVI', 'NDWI', 'NDVI_1', 'NDWI_1', 'NDVI_2', 'NDWI_2', 'NDVI_3', 'NDWI_3']
	DATA = numpy.loadtxt(tab_data, delimiter = ' ', skiprows = 1)

	print(DATA[0])
	DATA_T = numpy.transpose(DATA)
	print(DATA_T[0])
	for i, var in zip(range(len(NAME)), NAME):
		#print(i)
		#print(var)
		globals()[var] = DATA_T[i]
	n_class = int(numpy.amax(CLASS))
	OUT_DATA = numpy.zeros((n_class+1,12))

	for basin in numpy.arange(1,n_class+1,1):
		paser = numpy.logical_and(CLASS==basin, 1)
		#print(paser)
		OUT_DATA[basin][0] = basin
		OUT_DATA[basin][3] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI[paser]))
		OUT_DATA[basin][5] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI_1[paser]))
		OUT_DATA[basin][7] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI_2[paser]))
		OUT_DATA[basin][9] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI_3[paser]))

		OUT_DATA[basin][4] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI[paser]))
		OUT_DATA[basin][6] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI_1[paser]))
		OUT_DATA[basin][8] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI_2[paser]))
		OUT_DATA[basin][10] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI_3[paser]))
		OUT_DATA[basin][1] = numpy.argmax([OUT_DATA[basin][3], OUT_DATA[basin][5], OUT_DATA[basin][7], OUT_DATA[basin][9]])
		OUT_DATA[basin][2] = numpy.argmax([OUT_DATA[basin][4], OUT_DATA[basin][6], OUT_DATA[basin][8], OUT_DATA[basin][10]])
		
		print(OUT_DATA[basin])
	numpy.savetxt(folder + 'corrcoef.csv', OUT_DATA, header = "BASIN, NDVI_MAX, NDWI_MAX, NDVI_0, NDWI_0, NDVI_1, NDWI_1, NDVI_2, NDWI_2, NDVI_3, NDWI_3", delimiter=",",fmt='%.4f')
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

import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
import glob
import os.path
import shutil
from sys import argv
from utils import mapDict

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to output data
        -fi     folter to input data
        -ty     type input data (grb, tar, hrf, h5)
""" % argv[0]


yearStart  = 2015
yearEnd    = 2017
monthStart = 1
monthEnd   = 13
vet_month = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

class main():

    argDict = mapDict(argv, usage)

    if "-ty" in argDict and "-fo" in argDict and "-fi" in argDict:
        typeData = argDict["-ty"]
        dirData = argDict["-fi"]
        dirOutput = argDict["-fo"]

    else:
        exit(usage)     
    if not(typeData == 'grb' or typeData == 'tar' or typeData == 'hdf' or typeData == 'h5'):
        exit(usage)

    print('.Organize data in folder: ' + dirData)
    for year in range(yearStart, yearEnd):
        print('...........Organize data in year: ' + str(year))
        folder = dirOutput + '\\' + str(year)
        if not(os.path.isdir(folder)):
            os.makedirs(folder)
        for month in range(monthStart, monthEnd):
            print('......Organize data in month: ' + vet_month[month])
            folder = dirOutput + str(year) + '\\'  + vet_month[month]
            print('.Organize files in folder: ' + folder)
            files = glob.glob(dirData + '*-' + str(year) + vet_month[month] + '*.' + typeData)
            rows = len(files)
            if(rows == 0):
                print('. No files to folder: ' + folder)
                continue
            if not(os.path.isdir(folder)):
                os.makedirs(folder)
            id = len(dirData)
            for i in range (0, rows, 1):
                try:                         
                    shutil.move(files[i], folder + '\\' + files[i][id:]) 
                    print('.Move file: ' + files[i][id:]) 
                except:
                    print('.File exist in output folder:' + files[i][id:])
                    try:
                        os.remove(files[i])
                        print('....' + files[i])
                        print('.. Delete file: ' + files[i][id:])
                    except:
                        print('.. Problem to delete file: ' + files[i][id:])
import pandas as pd
import plotly.plotly as py
from plotly.graph_objs import *
import plotly.tools as tls

df = pd.read_csv(
    filepath_or_buffer='https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data', 
    header=None, 
    sep=',')

df.columns=['sepal_len', 'sepal_wid', 'petal_len', 'petal_wid', 'class']
df.dropna(how="all", inplace=True) # drops the empty line at file-end

df.tail()
# split data table into data X and class labels y

X = df.ix[:,0:4].values
y = df.ix[:,4].values

# plotting histograms

traces = []

legend = {0:False, 1:False, 2:False, 3:True}

colors = {'Iris-setosa': 'rgb(31, 119, 180)', 
          'Iris-versicolor': 'rgb(255, 127, 14)', 
          'Iris-virginica': 'rgb(44, 160, 44)'}

for col in range(4):
    for key in colors:
        traces.append(Histogram(x=X[y==key, col], 
                        opacity=0.75,
                        xaxis='x%s' %(col+1),
                        marker=Marker(color=colors[key]),
                        name=key,
                        showlegend=legend[col]))

data = Data(traces)

layout = Layout(barmode='overlay',
                xaxis=XAxis(domain=[0, 0.25], title='sepal length (cm)'),
                xaxis2=XAxis(domain=[0.3, 0.5], title='sepal width (cm)'),
                xaxis3=XAxis(domain=[0.55, 0.75], title='petal length (cm)'),
                xaxis4=XAxis(domain=[0.8, 1], title='petal width (cm)'),
                yaxis=YAxis(title='count'),
                title='Distribution of the different Iris flower features')

fig = Figure(data=data, layout=layout)
py.iplot(fig)

                

        
 

   


