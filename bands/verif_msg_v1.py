# Code to verify MSG data in hour interval. 
# Input:  123 reflectance
# Output: NDVI daily, samples daily and monthy
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
import itertools


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
	outputSubfolder = outputDir + 'DAILY\\'
	if not(os.path.isdir(outputSubfolder)):
		os.makedirs(outputSubfolder)

	HOURS = [1300, 1315, 1330, 1345, 1400]
	ds = gdal.Open(fileMask123)
	y_123 = ds.RasterYSize
	x_123 = ds.RasterXSize
	OUT_MONTH = numpy.zeros((y_123, x_123))
	# Process in range of data
	date = dateStart
	while(date < dateEnd):
		print('. Create data to: %s'%(date))
		#for hour in HOURS:
		OUT_DAY  = numpy.zeros((y_123, x_123))
		OUT_NDVI = numpy.zeros((y_123, x_123))
		# DATA[3][y_123][x_123] = [hour][y_123][x_123]
		DATA = numpy.zeros((len(HOURS), y_123, x_123))
		FILT = numpy.zeros((len(HOURS), y_123, x_123))
		for i in range(len(HOURS)):
			hour = HOURS[i]
			print('hour: %s'%(hour))
			# Read bands RED and NIR
			inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
			filename = inputSubDir.replace('___', '123_ref') + str(date) + str(hour) + '_123.tif'
			if os.path.isfile(filename) == 0:
				print('..Warning! Invalid data file: %s' %filename)
			else:
				RED = readFileBand(filename, 1)
				NIR = readFileBand(filename, 2)
				NDVI = numpy.zeros((y_123, x_123))
				for x in range(x_123):
					for y in range(y_123):
						try:

							if(RED[y,x] != 0 or NIR[y,x] != 0):
								NDVI[y,x] = (NIR[y,x] - RED[y,x])/(NIR[y,x] + RED[y,x])
							else:
								NDVI[y,x] = 0
						except:
							print('...Error process data file: %s'%(filename))
				# Read CLM data
				filename = inputSubDir.replace('___', 'CLM') + str(date) + str(hour) + '_CLM.tif'
				if os.path.isfile(filename) == 0:
					print('..Error! Invalid data file: %s' %filename)
				else:
					CLM = readFileBand(filename, 1)
					DATA[i] =   NDVI *(CLM <=1.05)
					DATA[i] = DATA[i]*(CLM > 0)
					FILT[i] = DATA[i]*(NDVI<1)
			OUT_DAY = ((FILT[0] + FILT[1] + FILT[2]) >= 1)
			for x in range(x_123):
				for y in range(y_123):
					for z in range(len(HOURS)):
						if(DATA[z][y,x] != 0):
							OUT_NDVI[y,x]  = DATA[z][y,x]
							break

			if(numpy.mean(OUT_DAY) > 0):
				outputSubfolder2 = outputSubfolder + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
				if not(os.path.isdir(outputSubfolder2)):
					os.makedirs(outputSubfolder2)
				output_filename = outputSubfolder2 + str(date) + '_NDVI.tif'
				saveRASTERfile(output_filename, fileMask123, OUT_NDVI)
				output_filename = outputSubfolder2 + str(date) + '_SAMP.tif'
				saveRASTERfile(output_filename, fileMask123, OUT_DAY)
		new_date = incDate(date)
		# If add month save data
		if(str(date)[4:6] != str(new_date)[4:6]):
			output_filename = outputSubfolder + '\\' + str(date)[0:4] + str(date)[4:6] + '_SAMP.tif'
			saveRASTERfile(output_filename, fileMask123, OUT_MONTH)			
			OUT_MONTH = numpy.zeros((y_123, x_123))
		else:
			OUT_MONTH = OUT_MONTH + OUT_DAY
		date = new_date