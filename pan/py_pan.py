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

