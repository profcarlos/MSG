import os
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import datetime
import time
import tarfile
from time import  strftime, localtime, sleep
import glob
from sys import argv
import msvcrt
from scipy.signal import savgol_filter
import sys

path = os.getcwd()

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
	   
def nearzero(value):
	if(value < 0.0001):
		return True
	else:
		return False
		
def calcNDVI(NIR, RED):
	x_len = 0
	try:
		x_len = len(NIR)

	except:
		if(RED > 0 and NIR > 0):
			NDVI = (NIR - RED)/(NIR + RED)
			return NDVI
		else:
			NDVI = 0
			return NDVI
	NDVI = numpy.zeros((x_len))
	for x in range(x_len):
		if(RED[x] > 0 or NIR[x] > 0):
			NDVI[x] = (NIR[x] - RED[x])/(NIR[x] + RED[x])
		else:
			NDVI[x] = 0       
	return NDVI

def wait():
	msvcrt.getch()

def calcDays(dateStart, dateEnd):
	dif_days = dateEnd - dateStart
	try:
		if(dif_days  < 31 ):
			num_days = dif_days + 1
		elif(dif_days < 1131):
			dif_days = str('%04d'%dif_days)
			num_days = int(dif_days[:2])*31 + int(dif_days[2:4]) + 1
		elif(dif_days < 31131):
			dif_days = str('%06d'%dif_days)
			num_days = int(dif_days[:2])*366 + int(dif_days[2:4])*31 + int(dif_days[4:6]) + 1
		else:
			return -1
		return num_days
	except:
		print('Verify data:')
		print(dateStart)
		print(dateEnd)
		print(dif_days)
		sys.exit()
	
def readPixel(file, band, posx, posy):
	DATA = readFileBand(file, band)
	if (len(DATA) != 1):
		if(DATA[posx, posy] == -9999):
			return False, 0
		else:
			return True, DATA[posx, posy]
	else:
		return False, 0

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

def readMaskGoias():
	try:
		vet = numpy.loadtxt(os.getcwd() + '\\shape\\maskGoias.txt')
	except:
		print('Error reading file MaskGoias')
		sys.exit(1)
	return vet

def createMask(inputFile, outputFile):
	
	DATA = readFileBand(inputFile, 1)
	if(DATA == []):
		print("Error reading input file: %s"%(inputFile))
		return False
	else:
		ds = gdal.Open(inputFile)
		y_123 = ds.RasterYSize
		x_123 = ds.RasterXSize
		DATA = numpy.zeros((2, y_123, x_123))
		for x in range(x_123):
			for y in range(y_123):
				DATA[0][y][x] = x
				DATA[1][y][x] = y
		saveRASTERfile(outputFile, inputFile, DATA)
	return True

def incHour(hour):
	# Increment hour and minuts in 15 minuts
	strHour = str(hour)
	tmpHour = int(strHour[:2]) 
	tmpMin = int(strHour[2:])
	if tmpMin < 45:
		tmpMin = tmpMin + 15
	else:
		tmpMin = 0
		tmpHour = tmpHour + 1
	strHour = str(tmpHour)
	if tmpMin == 0:
		strHour = strHour + '0' + str(tmpMin)
	else:
		strHour = strHour + str(tmpMin)
	hour = int(strHour)
	return hour

def inc_date_doy(date, int_doy):
	for inc in range(int_doy):
		date = incDate(date)
	return date

def incDate(date):

	dayStr = str(date)
	printlog(False, dayStr)
	day   = 10*int(dayStr[6]) + int(dayStr[7])
	month = 10*int(dayStr[4]) + int(dayStr[5])
	year  = 1000*int(dayStr[0]) + 100* int(dayStr[1]) + 10*int(dayStr[2]) + int(dayStr[3])
	day = day + 1
	if(((day == 30 and year == 2016) or (day == 29 and year != 2016) and month == 2)) or (day == 31 and (month == 4 or month == 6 or month == 9 or month == 11)) or (day == 32 and (month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12)):
		if(month == 12):
			year = year + 1
			month = 1
			day = 1
		else:
			month = month + 1
			day = 1
				
	dayStr = str(year)
	if month < 10:
		dayStr = dayStr + '0' + str(month)
	else:
		dayStr = dayStr + str(month)
	if day < 10:
		dayStr = dayStr + '0' + str(day)
	else:
		dayStr = dayStr + str(day)
	return int(dayStr)

def decDate(date):
	
	dayStr = str(date)
	day   = 10*int(dayStr[6]) + int(dayStr[7])
	month = 10*int(dayStr[4]) + int(dayStr[5])
	year  = 1000*int(dayStr[0]) + 100* int(dayStr[1]) + 10*int(dayStr[2]) + int(dayStr[3])
	day = day - 1
	if(day == 0): 
		if(month == 1):
			year = year - 1
			month = 12
			day = 31
		else:
			month = month - 1
			if (year != 2016 and month == 2):
				day = 28    
			else:
				day = 29       
			if(month == 4 or month == 6 or month == 9 or month == 11):
				day = 30    
			if(month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12):
				day = 31           
	dayStr = str(year)
	if month < 10:
		dayStr = dayStr + '0' + str(month)
	else:
		dayStr = dayStr + str(month)
	if day < 10:
		dayStr = dayStr + '0' + str(day)
	else:
		dayStr = dayStr + str(day)
	return int(dayStr)    
'''    
def calcDoy(year, month, day):
	doy = datetime.datetime(year, month, day).timetuple().tm_yday
	return doy
'''
# Reference: http://stackoverflow.com/questions/620305/convert-year-month-day-to-day-of-year-in-python
def calcDoy_date(date):
	""" given year, month, day return day of year
		Astronomical Algorithms, Jean Meeus, 2d ed, 1998, chap 7 """
	dayStr = str(date)  
	D  = 10*int(dayStr[6]) + int(dayStr[7])
	M  = 10*int(dayStr[4]) + int(dayStr[5])
	Y  = 1000*int(dayStr[0]) + 100* int(dayStr[1]) + 10*int(dayStr[2]) + int(dayStr[3])
	if is_leap_year(Y):
		K = 1
	else:
		K = 2
	N = int((275 * M) / 9.0) - K * int((M + 9) / 12.0) + D - 30
	return N

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

def delFile(file):
	#print( '..delFile: %s' % str(glob.glob(file)))
	#print('Remove files:')
	for x in glob.glob(file):
		try:
			os.remove(x)
			#print(x)
		except: 
			print('Error in del file: %s' % str(glob.glob(file)))

def extractTarfile(file, dataDir):
	if(file.find('.tar.gz') != -1):
		try:
			print( 'Try extract gz')
			tar = tarfile.open(file, 'r:gz')
			tar.extractall(dataDir)
			tar.close()
			return True
		except:
			print( 'Error in gz file')
			return False
	if(file.find('.tar') != -1):            
		try:
			print( 'Try extract tar')
			tar = tarfile.open(file, 'r:tar')
			tar.extractall(dataDir)
			tar.close()
			return True
		except:
			print( 'Error in tar file')   
			return False
	return False
		
def infoFile(file):
	if os.path.isfile(file) == 0:
		printlog(False, '..Error! Invalid File: ' + file)
	else:
		ds1 = gdal.Open(file)
		band1 = ds1.GetRasterBand(1)
		# Get table informations
		rows = ds1.RasterYSize
		cols = ds1.RasterXSize
		bands = ds1.RasterCount
		# get georeference info
		transform = ds1.GetGeoTransform()
		xOrigin = transform[0]
		yOrigin = transform[3]
		pixelWidth = transform[1]
		pixelHeight = transform[5]
		print('. rows: ' + str(rows))
		print('. cols :' + str(cols))
		print('. bands: ' + str(bands))
		print('. xOrigin: ' + str(xOrigin))
		print('. yOrigin: ' + str(yOrigin))
		print('. pixelWidth: ' + str(pixelWidth))
		print('. pixelHeight: ' + str(pixelHeight))

def printlog(flag, text):
	timeFile = strftime("%Y%m%d", localtime())
	time = strftime("%Y/%m/%d-%H:%M:%S", localtime())
	logFile = open(path + timeFile + '_logFile.txt', "a+")
	#print(logFile)
	logFile.write(str(text)+'\n')
	if(flag == True): 
		print(text)
	logFile.close()


def mapDict(argv, msg):
	argd = { }

	for i in range(len(argv)):
		if argv[i] == "-fo": # folder output data
			try:
				argd["-fo"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-fi": # folder input data
			try:
				argd["-fi"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-ds": # date start
			try:
				argd["-ds"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-de": # date end
			try:
				argd["-de"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-ty": # data type
			try:
				argd["-ty"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-py": # pixel y position
			try:
				argd["-py"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-px": # pixel x position
			try:
				argd["-px"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-class": # class data
			try:
				argd["-class"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-ps": # type process data
			try:
				argd["-ps"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-hs": # hour start
			try:
				argd["-hs"] = argv[i + 1]
			except:
				exit(msg)
		elif argv[i] == "-he": # hour end
			try:
				argd["-he"] = argv[i + 1]
			except:
				exit(msg)
	return argd 
	 
def sleepTime(time_hours):
	now = datetime.datetime.today()
	future = now + datetime.timedelta(hours=time_hours)
	print( 'Sleep Time: ' + str(int((future - now).seconds/(60*60))) + ' hours')     
	time.sleep((future - now).seconds)

def savgolFilter(DATA, type):
	# Apply Savtisk Golae filter
	lines = 0
	for lines in range(len(DATA)):
		if (numpy.sum(DATA[lines]) == 0):
			break
	lines = lines + 1        
	cols = len(DATA[0])        
	DATA_OUT = numpy.zeros((lines, cols))
	DATA_SAV = numpy.zeros((cols, lines))
	DATA_TMP = numpy.zeros((cols, lines))
	for i in range(cols):
		for j in range(lines):
			DATA_TMP[i][j] = DATA[j][i]
	DATA_SAV[0] = DATA_TMP[0]
	if(type == 'MOD13'):
		for i in range(1, cols):
			#print('DATA_TMP[%d]'%(i))
			#print(DATA_TMP[i])
			DATA_SAV[i] = savgol_filter(DATA_TMP[i], 5, 2, mode ='nearest')
			#print('DATA_SAV[%d]'%(i))
			#print(DATA_SAV[i])
	if(type == 'MSG'):
		DATA_SAV[1] = DATA_TMP[1]
		for i in range(2, cols):
			DATA_SAV[i] = savgol_filter(DATA_TMP[i], 17, 12, mode ='nearest')
	for i in range(lines):
		for j in range(cols):
			DATA_OUT[i][j] = DATA_SAV[j][i]
	return DATA_OUT
