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
