import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import sys
import subprocess
import os
import glob
from utils import *

outputDir = 't:\\OUTROS\\STAT_GAB\\'
stat_csv  = outputDir + 'stations_goias.csv'
stat_data = outputDir + 'estacaoes_par_efic.csv'
grid_stations ='C:\\Users\\carlos.silveira\\Dropbox\\newPython\\kriging\\grid_stations_gab.csv'
exeDir = 'C:\\osgeo4w\\bin\\'
modiscutfile = outputDir + 'mod_goias.shp'
fileMask  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
'''
Driver: GTiff/GeoTIFF
Files: T:\MODIS\FPAR_calc\2014337_FPAR.tif
Size is 3858, 3730
Coordinate System is:
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0],
    UNIT["degree",0.0174532925199433],
    AUTHORITY["EPSG","4326"]]
Origin = (-53.252403461635794,-12.393279067399060)
Pixel Size = (0.001905476965202,-0.001905476965202)
Metadata:
  AREA_OR_POINT=Area
Image Structure Metadata:
  INTERLEAVE=BAND
Corner Coordinates:
Upper Left  ( -53.2524035, -12.3932791) ( 53d15' 8.65"W, 12d23'35.80"S)
Lower Left  ( -53.2524035, -19.5007081) ( 53d15' 8.65"W, 19d30' 2.55"S)
Upper Right ( -45.9010733, -12.3932791) ( 45d54' 3.86"W, 12d23'35.80"S)
Lower Right ( -45.9010733, -19.5007081) ( 45d54' 3.86"W, 19d30' 2.55"S)
Center      ( -49.5767384, -15.9469936) ( 49d34'36.26"W, 15d56'49.18"S)
Band 1 Block=3858x1 Type=Float32, ColorInterp=Gray
  NoData Value=0
'''  
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
	#return D
	return M, D 

class main():
	#STATIONS = [NAME, LONG, LAT]
	STATIONS = numpy.loadtxt(stat_csv, delimiter = ',', skiprows = 1)
	n_stations = len(STATIONS)
	VARS = ['STAT','YEAR','MONTH','DAY','TMIN','VPD','PAR','EFIC']
	DATA = numpy.genfromtxt(stat_data, delimiter = ';', skip_header = 1,unpack=True)
	for i in range(len(VARS)):
		globals()[VARS[i]] = DATA[i]
	year = 2014
	for doy in range(1,365,16):
		print('...process doy: %d\n\n'%(doy))
		outfile = outputDir + 'PAR\\' + str(year) + str('%03d'%doy) + '_PAR.tif'#'_EFIC.tif' #
		print('outfile: %s'%(outfile))
		'''
		if(os.path.isfile(outfile)):
			#print('file exist: %s'%(outfile))
			TEST = readFileBand(outfile, 1)
			has_nan = numpy.any(numpy.isnan(TEST))
			if(has_nan):
				os.remove(outfile)
				print('...REMOVING [has nan] file %s'%(outfile))
		'''
		if(not os.path.isfile(outfile)):
			tmpfile1 = outfile.replace('.tif', '_tmpfile1.tif')
			tmpfile2 = outfile.replace('.tif', '_tmpfile2.tif')
			tmpfile3 = outfile.replace('.tif', '_tmpfile3.tif')

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
			file.write('LONG,LAT,DATA\n')
			if(os.path.isfile(outfile)):
				os.remove(outfile)
			n_error = 0
			for i, stat in zip(range(n_stations), numpy.transpose(STATIONS)[0]):
				#print('search in %d of %d %d %d'%(stat, year, month, day))
				#print(STAT[0])
				#print(YEAR[0])
				#print(MON[0])
				month, day = calcDay(year, doy)
				month_up, day_up = calcDay(year, doy+16)
				paser = numpy.logical_or(numpy.logical_and(numpy.logical_and(YEAR==year, MONTH<=month_up), DAY<=day_up), numpy.logical_and(numpy.logical_and(YEAR==year, MONTH>=month), DAY>=day))
				paser = numpy.logical_and(paser,  STAT==stat)
				paser = numpy.logical_and(paser, numpy.logical_and(PAR>0, PAR<300))#(EFIC>0, EFIC<=1))
				if(numpy.sum(paser) == 0):
					n_error = n_error + 1
					#print('error...')
				else:
					# Get data in array

					index = numpy.mean(PAR[paser])#(EFIC[paser])
					if(not (numpy.isnan(index) or index < 0)):
						stat_data = str('%.6f,%.6f,%.4f'%(STATIONS[i][1], STATIONS[i][2], index))
						print(stat_data)
						if(i != n_stations-1):
							stat_data = stat_data + '\n'
						file.write(stat_data)
					else:
						n_error = n_error + 1
			file.close()
			if(n_error > 0.1*n_stations):
				print('...error to generate file !! [n_error: %d]'%(n_error))
				continue
			print('...process create image')
			print(exeDir + 'gdal_grid -a invdist:power=2.0:smoothing=1.0 -zfield DATA -l grid_stations_gab  grid_stations_gab.vrt ' + ' ' + tmpfile1)
			subprocess.call(exeDir + 'gdal_grid -a invdist:power=2.0:smoothing=1.0 -zfield DATA -l grid_stations_gab grid_stations_gab.vrt '+ ' ' + tmpfile1)

			#print(exeDir + 'gdalwarp.exe -t_srs "+proj=latlong +datum=WGS84" ' + tmpfile1 + ' ' + tmpfile2)
			#subprocess.call(exeDir + 'gdalwarp.exe -t_srs "+proj=latlong +datum=WGS84" ' + tmpfile1 + ' ' + tmpfile2)
			print('...process tmp2')			
			print(exeDir + 'gdalwarp -tr 0.001905476965202 0.001905476965202 -r "cubic" -cutline ' + modiscutfile + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
			subprocess.call(exeDir + 'gdalwarp -tr 0.001905476965202 0.001905476965202 -r "bilinear" -cutline ' + modiscutfile + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
			print('...proces tmp3')
			print(exeDir + 'gdal_translate.exe -projwin -53.2524035 -12.3932791 -45.9010733 -19.5007081 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
			subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2524035 -12.3932791 -45.9010733 -19.5007081 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
			print('...process outfile')
			print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outfile)
			subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outfile)
			#sys.exit()
			delFile(tmpfile1)
			delFile(tmpfile2)
			delFile(tmpfile3)
			#sys.exit()
