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

class main():
	#STATIONS = [NAME, LONG, LAT]
	STATIONS = numpy.loadtxt(stat_csv, delimiter = ',', skiprows = 1)
	n_stations = len(STATIONS)
	file_station = 'z:\OUTROS\PLUVIOSIDADE\\BCK_ANO_2015.txt'
	print('. reading file: %s'%(file_station))
	with open(file_station, 'rb') as file_read:
		last_name_station = []
		for line in file_read:
			#print('line : %s'%line)
			line = str(line)
			line = line.replace('\\n', '')
			line = line.replace('b''', '')
			line = line.replace('\\r', '')
			#print('line : %s'%line)
			name_station = line[1:5]
			
			if (name_station != last_name_station):
				print('name_station: %s'%(name_station))
				if(not last_name_station == []):
					file_out.close()
				output_filename = inputDir + '\\BCK_ANO_2015\\'+ 'BCK_ANO_2015_EST_' + name_station + '.txt' 
				file_out = open(output_filename, 'a')

			file_out.write(line[1:len(line)-1] + '\n' )
			last_name_station = name_station
		else:
			print('EOF')
			sys.exit()




	