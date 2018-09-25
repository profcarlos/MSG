import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

from sys import argv
from utils import *
import subprocess
import os
import numpy
from scipy.special import factorial

csv_file   = 'C:\\Users\\carlos.silveira\\Dropbox\\Gabriel\\pontos_coord_v2.csv'
input_dir  = 't:\\modis\\2km\\'


def convert (x_input, y_input):

    x_orig = -53.250722730884448
    y_orig = -12.394785590574308

    x_size = 0.018000436173230
    y_size = -0.018030217830345

    x_sample = int(abs((x_orig - x_input)/x_size))
    y_sample = int(abs((y_orig - y_input)/y_size))

    return(y_sample, x_sample)

class main():
  

    COORD = numpy.loadtxt(csv_file, delimiter = ',', skiprows = 1)
    DOT = numpy.zeros((len(COORD), len(COORD[0])))
    for x in range(len(DOT)):
        DOT[x][0] = COORD[x][0]
        (DOT[x][1], DOT[x][2]) =  convert(COORD[x][1], COORD[x][2])
    numpy.savetxt(input_dir + 'dot_file.csv', DOT, fmt='%d', comments='', header = "ID X Y", delimiter = ',', newline='\n')

    DATA_OUT = datafile('', input_dir, '_ndvi_fpar' + '_.csv', "ID LAT LON YEAR DOY NDVI FPAR_1 FPAR_2", len(DOT)*(365*(2015-2000)/8), 8)
    for year in range(2000, 2016, 1):
        for doy in range(1, 366, 16):
            # Read MODIS data files
            file_MOD13 = input_dir + '\\MOD13\\' + str(year) + str('%03d'%doy) + '_MOD13.tif'
            if(os.path.isfile(file_MOD13)):
                NDVI = readFileBand(file_MOD13, 1)
            else:
                print('Error to read file: %s'%(file_MOD13))
                continue
            file_MOD15_1 = input_dir + '\\MOD15\\' + str(year) + str('%03d'%doy) + '_MOD15.tif'
            if(os.path.isfile(file_MOD15_1)):
                FPAR_1 = readFileBand(file_MOD15_1, 1)
            else:
                print('Error to read file: %s'%(file_MOD15_1))
                continue

            file_MOD15_2 = input_dir + '\\MOD15\\' + str(year) + str('%03d'%(doy-8)) + '_MOD15.tif'
            if(os.path.isfile(file_MOD15_2)):
                FPAR_2 = readFileBand(file_MOD15_2, 1)
            else:
                print('Error to read file: %s'%(file_MOD15_2))
                continue
            print('Save data file: %s'%(str(year) + str('%03d'%doy)))
            # Get data for all dot    
            for x in range(len(DOT)):
                data_out = (DOT[x][0], COORD[x][1], COORD[x][2], year, doy, NDVI[DOT[x][1], DOT[x][2]], FPAR_1[DOT[x][1], DOT[x][2]], FPAR_2[DOT[x][1], DOT[x][2]])
                #print(data_out)
                DATA_OUT.addAtrib(data_out)
                DATA_OUT.addSample()
        # Save data file
        DATA_OUT.save()
    	







