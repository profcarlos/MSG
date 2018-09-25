'''
Coordinate System is:
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0],
    UNIT["degree",0.0174532925199433],
    AUTHORITY["EPSG","4326"]]
Origin = (-73.991975501456054,5.273674044907336)
Pixel Size = (0.002091294669953,-0.002091294669953)

Corner Coordinates:
Upper Left  ( -73.9919755,   5.2736740) ( 73d59'31.11"W,  5d16'25.23"N)
Lower Left  ( -73.9919755, -33.7498845) ( 73d59'31.11"W, 33d44'59.58"S)
Upper Right ( -28.8367410,   5.2736740) ( 28d50'12.27"W,  5d16'25.23"N)
Lower Right ( -28.8367410, -33.7498845) ( 28d50'12.27"W, 33d44'59.58"S)
Center      ( -51.4143582, -14.2381052) ( 51d24'51.69"W, 14d14'17.18"S)
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

usage = """\
Usage: %s [OPTIONS]
        -ps     file to processe data (.csv). If not used then use my file in 'C:\\Users\\carlos.silveira\\Dropbox\\Gabriel\\bfast_results.csv'.
        -fi     folder input data (used y:\)
        -fo     folder output data
        -ds     date to start process (format: yyyy|doy)
        -de     date to end process (format: yyyy|doy)

""" % argv[0]

csv_file   = 'C:\\Users\\carlos.silveira\\Dropbox\\Gabriel\\bfast_results_v3.csv'
modis_file = ['pa_br_ndvi_250_lapig', 'pa_br_fpar_1000_lapig'] 
class main():
  
    argDict = mapDict(argv, usage)

    if("-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict):

        yearStart   = int(argDict["-ds"][0:4])
        yearFinish  = int(argDict["-de"][0:4])
        doyStart  = int(argDict["-ds"][4:6])
        doyFinish = int(argDict["-de"][4:6])
		inputDir = argDict["-fi"]
        outputDir = argDict["-fo"]
        dateStart = int(argDict["-ds"])
        dateEnd   = int(argDict["-de"])            
    else:
        sys.exit(usage)

    if("-ps" in argDict):
    	processData = str(argDict["-ps"])
    else:
    	processData = csv_file
    try:
    	DATA_IN  = numpy.loadtxt(processData, dtype='float32', delimiter = ';', skiprows = 1)
    	DATA_OUT = numpy.zeros((6))
    except:
    	print('...Error reading csv file.')
    	sys.exit()

    if(typeData == )


    if(typeData == 'CHANGE'):
        y_an = 0
        x_an = 0
        date_change = 0
        for row in range(len(DATA_IN[0])):
        	[y x lon lat date mag trend_before trend_after slop_before slop_after] = DATA_IN[row]
        	if(y =! y_an or x =! x_an):
                DATA_OUT = [y x lon lat date_change mag_change]
                date_change = 200001
        		y_an = y
        		x_an = x
            # Test if has change in ground use
            if(mag > 0.3):
                date_change = date
                mag_change = mag
        numpy.savetxt('c:\\atmfile_hist.csv', DATA_OUT, fmt='%1.4f', comments='', header = "Y X LON LAT DATE MAG", delimiter = ' ', newline='\n')



    	







