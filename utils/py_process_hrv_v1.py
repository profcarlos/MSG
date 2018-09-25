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
from Py6S import *
import psutil
import random

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
fileGoiasHRV = os.getcwd() + '\\shape\\201301011700_HRV.tif'


def process_hrv(inputDir, outputDir, date):
    
    row = 612
    col = 632

    row_clm = 204
    col_clm = 211
    HRV   = numpy.zeros((4, row, col), dtype=numpy.float)
    CLM     = numpy.zeros((4, row_clm, col_clm), dtype=numpy.float)
    OUT_HRV = numpy.zeros((row, col), dtype=numpy.float)

    hours = ['1700', '1715', '1730', '1745']
    nhours = 0
    printlog(True, '----- Reading files')
    # Read four files
    for x in range(4):
        filename = inputDir + str(date) + hours[x] + '_HRV.tif'
        if(os.path.isfile(filename) and os.path.isfile(filename.replace('HRV', 'CLM'))):
            #printlog(True, '... Files exist (HRV and CLM): ' + str(date))
            HRV[nhours] = readFileBand(filename, 1)
            CLM[nhours] = readFileBand(filename.replace('HRV', 'CLM'), 1)
            if(len(HRV[nhours]) != [] or len(CLM[nhours]) != [] ):
                nhours = nhours + 1
    VET_HRV = numpy.zeros((row, col, nhours), dtype=numpy.float)
    if(nhours == 0):
        printlog(True, '----- Error to acess data: ' + str(date))
        return False
    
    for i in range(row):
        for j in range(col):
            for x in range(nhours):
                if(CLM[x][int(i/3)][int(j/3)] > 1):
                    VET_HRV[i][j][x] = 0
                else:
                    VET_HRV[i][j][x] = HRV[x][i][j]
            
            VET_HRV[i][j].sort()
            OUT_HRV[i][j] = VET_HRV[i][j][0]
    # Save georeference file
    ds = gdal.Open(fileGoiasHRV, GA_ReadOnly)
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")
    filename = outputDir + str(date) + '_HRV.tif'
    dsOut  = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 1, band.DataType)
    CopyDatasetInfo(ds,dsOut)

    bandOut=dsOut.GetRasterBand(1)
    BandWriteArray(bandOut, OUT_HRV)


class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    gdal.UseExceptions()    
    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict:
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
    else:
        exit(usage)
    date = dateStart
    while(date < dateEnd):
        inputsubDir = inputDir + 'HRV\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        outputsubDir = outputDir  + 'HRV\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        if(process_hrv(inputsubDir, outputsubDir, date) == 0):
            printlog(True, "----- Error to generate HRV: %d"%(date))
        else:
            printlog(True, "----- Generate HRV: %d"%(date))
        date = incDate(date)

