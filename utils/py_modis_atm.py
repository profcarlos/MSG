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

