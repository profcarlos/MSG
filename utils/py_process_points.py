import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 't:\\carlos\\newPython\\utils')

import subprocess
import os
import glob
import datetime
import gc
import subprocess
import shutil
import copy
from utils import *
from sys import argv
import numpy
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import time
import pandas as pd


usage = """\
Usage: %s [OPTIONS]
        -ty     type data (MSG, MODIS) in file name
        -fi     folder input data
        -fo     folder output data
""" % argv[0]


class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    gdal.UseExceptions()    
    if ("-fo" in argDict and "-fi" in argDict):
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData  = argDict["-ty"] 
    else:
        exit(usage)
    dateStart = 20130101
    dateEnd = 20131231
    files = glob.glob(inputDir + '*' + typeData + '*.csv')
    print (files)
    if(len(files) == 0):
        print('Verify, zero files.')
        sys.exit()
    for file in files:
        DATA = pd.read_csv(file, sep = ",", index_col = False, usecols = ['DATE', 'NDVI'])
        DATA_1km = DATA[DATA['DATE'].str.contains("NDVI")].sort(ascending = False)
        DATA_3km = DATA[DATA['DATE'].str.contains("BRDF")].sort(ascending = False)

        date_1km = DATA_1km['DATE'].to_string().replace('_NDVI','')
        date_1km = date_1km.split(' ')
        ndvi_1km = DATA_1km['NDVI'].tolist()
        date_3km = DATA_3km['DATE'].to_string().replace('_BRDF','')
        date_3km = date_3km.split(' ') 
        ndvi_3km = DATA_3km['NDVI'].tolist()
        if(len(date_1km) != len(date_3km)):
            print('Verify, not equal len to 1km and 3km data.')
            print(date_1km)
            print(date_3km)
            sys.exit()
        samples = len(DATA_1km['DATE'].tolist())
        print('samples: %d'%(samples))
        SAMPLES = numpy.zeros((samples, 3))
        date = dateStart
        n = 0
        while(date < dateEnd):
            doy = calcDoy_date(date)
            str_date = str(date)
            if(str_date in date_1km or str_date in date_3km):
                SAMPLES[n][0] = doy
                if(str_date in date_1km):
                    id = date_1km.find(str_date)
                    print(id)
                    if(id != 0):
                        SAMPLES[n][1] = ndvi_1km[id]
                if(str_date in date_3km):
                    id = date_3km.find(str_date)
                    print(id)
                    if(id != 0):
                        SAMPLES[n][2] = ndvi_3km[id]
                        print(SAMPLES[n])
                n = n + 1
            date = incDate(date)
            if(date > 20130107):
                sys.exit()
