import sys
import subprocess
import os
import glob
import datetime
import subprocess
from sys import argv
import numpy
'''
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
'''
import time
from datetime import datetime
from threading import Timer
import schedule
from time import  strftime, localtime, sleep

# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileMask = fileBrasil

logfile = os.getcwd() + '\\logs\\logfile_goes.log'
#exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
exeDir = 'C:\\Program Files\\GDAL\\'
gdalIlwisDir = "C:\\Ilwis372\\Extensions\\Geonetcast-Toolbox\\GDAL\\bin\\" 
utilDir = 'C:\\ILWIS372\\Extensions\\GEONETCast-Toolbox\\util\\'
tmpDir    = 'c:\\tmp\\'
osgeoDir = 'c:\\osgeo4w\\bin\\'
outputDir = 'Z:\\product\\GOES\\RED\\'
inputDir  = 'Z:\\receive\\GOES-R-CMI-Imagery\\Band02\\'

def delFile(file):
    #print( '..delFile: %s' % str(glob.glob(file)))
    #print('Remove files:')
    for x in glob.glob(file):
        try:
            os.remove(x)
            #print(x)
        except: 
            print('Error in del file: %s' % str(glob.glob(file)))

def printlog(text):
    timeFile = strftime("%Y%m%d", localtime())
    time = strftime("%Y/%m/%d-%H:%M:%S", localtime())
    logFile = open(logfile, "a+")
    #print(logFile)
    logFile.write(str(text)+'\n')
    logFile.close()

def generateGOES(tmpDir, inputDir, outputDir, filename, date):
    printlog('Generate GOES Band02: %s'%(date))
    delFile(tmpDir + '*' + date + '*')
    tmpfile1 = tmpDir + date + '_tmp1.tif'
    tmpfile2 = tmpDir + date + '_tmp2.tif'
    tmpfile3 = tmpDir + date + '_tmp3.tif'
    tmpfile4 = tmpDir + date + '_tmp4.tif'
    outputfile = outputDir + date + '_GOES_RED.tif'

    # Copia apenas a banda do NDVI max
    print(exeDir +  'gdal_translate.exe netCDF:' + filename + '://CMI '  +  tmpfile1)
    subprocess.call(exeDir +  'gdal_translate.exe netCDF:' + filename + '://CMI '  +  tmpfile1)

    # Traduz imagem para coordenadas latitude e longitude
    print(exeDir + 'gdal_translate -a_srs  "+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027489803406172 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs" -a_ullr  -5434894.885056 5434894.885056 5434894.885056 -5434894.885056 ' + tmpfile1 + ' ' + tmpfile2)
    subprocess.call(exeDir + 'gdal_translate -a_srs  "+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027489803406172 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs" -a_ullr  -5434894.885056 5434894.885056 5434894.885056 -5434894.885056 ' + tmpfile1 + ' ' + tmpfile2)

    # Reprojeta imagem para WGS84
    print(exeDir + 'gdalwarp.exe  -s_srs "+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027489803406172 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs" -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"' + tmpfile2 + ' ' + tmpfile3)
    subprocess.call(exeDir + 'gdalwarp.exe  -s_srs "+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027489803406172 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs" -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"' + tmpfile2 + ' ' + tmpfile3)
    
    #Recorta para o Brasil
    print(exeDir + 'gdalwarp.exe  -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outputfile)
    subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outputfile)   
    
    # Remove temporal files
    sys.exit()
    delFile(tmpDir + '*' + date + '*')
    return True

def verifyGenerateFiles():
    printlog('Verify generate files')
    files = glob.glob(inputDir + 'OR_ABI-L*' + '_s201????1500*.nc')
    if(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir)
        sys.exit()
    print('... verify read files')
    for filename in files:
        #print(filename)
        id_date = filename.index('_c201')
        date = filename[id_date+2:id_date+13]
        if(not os.path.isfile(outputDir + date + '_GOES_RED.tif')):
            print('... generate NDVI date: %s'%date)
            if(not generateGOES(tmpDir, inputDir, outputDir, filename, date)):
                printlog('Error to generateGOES [tmpDir: %s, inputDir: %s, filename: %s, date: %s]'%(tmpDir, inputDir, outputDir, filename, date))
    return

class main:

    verifyGenerateFiles()
    schedule.every().day.at("15:00").do(verifyGenerateFiles)
    print('... waiting new files to process')
    while True:
        schedule.run_pending()
        time.sleep(60) # wait one minute
