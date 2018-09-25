import sys
import subprocess
import os
import glob
import datetime
import subprocess
from sys import argv
import numpy
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import time
from datetime import datetime
from threading import Timer
import schedule
from time import  strftime, localtime, sleep

# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileMask = fileBrasil

logfile = os.getcwd() + '\\logs\\logfile.log'
exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
gdalIlwisDir = "C:\\Ilwis372\\Extensions\\Geonetcast-Toolbox\\GDAL\\bin\\" 
utilDir = 'C:\\ILWIS372\\Extensions\\GEONETCast-Toolbox\\util\\'
tmpDir    = 'c:\\tmp\\'
osgeoDir = 'c:\\osgeo4w\\bin\\'
inputDir  = 'Z:\\receive\\MSG-0degree\\MetProducts\\'
outputDir = 'Z:\\product\\MSG\\NDVI\\'

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

def generateNDVI(tmpDir, inputDir, outputDir, filename, date):
    printlog('Generate NDVI of: %s'%(date))
    delFile(tmpDir + '*' + date + '*')
    tmpfile1 = filename[filename.index('L-000'):] + '.xyz'
    tmpfile2 = tmpDir + date + '_tmp2.hdf'
    tmpfile3 = tmpDir + date + '_tmp3.tif'
    tmpfile4 = tmpDir + date + '_tmp4.tif'
    tmpfile5 = tmpDir + date + '_tmp5.tif'
    tmpfile6 = tmpDir + date + '_tmp6.tif'
    outputfile = outputDir + date + '_NDVI.tif'
    # Junta dados dos arquivos (PRO e 000001)
    print(utilDir +  'joinmsg.exe ' + filename + ' '  +  tmpDir)
    subprocess.call(utilDir +  'joinmsg.exe ' + filename + ' '  +  tmpDir)
    # Renomeia arquivo hdf
    print('rename' +  ' ' + tmpDir + tmpfile1 + ' ' + tmpfile2)
    os.rename(tmpDir + tmpfile1, tmpfile2)
    # Copia apenas a banda do NDVI max
    print(gdalIlwisDir +  'gdal_translate.exe HDF5:' + tmpfile2 + '://NDVImax '  +  tmpfile3)
    subprocess.call(gdalIlwisDir +  'gdal_translate.exe HDF5:' + tmpfile2 + '://NDVImax '  +  tmpfile3)

    # Traduz imagem para coordenadas latitude e longitude
    print(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr  -5568000 5568000 5568000 -5568000 ' + tmpfile3 + ' ' + tmpfile4)
    subprocess.call(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr  -5568000 5568000 5568000 -5568000 ' + tmpfile3  + ' ' + tmpfile4)

    # Reprojeta imagem para WGS84
    print(exeDir + 'gdalwarp.exe -srcnodata 255 -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + tmpfile5)
    subprocess.call(exeDir + 'gdalwarp.exe  -srcnodata 255 -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + tmpfile5)

    # Converte para escala [0 1]
    print('python ' + osgeoDir + 'gdal_calc.py -A ' + tmpfile5 + ' --outfile=' + tmpfile6 + ' --type=Float32 --calc="0.01*A*(A<255)" --NoDataValue=0')
    subprocess.call('python ' + osgeoDir + 'gdal_calc.py -A ' + tmpfile5 + ' --outfile=' + tmpfile6 + ' --type=Float32 --calc="0.01*A*(A<255)" --NoDataValue=0') 
    
    #Recorta para o Brasil
    print(exeDir + 'gdalwarp.exe  -cutline ' + fileMask +  ' ' + tmpfile6 + ' ' + outputfile)
    subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile6 + ' ' + outputfile)   
    
    # Remove temporal files
    #sys.exit()
    delFile(tmpDir + '*' + date + '*')
    return True

def verifyGenerateFiles():
    printlog('Verify generate files')
    files = glob.glob(inputDir + '*L-000-MSG3__-MPEF________-NDVI_____-000001___-*')
    if(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir + file)
        sys.exit()
    print('... verify read files')
    for filename in files:
        print(filename)
        id_date = filename.index('-201')
        date = filename[id_date+1:id_date+9]
        if(not os.path.isfile(outputDir + date + '_NDVI.tif')):
            print('... generate NDVI date: %s'%date)
            if(not generateNDVI(tmpDir, inputDir, outputDir, filename, date)):
                printlog('Error to generateNDVI [tmpDir: %s, inputDir: %s, filename: %s, date: %s]'%(tmpDir, inputDir, outputDir, filename, date))
        else:
            print('... file exist date: %s'%date)
    return

class main:
    print('... verify files to generate')
    verifyGenerateFiles()
    schedule.every().day.at("23:00").do(verifyGenerateFiles)
    print('... waiting new files to process')
    while True:
        schedule.run_pending()
        time.sleep(60) # wait one minute
