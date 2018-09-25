import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import subprocess
from subprocess import Popen, PIPE
import os
import glob
import datetime
import numpy
import gc
import shutil
import sys
import copy
import gdal
from osgeo import gdal, osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from time import  strftime, localtime, sleep
from utils import printlog, sleepTime, mapDict, sleepTime, delFile, infoFile, incHour, incDate
from osgeo import gdal, osr
import math
import copy

from sys import argv

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to save output files (optional, default is \\files\\)
        -ds     date to start process (format: yyyymmdd)
        -de     date to end process   (format: yyyymmdd) 
        -hs     hour to start process (format: hhMM)
        -he     hour to end process   (format: hhMM)
""" % argv[0]

# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMask = fileGoias
angleDir = os.getcwd() +  '\\angle\\'
#exeDir = 'C:\\Python34\\Lib\\site-packages\osgeo\\'
exeDir = 'C:\\OSGeo4W\\bin\\'
outputDir = os.getcwd() + '\\file\\'
tmpDir = 'C:\\tmp\\'



def timeAngle(day, hour, form):
    day = str(day)
    hour = str(hour)
    minuto = 10*int(hour[2])+int(hour[3]) 
    if(minuto == 00):
        hour = hour[:2]+'.00'
    else:
        if(minuto <= 15):
            hour = hour[:2]+'.15'
        else:
            if(minuto <= 30):    
                hour = hour[:2]+'.30'
            else:
                hour = hour[:2]+'.45'   
    
    if(form == 'SPACE'):
       string = day[:4]+' '+day[4]+day[5]+' '+day[6]+day[7]+' '+ hour
    if(form == 'UNDER'):
       string = '_' + day[:4]+'_'+day[4]+day[5]+'_'+day[6]+day[7]+'_'+ hour
    return string

    
def array2raster(outFile,array):

    #Xorigin = -66.96422963471274   # Oeste 
    #Xend    = +66.96364537142968   # Leste
    #Yend    = -62.27591109399421   # South 
    #Yorigin = +62.27979992238434   # North
    #Xpixel  = 0.034840758326259734
    #Ypixel  = 0.034840758326259734

    Xpixel = 0.45
    Ypixel = 0.45
 
    cols = 378
    rows = 378

    Xorigin = -85 #- Xpixel/2
    Yorigin = -85 #- Ypixel/2

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outFile, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((Xorigin, Xpixel, 0, Yorigin, 0, Ypixel))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    # this assumes the projection is Geographic lat/lon WGS 84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outRaster.SetProjection(srs.ExportToWkt())
    outband.FlushCache()
    outRaster = None
    infoFile(outFile)

def transformAngle(outputDir, day, hour, vetAngle):

    timeprocess = str(day) + str(hour)
    tmpfile1 = tmpDir + timeprocess + '_tmp1.tif'
    tmpfile2 = tmpDir + timeprocess + '_tmp2.tif'
    tmpfile3 = tmpDir + timeprocess + '_tmp3.tif'
    tmpfile4 = tmpDir + timeprocess + '_tmp4.tif'
    tmpfile5 = tmpDir + timeprocess + '_tmp5.tif'
    fileIn   = vetAngle[0]
    fileOut  = vetAngle[1]
    minValue = vetAngle[2]
    maxValue = vetAngle[3]
    resampling_method = vetAngle[4]
       
    print(True,'transform Angle: ' + fileIn)

    arrayTemp = numpy.fromfile(fileIn, dtype='>f')
    matrixTemp = numpy.reshape(arrayTemp,(-1,378))
    rows = 378
    cols = 378
    matrix = numpy.zeros((rows, cols), dtype = float)
    for i in range (0, rows, 1):
        for j in range (0, cols, 1):
            temp = float("{0:.5f}".format(matrixTemp[i][j]))
            if ('sunzen' in fileIn) and math.isnan(temp):
                #print(False,'...none: ' + str(temp))
                temp = 90.0001
            if temp > maxValue:
                printlog(False, '...maior: ' + str(temp))
                if ('sunzen' in fileIn): 
                    temp = 90.0000
                else:
                    temp = maxValue
            if temp < minValue:
                printlog(False, '...menor:'  + str(temp))
                temp = minValue
            matrix[i][j]=temp
    del arrayTemp
    del matrixTemp
    array2raster(tmpfile1, matrix)
    
    print(exeDir + 'gdalwarp -tr 0.034840758326259734 0.034840758326259734 -r  bilinear ' + ' ' + tmpfile1 + ' ' + tmpfile2)
    subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326259734 0.034840758326259734 -r  bilinear ' + ' ' + tmpfile1 + ' ' + tmpfile2)

    print(exeDir + 'gdalwarp -cutline ' + fileDataRetriever + ' -crop_to_cutline -overwrite '  + tmpfile2 + ' ' + tmpfile3)
    subprocess.call(exeDir + 'gdalwarp -cutline ' + fileDataRetriever + ' -crop_to_cutline -overwrite '  + tmpfile2 + ' ' + tmpfile3)
    
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)

    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + fileOut)
    subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + fileOut)

    #print(fileOut)
    #sys.exit()
    delFile(tmpDir + timeprocess + '_tmp' + '*')


def create1file(inputFile, outputFile):
    # Read input and output data
    printlog(False, '...create1file:')
    printlog(False, '..inputFile: %s' % inputFile)
    try:
        dsIn = gdal.Open(inputFile[0])
    except(RuntimeError, e):
        raise RuntimeError('...Error: file not found!')
    bandIn = dsIn.GetRasterBand(1)
    Nbands = len(inputFile)
    printlog(False, '..outputFile: %s' % outputFile)
    driver = gdal.GetDriverByName("GTiff")
    delFile(outputFile)
    dsOut  = driver.Create(outputFile, dsIn.RasterXSize, dsIn.RasterYSize, Nbands, bandIn.DataType)
    CopyDatasetInfo(dsIn,dsOut)
    # Write output file
    for x in range(0, len(inputFile), 1):
        try:
            dsIn = gdal.Open(inputFile[x])
        except(RuntimeError, e):
            raise RuntimeError('...Error: file %s not found!' % inputFile[x])
        bandIn = dsIn.GetRasterBand(1)
        arrayIn = BandReadAsArray(bandIn)
        bandOut = dsOut.GetRasterBand(x+1)
        BandWriteArray(bandOut, arrayIn)
        arrayIn = None
        bandIn = None
        dsIn = None
        delFile(inputFile[x])
        
class main():

    argDict = mapDict(argv, usage)

    if "-ds" in argDict and "-de" in argDict and "-hs" in argDict and "-he" in argDict:

        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        hourStart = int(argDict["-hs"])
        hourEnd = int(argDict["-he"])
        if "-fo" in argDict:
            outputDir = argDict["-fo"]
    else:
        exit(usage)
    gc.collect()
    gdal.UseExceptions()
    #printlog(True, 'os.getcwd() : ' + os.getcwd())
    inputFile = ['test', 'test', 'test', 'test']
    arrayAngle =  [['sataz','msg_azres',0, 360, 'average'], ['sunaz','sol_azres',0, 360, 'average'], ['satzen', 'msg_zenres', 0, 99, 'average'], ['sunzen', 'sol_zenres', 0, 98.9998, 'average']]
    day = dateStart
    
    while (day <= dateEnd):
        year  = str(day)[0:4]
        month = str(day)[4:6]
        subfolder = year + '\\' + month + '\\'
        if not(os.path.isdir(outputDir + subfolder)):
            os.makedirs(outputDir + subfolder)
        hour = hourStart
        while(hour <= hourEnd):
            outputFile = outputDir + subfolder + str(day) + str(hour) + '_angles.tif'
            if os.path.isfile(outputFile):
                printlog(True,'File exist: ' + str(day) + str(hour) + '_angles.tif')
                hour = incHour(hour)
                continue
            arrayAngle =  [['sataz','msg_azres',0, 360, 'average'], ['sunaz','sol_azres',0, 360, 'average'], ['satzen', 'msg_zenres', 0, 99, 'average'], ['sunzen', 'sol_zenres', 0, 98.9998, 'average']]
            matAngle = copy.copy(arrayAngle)
            #printlog(False, 'matAngle:')
            #printlog(False, matAngle)
            #printlog(False, 'arrayAngle:')
            #printlog(False, arrayAngle)

            for x in range(0, len(matAngle), 1):
                delFile(os.getcwd()  + '\\' + matAngle[x][0] + '_*')
                matAngle[x][0] += timeAngle(day, hour,'UNDER')
                matAngle[x][1] = str(day) + str(hour) + '_' + matAngle[x][1] + '.tif'
                inputFile[x] = os.getcwd() + '\\' + matAngle[x][1]
            
            #ans = check_output('java.exe -cp .;C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle\\operation.jar;C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle\\Jama-1.0.1.jar;C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle AngleMaps '+ timeAngle(day, hour,'SPACE'))   
            printlog(False, 'java -cp .;' + angleDir + 'operation.jar;' + angleDir  + 'Jama-1.0.1.jar;' + os.getcwd() +'\\angle AngleMaps ' +  timeAngle(day, hour,'SPACE'))
            try:
                subprocess.call('java -cp .;' + angleDir + 'operation.jar;' + angleDir  + 'Jama-1.0.1.jar;' + os.getcwd() +'\\angle AngleMaps ' +  timeAngle(day, hour,'SPACE'))
            except:
                printlog(True, '..Fail Generate Azimuth Files')
                sys.exit()
            for x in range(0, len(matAngle), 1):
                transformAngle(outputDir , day, hour, matAngle[x])
            printlog(False, 'inputFile: %s'  % inputFile)
            printlog(False, 'outputFile: %s' % outputFile)
            create1file(inputFile, outputFile)
            #infoFile(outputFile)
            for x in range(0, len(matAngle), 1):
                delFile(inputFile[x])
            hour = incHour(hour)
        day = incDate(day)


