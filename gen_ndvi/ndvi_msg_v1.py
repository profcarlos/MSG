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

# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileMask = fileBrasil

filelog = os.getcwd() + '\\logs\\logfile.log'
exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
gdalIlwisDir = "C:\\Ilwis372\\Extensions\\Geonetcast-Toolbox\\GDAL\\bin\\" 
utilDir = 'C:\\ILWIS372\\Extensions\\GEONETCast-Toolbox\\util\\'
ilwDir = 'C:\\ilwis372\\'
tmpDir    = 'c:\\tmp\\'
inputDir  = 'Z:\\receive\\MSG-0degree\\MetProducts\\'
outputDir = 'Z:\\product\\NDVI\\'

def delFile(file):
    #print( '..delFile: %s' % str(glob.glob(file)))
    #print('Remove files:')
    for x in glob.glob(file):
        try:
            os.remove(x)
            #print(x)
        except: 
            print('Error in del file: %s' % str(glob.glob(file)))

def generateNDVI(tmpDir, inputDir, outputDir, filename, date):

    delFile(tmpDir + '*' + date + '*')
    tmpfile1 = filename[filename.index('L-000'):] + '.xyz'
    tmpfile2 = tmpDir + date + '_tmp2.hdf'
    tmpfile3 = tmpDir + date + '_tmp3'
    tmpfile4 = tmpDir + date + '_tmp4.mpr'
    tmpfile5 = tmpDir + date + '_tmp5.tif'
    tmpfile6 = tmpDir + date + '_tmp6.tif'
    outputfile = outputDir + date + '_NDVI.tif'

    print(utilDir +  'joinmsg.exe ' + filename + ' '  +  tmpDir)
    subprocess.call(utilDir +  'joinmsg.exe ' + filename + ' '  +  tmpDir)

    print('rename' +  ' ' + tmpDir + tmpfile1 + ' ' + tmpfile2)
    os.rename(tmpDir + tmpfile1, tmpfile2)

    print(gdalIlwisDir +  'gdal_translate.exe -of ILWIS HDF5:' + tmpfile2 + '://NDVImax '  +  tmpfile3)
    subprocess.call(gdalIlwisDir +  'gdal_translate.exe -of ILWIS HDF5:' + tmpfile2 + '://NDVImax '  +  tmpfile3)

    print(ilwDir + 'ilwis.exe -C ' + tmpfile4 + ':=iff(' + tmpfile3 +' le 100,' + tmpfile3 + '/100,?);')
    subprocess.call(ilwDir + 'ilwis.exe -C ' + tmpfile4 +':=iff(' + tmpfile3 +' le 100,' + tmpfile3 + '/100,?);')

    # Config data file
    print(ilwDir + 'ilwis.exe -C setgrf ' + tmpfile4 + ' ' + utilDir +'lritmsg;')
    subprocess.call(ilwDir + 'ilwis.exe -C setgrf ' + tmpfile4 + ' ' + utilDir +'lritmsg;')

    print(exeDir + 'gdal_translate.exe HDF5:' + tmpfile2 + '://NDVImax '  +  tmpfile3)
    subprocess.call(exeDir + 'gdal_translate.exe HDF5:' + tmpfile2 + '://NDVImax '  +  tmpfile3)

    # Traduz imagem para coordenadas latitude e longitude
    print(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5568000 5568000 5568000 -5568000 ' + tmpfile4 + ' ' + tmpfile5)
    subprocess.call(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5568000 5568000 5568000 -5568000 ' + tmpfile4  + ' ' + tmpfile5)

    # Reprojeta imagem para WGS84 near/bilinear/cubic
    print(exeDir + 'gdalwarp.exe -dstnodata -3000 -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask +  ' ' + tmpfile5 + ' ' + outputfile)
    subprocess.call(exeDir + 'gdalwarp.exe -dstnodata -3000 -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask +  ' ' + tmpfile5 + ' ' + outputfile)
    
    # Remove temporal files
    sys.exit()
    delFile(tmpDir + '*' + date + '*')
    return True

class main:
    files = glob.glob(inputDir + '*L-000-MSG3__-MPEF________-NDVI_____-000001___-*')
    if(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir + file)
        sys.exit()
    print('... verify read files')
    for filename in files:
        print(filename)
        id_date = filename.index('-201')
        date = filename[id_date+1:id_date+9]
        if(not os.path.isfile(outputDir + date + '.tif')):
            print('... generate NDVI date: %s'%date)
            generateNDVI(tmpDir, inputDir, outputDir, filename, date)
            sys.exit()


