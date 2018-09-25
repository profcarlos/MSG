import subprocess
from subprocess import Popen, PIPE
import os
import glob
import datetime
import numpy as np
import gc
from subprocess import check_output
from osgeo import gdal, osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import shutil
import sys
import copy
from time import  strftime, localtime, sleep

# Constantes
goiasShape = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'
osgeo4WDir = 'C:\OSGeo4W\\bin\\'
pyDir = "C:\\python34\\"
ilwisAngleDir = "C:\\ILWIS38\\Extensions\\Geonetcast-Toolbox\\toolbox_startscript\\Angle\\"

def command(cmd):
    exec = False
    while(exec != True):
        try:
            subprocess.call(cmd)
            exec = True
        except:
            time.sleep(randint(5, 50))
            exec = False

def delFile(file):
    #printlog('..delFile: %s' % str(glob.glob(file)))
    for x in glob.glob(file):
        os.remove(x)
        
class main():
    
    gc.collect()
    gdal.UseExceptions()

    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')

    dataDir = 'D:\\Dados\\MODIS\\'
    outputDir = os.getcwd() + '\\files\\'
    #file = 'MOD09CMG.A2014234.005.2014236050225.hdf'
    file = 'MOD08_D3.A2015058.006.2015061092031.hdf'
    tempFile = outputDir + 'tempFile.tif'
    outputFile = outputDir + 'modis_teste.tif'
    delFile(tempFile)
    delFile(outputFile)
    filename = dataDir + file
    #band1_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Surface Reflectance Band 1'
    #band2_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Surface Reflectance Band 2'
    #band3_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Ozone'

    band1_sub = 'HDF4_EOS:EOS_GRID:"D:\Dados\MODIS\MOD08_D3.A2015058.006.2015061092031.hdf":mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57
    #band2_sub = 'HDF4_EOS:EOS_GRID:"D:\Dados\MODIS\MOD08_D3.A2015058.006.2015061092031.hdf":mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean'

    band2_sub = 'HDF4_EOS:EOS_GRID:"D:\Dados\MODIS\MOD08_D3.A2015058.006.2015061092031.hdf":mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37
    band3_sub = 'HDF4_EOS:EOS_GRID:"D:\Dados\MODIS\MOD08_D3.A2015058.006.2015061092031.hdf":mod08:Total_Ozone_Mean' # Band 829
    band4_sub = 'HDF4_EOS:EOS_GRID:"D:\Dados\MODIS\MOD08_D3.A2015058.006.2015061092031.hdf":mod08:Atmospheric_Water_Vapor_Mean' # Band 869
        
    Xorigin = -180 
    Yorigin = -90
    Xpixel = 1
    Ypixel = 1
    cols = 180
    rows = 360
    
    # Read data
    ds1 = gdal.Open(filename) 
    subdata = ds1.GetSubDatasets()
    '''
    x = len(subdata)
    y = len(subdata[1])
    for i in range (0, x, 1):
        for j in range (0, y, 1):
            print(subdata[i][j])
    '''
    band1 = gdal.Open(band1_sub)
    band2 = gdal.Open(band2_sub)
    band3 = gdal.Open(band3_sub)
    band4 = gdal.Open(band4_sub)
    #BAND1 = np.int16(band1.ReadAsArray())
    #BAND2 = np.int16(band2.ReadAsArray())
    #BAND3 = np.int16(band3.ReadAsArray())
    #BAND4 = np.int16(band4.ReadAsArray())
    BAND1 = band1.ReadAsArray()
    BAND2 = band2.ReadAsArray()
    BAND3 = band3.ReadAsArray()
    BAND4 = band4.ReadAsArray()
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(tempFile, rows, cols, 5, gdal.GDT_Int16)
    outRaster.SetGeoTransform((Xorigin, Xpixel, 0, Yorigin, 0, Ypixel))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outRaster.SetProjection(srs.ExportToWkt())
    bandOut1 = outRaster.GetRasterBand(1)
    bandOut2 = outRaster.GetRasterBand(2)
    bandOut3 = outRaster.GetRasterBand(3)
    bandOut4 = outRaster.GetRasterBand(4)
    bandOut5 = outRaster.GetRasterBand(5)
    
    bandOut1.WriteArray(BAND1)
    bandOut2.WriteArray(BAND2[1])
    bandOut3.WriteArray(BAND2[2])
    bandOut4.WriteArray(BAND3)
    bandOut5.WriteArray(BAND4)
    
    bandOut1.FlushCache()
    bandOut2.FlushCache()
    bandOut3.FlushCache()
    bandOut4.FlushCache()
    bandOut4.FlushCache()
    outRaster = None
    
    print  (osgeo4WDir + 'gdalwarp -cutline ' + goiasShape + ' -crop_to_cutline ' + tempFile + ' ' + outputFile)
    command(osgeo4WDir + 'gdalwarp -cutline ' + goiasShape + ' -crop_to_cutline ' + tempFile + ' ' + outputFile)
