import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import subprocess
import os
import glob
import datetime
import gc
import subprocess
import shutil
import copy
import tarfile
from utils import *
from sys import argv
import numpy
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import threading
import queue
import time
import shutil

usage = """\
Usage: %s [OPTIONS]
        -ty     type of data to process (123, CLM, NDVI, MOD08, MYD08, MOD09, HRV, MOD13)
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)

""" % argv[0]


# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoiasGeo = os.getcwd() + '\\shapes\\201401011300_123.tif'
fileMaskRect   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMask = fileGoias

pathlog = os.getcwd() + '\\logs\\'
MRTDir =  os.getcwd() + '\\MRT\\'
dataRetrieve = os.getcwd() + '\\MSGDataRetriever\\'
pyModisDir = os.getcwd() + '\\\pyModis-2.0.5\\scripts'
exeDir = 'C:\\Python34\\Lib\\site-packages\\osgeo\\'
gdalIlwisDir = "C:\\Ilwis372\\Extensions\\Geonetcast-Toolbox\\GDAL\\bin\\" 
utilDir = 'C:\\ILWIS372\\Extensions\\GEONETCast-Toolbox\\util\\'
ilwDir = 'C:\\ilwis372\\'
osgeoDir = 'C:\\OSGeo4W\\bin\\'

tmpDir    = 'c:\\tmp\\'
logTarError = tmpDir + 'logTarError.txt'

threadList = ["Thread-1", "Thread-2", "Thread-3", "Thread-4", "Thread-5", "Thread-6"]
queueLock = threading.Lock()
workQueue = queue.Queue(12)
threads = []
exitFlag = 0

class myThread (threading.Thread):
    def __init__(self, threadID, name, q, typeData, tmpDir):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        self.typeData = typeData
        self.tmpDir = tmpDir
    def run(self):
        print("Starting " + self.name)
        #print( "q: " + str(self.q) + " typeData: " + self.typeData + " donwloadDir: " + self.downloadDir)
        process_data(self.name, self.q, self.typeData, self.tmpDir)
        print("Exiting " + self.name)

def process_data(threadName, q, typeData, tmpDir):
    global queueLock
    global workQueue
    global exitFlag
    #print( "Enter in process_data. exitFlag: " + str(exitFlag))

    while not exitFlag:
        #print( "exitFlag is False")
        queueLock.acquire()
        #print( "\nqueueLock acquire")
        if not workQueue.empty():
            #print( "Queue is not empty")
            inputSubdir = q.get()
            outputSubdir = q.get()
            tiffile = q.get()
            queueLock.release()
            #print( "\nqueueLock release")
            print("\n%s processing %s" % (threadName, tiffile))
            # Generate file
            flagGenFile = False
            if(typeData == '123' or typeData == 'HRV'):
                flagGenFile = generateBands(tmpDir, inputSubdir, outputSubdir, tiffile, typeData)
            if(typeData == 'CLM'):        
                [flagGenFile, error] = generateCloudMask(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'NDVI'):
                flagGenFile = generateNDVI(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD08' or typeData == 'MYD08'):
                flagGenFile = generateMOD08(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD09' or typeData == 'MYD09'):
                flagGenFile = generateMOD09(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD13'):
                flagGenFile = generateMOD13(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD15'):
                flagGenFile = generateMOD15(tmpDir, inputSubdir, outputSubdir, tiffile)
            if(typeData == 'MOD17'):
                flagGenFile = generateMOD17(tmpDir, inputSubdir, outputSubdir, tiffile)

            if(flagGenFile == True):    
                print( '... File generate type: ' + outputSubdir + tiffile)
            else:
                print( '... Error to generate file: ' + tiffile)
                print("FileExtractError: %s error: %d \n" % (tiffile, flagGenFile))
                logFileError = os.getcwd() + '\\logExtractError.txt'
                logfile = open(logFileError, "a+")
                logfile.write("FileExtractError: %s error: %d \n" % (tiffile, flagGenFile))
                logfile.close()
                sys.exit()
                if(typeData == 'CLM'):
                    print('Error ' + str(error) + ' in ' + typeData + 'file:'  + tiffile)        
                # Here is necessary start another process to get this data
                #sys.exit()
            
        else:
            #print( "Queue is empty")
            try:
                queueLock.release()
                #print( "\nqueueLock release")
            except:
                print( 'Fail to release queueLock !')
        time.sleep(1)
    #print( "exitFlag is True")

def generateBands(tmpDir, inputDir, outputDir, file, typeData):

    hour = file[8:12]
    if(hour[2:] == '45'): hour2 = hour[:2] + '57'
    if(hour[2:] == '30'): hour2 = hour[:2] + '42'
    if(hour[2:] == '15'): hour2 = hour[:2] + '27'
    if(hour[2:] == '00'): hour2 = hour[:2] + '12'
    filename = file[:file.rfind(hour)] + hour2

    #print('file: ' + file)
    #print('filename: ' + filename)
    # Find files
    files = glob.glob(inputDir + '*' + filename + '*.*')
    if(len(files) == 0):
        print( 'Not have this file in subfolder: ' + inputDir)
        return False
    filename = ''
    for x in range(len(files)):
        filename = files[x]
        # extract tar file
        if(tarfile.is_tarfile(filename)):
            print('. Trying extract file: ' + filename)
            if(extractTarfile(filename, tmpDir) == False):
                print('... File extract Error: ' + filename)
                print('save log in: %s' %logTarError)
                logfile = open(logTarError, "a+")
                logfile.write("FileExtractError: %s\n" % filename)
                logfile.close()
                if(x == (len(files)-1)):
                    return False
            else:
                break
        else: 
            print('... Not is tar file: ' + filename)
            logfile = open(logTarError, "a+")
            logfile.write("NotIsTarFile: %s\n" % filename)
            logfile.close()
            if(x == (len(files)-1)):
                return False
    # Process file
    tmpfile1   = tmpDir + 'temp1' + file + '.tif'
    tmpfile2   = tmpDir + 'temp2' + file + '.tif'
    delFile(tmpfile1)
    delFile(tmpfile2)    
    outputFile = outputDir + file + '_' + typeData + '.tif'
    print('. GenerateFile outputFile: ' + outputFile)
    if(typeData == '123'):

        # Radiance without cut
        # Use option L to radiance and T to reflectance 
        print(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',(1,2,3,4,5,9),Y,T,1,1) ' + tmpfile1)
        try:
            subprocess.call(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',(1,2,3,4,5,9),Y,T,1,1) ' + tmpfile1)
        except:
            return False

        print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        try:
            subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        except:
            return False

        print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        except:
            return False
   

    else:
        # Radiance without cut
        print(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',12,Y,L,1,1) ' + tmpfile1)
        try:
            subprocess.call(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',12,Y,L,1,1) ' + tmpfile1)
        except:
            return False

        print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        try:
            subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile1 + ' ' + tmpfile2)
        except:
            return False

        print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        except:
            return False       
    delFile(tmpfile1)
    delFile(tmpfile2)
    delFile(tmpDir + '*' + file + '*')
    #sys.exit()
    return True

def generateCloudMask(tmpDir, inputDir, outputDir, file): # reviewed

    tmpfile1   = tmpDir + file + 'temp1_clm.tif'
    tmpfile2   = tmpDir + file + 'temp2_clm.tif'
    tmpfile3   = tmpDir + file + 'temp3_clm.tif'
    delFile(tmpfile1)
    delFile(tmpfile2)
    outputFile = outputDir + file + '_CLM.tif'
    files = glob.glob(inputDir + '*CLM*' + file + '*')
    print(inputDir + '*CLM*' + file + '*')
    if(len(files) > 1):
        print('Many files to time: ' + file)
    elif(len(files) == 0):
        print('Not have this file in subfolder: ' + inputDir)
        return False
    for x_files in range(len(files)):
        filename = files[0]      
        # Reproject image to WGS84
        print(exeDir + 'gdalwarp  -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" ' + filename + ' ' + tmpfile1)  
        try:
            subprocess.call(exeDir + 'gdalwarp -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" ' + filename + ' ' + tmpfile1)
        except:
            if(x_files < range(len(files))):
                continue
            else:
                return [False, -1]
        
        # Cut to dataRetrieve shape
        print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "bilinear" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
        try:
            subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r "bilinear" -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
        except:
            if(x_files < range(len(files))):
                delFile(filename)
                continue
            else:
                return [False, -2]
        # Translate data to rectangle Goias shape file                                      
        print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
        try:
            subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
        except:
            if(x_files < range(len(files))):
                continue
            else:
                return [False, -3]

        # Data cut using mask file
        print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outputFile)
        except:
            if(x_files < range(len(files))):
                continue
            else:
                return [False, -4]
    # Remove temp files
    delFile(tmpfile1)
    delFile(tmpfile2)
    delFile(tmpfile3)
    #sys.exit()
    return [True, 0]

def generateNDVI(tmpDir, inputDir, outputDir, file):
    files = glob.glob(inputDir + '*MSG3-SEVI-MSGNDVE*' + file + '*')
    if(len(files) > 1):
        print('Many files to time: ' + file)
    elif(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir + file)
        return False
    fileName = files[0]
    print('file: %s' % file)
    delFile(tmpDir + '*' + file + '*')
    print('fileName:' + fileName)
   
    print(gdalIlwisDir +  'gdal_translate.exe -of ILWIS HDF5:' + fileName + '://NDVImax '  +  tmpDir +  'tNDVImax'  +  file)
    try:
        subprocess.call(gdalIlwisDir +  'gdal_translate.exe -of ILWIS HDF5:' + fileName + '://NDVImax '  +  tmpDir +  'tNDVImax'  +  file)
    except:
        return False
    print(ilwDir + 'ilwis.exe -C ' + tmpDir + 'NDVImax_'   + file +'.mpr:=iff(' + tmpDir + 'tNDVImax'  + file +' le 100,' + tmpDir + 'tNDVImax'  + file + '/100,?);')
    try:
        subprocess.call(ilwDir + 'ilwis.exe -C ' + tmpDir + 'NDVImax_'   + file +'.mpr:=iff(' + tmpDir + 'tNDVImax'  + file +' le 100,' + tmpDir + 'tNDVImax'  + file + '/100,?);')
    except:
        return False
    # Convert study and define..
    print(ilwDir + 'ilwis.exe -C setgrf ' + tmpDir + 'NDVImax_'  + file + '.mpr ' + utilDir +'lritmsg;')
    try:
        subprocess.call(ilwDir + 'ilwis.exe -C setgrf ' + tmpDir + 'NDVImax_'  + file + '.mpr ' + utilDir +'lritmsg;')
    except:
        return False
    # Traduz imagem para coordenadas latitude e longitude
    print(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5568000 5568000 5568000 -5568000 ' + tmpDir + 'ndvimax_'  + file +'.mpr ' + tmpDir + '\\ndvimax_'+ file + '.tif')
    try:
        subprocess.call(exeDir + 'gdal_translate -a_srs  "+proj=geos +a=6378169 +b=6356583.8 +lon_0=0 +h=35785831" -a_ullr -5568000 5568000 5568000 -5568000 ' + tmpDir + 'ndvimax_'  + file +'.mpr ' + tmpDir + file + '_tNDVImax.tif')
    except:
        return False

    # Reprojeta imagem para WGS84 near/bilinear/cubic
    print(exeDir  + 'gdalwarp  -tr 0.0348420636891 0.0348426244628 -r near -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask + ' -crop_to_cutline -dstalpha ' + tmpDir + file +'_tNDVImax.tif ' + outputDir +  file + '_NDVI.tif')
    try:
        subprocess.call(exeDir + 'gdalwarp  -tr 0.0348420636891 0.0348426244628 -r near -s_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8" -t_srs "+proj=latlong +datum=WGS84" -cutline ' + fileMask + ' -crop_to_cutline -dstalpha ' + tmpDir + file +'_tNDVImax.tif ' + outputDir +  file + '_NDVI.tif')
    except:
        return False


    # Remove temporal files
    delFile(tmpDir + '*' + file + '*')
    return True

def generateMOD08(tmpDir, inputDir, outputDir, file):

    tmpfile1 = tmpDir + file[:file.rfind('.')] + '_tmp1.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    outputFile = outputDir + file

    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')
    fileMODIS = file[0:4] + str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    
    files = glob.glob(inputDir + 'M?D08_D3.A' + fileMODIS + '*.hdf')

    if(len(files) > 1):
        print('Many files to time: ' + file)
    elif(len(files) == 0):
        print('No have this file in subfolder: ' + inputDir)
        #print('findfile: '+ inputDir + '*' + file + '*')
        #print(files)
        #print(file)
        #sys.exit()
        return False
    filename = files[0]
    '''
    if('MOD08' in fileName):
        typeFile = '_MOD08.tif'
    else:
        typeFile = '_MYD08.tif'
    '''    
    #band1_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Surface Reflectance Band 1'
    #band2_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Surface Reflectance Band 2'
    #band3_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":MOD09CMG:Coarse Resolution Ozone'

    band1_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57
    band2_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37
    band3_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Total_Ozone_Mean' # Band 829
    band4_sub = 'HDF4_EOS:EOS_GRID:"'+ filename +'":mod08:Atmospheric_Water_Vapor_Mean' # Band 853

    #print('-------')
    #print(band1_sub)
    #print('-------')
    #print(band2_sub)
    #print(fileName)
    #sys.exit(1)
        
    Xorigin = -180 
    Yorigin = 90
    Xpixel = 1
    Ypixel = -1
    cols = 180
    rows = 360
    
    # Read data
    try:
        ds1 = gdal.Open(filename)
    except:
        print('Error to open file %s' %(filename))
        return False
    subdata = ds1.GetSubDatasets()

    band1 = gdal.Open(band1_sub)
    band2 = gdal.Open(band2_sub)
    band3 = gdal.Open(band3_sub)
    band4 = gdal.Open(band4_sub)

    BAND1 = band1.ReadAsArray()
    BAND2 = band2.ReadAsArray()
    BAND3 = band3.ReadAsArray()
    BAND4 = band4.ReadAsArray()
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(tmpfile1, rows, cols, 4, gdal.GDT_Int16)
    outRaster.SetGeoTransform((Xorigin, Xpixel, 0, Yorigin, 0, Ypixel))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outRaster.SetProjection(srs.ExportToWkt())
    bandOut1 = outRaster.GetRasterBand(1)
    bandOut2 = outRaster.GetRasterBand(2)
    bandOut3 = outRaster.GetRasterBand(3)
    bandOut4 = outRaster.GetRasterBand(4)
    
    bandOut1.WriteArray(BAND1)
    bandOut2.WriteArray(BAND2[1])
    bandOut3.WriteArray(BAND3)
    bandOut4.WriteArray(BAND4)
    
    bandOut1.FlushCache()
    bandOut2.FlushCache()
    bandOut3.FlushCache()
    bandOut4.FlushCache()
    outRaster = None
    # Cut to dataRetrieve shape
    print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r near -srcnodata [-9999] -dstnodata [-9999] -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
    try:
        subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r near -srcnodata [-9999] -dstnodata [-9999] -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + tmpfile1 + ' ' + tmpfile2) 
    except:                                  # Dimension in 123 file: 0.034840758326260 0.034840758326260
        return False
    # Translate data to rectangle Goias shape file                                      
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile2 + ' ' + tmpfile3)
    except:
        return False

    # Data cut using mask file
    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile2 + ' ' + outputFile)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile3 + ' ' + outputFile)
    except:
        return False
    # Remove temp files
    delFile(tmpfile1)
    delFile(tmpfile2)
    delFile(tmpfile3)
    #sys.exit()
    return True

def generateMOD09(tmpDir, inputDir, outputDir, file):

    tmpfile1a = tmpDir + file[:file.rfind('.')] + '_tmp1a.tif'
    tmpfile1b = tmpDir + file[:file.rfind('.')] + '_tmp1b.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    tmpfile4 = tmpDir + file[:file.rfind('.')] + '_tmp4.tif'
    tmpfile5 = tmpDir + file[:file.rfind('.')] + '_tmp5.tif'
    outputFile = outputDir + file
    if (os.path.isfile(outputFile)):
        return True
    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')
    doy = str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    fileMODIS = file[0:4] + doy
    find_files = inputDir + file[0:4] + '\\' + doy + '\\M?D09GA.A' + fileMODIS + '*.hdf'
    files = glob.glob(find_files)
    print('--------------- MOD09')
    print(find_files)
    print(files)

    if(len(files) == 2):
        if(not 'h13v10' in files[0] and not 'h12v10' in files[1]):
            file_tmp = files[0]
            files[0] = files[1]
            files[1] = file_tmp
    if(len(files) == 0):
        print('...Error not have files to process')
        return(0)

    band = ['":MODIS_Grid_500m_2D:sur_refl_b01_1', '":MODIS_Grid_500m_2D:sur_refl_b02_1', '":MODIS_Grid_1km_2D:state_1km_1', '":MODIS_Grid_500m_2D:QC_500m_1']
    data_type = ['Int16', 'Int16', 'Int32', 'Int32']
    logModisError = os.getcwd() + '\\logModError.txt'
    
    delFile(tmpfile3)
    delFile(tmpfile4)
    delFile(tmpfile5)  
    for i in range(len(band)):
        delFile(tmpfile1a)
        delFile(tmpfile1b)
        delFile(tmpfile2)
        print('--------------- warp')
        print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] + ' HDF4_EOS:EOS_GRID:"' + files[0] + band[i] + ' ' + tmpfile1a)
        try:
            subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] + ' HDF4_EOS:EOS_GRID:"' + files[0] + band[i] + ' ' + tmpfile1a)
        except:
            return (-1)
        if(len(files) == 2):
            print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] +  ' HDF4_EOS:EOS_GRID:"' + files[1] + band[i] + ' ' + tmpfile1b)
            try:
                subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type[i] +  ' HDF4_EOS:EOS_GRID:"' + files[1] + band[i] + ' ' + tmpfile1b)
            except:
                return (-1)    
        
        if(len(files) == 2): 
            print('--------------- merge')
            print('python ' + osgeoDir + 'gdal_merge.py  -n -28672  -ot ' + data_type[i] + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py  -n -28672 -ot ' + data_type[i] + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            except:
                return (-2)
        else:
            os.rename(tmpfile1a, tmpfile2)

        if('MODIS_Grid_1km_2D' in band[i]):
            print('------------- reproject')
            BAND_1km = readFileBand(tmpfile2, 1)
            ds = gdal.Open(tmpfile3)
            y_len = ds.RasterYSize
            x_len = ds.RasterXSize
            ds = None
            ds = gdal.Open(tmpfile2)
            y_len_1km = ds.RasterYSize
            x_len_1km = ds.RasterXSize
            ds = None
            print('tmptile2 [%d %d] tmpfile3 [%d %d]'%(y_len_1km, x_len_1km, y_len, x_len))
            BAND_1km_rep = numpy.zeros((y_len, x_len))
            for x in range(x_len):
                for y in range(y_len):
                    try:
                        BAND_1km_rep[y][x] = BAND_1km[int(y/2)][int(x/2)]
                    except:
                        return (-3)  
            delFile(tmpfile2)
            saveRASTERfile(tmpfile2, tmpfile3, BAND_1km_rep)
        if(not os.path.isfile(tmpfile3)):
            os.rename(tmpfile2, tmpfile3)
        else:
            print('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type[i] + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type[i] + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            except:
                return (-4)
            try:
                delFile(tmpfile3)
                os.rename(tmpfile4, tmpfile3)
            except:
                return (-5)
    # Translate data to rectangle Goias shape file                                      
    print('--------------- translate')
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856  -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    except:
        return (-6)
    print('--------------- warp')
    # Data cut using mask file
    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + tmpfile5)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + tmpfile5)
    except:
        return (-7)
    print('--------------- calc NDVI')
    RED = numpy.float32(readFileBand(tmpfile5, 1))*0.0001
    NIR = numpy.float32(readFileBand(tmpfile5, 2))*0.0001
    QUAL1 = numpy.int32(readFileBand(tmpfile5, 3))
    QUAL2 = numpy.int32(readFileBand(tmpfile5, 4))

    ds = gdal.Open(tmpfile5)
    y_len = ds.RasterYSize
    x_len = ds.RasterXSize
    ds = None
    NDVI = numpy.zeros((y_len, x_len))
    # Filter data, verify in Table 2: 500-meter QA Descriptions (32-bit)
    # https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod09ga
    for x in range(x_len):
        for y in range(y_len):
            if(RED[y][x] > 0 and NIR[y][x] > 0):
                NDVI[y][x] = (NIR[y][x] - RED[y][x])/(NIR[y][x] + RED[y][x])
                if(NDVI[y][x] > 1 or NDVI[y][x] < 0):
                    NDVI[y][x] = 0
            else:
                NDVI[y][x] = 0

    DATA  = [NDVI, RED, NIR, QUAL1, QUAL2] 
    saved = saveRASTERfile(outputFile, tmpfile5, DATA)
    if(not saved):
        return (-8)
    # Remove temp files
    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    delFile(tmpfile3)
    delFile(tmpfile4)
    delFile(tmpfile5)

    return True
def generateMOD17(tmpDir, inputDir, outputDir, file):

    tmpfile1a = tmpDir + file[:file.rfind('.')] + '_tmp1a.tif'
    tmpfile1b = tmpDir + file[:file.rfind('.')] + '_tmp1b.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    tmpfile4 = tmpDir + file[:file.rfind('.')] + '_tmp4.tif'


    doy = str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    fileMODIS = file[0:4] + doy
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_GPP.tif'

    if (os.path.isfile(outputFile)):
        return True
    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')

    find_files = inputDir + 'MOD17A2H.A' + fileMODIS + '*.hdf'

    files = glob.glob(find_files)
    print('--------------- MOD17AH2')
    print(find_files)
    print(files)

    if(len(files) == 2):
        if(not 'h13v10' in files[0] and not 'h12v10' in files[1]):
            file_tmp = files[0]
            files[0] = files[1]
            files[1] = file_tmp
    if(len(files) == 0):
        print('...Error not have files to process')
        return(0)


    # HDF4_EOS:EOS_GRID:"MOD17A2H.A2014001.h12v10.006.2015272092218":MOD_Grid_MOD17A2H:Gpp_500m
    # HDF4_EOS:EOS_GRID:"MOD17A2H.A2014001.h12v10.006.2015272092218":MOD_Grid_MOD17A2H:Psn_QC_500m
    bands = ['":MOD_Grid_MOD17A2H:Gpp_500m', '"::MOD_Grid_MOD17A2H:Psn_QC_500m']
    data_type = 'Int16'
    logModisError = os.getcwd() + '\\logModError.txt'
    
    delFile(tmpfile3)
    delFile(tmpfile4)
    for i in range(len(bands)):
        delFile(tmpfile1a)
        delFile(tmpfile1b)
        delFile(tmpfile2)
        print('--------------- warp')
        print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + bands[i] + ' ' + tmpfile1a)
        try:
            subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + bands[i] + ' ' + tmpfile1a)
        except:
            return (-1)
        if(len(files) == 2):
            print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + bands[i]  + ' ' + tmpfile1b)
            try:
                subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + bands[i]  + ' ' + tmpfile1b)
            except:
                return (-1)    
        
        if(len(files) == 2): 
            print('--------------- merge')
            print('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -ot ' + data_type + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -ot ' + data_type + ' -o  ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            except:
                return (-2)
        else:
            os.rename(tmpfile1a, tmpfile2)

        if(not os.path.isfile(tmpfile3)):
            os.rename(tmpfile2, tmpfile3)
        else:
            print('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type + ' -o ' +  tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -separate  -ot ' + data_type + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            except:
                return (-3)
            try:
                delFile(tmpfile3)
                os.rename(tmpfile4, tmpfile3)
            except:
                return (-4)

    # Data cut using mask file
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    except:
        return False
    # Cut using fileMask without change pixel position
    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + outputFile)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + outputFile)
    except:
        return False

    # Remove temp files
    #sys.exit()
    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    delFile(tmpfile3)
    delFile(tmpfile4)

    return True
def generateMOD15(tmpDir, inputDir, outputDir, file):

    tmpfile1a = tmpDir + file[:file.rfind('.')] + '_tmp1a.tif'
    tmpfile1b = tmpDir + file[:file.rfind('.')] + '_tmp1b.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    tmpfile4 = tmpDir + file[:file.rfind('.')] + '_tmp4.tif'


    doy = str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    fileMODIS = file[0:4] + doy
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_LAI.tif'

    if (os.path.isfile(outputFile)):
        return True
    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')

    find_files = inputDir + 'MOD15A2H.A' + fileMODIS + '*.hdf'

    files = glob.glob(find_files)
    print('--------------- MOD15AH2')
    print(find_files)
    print(files)

    if(len(files) == 2):
        if(not 'h13v10' in files[0] and not 'h12v10' in files[1]):
            file_tmp = files[0]
            files[0] = files[1]
            files[1] = file_tmp
    if(len(files) == 0):
        print('...Error not have files to process')
        return(0)

    # Data name:HDF4_EOS:EOS_GRID:"MCD15A2H.A2014001.h12v10.006.2015273193555"

    bands = ['":MOD_Grid_MOD15A2H:Lai_500m', '":MOD_Grid_MOD15A2H:FparLai_QC', '":MOD_Grid_MOD15A2H:LaiStdDev_500m']
    data_type = 'Int16'
    logModisError = os.getcwd() + '\\logModError.txt'
    
    delFile(tmpfile3)
    delFile(tmpfile4)
    for i in range(len(bands)):
        delFile(tmpfile1a)
        delFile(tmpfile1b)
        delFile(tmpfile2)
        print('--------------- warp')
        print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + bands[i] + ' ' + tmpfile1a)
        try:
            subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + bands[i] + ' ' + tmpfile1a)
        except:
            return (-1)
        if(len(files) == 2):
            print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + bands[i]  + ' ' + tmpfile1b)
            try:
                subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + bands[i]  + ' ' + tmpfile1b)
            except:
                return (-1)    
        
        if(len(files) == 2): 
            print('--------------- merge')
            print('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -ot ' + data_type + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -ot ' + data_type + ' -o  ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            except:
                return (-2)
        else:
            os.rename(tmpfile1a, tmpfile2)

        if(not os.path.isfile(tmpfile3)):
            os.rename(tmpfile2, tmpfile3)
        else:
            print('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type + ' -o ' +  tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -separate  -ot ' + data_type + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            except:
                return (-3)
            try:
                delFile(tmpfile3)
                os.rename(tmpfile4, tmpfile3)
            except:
                return (-4)

    # Data cut using mask file
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    except:
        return False
    # Cut using fileMask without change pixel position
    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + outputFile)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + outputFile)
    except:
        return False

    # Remove temp files
    #sys.exit()
    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    delFile(tmpfile3)
    delFile(tmpfile4)

    return True

def generateMOD13(tmpDir, inputDir, outputDir, file):

    tmpfile1a = tmpDir + file[:file.rfind('.')] + '_tmp1a.tif'
    tmpfile1b = tmpDir + file[:file.rfind('.')] + '_tmp1b.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    tmpfile4 = tmpDir + file[:file.rfind('.')] + '_tmp4.tif'


    doy = str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    fileMODIS = file[0:4] + doy
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_MOD13.tif'

    if (os.path.isfile(outputFile)):
        return True
    if not gdal.GetDriverByName('HDF5'):
        raise Exception('HDF5 driver is not available')

    find_files = inputDir + 'MOD13Q1.A' + fileMODIS + '*.hdf'

    files = glob.glob(find_files)
    print('--------------- MOD13AH2')
    print(find_files)
    print(files)

    if(len(files) == 2):
        if(not 'h13v10' in files[0] and not 'h12v10' in files[1]):
            file_tmp = files[0]
            files[0] = files[1]
            files[1] = file_tmp
    if(len(files) == 0):
        print('...Error not have files to process')
        return(0)
    '":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days VI Quality'
    bands = ['":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days NDVI"', '":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days composite day of the year"','":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days pixel reliability"',  '":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days VI Quality"']
    # Data name:HDF4_EOS:EOS_GRID:"MCD15A2H.A2014001.h12v10.006.2015273193555"
    #bands = ['":MOD_Grid_MOD15A2H:Fpar_500m', '":MOD_Grid_MOD15A2H:FparExtra_QC', '":MOD_Grid_MOD15A2H:FparLai_QC', '":MOD_Grid_MOD15A2H:FparStdDev_500m']
    #bands = ['":MODIS_Grid_16DAY_1km_VI:"1 km 16 days NDVI"', '":MODIS_Grid_16DAY_1km_VI:"1 km 16 days composite day of the year"', '":MODIS_Grid_16DAY_1km_VI:"1 km 16 days pixel reliability"', '":MODIS_Grid_16DAY_1km_VI:"1 km 16 days VI Quality"']
    data_type = 'Int16'
    logModisError = os.getcwd() + '\\logModError.txt'
    
    delFile(tmpfile3)
    delFile(tmpfile4)
    for i in range(len(bands)):
        delFile(tmpfile1a)
        delFile(tmpfile1b)
        delFile(tmpfile2)
        print('--------------- warp')
        print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + bands[i] + ' ' + tmpfile1a)
        try:
            subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type + ' HDF4_EOS:EOS_GRID:"' + files[0] + bands[i] + ' ' + tmpfile1a)
        except:
            return (-1)
        if(len(files) == 2):
            print(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + bands[i]  + ' ' + tmpfile1b)
            try:
                subprocess.call(osgeoDir + 'gdalwarp.exe -srcnodata -28672 -t_srs "+proj=latlong +datum=WGS84" -ot ' + data_type +  ' HDF4_EOS:EOS_GRID:"' + files[1] + bands[i]  + ' ' + tmpfile1b)
            except:
                return (-1)    
        
        if(len(files) == 2): 
            print('--------------- merge')
            print('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -ot ' + data_type + ' -o ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -n -28672 ' + ' -ot ' + data_type + ' -o  ' + tmpfile2 + ' ' + tmpfile1a + ' ' + tmpfile1b)
            except:
                return (-2)
        else:
            os.rename(tmpfile1a, tmpfile2)

        if(not os.path.isfile(tmpfile3)):
            os.rename(tmpfile2, tmpfile3)
        else:
            print('python ' + osgeoDir + 'gdal_merge.py -separate -ot ' + data_type + ' -o ' +  tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            try:
                subprocess.call('python ' + osgeoDir + 'gdal_merge.py -separate  -ot ' + data_type + ' -o ' + tmpfile4 + ' ' + tmpfile3 + ' ' + tmpfile2)
            except:
                return (-3)
            try:
                delFile(tmpfile3)
                os.rename(tmpfile4, tmpfile3)
            except:
                return (-4)

    # Data cut using mask file
    print(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    try:
        subprocess.call(exeDir + 'gdal_translate.exe -projwin -53.2507227 -12.3947856 -45.8993227 -19.5023003 -of GTiff ' + tmpfile3 + ' ' + tmpfile4)
    except:
        return False

    print(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + outputFile)
    try:
        subprocess.call(exeDir + 'gdalwarp.exe -cutline ' + fileMask +  ' ' + tmpfile4 + ' ' + outputFile)
    except:
        return False

    # Remove temp files
    #print('Stop process...')
    #sys.exit()
    delFile(tmpfile1a)
    delFile(tmpfile1b)
    delFile(tmpfile2)
    delFile(tmpfile3)
    delFile(tmpfile4)

    return True

def generateMOD15_lapig(tmpDir, inputDir, outputDir, file):

    tmpfile1 = tmpDir + file[:file.rfind('.')] + '_tmp1.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    band_names = ['fpar']
    fileMODIS = file[0:4] + str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    #outputFile = outputDir + file
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_MOD15.tif'


    num_band  = 1
    for pos_band in range(len(band_names)):
        subname    = inputDir + 'pa_br_'+ band_names[pos_band] +'_1000_lapig'
        filesearch = inputDir + 'pa_br_'+ band_names[pos_band] +'_1000_lapig' + '\\' + 'pa_br_'+ band_names[pos_band] +'_1000_'+ fileMODIS + '_lapig.tif'
        files = glob.glob(filesearch)
        print(filesearch)
        print(files)
        if(len(files) > 1):
            print('Many files to time: ' + file)
        elif(len(files) == 0):
            print('No have this file in subfolder: ' + filesearch)
            if(pos_band == 0):
                return False

        else:
            filename = files[0]


            print(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
            try:
                subprocess.call(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
            except:
                return False

            # Data cut using mask file
            print(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
            try:
                subprocess.call(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
            except:
                return False


            VET_DATA = readFileBand(tmpfile3, 1)

            if(pos_band == 0):
               # Save georeference file
                #print('------------------')
                #print(num_band)
                #print(outputDir + file)
                ds = gdal.Open(tmpfile3)
                band = ds.GetRasterBand(1)
                driver = gdal.GetDriverByName("GTiff")
                dsOut  = driver.Create(outputFile, ds.RasterXSize, ds.RasterYSize, len(band_names), band.DataType)
                CopyDatasetInfo(ds,dsOut)
                bandOut = dsOut.GetRasterBand(num_band)
                BandWriteArray(bandOut, VET_DATA)
            else:
                num_band = num_band + 1
                #print('------------------')
                #print(num_band)
                bandOut = dsOut.GetRasterBand(num_band)
                BandWriteArray(bandOut, VET_DATA) 
            #print('---- go out to test...')
            delFile(tmpDir + file[:file.rfind('.')] + '_tmp?.tif')
    return True

def generateMOD13_lapig(tmpDir, inputDir, outputDir, file):

    tmpfile1 = tmpDir + file[:file.rfind('.')] + '_tmp1.tif'
    tmpfile2 = tmpDir + file[:file.rfind('.')] + '_tmp2.tif'
    tmpfile3 = tmpDir + file[:file.rfind('.')] + '_tmp3.tif'
    #band_names = ['NDVI', 'RED', 'NIR', 'MIR', 'composite_day', 'pixel_reliability']
    band_names = ['NDVI']
    fileMODIS = file[0:4] + str('%03d'%calcDoy (int(file[0:4]), int(file[4:6]), int(file[6:8])))
    #outputFile = outputDir + file
    outputFile = outputDir.replace(file[0:4]+'\\'+file[4:6], '') + fileMODIS + '_MOD13.tif'


    num_band  = 1
    for pos_band in range(len(band_names)):
        subname    = inputDir + 'pa_br_'+ band_names[pos_band] +'_250_lapig'
        filesearch = inputDir + 'pa_br_'+ band_names[pos_band] +'_250_lapig' + '\\' + 'pa_br_'+ band_names[pos_band] +'_250_'+ fileMODIS + '_lapig.tif'
        files = glob.glob(filesearch)
        print(filesearch)
        print(files)
        if(len(files) > 1):
            print('Many files to time: ' + file)
        elif(len(files) == 0):
            print('No have this file in subfolder: ' + inputDir + subname)
            #print('findfile: '+ inputDir + '*' + file + '*')
            #print(files)
            #print(file)
            #sys.exit()
            return False
        filename = files[0]

        # don't transform spatial resolution 
        #tmpfile1 = filename

        '''
        # Cut to dataRetrieve shape and transform 3 x 3 km
        print(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r cubic -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + filename + ' ' + tmpfile1) 
        try:
            subprocess.call(exeDir + 'gdalwarp -tr 0.034840758326260 0.034840758326260 -r cubic -cutline ' + fileDataRetriever + ' -crop_to_cutline ' + filename + ' ' + tmpfile1) 
        except:                                  # Dimension in 123 file: 0.034840758326260 0.034840758326260
            return False
        '''
        print(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -tr 0.0180110966371816 0.0180110966371816 -r cubic -cutline ' + fileMaskRect + ' -crop_to_cutline -of GTiff ' + filename + ' ' + tmpfile1)
        except:
            return False

        # Data cut using mask file
        print(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
        try:
            subprocess.call(exeDir + 'gdalwarp.exe -crop_to_cutline -cutline ' + fileMask +  ' ' + tmpfile1 + ' ' + tmpfile3)
        except:
            return False

        VET_DATA = readFileBand(tmpfile3, 1)

        if(pos_band == 0):
           # Save georeference file
            #print('------------------')
            #print(num_band)
            #print(outputDir + file)
            ds = gdal.Open(tmpfile3)
            band = ds.GetRasterBand(1)
            driver = gdal.GetDriverByName("GTiff")
            dsOut  = driver.Create(outputFile, ds.RasterXSize, ds.RasterYSize, len(band_names), band.DataType)
            CopyDatasetInfo(ds,dsOut)
            bandOut = dsOut.GetRasterBand(num_band)
            BandWriteArray(bandOut, VET_DATA)
        else:
            num_band = num_band + 1
            #print('------------------')
            #print(num_band)
            bandOut = dsOut.GetRasterBand(num_band)
            BandWriteArray(bandOut, VET_DATA) 
        #print('---- go out to test...')
        delFile(tmpDir + file[:file.rfind('.')] + '_tmp?.tif')
    return True

class main():
    global exitFlag
 
    argDict = mapDict(argv, usage)

    if "-ty" in argDict and "-ds" in argDict and "-de" in argDict:
        typeData = str(argDict["-ty"])
        yearStart = int(argDict["-ds"][0:4])
        yearFinish = int(argDict["-de"][0:4])
        if('MOD' in typeData or 'MYD' in typeData):
            monthStart = int(argDict["-ds"][4:7])
            monthFinish = int(argDict["-de"][4:7])
        else:    
            monthStart = int(argDict["-ds"][4:6])
            monthFinish = int(argDict["-de"][4:6])
    else:
        exit(usage)
    if "-fo" in argDict and "-fi" in argDict:
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
    else:
        inputDir  = 'f:\\DADOS\\'
        outputDir = inputDir

    if (typeData not in ['CLM', '123', 'NDVI', 'MOD08', 'MYD08', 'MOD09', 'MYD09', 'HRV', 'MOD13', 'MOD15', 'MOD17']):
        print('Error in typeData: %s'%(typeData))
        exit(usage)

    if(typeData == '123'):
        inputDir = inputDir + 'HRIT\\'
        outputDir = outputDir + '123\\'
    elif(typeData == 'HRV'):
        inputDir = inputDir + 'HRIT\\'
        outputDir = outputDir + 'HRV\\'        
    elif(typeData == 'CLM'):
        inputDir = inputDir + 'GRB\\'
        outputDir = outputDir + 'CLM\\'
    elif(typeData == 'NDVI'):
        inputDir = inputDir + 'HDF\\'
        outputDir = outputDir + 'NDVI\\'
    elif(typeData == 'MOD08'):
        inputDir = inputDir + 'MOD08\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MYD08'):
        inputDir = inputDir + 'MYD08\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MOD09'):
        inputDir = inputDir + 'MOD09\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MYD09'):
        inputDir = inputDir + 'MYD09\\'
        outputDir = outputDir + 'MODIS\\' 
    elif(typeData == 'MOD13'):
        inputDir = inputDir + 'MOD13Q1\\'      
        outputDir = outputDir + 'MOD13_250m\\' 
    elif(typeData == 'MOD15'):
        inputDir = inputDir + 'MOD15A2H\\'      
        outputDir = outputDir + 'MOD15\\' 
    elif(typeData == 'MOD15'):
        inputDir = inputDir + 'MOD17A2H\\'      
        outputDir = outputDir + 'MOD17\\' 
    gc.collect()
    gdal.UseExceptions()
    lst = []
    print('....................... Generate data in folder: ' + inputDir)
    for year in range(yearStart, yearFinish+1):
        print('...........Generate data in year: ' + str(year))
        folder = outputDir + '\\' + str(year)
        if not(os.path.isdir(folder)):
            os.makedirs(folder)
        if('MOD' in typeData or 'MYD' in typeData):
            if(yearStart == yearFinish):
                monthInic = monthStart
                monthEnd  = monthFinish                    
            elif(year == yearStart):
                monthInic = monthStart
                monthEnd  = 365
            elif(year == yearFinish):
                monthInic = 1
                monthEnd  = monthFinish
            else:
                monthInic = 1
                monthEnd  = 365

        else:
            if(yearStart == yearFinish):
                monthInic = monthStart
                monthEnd  = monthFinish  
            elif(year == yearStart):
                monthInic = monthStart
                monthEnd  = 12
            elif(year == yearFinish):
                monthInic = 1
                monthEnd  = monthFinish
            else:
                monthInic = 1
                monthEnd  = 12                   

        for month in range(monthInic, monthEnd+1):
            str_month = ''
            if('MOD' in typeData or 'MYD' in typeData):
                str_month = str('%03d'%month) 
                outputSubdir = outputDir + str(year) + '\\' + str('%02d'%calcMonth(year, month)) + '\\'
                if(typeData == 'MOD08' or typeData == 'MYD08'):
                    inputSubdir  = inputDir + str(year) + '\\' + str_month + '\\'
                else:
                    inputSubdir  = inputDir
            else:
                str_month = str('%02d'%month)
                outputSubdir = outputDir + str(year) + '\\'  + str_month + '\\'
                inputSubdir  = inputDir + str(year) + '\\'  + str_month + '\\'
            #print('str_month: ' + str_month)
            #print('typeData: ' + typeData)
            #print('month: ' + str(month))
            #sys.exit()
            print('....... Generate data in month/doy: ' + str_month)

            print('... Generate files in folder: ' + outputSubdir)
            findfile = inputDir+ str(year) + '\\' + str_month + '\\' + '*-' + str(year) + str_month + '*.*'
            if(typeData == 'HRV'):
                findfile = findfile.replace('*.*', '??19*.*')
            if(typeData == 'CLM'):
                findfile = findfile.replace('*.*', '*.grb')
            if(typeData == 'NDVI'):
                findfile = findfile.replace('*.*', '*.h5')
            if(typeData == 'MOD08' or typeData == 'MYD08'):
                findfile = inputDir + str(year) + '\\' + str_month + '\\' + '*A' + str(year) + str_month + '*.hdf'
            if(typeData == 'MOD09' or typeData == 'MYD09'):
                findfile = inputDir + str(year) + '\\' + str_month + '\\' + '*A' + str(year) + str_month + '*h13v10*.hdf'
            if(typeData == 'MOD13'):
                #findfile = inputDir +'pa_br_ndvi_250_lapig\\pa_br_ndvi_250_' + str(year) + str_month + '_lapig.tif'
                findfile = inputDir + 'MOD13Q1.A' + str(year) + str_month + '*h13v10*.hdf'
            if(typeData == 'MOD15'):
                #findfile = inputDir +'pa_br_fpar_1000_lapig\\pa_br_fpar_1000_' + str(year) + str_month + '_lapig.tif'
                findfile = inputDir + 'MOD15A2H.A' + str(year) + str_month + '*h13v10*.hdf'
            if(typeData == 'MOD17'):
                #findfile = inputDir +'pa_br_fpar_1000_lapig\\pa_br_fpar_1000_' + str(year) + str_month + '_lapig.tif'
                findfile = inputDir + 'MOD17A2H.A' + str(year) + str_month + '*h13v10*.hdf'              
            print('... Finding file in folder: ' + findfile)
            files = glob.glob(findfile)
            print(findfile)
            print(files)
            #print('findfile: ' + findfile)
            rows = len(files)
            if(rows == 0):
                print('... No files to folder: ' + outputSubdir)
                continue
            if not(os.path.isdir(outputSubdir)):
                os.makedirs(outputSubdir)
            id = len(inputDir)
            for i in range (0, rows, 1):
                print('... file files['+str(i)+']: '+ files[i])
                # Verify if tar file has files to this year and month
                if not ('MOD' in typeData or 'MYD' in typeData):
                    id = files[i].index('-201')
                    day  = files[i][id+7:id+9]
                    hour = files[i][id+9:id+13]
                    if(typeData == '123' or typeData == 'HRV'):
                        if(hour[2:] == '57'): hour = hour[:2] + '45'
                        if(hour[2:] == '42'): hour = hour[:2] + '30'
                        if(hour[2:] == '27'): hour = hour[:2] + '15'
                        if(hour[2:] == '12'): hour = hour[:2] + '00'

                print('. Verify file: ' + files[i])
                
                if('MOD' in typeData or 'MYD' in typeData):
                    tiffile = str(year) + str('%02d'%calcMonth(year, month)) + str('%02d'%calcDay(year, month))
                    extension = '_' + typeData + '.tif'
                    tiffile = tiffile + extension
                    testfile = outputSubdir + tiffile 
                else:
                    tiffile = str(year) + str_month + str(day) + str(hour)
                    extension = '_' + typeData + '.tif'
                    testfile = outputSubdir + tiffile + extension
                if(os.path.isfile(testfile)):
                    print('... File exist in: ' + testfile)
                    if(typeData == 'MOD09'):
                        if(numpy.sum(readFileBand(testfile, 5)) == 0 or numpy.sum(readFileBand(testfile, 4)) == 0):
                            print('... Error in file band, I''m sorry, I delete file!')
                            delFile(testfile)
                            lst.append(inputSubdir)
                            lst.append(outputSubdir)
                            lst.append(tiffile)
                    continue
                else:
                    print('... File not exist in: ' + testfile)
                    lst.append(inputSubdir)
                    lst.append(outputSubdir)
                    lst.append(tiffile)

    if (lst == []):
        print('Waiting new files')
    else:
        print( "----- Create new threads")
        threadID = 0
        for tName in threadList:
            thread = myThread(threadID, tName, workQueue, typeData, tmpDir)
            thread.start()
            threads.append(thread)
            threadID += 1

        print( "----- Fill the queue and wait case many files")
        
        x = 0
        print('----- lst')
        for i in range(len(lst)):
            while(x < len(lst)):
                queueLock.acquire()
                if not workQueue.full():
                    word =  lst[x]
                    #print( 'x: ' + str(x) + '\tfilename: ' + word)
                    workQueue.put(word)
                    x = x + 1
                    queueLock.release()
                    
                else:
                    queueLock.release()
                    time.sleep(5)

        print( "----- Wait for queue to empty")
        while not workQueue.empty():
            #print('Files in Queue: ' + str(workQueue.qsize()))
            #time.sleep(5)
            pass
        print( '----- Queue is empty')

        # Notify threads it's time to exit
        exitFlag = 1

        # Wait for all threads to complete
        for t in threads:
            t.join()
        print( "Exiting Main Thread")

   
