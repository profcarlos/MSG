import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import subprocess
import os
import glob
import subprocess
import numpy
import logging

# Constantes
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'

fileMask = fileBrasil

logfile = os.getcwd() + '\\logs\\log_gen_ndvi_msg.log'
dataRetrieve = os.getcwd() + '\\MSGDataRetriever\\'
exeFolder  = 'c:\\Program Files\\GDAL\\'
tmpFolder  = 'c:\\tmp\\'
path_rec = 'z:\\receive\\'
path_MSG = '\\MSG-0degree\\IMG-3h\\'

inputFolder = path_rec + path_MSG
outputFolder = 'z:\\produtc\\NDVI\\' 

def generateBands(tmpDir, inputDir, outputDir, file, typeData):

    filename = ''
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
        print(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',(1,2,3),Y,L,1,1) ' + tmpfile1)
        try:
            subprocess.call(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',(1,2,3),Y,L,1,1) ' + tmpfile1)
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
        print(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',12,Y,T,1,1) ' + tmpfile1)
        try:
            subprocess.call(dataRetrieve + 'gdalwarp.exe --config GDAL_CACHEMAX 1000 -t_srs "+proj=latlong +datum=WGS84"  -of GTiff MSG(' + tmpDir + ',' + file + ',12,Y,T,1,1) ' + tmpfile1)
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

class main:

    logging.basicConfig(filename=logfile,level=logging.DEBUG)
    logging.debug('Start process to generate NDVI MSG')

    filename = 'H-000-MSG3__-MSG3________-_________-PRO______-201*1200-__'
    inputDir = 'z:\\receive\\' + path_MSG
    files = glob.glob(inputDir + '*' + filename + '*.*')
    if(len(files) == 0):
        print('Error. Not have this file in subfolder: ' + inputDir)
        logging.debug ('Error. Not have this file in subfolder: ' + inputDir)
        return False
        sys.exit()
    else:
        for file in files:
            id_date = file.index('-201') - 3
            date_pro_file = file[id_date:id_date+12]
            print(id_date)
            print(date_pro_file)
            filename = date_pro_file + '_123.tif'
        if(not os.path.isfile(outputFolder + filename)):
            if(not generateBands(tmpFolder, inputFolder, outputFolder, date_pro_file, typeData='123'))
                print('Error to generate file: ' + filename)
                logging.debug('Error to generate file: ' + filename)
                sys.exit()
            else:
                print('Generate file: ' + filename)
                logging.debug('Generate file: ' + filename)