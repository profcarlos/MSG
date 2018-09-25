import sys
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\SMAC')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\SMAC')

import numpy
import matplotlib.pyplot as plt
import itertools
import kernels
import os
import gc
import glob
import sys
from utils import * 
from sys import argv
from kernels import Kernels # Import ther kernels
from Py6S import *
from SMAC_CESBIO import *
import time

otbDir = 'E:\\Install\\OTB-5.4.0-win64\\OTB-5.4.0-win64\\bin\\'
exeDir = 'C:\\OSGeo4W\\bin\\'
fileBrasil = os.getcwd() + '\\shapes\\pa_br_bioma_5000_2004_IBGE.shp'
fileGoias  = os.getcwd() + '\\shapes\\limite_go_df_WGS84.shp'
fileGoias123 = os.getcwd() + '\\shapes\\201401011300_123.tif'
fileGoiasHRV = os.getcwd() + '\\shapes\\201401011700_HRV??.tif' # Colocar arquivo atualizado
fileMask   = os.getcwd() + '\\shapes\\maskVectorGoias.shp'
fileDataRetriever = os.getcwd() + '\\shapes\\dataRetrieverShape.shp'

fileMaskSHP = fileGoias
fileMaskHRV = fileGoiasHRV
fileMask123 = fileGoias123

usage = """\
Usage: %s [OPTIONS]
        -ty     type data (1day, 2days)
        -ps     type input data (123, 6S)
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)

""" % argv[0]

def importATMdata(inputDir, date):
    #print('ATMclass read')
    # Read ATM file
    # band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
    # band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
    # band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
    # band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)
    # ATM to 6S: Vector with atmospheric data [WV O AOT]

    ds = gdal.Open(fileMask123)
    y_123 = ds.RasterYSize
    x_123 = ds.RasterXSize
    # DATA[3][y_123][x_123] = [WV, OZONE, AOT550]
    DATA = numpy.zeros((3, y_123, x_123))

    inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
    filename = inputSubDir.replace('___', 'ATM') + str(date) + '_ATM.tif'
    if os.path.isfile(filename) == 0:
        #print('Error to read ATM file: %s' % (filename))
        return []
    else:
        #print('Read ATM file: %s' % (filename))
        '''
        # get Aerosol Data Band
        vet = readFileBand(filename, 4)
        if(not len(vet)):
            return [0, 0, 0]
        else:
            DATA[2] = vet[i][j]
        '''
        # get Water_Vapour  
        vet = readFileBand(filename, 10)
        if(not len(vet)):
            return []
        else:
            DATA[0] = vet

        # get Total Ozone Data Band
        vet = readFileBand(filename, 7)
        if(not len(vet)):
            return []
        else:
            DATA[1] = vet/1000 # Convert Dobson Unit to cm-atm
        # get Deep_Blue Data Band
        vet = readFileBand(filename, 1)
        if(not len(vet)):
            return []
        else:
            DATA[2] = vet
    return DATA


def importBRDFdata(inputDir, date, type):
    startHour = 1300
    endHour   = 1500    

    ds = gdal.Open(fileMask123)
    y_123 = ds.RasterYSize
    x_123 = ds.RasterXSize
    # DATA[17][8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
    DATA = numpy.zeros((17, 8, y_123, x_123))
    inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
    hour = startHour
    h = 0

    while(hour <= endHour):
        fileError = False
        # --------- Process ANGLE file
        filename = inputSubDir.replace('___', 'ANGLE') + str(date) + str(hour) + '_ANGLES.tif'
        if os.path.isfile(filename) == 0:
            print('..Error! Invalid data file: %s' %filename)
            fileError = True
        else:
            #print('.Read ANGLE file:%s'%(filename))
            # get VZA View Zenith Angle (msg_zenres.tif)
            vet = readFileBand(filename, 3)
            if(not len(vet)):
                fileError = True
            else:
                DATA[h][0] = vet
            # get VAZ View Azimuth Angle (msg_azres.tif)
            vet = readFileBand(filename, 1)
            if(not len(vet)):
                fileError = True
            else:
                DATA[h][1] = vet
            # get SZA Sun Zenith Angle (sun_zenres.tif)
            vet = readFileBand(filename, 4)
            if(not len(vet)):
                fileError = True
            else:
                DATA[h][2] = vet
            # get SAZ Sun Azimuth Angle (sol_azres.tif)
            vet = readFileBand(filename, 2)
            if(not len(vet)):
                fileError = True
            else:
                DATA[h][3] = vet
        # --------- Process 123 file
        if(fileError == False):
            #print('Process 123 file')
            filename = inputSubDir.replace('___', type) + str(date) + str(hour) + '_' + type + '.tif'
            if os.path.isfile(filename) == 0:
                print('..Error! Invalid data file: %s' %filename)
                fileError = True
            else:
                #print('.Read 123 file:%s'%(filename))
                # get RED Data Band
                vet = readFileBand(filename, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    DATA[h][4] = vet
                # get NIR Data Band
                vet = readFileBand(filename, 2)
                if(not len(vet)):
                    fileError = True
                else:
                    DATA[h][5] = vet
                # get WIR Data Band
                vet = readFileBand(filename, 3)
                if(not len(vet)):
                    fileError = True
                else:
                    DATA[h][6] = vet

        # --------- Process CLM file
        if(fileError == False):
            #print('Process CLM file')
            filename = inputSubDir.replace('___', 'CLM') + str(date) + str(hour) + '_CLM.tif'
            if os.path.isfile(filename) == 0:
                print('..Error! Invalid data file: %s' %filename)
                fileError = True
            else:
                #print('.Read CLM file:%s'%(filename))
                vet = readFileBand(filename, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    DATA[h][7] = vet
        # Verify if error clear data, else inc position
        if(fileError == False):
            h = h + 1
        hour = incHour(hour)
    return DATA[:h]

def calc_lstsq(BAND, kern):
    # Calc to RED and NIR
    samples = len(kern)
    len_data= len(BAND)
    NDVI = numpy.zeros((3, samples))
    COV  = [0,0]
    # Receive n_sample to process NDVI
    s_sample = 99
    for band in range(2):
        #print('passer: %d  kern: %d'%(len(passer), len(kern)))
        tmp_obs = BAND
        obs_col = numpy.zeros((len(tmp_obs),1))
        obs_lin = numpy.zeros((len(tmp_obs)))
        for x_obs in range(len(tmp_obs)):
            obs_col[x_obs] = tmp_obs[x_obs][band]
            obs_lin[x_obs] = tmp_obs[x_obs][band]

        if(samples > len_data):
            K = kern[len_data, :]
        else:
            K = kern
        r_sample = [samples]
        for x in range(samples-2, 4, -2):
            if(x >= 5):
                r_sample.append(x)
            else:
                break     
        
        #print('r_sample:')
        #print(r_sample)
        for n_sample in r_sample:
            (f, rmse, rank, svals ) = numpy.linalg.lstsq(K[:n_sample], obs_col[:n_sample])
            '''
            print('n_samples: %d'%(n_sample))
            print('K:')
            print(K[:n_sample])
            print('obs_col:')
            print(obs_col[:n_sample])

            print('----- Results Band[%d]' %(band))
            print ("%-20s %20s" % ( "Kernel", "Value"))
            for i, k in enumerate( ["Isotropic", "Ross-Thick", "Li-Sparce"] ):
                print ("%-20s %20f" % ( k, f[i] ))
            '''
            if(f[0] > 0):
                break

        # Error if k_iso <= 0
        if(f[0] <= 0):
            #print('...return [0] in else')       
            return [0]
        else:
            if(s_sample > n_sample):
                s_sample = n_sample

        fwd = K.dot(f)
        fwd_lin = numpy.zeros(len(fwd))
        for x_fwd in range(len(fwd)):
            fwd_lin[x_fwd] = fwd[x_fwd]
        COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
        NDVI[band] = fwd_lin
    # Calc NDVI = (NIR - RED)/(NIR + RED)
    NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
    #print('...return NDVI (s_sample: %d of %d):'%(s_sample, len(NDVI[0])))
    NDVI_filt = numpy.zeros((3,s_sample))
    for x in range(3):
        NDVI_filt[x] = NDVI[x][:s_sample]
    #print(NDVI_filt)
    return NDVI_filt

class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    
    if ("-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict and "-ps" in argDict):
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData = argDict["-ty"]
        processData = argDict["-ps"]
    else:
        exit(usage)
    if(typeData != '1day' and typeData != '2days'):
        print('Error in typeData, verify!')
        exit(usage)

    if(processData != '123' and processData != '6S'):
        print('Error in processData, verify!')
        exit(usage)


    # Go out if folder not exist
    if(os.path.isdir(inputDir) == False or os.path.isdir(outputDir) == False):
        print('... Error in folder adress!')
        sys.exit()
    # Create output data matrix
    ds = gdal.Open(fileMask123)
    y_123 = ds.RasterYSize
    x_123 = ds.RasterXSize
    KERNEL = numpy.zeros((y_123, x_123, 5, 5))
    ds = None   
    dateSaveKernel = dateStart
    # Process in range of data
    date = dateStart
    while(date < dateEnd):
        # Output file name   
        outputFile = outputDir + 'BRDF_' + processData + '\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        # Verify and create folders
        if(not os.path.isdir(outputFile)):
            os.makedirs(outputFile)
        kernelFile = outputFile + str(date) + '_KERN.npy'
        if(date == dateStart):
            if(os.path.isfile(kernelFile)):
                KERNEL = numpy.load(kernelFile)
                print('Reading Kernel file data: %s'%(kernelFile))
        outputFileBRDF = outputFile + str(date) + '_BRDF.tif'    

        # Verify if files exist
        if(os.path.isfile(outputFileBRDF)):
            print('...File exist in: %s'%(outputFileBRDF))
            date = incDate(date)
            continue

        NDVI_BRDF = numpy.zeros((5, y_123, x_123))

        # DATA[17][8][y_123][x_123] = [hour][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][y_123][x_123]
        DATA = importBRDFdata(inputDir, date, processData)
        if(len(DATA) == 0):
            print('..Error! Don\'t is possible create %s_BRDF.tif' %(date))
            date = incDate(date)
            continue
        print('...Create BRDF data: %s'%(date))
        tac = time.clock()
        hours = len(DATA)
        #print('hours:%d [%d, %d]'%(hours, y_123, x_123))
        for x in range(x_123):
            for y in range(y_123):
                if(DATA[0][0][y][x] == 0):
                    continue

                # Copy data to vector process
                VZA = numpy.zeros((hours,1))
                VAZ = numpy.zeros((hours,1))
                SZA = numpy.zeros((hours,1))
                SAZ = numpy.zeros((hours,1))
                RAA = numpy.zeros((hours,1))
                RED = numpy.zeros((hours,1))
                NIR = numpy.zeros((hours,1))
                CLM = 2*numpy.ones((hours,1))
                BAND = numpy.zeros((hours,2))
                HOUR = numpy.zeros((hours,1)) 
                NDVI = numpy.zeros((hours,1))
                hour = 1300

                for h in range (hours):
                    if(DATA[h][0][y][x] != 0):
                        VZA[h]  = DATA[h][0][y][x]
                        VAZ[h]  = DATA[h][1][y][x]
                        SZA[h]  = DATA[h][2][y][x]
                        SAZ[h]  = DATA[h][3][y][x]
                        RED[h]  = DATA[h][4][y][x]
                        NIR[h]  = DATA[h][5][y][x]
                        CLM[h] = DATA[h][7][y][x]
                        HOUR[h] = hour
                        BAND[h,0] = RED[h]
                        BAND[h,1] = NIR[h]
                        NDVI[h] = calcNDVI(NIR[h], RED[h])
                        # WIR = DATA[h][0][y][x]
                        hour = incHour(hour)
                RAA  = SAZ - VZA
                # Verify number of samples to process BRDF
                if(processData == '123'):
                    passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(RED > 0, NIR > 0))
                else:
                    passer = numpy.logical_and(RED > 0, NIR > 0)
                samples = len(passer[passer == True])
                if(samples == 0):
                    continue
                '''    
                NDVI_mean = numpy.mean(NDVI)
                SAZ_mean  = numpy.mean(SAZ[passer1])
                print('NDVI_mean: %.4f SAZ_mean: %.4f'%(NDVI_mean, SAZ_mean)) 
                if(processData == '123'): 
                    passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(RED > 0, numpy.logical_and(NIR> 0, \
                         numpy.logical_and (NDVI > 0.8*NDVI_mean, numpy.logical_and(NDVI < 1.2*NDVI_mean, \
                         numpy.logical_and (SAZ > 0.2*SAZ_mean, SAZ < 1.6*SAZ_mean))))))
                else:
                    passer = numpy.logical_and(RED > 0, numpy.logical_and(NIR> 0, \
                         numpy.logical_and (NDVI > 0.8*NDVI_mean, numpy.logical_and(NDVI < 1.2*NDVI_mean, \
                         numpy.logical_and (SAZ > 0.2*SAZ_mean, SAZ < 1.6*SAZ_mean)))))
                samples = len(passer[passer == True])
                if(samples == 0):
                    print('...continue in samples == 0 second test')
                    date = incDate(date)
                    continue
                '''
                VZA = VZA[passer]
                SZA = SZA[passer]
                RAA = RAA[passer]
                VAZ = VAZ[passer]
                SAZ = SAZ[passer]
                RED = RED[passer]
                NIR = NIR[passer]
                CLM = CLM[passer]
                NDVI = NDVI[passer]
                HOUR = HOUR[passer]
                if (samples >= 5):                  
                    BAND = numpy.zeros((len(RED),2))
                    BAND[:,0] = RED
                    BAND[:,1] = NIR
                else:
                    BAND = numpy.zeros((5,2))
                    for i in range(len(RED)):
                        BAND[i,0] = RED[i]
                        BAND[i,1] = NIR[i]                              

                id_brdf   = 0 
                red_brdf  = 0
                nir_brdf  = 0
                ndvi_brdf = 0
                hour_brdf = 0

                geo_kernel = 'Sparse'
                vol_kernel = 'Thick'
                flag_kernel = False
                if(samples >= 5):
                    #print('--------- Calc NDVI using Ambrals (2000)')
                        
                    # Generate the semiempirical kernels
                    K_obs =  Kernels( VZA, SZA, RAA, LiType=geo_kernel, doIntegrals=False, \
                        normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, RossType= vol_kernel )
                    kern = numpy.ones (( numpy.sum(passer==True), 3 )) # Store the kernels in an array
                    #print(K_obs)
                    kern[ :, 1 ] = K_obs.Ross
                    kern[ :, 2 ] = K_obs.Li
                    NDVI = calc_lstsq(BAND, kern)
                    # If kernel is negative not use in the calc
                    if(len(NDVI) > 1):
                        # Did calc BRDF in this pixel
                        flag_kernel = True
                        samples = len(NDVI[0])
                        # Calc best NDVI
                        id_brdf   = numpy.argmax(NDVI[2]) 
                        red_brdf  = NDVI[0, id_brdf]
                        nir_brdf  = NDVI[1, id_brdf]
                        ndvi_brdf = NDVI[2, id_brdf]
                        hour_brdf = HOUR[id_brdf]
                        qual_brdf = 0
                        # save BRDF KERNEL data to process future day
                        id_ord = numpy.argsort(NDVI[2])
                        id_ord = id_ord[::-1]
                        n = 0
                        for i in id_ord[:5]:
                            KERNEL[y,x,n] = [kern[i, 0], kern[i, 1], kern[i,2], RED[i], NIR[i]]
                            n = n + 1
                    else:
                        continue
                if(flag_kernel == False):
                    if(samples > 5):
                        samples = 3
                    if(numpy.mean(KERNEL[y,x,0]) == 0 or numpy.mean(KERNEL[y,x,1]) == 0):
                        continue
                    # Generate the semiempirical kernels
                    K_obs =  Kernels( VZA, SZA, RAA, LiType=geo_kernel, doIntegrals=False, \
                        normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, RossType= vol_kernel )
                    kern = numpy.ones ((5, 3)) # Store the kernels in an array
                    #print(K_obs)
                    for i in range(len(RED[:samples])):
                        kern[i,1] = K_obs.Ross[i]
                        kern[i,2] = K_obs.Li[i]
                    for i in range(5):
                        if(samples <= i):
                            kern[i] = KERNEL[y,x,i,:3]
                            BAND[i] = KERNEL[y,x,i,3:]
                    NDVI = calc_lstsq(BAND, kern)
                     # If kernel is negative not use in the calc
                    if(len(NDVI) > 1):
                        #samples = len(NDVI[0])
                        # Calc best NDVI in this day
                        id_brdf   = numpy.argmax(NDVI[2, :samples])
              
                        red_brdf  = NDVI[0, id_brdf]
                        nir_brdf  = NDVI[1, id_brdf]
                        ndvi_brdf = NDVI[2, id_brdf]
                        hour_brdf = HOUR[id_brdf]
                        qual_brdf = 1
                    else:
                        continue
                '''        
                print('samples: %d id_brdf: %d'%(samples, id_brdf))
                print('y: %d x: %d date %s'%(y, x, date))
                print('BAND')
                print(BAND[:samples])
                print('kern:')
                print(kern[:samples]) 
                print('samples:')
                print(samples)
                print('NDVI')
                print(NDVI)
                print('HOUR:')
                print(HOUR[:samples])
                '''
                NDVI_BRDF[0][y][x] = ndvi_brdf
                NDVI_BRDF[1][y][x] = red_brdf
                NDVI_BRDF[2][y][x] = nir_brdf
                NDVI_BRDF[3][y][x] = hour_brdf
                NDVI_BRDF[4][y][x] = qual_brdf
        #print('.Calc data time: %.2f s' %(time.clock() - tac))
        saveRASTERfile(outputFileBRDF, fileMask123, NDVI_BRDF)
        # Save kernel file data in 5 to 5 days
        if(calcDays(dateSaveKernel, date) > 5):
            numpy.save(kernelFile, KERNEL)
            dateSaveKernel = date
        #print('.Save data time: %.2f s' %(time.clock() - tac))
        #sys.exit()
        date = incDate(date)


