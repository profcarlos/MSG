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
import psutil
from Py6S import *
from SMAC_CESBIO import *

numpy.set_printoptions(threshold=np.inf)
usage = """\
Usage: %s [OPTIONS]
        -ps     process (BRDF or 6S)
        -ty     type data (1day, 2days)
        -fi     folder input data
        -fo     folder output data
        -py     pixel y to process
        -px     pixel x to process
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
        -class  data class (culture or pasture)

Points to verify data:
        1. [143, 48] [-51.557, -17.376] [Rio Verde] Agricultura Anual  (Terraclass 1)
        2. [178, 14] [-52.751, -18.624] [Chapadão do Céu]
        3. [155, 32] [-52.105, -17.820] [Jatai]
        4. [144, 38] [-52.084, -17.668] [Jatai]
        5. [146, 33] [-52.068, -17.510] [Perolandia]
        6. [146, 58] [-51.202, -17.504] [Montividiu]
        7. [139, 57] [-51,256, -17.265] [Montividiu]
        8. [144, 67] [-50.916, -17.439] [Rio Verde]
        9. [134, 79] [-50.501, -17.089] [Parauna]
        10. [130, 63] [-51.028, -16.949] [Parauna]
        [84, 122] [-48.990, -15.330] Agricultura Perene (Terraclass 2)
        [16, 172] Natural Campestre  (Terraclass 5)
        [32, 159] Natural Florestal  (Terraclass 6)
        11. [35,  81] [-50.427, -13.634] [Novo Mundo] Pastagem           (Terraclass 11)
        12. [8 ,  83] [-50.361, -12.707] [São Miguel do Araguaia] 
        13. [32, 118] [-49.122, -13.531] [Porangatu]
        14. [34, 123] [-48.969, -13.607] [Santa Tereza de Goias]
        15. [50, 83]  [-50.352, -14.161] [Nova Crixas]
        16. [123, 81] [-50.426, -16.701] [Aurilandia]
        17. [130, 88] [-50.163, -16.923] [Jandaia]
        18. [156, 108][-49.460, -17.830] [Joviania]
        19. [15, 80]  [-50.436, -12.941] [Sao Miguel do Araguaia]
        20. [44, 188] [-46,698, -13,956] [Iaciara]
        [43, 177] Natural Savânica   (Terraclass 12) *** EXCLUIR DAS AMOSTRAS DE PASTAGEM ***
""" % argv[0]

CULTURE = [[143,48],[178,14],[155,32],[144,38],[146,33],[146,58],[139,57],[144,67],[134,79],[130,63]]
PASTURE = [[35,81],[8,83],[32,118],[34,123],[50,83],[123,81],[130,88],[156,108], [15, 80], [44,188]]


def process6Sdata(RED, NIR, ATM, GEOM):

    if(ATM[0] != -1 and ATM[1] != -1 and ATM[2] != -1):
        # process 6S data
        print('WV: %.4f OZONE: %.4f AOT: %.2f'% (ATM[0], ATM[1], ATM[2]))
        print('SZA: %.2f SAZ %.2f VZA: %.2f VAZ: %.2f'%(GEOM[0], GEOM[1], GEOM[2], GEOM[3]))
        OUT_6S = sixsv1(RED, NIR, ATM, GEOM)
        if(OUT_6S[0] != -1 and OUT_6S[1] != -1):
            return OUT_6S
        else:
            print('...Error in OUT_6S data')
            sys.exit(1)
    else:
        print('...Error in ATM read data')
        sys.exit(1)

def sixsv1(RED, NIR, ATM, GEOM):
    
    # Config 6SV1
    [WV, OZONE, AOT] = ATM
    [SZA, SAZ, VZA, VAZ] = GEOM
    try:
        s = SixS(os.getcwd() + "\\sixsV1_1_lab.exe")
    except:
        printlog(True, ' Except in sixv1 path.')
        s.produce_debug_report()

    s.altitudes.set_sensor_satellite_level()
    s.altitudes.set_target_sea_level()
    s.geometry = Geometry.User()
    s.geometry.solar_z = SZA
    s.geometry.solar_a = SAZ
    s.geometry.view_z = VZA
    s.geometry.view_a = VAZ
    s.geometry.day = 1
    s.geometry.month = 1
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(WV, OZONE)
    s.aot550 = AOT


    # Calc to RED band
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(RED) 
    s.wavelength = Wavelength (0.56, 0.71)  

    try:
        s.run()
    except:
        printlog(True, ' Error in Py6S (RED)')
        OUT_6S =  [-1, 0, 0]
        return OUT_6S
    
    RED = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out RED: %.2f'%(RED))

    # Calc to NEAR INFRARED band    
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(NIR) 
    s.wavelength = Wavelength (0.74, 0.88)
    try:
        s.run()
    except:
        printlog(True, ' Error in Py6S (NIR)')
        OUT_6S =  [0, -1, 0]
        return OUT_6S
    #s.produce_debug_report()
    #time.sleep(random.randint(1,10))
    NIR = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out NIR: %.2f'%(NIR))
    NDVI = calcNDVI(NIR, RED)
    OUT_6S = [RED, NIR, NDVI]
    #printlog(True, OUT_6S)
    return OUT_6S

class ATMclass:

    def __init__(self):
        #print('ATMclass init')
        self.date = 0
        self.i = 0
        self.j = 0
        self.inputDir = ''
        self.DATA = [0, 0, 0]

    def read(self, inputDir, date, i, j):
        #print('ATMclass read')
        # Read ATM file
        # band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
        # band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
        # band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
        # band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)
        # ATM to 6S: Vector with atmospheric data [WV O AOT]
        if(self.date == date and self.i == i and self.j == j):
            print('----------- return self.DATA')
            return self.DATA
        
        self.inputDir = inputDir
        self.date = date
        self.i = i
        self.j = j
        inputSubDir = self.inputDir + '___\\' + str(self.date)[0:4] + '\\' + str(self.date)[4:6] + '\\'
        strFile = inputSubDir.replace('___', 'ATM') + str(self.date) + '_ATM.tif'
        DATA = numpy.zeros((3))
        if os.path.isfile(strFile) == 0:
            print('Error to read ATM file: %s' % (strFile))
            DATA = [-1, -1, -1]
        else:
            #print('Read ATM file: %s' % (strFile))
            '''
            # get Aerosol Data Band
            vet = readFileBand(strFile, 4)
            if(not len(vet)):
                return [0, 0, 0]
            else:
                DATA[2] = vet[i][j]
            '''
            # get Water_Vapour  
            vet = readFileBand(strFile, 10)
            if(not len(vet)):
                DATA[0] = -1
            else:
                DATA[0] = vet[i][j] 

            # get Total Ozone Data Band
            vet = readFileBand(strFile, 7)
            if(not len(vet)):
                DATA[1] = -1
            else:
                DATA[1] = vet[i][j]
                DATA[1] = DATA[1]/1000 # Convert Dobson Unit to cm-atm
            # get Deep_Blue Data Band
            vet = readFileBand(strFile, 1)
            if(not len(vet)):
                DATA[2] = -1
            else:
                DATA[2] = vet[i][j]
        self.DATA = DATA
        return self.DATA


def importBRDFdata(inputDir, outputDir, dateStart, dateEnd, i, j):
    startHour = 1300
    endHour   = 1700    
    date = dateStart
    num_days = calcDays(dateStart, dateEnd)
    if(num_days == -1):
        print(' Error in dates, please verify.')
        return []
    data = numpy.zeros((int(num_days*(endHour - startHour)/20) , 11))
    h = 0

    while(date <= dateEnd):

        print(' Reading data files: %s' %(date))
        inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        hour = startHour
        while(hour <= endHour):
            # Flag to read error file
            fileError = False
            #print('Reading hour: %d' %hour)
            data[h][0] = date
            data[h][1] = hour
            
            # --------- Process ANGLE file
            #print('Process ANGLE file')
            strFile = inputSubDir.replace('___', 'ANGLE') + str(date) + str(hour) + '_ANGLES.tif'
            if os.path.isfile(strFile) == 0:
                #print('..Error! Invalid data file: %s' %strFile)
                fileError = True
            else:
                # get VAZ View Azimuth Angle (msg_azres.tif)
                vet = readFileBand(strFile, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][4] = vet[i][j]
                # get SAZ Sun Azimuth Angle (sol_azres.tif)
                vet = readFileBand(strFile, 2)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][6] = vet[i][j]
                # get VZA View Zenith Angle (msg_zenres.tif)
                vet = readFileBand(strFile, 3)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][3] = vet[i][j]
                # get SZA Sun Zenith Angle (sun_zenres.tif)
                vet = readFileBand(strFile, 4)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][5] = vet[i][j]

            # --------- Process 123 file
            if fileError == False:
                #print('Process 123 file')
                strFile = inputSubDir.replace('___', '123') + str(date) + str(hour) + '_123.tif'
                if os.path.isfile(strFile) == 0:
                    #print('..Error! Invalid data file: %s' %strFile)
                    fileError = True
                else:
                    # get RED Data Band
                    vet = readFileBand(strFile, 1)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][7] = vet[i][j]
                    # get NIR Data Band
                    vet = readFileBand(strFile, 2)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][8] = vet[i][j]
                    # get WIR Data Band
                    vet = readFileBand(strFile, 3)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][9] = vet[i][j]

            # --------- Process CLM file
            if fileError == False:
                #print('Process CLM file')
                strFile = inputSubDir.replace('___', 'CLM') + str(date) + str(hour) + '_CLM.tif'
                if os.path.isfile(strFile) == 0:
                    #print('..Error! Invalid data file: %s' %strFile)
                    fileError = True
                else:
                    vet = readFileBand(strFile, 1)
                    if(not len(vet)):
                        fileError = True
                    else:
                        data[h][2] = vet[i][j] 
            # Verify if error clear data, else inc position
            if(fileError == False):
                h = h + 1
            hour = incHour(hour)
        date = incDate(date)
    # data[1:16] = [VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM]
    filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(i) +'_'+ str(j) +  '_BRDF.csv'
    print(' Save output file: %s'%(outputDir + filename))
    numpy.savetxt(outputDir + filename, data, fmt='%1.4f', comments = '', header = 'DATE HOUR CLM VZA VAZ SZA SAZ RED NIR WIR', delimiter = ' ', newline='\n')
    return data

def importMOD13data(inputDir, outputDir, dateStart, dateEnd, i, j):
    
    date = dateStart
    num_days = calcDays(dateStart, dateEnd)
    if(num_days == -1):
        print(' Error in dates, please verify.')
        return []
    data = numpy.zeros((int(num_days/16+1) , 4))
    h = 0
    # MOD13 band_names = ['NDVI', 'RED', 'NIR', 'MIR', 'day_of_the_year', 'pixel_reliability']
    while(date <= dateEnd):
        inputSubDir = inputDir + '___\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        strFile = inputSubDir.replace('___', 'MOD13') + str(date) + '_MOD13.tif'
        if os.path.isfile(strFile):
            print(' Reading data files: %s' %(date))
            # get DOY data band
            vet = readFileBand(strFile, 5)
            if(not len(vet)):
                return []
            else:
                data[h][0] = vet[i][j]
            # get RED data band
            vet = readFileBand(strFile, 2)
            if(not len(vet)):
                return []
            else:
                data[h][2] = vet[i][j]
            # get NIR data band  
            vet = readFileBand(strFile, 3)
            if(not len(vet)):
                return []
            else:
                data[h][3] = vet[i][j]
            # get NDVI data band
            vet = readFileBand(strFile, 1)
            if(not len(vet)):
                return []
            else:
                if(vet[i][j] != 0):
                    data[h][1] = vet[i][j]
                # if NDVI = 0
                else:
                    if(data[h][3] + data[h][2] != 0):
                        data[h][1] = (data[h][3] - data[h][2])/(data[h][3] + data[h][2])
                    else:
                       data[h][1] = 0 
            h = h + 1
            #print(data)
        #else:
            #print('...Error to read: %s'%(strFile))
        date = incDate(date)
    # data = [DOY, NDVI, RED, NIR]
    filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(i) +'_'+ str(j) +  '_MOD13.csv'
    print(' Save output file: %s'%(outputDir + filename))
    numpy.savetxt(outputDir + filename, data, fmt='%1.4f', comments = '', header = 'DOY NDVI RED NIR', delimiter = ' ', newline='\n')
    return data


class main():

    argDict = mapDict(argv, usage)
    gc.collect()
    
    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict and (("-px" in argDict and "-py" in argDict) or "-class" in argDict):
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData = argDict["-ty"]
        if("-class" in argDict):
            classData = argDict["-class"]
        else:            
            MSG_pix_x = int(argDict["-px"])
            MSG_pix_y  = int(argDict["-py"])
            classData = None
    else:
        exit(usage)
    if(typeData != '1day' and typeData != '2days'):
        print('Error in typeData, verify!')
        exit(usage)

    if(classData == 'culture'):
        CLASS = CULTURE
    elif(classData == 'pasture'):
        CLASS = PASTURE
    elif(classData == None):
        CLASS = [[MSG_pix_x, MSG_pix_y]]
    else:
        exit(usage)

    # Go out if folder not exist
    if(os.path.isdir(inputDir) == False or os.path.isdir(outputDir) == False):
        print('... Error in folder adress!')
        sys.exit()
    # Create class ATM to get atmospheric data
    ATMdata = ATMclass()
    num_days = calcDays(dateStart, dateEnd)
    filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) +  '_class_' + str(classData) + '_typ_' + str(typeData) + '_NDVI_ALL_.csv'
    header = 'DATE DOY HOUR_1330 MAX_DAY MIN_DAY MEAN_DAY MAX_SPARSE_THICK MIN_SPARSE_THICK MEAN_SPARSE_THICK KERN_FAIL_SPARCE_THICK 6S_SPARCE_THICK MAX_GAO MIN_GAO MEAN_GAO KERN_FAIL_GAO 6S_GAO'
    data_ALL = datafile('BRDF', outputDir, filename, header, num_days*10, 16)

    for pos_class in range(len(CLASS)):
        MSG_pix_x = CLASS[pos_class][0]
        MSG_pix_y = CLASS[pos_class][1]
        
        # Create BRDF data base
        filenameBRDF = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) +  '_typ_' + str(typeData) + '_NDVI_.csv'
        data_BRDF = datafile('BRDF', outputDir, filenameBRDF, header, num_days, 16)
 
        # Verify and/or process BRDF file
        filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_BRDF.csv'
        if(os.path.isfile(outputDir + filename)):
            print('... Load data file: %s'%(filename))
            DATA = numpy.loadtxt(outputDir + filename, delimiter = ' ', skiprows = 1)
        else:
            print('... Create data file: %s'%(filename))
            DATA = importBRDFdata(inputDir, outputDir, dateStart, dateEnd, MSG_pix_x, MSG_pix_y)
        
        # Read or create MOD13 data file
        filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_MOD13.csv'
        if(os.path.isfile(outputDir + filename)):
            print('... Load data file: %s'%(filename))
            MOD13 = numpy.loadtxt(outputDir + filename, delimiter = ' ', skiprows = 1)
        else:
            print('... Create data file: %s'%(filename))
            MOD13 = importMOD13data(inputDir, outputDir, dateStart, dateEnd, MSG_pix_x, MSG_pix_y)
            # To calc only MOD13 data
        #continue

        # Verify if files exist
        if(os.path.isfile(outputDir + filename.replace('.csv', 'BRDF.csv'))):
            print('... Data files NDVI exist! Verify!')
            continue

        DATE = DATA[:,0]
        HOUR = DATA[:,1]
        CLM  = DATA[:,2]
        VZA  = DATA[:,3]
        SZA  = DATA[:,5]
        RAA  = DATA[:,4] - DATA[:,6]
        VAZ  = DATA[:,4]
        SAZ  = DATA[:,6]
        BAND = DATA[:,7:10]
        RED  = DATA[:,7]
        NIR  = DATA[:,8]

        date = dateStart
        while(date <= dateEnd):
            doy_today = calcDoy(int(str(date)[0:4]), int(str(date)[4:6]), int(str(date)[6:8]))
            # Select number of days to create NDVI
            date_up   = date    
            if(typeData == '2days'):
                date_up   = incDate(date)
            passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))
            passer2= numpy.logical_and (CLM < 1.05, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(HOUR == 1330, numpy.logical_and(RED > 0, NIR > 0))))
            samples = len(passer[passer == True])
            if(samples < 5):
                date = incDate(date)
                continue
            VZA_filt = VZA[passer]
            SZA_filt = SZA[passer]
            RAA_filt = RAA[passer]
            VAZ_filt = VAZ[passer]
            SAZ_filt = SAZ[passer]
            RED_filt = RED[passer]
            NIR_filt = NIR[passer]
            BAND_filt = BAND[passer]
            print('............ doy: %d' %(doy_today))
            # Mount RESULT matrix
            data_BRDF.addAtrib([date, doy_today])
            data_ALL.addAtrib([date, doy_today])
            ATM = ATMdata.read(inputDir, date, MSG_pix_x, MSG_pix_y);
            # NDVI Sample 1330
            sample_1330 = len(passer2[passer2 == True])
            if(sample_1330 == 1):
                NDVI_1330 = calcNDVI(NIR[passer2], RED[passer2])
                data_BRDF.addAtrib(NDVI_1330)
                data_ALL.addAtrib(NDVI_1330)
                '''
                print('NDVI_1330')
                print(NDVI_1330)
                print(NIR[passer2])
                print(RED[passer2])
                input()
                '''
            else:
                data_BRDF.addAtrib(0)
                data_ALL.addAtrib(0)           
            # NDVI Mean day
            NDVI = numpy.zeros((3, samples))    #[RED NIR NDVI]
            NDVI[0] = RED_filt
            NDVI[1] = NIR_filt
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])

            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])])
            '''
            print('NDVI: %.4f %.4f %.4f'%(numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])))  
            print('... Array NDVI')
            print(NDVI)
            input()
            '''
            print('... Ross-Li Sparce-Thick Normal') 
            geo_kernel = 'Sparse'
            vol_kernel = 'Thick'
            # Generate the kernels, only bother for obs where the QA is OK
            K_obs =  Kernels( VZA_filt, SZA_filt, RAA_filt, \
                LiType=geo_kernel, doIntegrals=False, \
                normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                RossType= vol_kernel )
            kern = numpy.ones (( numpy.sum(passer==True), 3 )) # Store the kernels in an array
            kern[ :, 1 ] = K_obs.Ross
            kern[ :, 2 ] = K_obs.Li
            '''
            print ("%s\t%s\t%s\t%s\t%s" % ( "NUM", "PASSER", "CLM", "DATE", "KERN"))
            for x in range(len(passer[passer == True])):
                if(passer[x] == True):
                    print ("%2d\t%s\t%d\t%d\t%6f\t%6f\t%6f" % (x, passer[x], CLM[x], DATE[x], kern[x][0], kern[x][1], kern[x][2]))
            '''    
            # Calc to RED and NIR
            NDVI = numpy.zeros((3, samples))
            COV  = [0,0]
            kern_fails = 0
            for band in range(2):
                #print('passer: %d  kern: %d'%(len(passer), len(kern)))
                tmp_obs = BAND_filt
                obs_col = numpy.zeros((len(tmp_obs),1))
                obs_lin = numpy.zeros((len(tmp_obs)))
                for x in range(len(tmp_obs)):
                    obs_col[x] = tmp_obs[x][band]
                    obs_lin[x] = tmp_obs[x][band]

                if(len(kern) >= len(passer)):
                    K = kern[passer, :]
                else:
                    K = kern
                (f, rmse, rank, svals ) = numpy.linalg.lstsq( K, obs_col )
                print('----- Results Band[%d]' %(band))
                print ("%-20s %20s" % ( "Kernel", "Value"))
                for i, k in enumerate( ["Isotropic", "Ross-" + vol_kernel, "Li-" + geo_kernel] ):
                    print ("%-20s %20f" % ( k, f[i] ))
                if(f[0] < 0):
                    kern_fails = (band + 1) + kern_fails
                fwd = K.dot(f)
                fwd_lin = numpy.zeros(len(fwd))
                for x in range(len(fwd)):
                    fwd_lin[x] = fwd[x]
                COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
                NDVI[band] = fwd_lin
            print('kern_fails: %d'%(kern_fails))
            # Calc NDVI = (NIR - RED)/(NIR + RED)
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_BRDF.addAtrib(kern_fails)
            data_ALL.addAtrib(kern_fails)
            error = rmse
            '''
            print('--- NDVI BRDF')
            print(NDVI)
            '''
            print(' error: \t%f'%(error))
            # Process 6S correction in BRDF data
            id = numpy.argmax(NDVI[2]) 
            RED_BRDF  = NDVI[0][id]
            NIR_BRDF  = NDVI[1][id]
            NDVI_BRDF = NDVI[2][id]

            GEOM = [SZA_filt[id], SAZ_filt[id], VZA_filt[id], VAZ_filt[id]]
            [RED_6S, NIR_6S, NDVI_6S] = [0, 0, 0] #process6Sdata(RED_BRDF, NIR_BRDF, ATM, GEOM)
            data_BRDF.addAtrib(NDVI_6S)
            data_ALL.addAtrib(NDVI_6S)

            print('... Ross-Li Sparce-Thick Filter Data') 
            print('all samples: %d'%samples)
            
            # passer to +/- 20% NDVI variation
            NDVI = numpy.zeros((3, samples))    #[RED NIR NDVI]
            NDVI[0] = RED_filt
            NDVI[1] = NIR_filt
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
            NDVI = NDVI[2]
            NDVI_temp = np.partition(-NDVI, 4)
            NDVI_max = -NDVI_temp[:4]
            NDVI_mean = np.mean(NDVI_max)
            print('NDVI:')
            print(NDVI)
            print('NDVI_max:')
            print(NDVI_max)
            print('NDVI_mean: %.4f' %NDVI_mean) 
            
            NDVI = calcNDVI(NIR, RED)
            #print('NDVI_all')
            #print(NDVI)
            passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(NDVI > 0.2*NDVI_mean, numpy.logical_and(NDVI < 1.8*NDVI_mean, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))))

            SAZ_mean = numpy.mean(SAZ_filt)
            #passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(SAZ > 0.4*SAZ_mean, numpy.logical_and(SAZ < 1.6*SAZ_mean, numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))))
            passer = numpy.logical_and (CLM < 1.05, numpy.logical_and(NDVI > 0.8*NDVI_mean, numpy.logical_and(NDVI < 1.2*NDVI_mean, numpy.logical_and(SAZ > 0.2*SAZ_mean, numpy.logical_and(SAZ < 1.8*SAZ_mean,numpy.logical_and(numpy.logical_or(DATE==date, DATE == date_up), numpy.logical_and(RED > 0, NIR > 0)))))))

            samples = len(passer[passer == True])
            if(samples < 5):
                print('samples: %d'%(samples))
                date = incDate(date)
                data_BRDF.addAtrib([0,0,0,-1,0]) 
                data_ALL.addAtrib([0,0,0,-1,0]) 
                data_BRDF.addSample() 
                data_ALL.addSample()
                #input()
                continue
            VZA_filt = VZA[passer]
            SZA_filt = SZA[passer]
            RAA_filt = RAA[passer]
            VAZ_filt = VAZ[passer]
            SAZ_filt = SAZ[passer]
            RED_filt = RED[passer]
            NIR_filt = NIR[passer]
            BAND_filt = BAND[passer]

            '''
            geo_kernel = 'Sparse'
            vol_kernel = 'Thick'
            # Generate the kernels, only bother for obs where the QA is OK
            K_obs =  Kernels( VZA_filt, SZA_filt, RAA_filt, \
                LiType=geo_kernel, doIntegrals=False, \
                normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                RossType= vol_kernel )
            kern = numpy.ones (( numpy.sum(passer==True), 3 )) # Store the kernels in an array
            kern[ :, 1 ] = K_obs.Ross
            kern[ :, 2 ] = K_obs.Li

            print ("%s\t%s\t%s\t%s\t%s" % ( "NUM", "PASSER", "CLM", "DATE", "KERN"))
            for x in range(len(passer[passer == True])):
                if(passer[x] == True):
                    print ("%2d\t%s\t%d\t%d\t%6f\t%6f\t%6f" % (x, passer[x], CLM[x], DATE[x], kern[x][0], kern[x][1], kern[x][2]))
   
            # Calc to RED and NIR
            NDVI = numpy.zeros((3, samples))
            COV  = [0,0]
            kern_fails = 0
            for band in range(2):
                #print('passer: %d  kern: %d'%(len(passer), len(kern)))
                tmp_obs = BAND_filt
                obs_col = numpy.zeros((len(tmp_obs),1))
                obs_lin = numpy.zeros((len(tmp_obs)))
                for x in range(len(tmp_obs)):
                    obs_col[x] = tmp_obs[x][band]
                    obs_lin[x] = tmp_obs[x][band]

                if(len(kern) >= len(passer)):
                    K = kern[passer, :]
                else:
                    K = kern
                (f, rmse, rank, svals ) = numpy.linalg.lstsq( K, obs_col )
                print('----- Results Band[%d]' %(band))
                print ("%-20s %20s" % ( "Kernel", "Value"))
                for i, k in enumerate( ["Isotropic", "Ross-" + vol_kernel, "Li-" + geo_kernel] ):
                    print ("%-20s %20f" % ( k, f[i] ))
                if(f[0] < 0):
                    kern_fails = (band + 1) + kern_fails
                fwd = K.dot(f)
                fwd_lin = numpy.zeros(len(fwd))
                for x in range(len(fwd)):
                    fwd_lin[x] = fwd[x]
                COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
                NDVI[band] = fwd_lin
            print('kern_fails: %d'%(kern_fails))
            # Calc NDVI = (NIR - RED)/(NIR + RED)
            NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])]) 
            data_BRDF.addAtrib(kern_fails)
            data_ALL.addAtrib(kern_fails)
            error = rmse

            print('--- NDVI BRDF')
            print(NDVI)

            print(' error: \t%f'%(error))
            # Process 6S correction in BRDF data
            id = numpy.argmax(NDVI[2]) 
            RED_BRDF  = NDVI[0][id]
            NIR_BRDF  = NDVI[1][id]
            NDVI_BRDF = NDVI[2][id]

            '''
            # Calc NDVI using Gao et al. (2002)
            #print('---------------- Calc NDVI using Gao et al. (2002)')
            print('... Gao (2002)')
            NDVI_calc = numpy.zeros((3, samples))
            NDVI_calc[0] = RED[passer]
            NDVI_calc[1] = NIR[passer]
            NDVI_calc[2] = (NDVI_calc[1] - NDVI_calc[0])/(NDVI_calc[1] + NDVI_calc[0])
            NDVI_mean = numpy.mean(NDVI_calc[2])
            W_ini   = (NDVI_calc[2]/NDVI_mean)**2 
            geo_kernel = 'Sparce'
            vol_kernel = 'Thick'
            error = 1
            loop = 0
            while(error > 0.0001 and loop <= 5):

                K_obs =  Kernels( VZA_filt, SZA_filt, RAA_filt, \
                    LiType=geo_kernel, doIntegrals=False, \
                    normalise=1, RecipFlag=True, RossHS=False, MODISSPARSE=True, \
                    RossType= vol_kernel )
                # Store the kernels in an array
                kern = numpy.ones (( numpy.sum(passer==True), 3 )) 
                kern[ :, 1 ] = K_obs.Ross
                kern[ :, 2 ] = K_obs.Li
                # Calc to RED and NIR
                NDVI = numpy.zeros((3, samples))
                COV  = [0,0]
                kern_fails = 0
                for band in range(2):
                    #print('passer: %d  kern: %d'%(len(passer), len(kern)))
                    tmp_obs = BAND_filt
                    obs_col = numpy.zeros((len(tmp_obs),1))
                    obs_lin = numpy.zeros((len(tmp_obs)))
                    for x in range(len(tmp_obs)):
                        obs_col[x] = tmp_obs[x][band]
                        obs_lin[x] = tmp_obs[x][band]

                    if(len(kern) >= len(passer)):
                        K = kern[passer, :]
                    else:
                        K = kern
                    # Using Weight Inversion
                    K = K *numpy.sqrt(W_ini)[:,None]
                    #print('--- obs_col')
                    #print(obs_col)
                    #print('--- W_ini')
                    #print(W_ini)
                    #print('--- K')
                    #print(K)
                    obs_col = obs_col*numpy.sqrt(W_ini)[:,None]
                    #print('--- obs_col')
                    #print(obs_col)
                    (f, rmse, rank, svals ) = numpy.linalg.lstsq( K, obs_col )
                    print('----- Results Band[%d]' %(band))
                    print ("%-20s %20s" % ( "Kernel", "Value"))
                    for i, k in enumerate( ["Isotropic", "Ross-" + vol_kernel, "Li-" + geo_kernel] ):
                        print ("%-20s %20f" % ( k, f[i] ))
                    fwd = K.dot(f)
                    fwd_lin = numpy.zeros(len(fwd))
                    for x in range(len(fwd)):
                        fwd_lin[x] = fwd[x]
                    COV[band]  = numpy.corrcoef (obs_lin, fwd_lin)[1,0]
                    NDVI[band] = fwd_lin
                    if(f[0] < 0):
                        kern_fails = (1 + band) + kern_fails
                # Calc NDVI = (NIR - RED)/(NIR + RED)
                NDVI[2] = calcNDVI(NDVI[1], NDVI[0])
                W_obs = (NDVI[2]/NDVI_calc[2])**2
                error = 0
                #print('---- W_ini')
                #print(W_ini)
                #print('---- W_obs')
                #print(W_obs)
                for i in range(len(W_ini)):
                    error = error + (W_obs[i] - W_ini[i])**2
                W_ini = W_obs
                print(' error: %f in loop %d'%(error, loop))
                loop = loop + 1
            print('kern_fails: %d'%(kern_fails))
            data_BRDF.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])])
            data_ALL.addAtrib([numpy.amax(NDVI[2]), numpy.amin(NDVI[2]), numpy.mean(NDVI[2])])
            data_BRDF.addAtrib(kern_fails)
            data_ALL.addAtrib(kern_fails)
            # Process 6S correction in BRDF data
            id = numpy.argmax(NDVI[2]) 
            RED_BRDF  = NDVI[0][id]
            NIR_BRDF  = NDVI[1][id]
            NDVI_BRDF = NDVI[2][id]
            GEOM = [SZA_filt[id], SAZ_filt[id], VZA_filt[id], VAZ_filt[id]]
            [RED_6S, NIR_6S, NDVI_6S] = [0, 0, 0] #process6Sdata(RED_BRDF, NIR_BRDF, ATM, GEOM)
            data_BRDF.addAtrib(NDVI_6S)
            data_ALL.addAtrib(NDVI_6S)

            data_BRDF.addSample() 
            data_ALL.addSample()
            #input()   
            date = incDate(date)
            if(type == '2days'):
                date = incDate(date)
        data_BRDF.save()
    data_ALL.save()
    