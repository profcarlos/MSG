import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 't:\\carlos\\newpython\\utils')

import numpy
import math
import os
from time import  strftime, localtime, sleep
import sys
import gc
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from time import strftime, localtime
from sys import argv
from utils import mapDict, incHour, incDate, infoFile, printlog,  readFileBand, readMaskGoias
import time
import threading
import queue
from Py6S import *
import psutil
import random


#def __debug__

usage = """\
Usage: %s [OPTIONS]
        -ty     type of process BRDF (to create...)
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
""" % argv[0]

# Constantes
# To print lib version
# lib.__version__
fileGoiasGeo = os.getcwd() + '\\shape\\201408221300.tif'

# Variables of threads
NUM_THREADS = 4
queueLock = threading.Lock()
workQueue = queue.Queue(4*NUM_THREADS)
threads = []
exitFlag = 0

# Process time data
startHour = 1300
endHour   = 1700

# Constants of BRDF Process
br = 1       # b/r used in K_geo, to calc SZA_line and VZA
hb = 1.5       # h/b used in K_geo, to calc cos(t)       


def sixsv1(i,j,VET_BRDF, VET_ATM):
    # VET: Matrix with data [REF VZA VAZ SZA SAZ] in 3 bands data
    # ATM: Vector with atmospheric data [WV O AOT]
    
    # Config 6SV1
    try:
        s = SixS(os.getcwd() + "\\sixsV1_1_dell.exe")
        #s.test()
    except:
        s = SixS(os.getcwd() + "\\sixsV1_1_lab.exe")
        printlog(True, '----- Except in sixv1 path')
        #s.test()
    #s.produce_debug_report()

    # Set the atmosphere profile to be base
    # Set the wavelength to be that of
    # Data file Py6S_artigo_ex3.py
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(VET_ATM[0], VET_ATM[1])
    s.aot550 = VET_ATM[2]

    # Calc to RED band
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(VET_BRDF[0,0]) 
    s.wavelength = Wavelength (0.56, 0.71)  
    s.geometry = Geometry.User()
    s.geometry.view_z  = VET_BRDF[0,1]
    s.geometry.view_a  = VET_BRDF[0,2]
    s.geometry.solar_z = VET_BRDF[0,3]
    s.geometry.solar_a = VET_BRDF[0,4]
    try:
        s.run()
    except:
        printlog(True, '----- Error in Py6S (RED)')
        printlog(True, 'Using Memory: %.2f CPU: %.2f' %(psutil.virtual_memory().percent, psutil.cpu_percent()))
        DATA_BRDF =  [ 0, -1, 0, 0]
        return DATA_BRDF
    
    RED = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out RED: %.2f'%(RED))
    # Calc to NEAR INFRARED band    
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(VET_BRDF[1,0]) 
    s.wavelength = Wavelength (0.74, 0.88)
    s.geometry = Geometry.User()
    s.geometry.view_z  = VET_BRDF[1,1]
    s.geometry.view_a  = VET_BRDF[1,2]
    s.geometry.solar_z = VET_BRDF[1,3]
    s.geometry.solar_a = VET_BRDF[1,4]
    try:
        s.run()
    except:
        printlog(True, '----- Error in Py6S (NIR)')
        printlog(True, 'Using Memory: %.2f CPU: %.2f' %(psutil.virtual_memory().percent, psutil.cpu_percent()))
        DATA_BRDF =  [ 0, 0, -1, 0]
        return DATA_BRDF
    #s.produce_debug_report()
    #time.sleep(random.randint(1,10))
    NIR = float(s.outputs.pixel_reflectance)
    #printlog(True, '----- Out NIR: %.2f'%(NIR))
    # Calc to SHORT INFRARED band    
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(VET_BRDF[2,0]) 
    s.wavelength = Wavelength (1.50, 1.78)
    s.geometry = Geometry.User()
    s.geometry.view_z  = VET_BRDF[2,1]
    s.geometry.view_a  = VET_BRDF[2,2]
    s.geometry.solar_z = VET_BRDF[2,3]
    s.geometry.solar_a = VET_BRDF[2,4]
    try:
        s.run()
    except:
        printlog(True, '----- Error in Py6S (WIR)')
        printlog(True, 'Using Memory: %.2f CPU: %.2f' %(psutil.virtual_memory().percent, psutil.cpu_percent()))
        DATA_BRDF =  [ 0, 0, 0, -1]
        return DATA_BRDF
    #s.produce_debug_report()
    #time.sleep(random.randint(1,10))
    WIR = float(s.outputs.pixel_reflectance)
    if((NIR + RED) < 0.0001 and (NIR + RED) > -0.0001):
        printlog(True,'... Zero division warning')
        NDVI = 9999
        printlog(True, [i, j, NDVI, RED, NIR, WIR])
    else:
        NDVI = float((NIR - RED)/(NIR + RED))

    #printlog(True, '----- Out WIR: %.2f'%(WIR))
    '''
    printlog(False, 'RED: %.4f' %(RED))
    printlog(False, 'NIR: %.4f' %(NIR))
    printlog(False, 'NIR: %.4f' %(WIR))
    printlog(False, 'NDVI: %.4f' %(NDVI))
    '''
    DATA_BRDF = [NDVI, RED, NIR, WIR]
    #printlog(True, DATA_BRDF)
    return DATA_BRDF

class myThread (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        printlog (False, "Starting " + self.name)
        process_data(self.name, self.q)
        printlog (False, "Exiting " + self.name)

def process_data(threadName, q):
    global queueLock
    global workQueue
    global exitFlag
    global DATA_6S
    #printlog(True, "Enter in process_data. exitFlag: " + str(exitFlag))
    if __debug__: print('#')
    
    while not exitFlag:
        #printlog(True, "exitFlag is False")
        queueLock.acquire()
        #printlog(True, "\nqueueLock acquire")
        if not workQueue.empty():
            #printlog(True, "Queue is not empty")
            i = q.get()
            j = q.get()
            VET_BRDF = q.get()
            VET_ATM  = q.get()
            queueLock.release()
            DATA_6S[0,i,j] = -1
            if __debug__: print('l')
            while(DATA_6S[0,i,j] == -1):
                # Return NDVI, RED, NIR, WIR
                [DATA_6S[0,i,j], DATA_6S[1,i,j], DATA_6S[2,i,j], DATA_6S[3,i,j]]  = sixsv1(i, j, VET_BRDF, VET_ATM) 
                if __debug__: print('!')
        else:
            queueLock.release()
        time.sleep(1)
    
def calc_K_vol(VZA, SZA, RAA):
    # Pre-process variable
    EPS = math.acos(math.cos(SZA)*math.cos(VZA) + math.sin(SZA)*math.sin(VZA)*math.cos(RAA))
    # K_vol calc
    #printlog(False,'calc_k_vol:')
    #printlog(False,'SZA: %.2f\tVZA: %.2f\tRAA: %.2f\tEPS: %.2f' %(SZA, VZA, RAA, EPS))
    K_tmp = (math.pi/2 - EPS)*math.cos(EPS) + math.sin(EPS)
    K_vol  = K_tmp/(math.cos(SZA) + math.cos(VZA)) - math.pi/4
    #printlog(False,'K_tmp: %.2f\tK_vol: %.2f' %(K_tmp, K_vol))    
    #print '... BRDF >> K_vol = ' + str(K_vol)
    return K_vol
                         
def calc_K_geo(VZA, SZA, RAA):
    # Pre-process variable
    SZA = math.atan(br*math.tan(SZA))
    VZA = math.atan(br*math.tan(VZA))
    EPS = math.acos(math.cos(SZA)*math.cos(VZA) + math.sin(SZA)*math.sin(VZA)*math.cos(RAA))
    #printlog(False,'calc_k_geo:')
    #printlog(False,'SZA: %.2f\tVZA: %.2f\tRAA: %.2f\tEPS: %.2f' %(SZA, VZA, RAA, EPS))
    D = math.sqrt(math.pow(math.tan(SZA), 2) + math.pow(math.tan(VZA), 2) - 2*math.tan(SZA)*math.tan(VZA)*math.cos(RAA))
    t_tmp = hb*math.sqrt(math.pow(D, 2) + math.pow((math.tan(SZA)*math.tan(VZA)*math.sin(RAA)), 2))
    cost = t_tmp/(1/math.cos(SZA) + 1/math.cos(VZA))
    if(cost > 1 or cost < -1):
        #printlog(False, '---- Warning! cost not in range [-1, 1]. SZA: %.2f VZA: %.2f D: %.2f t_tmp: %.2f cost: %.2f' %(SZA, VZA, D, t_tmp, cost))
        if(cost > 1):
            t = math.acos(1)
            cost = 1
        else:
            t = math.acos(-1)
            cost = -1
    else:    
        t = math.acos(cost)
    O = (1/math.pi)*(t - math.sin(t)*cost)*(1/math.cos(SZA) + 1/math.cos(VZA))
    #printlog(False,'D:  %.2f\tO:  %.2f\tt_tmp: %.2f\tt: %.2f' %(D, O, t_tmp, t))

    # K_geo calc                     
    K_geo = O - 1/math.cos(SZA) - 1/math.cos(VZA) + 1/2*(1 + math.cos(EPS))*(1/math.cos(SZA)*1/math.cos(VZA)) 
    #print '... BRDF >> K_geo = ' + str(K_geo)
    return K_geo
                                    
def calc_BRDF(vet):
    # vet = [VZA VAZ SZA SAZ RED NIR WIR NDVI]
    row = len(vet)
    if (row < 5):
        printlog(True,'Error! Little samples to calc BRDF')
        return 0
    K = numpy.zeros((row,3))
    R = numpy.zeros((4,row))
    BRDF = numpy.zeros((3,row,3))
    reflect = numpy.zeros((3,row))
    COEFS = numpy.zeros((3,3))
    VZA = numpy.zeros((row))
    SZA = numpy.zeros((row))
    RAA = numpy.zeros((row))
    for i in range (0, row, 1):
        #printlog(False,'------- LOOP %d' %(i))
        VZA[i] = math.radians(vet[i][0])
        SZA[i] = math.radians(vet[i][2])
        RAA[i] = math.radians(vet[i][1] - vet[i][3])
        #printlog(True,'VZA: %.2fº (%.2f)\tSZA: %.2fº (%.2f)\tRAA: %.2fº (%.2f)' %(vet[i][0], VZA[i], vet[i][2], SZA[i], vet[i][1] - vet[i][3], RAA[i]))

        R[0][i] = vet[i][4]
        R[1][i] = vet[i][5]
        R[2][i] = vet[i][6]
        R[3][i] = vet[i][7]
        K[i][0] = 1
        K[i][1] = calc_K_vol(VZA[i], SZA[i], RAA[i])
        K[i][2] = calc_K_geo(VZA[i], SZA[i], RAA[i])
    #printlog(True, 'VZA_var: %.2f\t\tSZA_var: %.2f\t\tRAA_var: %.2f' %(numpy.var(VZA), numpy.var(SZA), numpy.var(RAA)))

    W_inic = numpy.zeros((row))
    neg_x = 0
    W_SZA = 1
    W_RAA = 1

    for i in range(0, row, 1):
        W_inic[i] = 1.6 - W_SZA*math.sin(abs(SZA[i])) - W_RAA*math.sin(abs(RAA[i]))
    #printlog(True, 'W_inic:')
    #printlog(True, W_inic)

    for loop in range (0, 4, 1):
        for i in range(0, row, 1):
            K[i][0] = W_inic[i]
        '''
        printlog(False,'----- K [K_iso K_vol K_geo]:')
        printlog(False, K)
        printlog(False,'----- R [R_vis R_nir R_wir]:')
        printlog(False, R)
        printlog(False,'----- lstsq:')
        printlog(False,'x[n]:')
        '''
        for i in range(0, 3, 1):
            x,resid,rank,s = numpy.linalg.lstsq(K, R[i])
            COEFS[i] = x
            '''    
            printlog(False,'----- calc to K and R[%d]:'%(i))   
            printlog(False, x)
            printlog(False,'resid:')
            printlog(False, resid)
            printlog(False,'rank:')
            printlog(False, rank)
            printlog(False,'s:')
            printlog(False, s)
            '''
            #printlog(False,'----- BRDF in R[%d]:' %(i))
            
            if(x[0] < 0):
                neg_x = neg_x + 1
            for j in range(0, row, 1):
                BRDF[i,j] = K[j]*x
                reflect[i,j] = numpy.sum(BRDF[i,j])
                #printlog(False,'BRDF: [%.4f %.4f %.4f]\treflect: [%.4f]' %(BRDF[i,j,0], BRDF[i,j,1], BRDF[i,j,2], reflect[i,j]))
            #printlog(False,'BRDF MEAN: [%.4f %.4f %.4f]\treflect: [%.4f]' %(numpy.mean(BRDF[i,:,0]), numpy.mean(BRDF[i,:,1]), numpy.mean(BRDF[i,:,2]), numpy.mean(reflect[i,:])))
        if(neg_x != 0):
            if(numpy.var(SZA) > numpy.var(RAA)):
                W_SZA = W_SZA - 0.2
                W_RAA = W_RAA + 0.2
                
            else:
                W_SZA = W_SZA + 0.2
                W_RAA = W_RAA - 0.2

            if(W_SZA < 0): W_SZA = 0
            if(W_RAA < 0): W_RAA = 0

            for i in range(0, row, 1):
                W_inic[i] = 1.6 - W_SZA*math.sin(abs(SZA[i])) - W_RAA*math.sin(abs(RAA[i]))
            '''
            printlog(False, 'W_SZA: %.2f\t\t W_RAA: %.2f neg_x: %d' %(W_SZA, W_RAA, neg_x)) 
            printlog(False, 'W_inic:')
            printlog(False, W_inic)
            '''
            neg_x = 0
        else:
            break
    # Calc NDVI to return result
    '''
    # Calc NDVI using 'num1' values with 10 percent variation of NDVI 
    NDVI_calc = numpy.zeros((row))
    NDVI_mean  = 0
    NDVI_mean2 = 0
    for i in range (0, row, 1):
        NDVI_calc[i] = (reflect[1,i] - reflect[0,i])/(reflect[1,i] + reflect[0,i])

    printlog(False, 'NDVI_calc:')
    printlog(False, NDVI_calc)
    printlog(False, 'NDVI_mean : %.4f' %(NDVI_mean))
    '''
    '''
    num1 = 0
    for i in range (0, row, 1):
        if NDVI_calc[i] < 1.1*NDVI_mean and NDVI_calc[i] > 0.9*NDVI_mean:
            NDVI_mean2 = NDVI_mean2 + NDVI_calc[i]
            num1 = num1 + 1
    NDVI_mean2 = NDVI_mean2/num1
    printlog(False, 'NDVI_mean2: %.4f (samples: %d)' %(NDVI_mean2, num1))
    '''
    # mat = [VZA VAZ SZA SAZ RED NIR WIR NDVI]
    reflect_mean = numpy.zeros(3)
    num = numpy.zeros(3)
    # vet = [VZA VAZ SZA SAZ RED NIR WIR NDVI]
    # Matrix to return data to 6S algorithm [REF VZA VAZ SZA SAZ] in 3 bands data
    # printlog(False, '----- MAT_BRDF:')
    MAT_BRDF = numpy.zeros((3,5))
    for i in range (0, 3, 1):
        reflect_mean[i] = numpy.mean(reflect[i,:])
        for j in range (0, row, 1):
            if(reflect[i,j] < 1.15*reflect_mean[i] and reflect[i,j] > 0.85*reflect_mean[i]):
                MAT_BRDF[i,0] = MAT_BRDF[i,0] + reflect[i,j]
                MAT_BRDF[i,1] = MAT_BRDF[i,1] + vet[i][0] # VZA
                MAT_BRDF[i,2] = MAT_BRDF[i,2] + vet[i][1] # VAZ
                MAT_BRDF[i,3] = MAT_BRDF[i,3] + vet[i][2] # SZA
                MAT_BRDF[i,4] = MAT_BRDF[i,4] + vet[i][3] # SAZ
                num[i] = num[i] + 1

        # Try against case data problem
        if(num[i] == 0):
            '''
            printlog(True, '. Problem (%d, %d):'%(i,j))
            printlog(True, 'reflect[i]:')
            printlog(True,  reflect[i])
            printlog(True, 'reflect_mean[i]:')
            printlog(True,  reflect_mean[i])
            printlog(True, 'vet[i]:')
            printlog(True,  vet[i])
            '''
            MAT_BRDF[i,:] = [0, 0, -1, -1, -1]
        else:
            MAT_BRDF[i,:] = MAT_BRDF[i,:]/num[i]
    
    RED  = MAT_BRDF[0,0]
    NIR  = MAT_BRDF[1,0]
    WIR  = MAT_BRDF[2,0]

    if(NIR + RED < 0.0001 and NIR + RED > -0.0001):
        NDVI = 9999 
        printlog(True, '...RuntimeWarning. RED: %.2f NIR: %.2f '%(RED, NIR))
    else:
        NDVI = (NIR - RED)/(NIR + RED)

    '''
    printlog(False, 'reflect_mean:')
    printlog(False, reflect_mean)
    printlog(False, 'num:')
    printlog(False, num)
    printlog(False, 'NDVI_mean_return: %.4f' %(NDVI_return))
    printlog(False, 'MAT_BRDF/num: ')
    printlog(False, MAT_BRDF)
    '''    
    
    '''
    printlog(True, 'COEFS:')
    printlog(True, COEFS)
    # wait to test routine
    wait()
    # return NDVI median of all values
    '''
    return NDVI, RED, NIR, WIR, COEFS, MAT_BRDF

def createMaskGoias(filename, band):
    # command to create maskGoias file: filename and band to read
    #createMaskGoias(inputDir + '123\\' + str(day) + str(hour) + '.tif', 4)
    if os.path.isfile(filename) == 0:
        printlog(True,'..Erro! Arquivo Invalido: %s' %strFile)
        return
    infoFile(filename)
    vetMask = readFileBand(filename, 4)
    vetMask = vetMask/255
    numpy.savetxt(os.getcwd() + '\\shape\\maskGoias.txt', vetMask, fmt='%1.0d')
    
def geraMascara(inputFile, outputFile):
    printlog(True,'...geraMascara [x][y]')
    ds1 = gdal.Open(inputFile, GA_ReadOnly )
    printlog(True,'inputFile: ' + inputFile)
    band1 = ds1.GetRasterBand(1)
    array = BandReadAsArray(band1)
    rows = len(array)
    cols = len(array[0])
    #printlog(False,'rows: ' + str(rows))
    #printlog(False,'cols: ' + str(cols))
    MASK1 = numpy.zeros((rows, cols))
    MASK2 = numpy.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            MASK1[i,j] = i
            MASK2[i,j] = j
    # Armazena em arquivo o resultado de selecao dos melhores pixels do dia
    driver = gdal.GetDriverByName("GTiff")
    dsOutMASK = driver.Create(outputFile, ds1.RasterXSize, ds1.RasterYSize, 2, band1.DataType)
    CopyDatasetInfo(ds1,dsOutMASK)
    bandOut1=dsOutMASK.GetRasterBand(1)
    bandOut2=dsOutMASK.GetRasterBand(2)
    BandWriteArray(bandOut1, MASK1)
    BandWriteArray(bandOut2, MASK2)
    #printlog(False,'outputFile: ' + outputFile)
      
def saveCSVfile(matrix, filename):
    timeFile = strftime("%Y%m%d%H%M%S", localtime()) 
    saveFile = os.getcwd()+ '\\logs\\' + filename + '_' + timeFile + '.csv'
    numpy.savetxt(saveFile, matrix, fmt='%.4f',delimiter=";")
    f = open(saveFile,'r')
    filedata = f.read()
    newdata = filedata.replace(".",",")
    f.close()
    f = open(saveFile,'w')
    f.write(newdata)
    f.close()
    
def generateBRDF(inputDir, outputDir, day):
    hour = startHour
    row = 204
    col = 211
    # data = [HOUR][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][X, Y]
    data = numpy.zeros((20, 8, 204, 211), dtype=numpy.float)
    DATA_BRDF = numpy.zeros((4, row, col), dtype=numpy.float)
    global DATA_6S 
    DATA_6S = numpy.zeros((4, row, col), dtype=numpy.float)
    COEFS = numpy.zeros((row, col, 3, 3), dtype=numpy)
    global workQueue
    global exitFlag

    lst = []
    time_brdf = 0
    cont_time_brdf = 0
    time_6S = 0
    cont_time_6S = 0
    start_time = 0
    #printlog(True,data.shape)
    # Read maskGoias
    vetMask = readMaskGoias()
    #printlog(True, 'generate BRDF: ' + str(day))
    #printlog(True,'vetMask: [%d] [%d] ' %(len(vetMask), len(vetMask[1])))
    h = 0 # mark number file in time 15 min
    errors = 0 # Get number errors to read data files
    printlog(True,'----- Reading data files')
    while(hour <= endHour):
        # Flag to read error file
        fileError = False
        printlog(False,'Reading hour: %d' %hour)
        
        # --------- Process ANGLE file
        printlog(False,'Process ANGLE file')
        strFile = inputDir.replace('___', 'ANGLE') + str(day) + str(hour) + '_ANGLES.tif'
        if os.path.isfile(strFile) == 0:
            printlog(False,'..Error! Invalid data file: %s' %strFile)
            fileError = True
        else:
            # get VAZ View Azimuth Angle (msg_azres.tif)
            vet = readFileBand(strFile, 1)
            if(not len(vet)):
                fileError = True
            else:
                data[h][1][:][:] = vet[:][:]
            # get SAZ Sun Azimuth Angle (sol_azres.tif)
            vet = readFileBand(strFile, 2)
            if(not len(vet)):
                fileError = True
            else:
                data[h][3][:][:] = vet[:][:]
            # get VZA View Zenith Angle (msg_zenres.tif)
            vet = readFileBand(strFile, 3)
            if(not len(vet)):
                fileError = True
            else:
                data[h][0][:][:] = vet[:][:]
            # get SZA Sun Zenith Angle (sun_zenres.tif)
            vet = readFileBand(strFile, 4)
            if(not len(vet)):
                fileError = True
            else:
                data[h][2][:][:] = vet[:][:]

        # --------- Process 123 file
        if fileError == False:
            printlog(False,'Process 123 file')
            strFile = inputDir.replace('___', '123') + str(day) + str(hour) + '_123.tif'
            if os.path.isfile(strFile) == 0:
                printlog(False,'..Error! Invalid data file: %s' %strFile)
                fileError = True
            else:
                # get RED Data Band
                vet = readFileBand(strFile, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][4][:][:] = vet[:][:]
                # get NIR Data Band
                vet = readFileBand(strFile, 2)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][5][:][:] = vet[:][:]
                # get WIR Data Band
                vet = readFileBand(strFile, 3)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][6][:][:] = vet[:][:]

        # --------- Process CLM file
        if fileError == False:
            printlog(False,'Process CLM file')
            strFile = inputDir.replace('___', 'CLM') + str(day) + str(hour) + '_CLM.tif'
            if os.path.isfile(strFile) == 0:
                printlog(False,'..Error! Invalid data file: %s' %strFile)
                fileError = True
            else:
                vet = readFileBand(strFile, 1)
                if(not len(vet)):
                    fileError = True
                else:
                    data[h][7][:][:] = vet[:][:] 
        # For next loop    
        h = h + 1
        hour = incHour(hour)
        if(fileError == True):
            errors = errors + 1
    if(errors > 10):
        printlog(False,'..Error! Invalid data file: %s' %strFile)
        return False
    else:
        printlog(True, '----- Now generate BRDF %s'%(day))
    # Read ATM file
    # band1 is mod08:Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean' # Band 57 (scale_factor = 0.001, range = 0 - 5000)
    # band2 is mod08:Aerosol_Optical_Depth_Land_Mean' # Band 37 (scale_factor = 0.001, range = 0 - 5000)
    # band3 is mod08:Total_Ozone_Mean' # Band 829 (scale_factor = 0.1, range = 0 - 5000)
    # band4 is mod08:Atmospheric_Water_Vapor_Mean' # Band 853 (scale_factor = 0.001, range = 0 - 20000)
    printlog(False,'Process ATM file')
    VET_ATM = numpy.zeros((4, row, col), dtype=numpy)
    strFile = inputDir.replace('___', 'ATM') + str(day) + '_ATM.tif'
    if os.path.isfile(strFile) == 0:
        printlog(False,'..Error! Invalid data file: %s' %strFile)
        return False
    else:
        # get Deep_Blue Data Band
        vet = readFileBand(strFile, 1)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[0][:][:] = vet[:][:]
        # get Aerosol Data Band
        vet = readFileBand(strFile, 2)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[1][:][:] = vet[:][:]
        # get Total Ozone Data Band
        vet = readFileBand(strFile, 3)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[2][:][:] = vet[:][:]
        # get Water_Vapour  
        vet = readFileBand(strFile, 3)
        if(not len(vet)):
            fileError = True
        else:
            VET_ATM[3][:][:] = vet[:][:]  
    # data = [HOUR][VZA, VAZ, SZA, SAZ, RED, NIR, WIR, CLM][X, Y]
    # mat = [VZA VAZ SZA SAZ RED NIR WIR]
    percent = 0
    start_time = time.time()
    for i in range(0, row, 1):
        
        if (float(100*i)/float(row) - 5) > percent :
            percent = float(100*i)/float(row)
            printlog(True, "... Process %.2f percent in %.2f minutes"%(percent, (time.time() - start_time)/60) )
            '''
            try:
                printlog(True, 'Time BRDF: %.4f \tTime 6S: %.4f' %(time_brdf/cont_time_brdf, time_6S/cont_time_6S))
            except:
                print('')
            '''
        for j in range(0, col, 1):
            # Verify if pixel is inside Goias area
            if(vetMask[i][j] == 1 and data[0,0,i,j]*data[0,4,i,j] !=0):
                # Verify cloud mask in pixel data
                cont_cloud_free = 0
                cont_images = 0
                #printlog(True, '------- data:')
                for x in range(0, 20, 1):
                    # Verify if pixel had cloud
                    # Cloud Mask: 0 - Clear sky over water, 1 - Clear sky over land, 2 - Cloud, 3 - No data
                    #printlog(True, data[x,0:8,i,j])
                    if(data[x,7,i,j] <= 1 and abs(data[x,1,i,j] - data[x,3,i,j]) < 100 and data[x,0,i,j] != 0):
                        cont_cloud_free = cont_cloud_free + 1
                    else:
                        cont_images = cont_images + 1
                if(cont_cloud_free > 5):
                    #printlog(True,'[%d, %d] Pixel generate BRDF (%d/%d)'%(i, j, cont_cloud_free, cont_cloud_free+cont_images))
                    mat = numpy.zeros((cont_cloud_free, 8))
                    ccf = 0 # count cloud free to loop
                    for x in range(0, 20, 1):
                        # Verify if pixel had cloud
                        if(data[x,7,i,j] <= 1 and abs(data[x,1,i,j] - data[x,3,i,j]) < 100  and data[x,0,i,j] != 0):
                            # To put NDVI in data vector for analyse
                            try:
                                data[x,7,i,j] = (data[x,5,i,j]-data[x,4,i,j])/(data[x,5,i,j]+data[x,4,i,j])
                            except ZeroDivisionError:
                                data[x,7,i,j] = 0
                            mat[ccf] = data[x,0:8,i,j]
                            ccf = ccf + 1
                    '''
                    printlog(False,'-------- mat:')
                    printlog(False, '[VZA, VAZ, SZA, SAZ, RED, NIR, WIR, NDVI]')
                    printlog(False, mat)
                    '''
                    #start_time = time.time()
                    DATA_BRDF[0,i,j], DATA_BRDF[1,i,j], DATA_BRDF[2,i,j],DATA_BRDF[3,i,j],COEFS[i,j,:], MAT_BRDF = calc_BRDF(mat)
                    #print(".time do BRDF: %.4f s" %(time.time() - start_time))

                    if (DATA_BRDF[0,i,j] + DATA_BRDF[1,i,j] + DATA_BRDF[2,i,j] + DATA_BRDF[3,i,j])!= 0:
                        #time_brdf = time_brdf + (time.time() - start_time)
                        #cont_time_brdf = cont_time_brdf + 1
                        #start_time = time.time()
                        #DATA_6S[i,j]= sixsv1(VET, 'angles')
                        #print(".time do 6SV1: %.4f s" %(time.time() - start_time))
                        #time_6S = time_6S + (time.time() - start_time)
                        #cont_time_6S = cont_time_6S + 1
                        if(workQueue.full()):
                            if __debug__: print('f')
                            while(workQueue.full()):
                                pass
                        #Send VET_ATM = [Water_Vapour, Total Ozone, Deep_Blue, Aerosol]
                        ATM = [VET_ATM[3,i,j], VET_ATM[2,i,j], VET_ATM[0,i,j], VET_ATM[1,i,j]]
                        queueLock.acquire()
                        workQueue.put(i)
                        workQueue.put(j)
                        workQueue.put(MAT_BRDF)
                        workQueue.put(ATM)
                        queueLock.release()
                        if __debug__: print('.')
                '''
                else:
                    printlog(True,'[%d, %d] Pixel insuficient observations'%(i,j))
            else:
                if(vetMask[i][j] == 0):
                    printlog(True,'[%d, %d] Pixel out area'%(i,j))
                else:
                    printlog(True,'[%d, %d] Pixel border area'%(i,j))
            '''
    printlog(True, "----- Wait for queue to empty")
    while not workQueue.empty():
        #print('Files in Queue: ' + str(workQueue.qsize()))
        #time.sleep(5)
        pass
    printlog(True, '----- Queue is empty')

    # Notify threads it's time to exit
    exitFlag = 1

    # Wait for all threads to complete
    for t in threads:
        if (t.isAlive()):
            t.join(5.0)
                
    printlog(True, "----- Exiting Main Thread")
    exitFlag = 0

    # Save DATA BRDF in georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly )
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")    
    filename = outputDir + str(day) + '_BRDF.tif'
    dsOut = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 4, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    for i in range(0, 4, 1):
        bandOut = dsOut.GetRasterBand(i+1)
        BandWriteArray(bandOut, DATA_BRDF[i])
   
    # Save DATA 6S in georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly )
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")    
    filename = outputDir + str(day) + '_6S.tif'
    dsOut = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 4, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    for i in range(0, 4, 1):
        bandOut = dsOut.GetRasterBand(i+1)
        BandWriteArray(bandOut, DATA_6S[i])

    # Save COEFs in georeference file
    ds = gdal.Open(fileGoiasGeo, GA_ReadOnly )
    band = ds.GetRasterBand(1)
    driver = gdal.GetDriverByName("GTiff")    
    filename = outputDir + str(day) + '_COEFS.tif'
    dsOut = driver.Create(filename, ds.RasterXSize, ds.RasterYSize, 9, band.DataType)
    CopyDatasetInfo(ds,dsOut)
    for i in range(0, 3, 1):
        for j in range(0, 3, 1):
            bandOut = dsOut.GetRasterBand(3*i+j+1)
            BandWriteArray(bandOut, COEFS[:,:,i,j])
    return True

class main():

    argDict = mapDict(argv, usage)
    lst = []
    global exitFlag
    global threads
    global workQueue

    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict:
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
    else:
        exit(usage)

    date = dateStart
    gc.collect()
    gdal.UseExceptions()
    start_time = time.time()
    exitFlag = 0

    numpy.set_printoptions(formatter={'float': '{: 0.4f}'.format})  
    date = dateStart
    while (date <= dateEnd):
        datestr = str(date)
        printlog(True, '----- Verify BRDF: %s' %(date))
        outputSubdir = outputDir + 'BRDF\\' + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        inputSubdir  = inputDir + '___\\'   + str(date)[0:4] + '\\' + str(date)[4:6] + '\\'
        if(os.path.isfile(outputSubdir + str(date) + '_BRDF.tif')):
            printlog(True, '... File exist in: %s' %(outputSubdir))
            date = incDate(date)
            continue
        else:
            printlog(True, '... Generate file: %s' %(date))
        if not(os.path.isdir(outputSubdir)):
            os.makedirs(outputSubdir)
        printlog(True, "----- Create new threads")
        threadID = 0
        for num_thread in range(0,NUM_THREADS,1):
            thread = myThread(threadID, 'T-' + str(num_thread), workQueue)
            thread.start()
            threads.append(thread)
            threadID += 1
        if(generateBRDF(inputSubdir, outputSubdir, date) == False):
            printlog(True, '----- Error to generate BRDF %s'%(date))
        date = incDate(date)
    printlog(True, '...Finish process')
    # Notify threads it's time to exit
    exitFlag = 1
    