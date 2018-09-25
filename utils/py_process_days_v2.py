import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 't:\\carlos\\newPython\\utils')

import subprocess
import os
import glob
import datetime
import gc
import subprocess
import shutil
import sys
import copy
from utils import * #printlog, sleepTime, delFile, mapDict, extractTarfile, calcMonth, calcDay, calcDoy, calc,Days, readFileBand, incDate
from sys import argv
import numpy
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import time
from scipy.signal import savgol_filter

usage = """\
Usage: %s [OPTIONS]
		-ty     type data (NDVI, BRDF, 6S, MOD13, ATM)
        -fi     folder input data
        -fo     folder output data
        -ds     date to start process (format: yyyyMMdd)
        -de     date to end process (format: yyyyMMdd)
""" % argv[0]

class main():


    argDict = mapDict(argv, usage)

    gc.collect()
    gdal.UseExceptions()    
    if  "-ds" in argDict and "-de" in argDict and "-fo" in argDict and "-fi" in argDict and "-ty" in argDict:
        dateStart = int(argDict["-ds"])
        dateEnd = int(argDict["-de"])
        outputDir = argDict["-fo"]
        inputDir  = argDict["-fi"]
        typeData = argDict["-ty"]
    else:
        exit(usage)

    pontos = numpy.matrix([[164, 9, 'Unid. Cons.', 'green'], [150, 65, 'Sul Agric.', 'red'], [50, 119, 'Norte Agric', 'blue'], [25, 117, 'Norte Pas', 'magenta'], [177, 58, 'Sul Pas', 'orange']])
    
    num_days = calcDays(dateStart, dateEnd)

    dados = numpy.zeros((len(pontos), num_days))
    dados_zeros  = numpy.zeros((len(pontos), num_days))
    dados_savgol = numpy.zeros((len(pontos), num_days))
    pontosX = numpy.zeros((num_days))

    num_days = -1
    date = dateStart
    while(date <= dateEnd):
        num_days = num_days + 1
        pontosX[num_days] = num_days + 1
        if(typeData == 'NDVI'):
            filename = inputDir + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date) + '1200_' + typeData + '.tif'
        else:
            filename = inputDir + str(date)[0:4] + '\\' + str(date)[4:6] + '\\' + str(date) + '_' + typeData + '.tif'
        if(os.path.isfile(filename) == True):
            print('Reading file: %s' %(filename))
            DATA_FILE = readFileBand(filename, 1)
            for pos in range(len(pontos)):
                dados[pos, num_days] = DATA_FILE[pontos[pos,0],pontos[pos,1]]
                #print('dados[%d, %d] = DATA_FILE[pontos[%d,0](%d),pontos[%d,1](%d)] = %.4f (len DATA_FILE: %d)'%(pos, num_days, pos, int(pontos[pos, 0]), pos, int(pontos[pos, 1]), DATA_FILE[pontos[pos,0],pontos[pos,1]], len(DATA_FILE)))
        #else:
            #print('Error to read file: %s' %(filename))
        date = incDate(date)
    filename = 'data_ds' + str(dateStart) + '_de' + str(dateEnd) + '_' + typeData + '.csv'
    print('--- Save output file: %s'%(filename))
    numpy.savetxt(outputDir + filename, dados, fmt='%1.3f', delimiter = ' ', newline='\n')

    # Complete vector in positio without data
    
    for i in range(len(dados)):
        for j in range(num_days+1):
            # Filter extremity data
            if(dados[i,j] >= 0.9 or dados[i,j] <= -0.9):
                dados[i,j] == 0 
            # Filter nan and inf data
            if(numpy.isnan(dados[i, j]) == True or numpy.isinf(dados[i,j]) == True or dados[i, j] == 0):
                #Search nearst data
                dado_down = 0
                dado_up = 0

                for k_down in range(j):
                    if(numpy.isnan(dados[i,j-k_down]) == False and numpy.isinf(dados[i,j-k_down]) == False and dados[i,j-k_down] != 0 and dados[i,j-k_down] < 1 and dados[i,j-k_down] > -1):
                        dado_down =  dados[i,j-k_down]
                        break
                for k_up in range(num_days - j):
                    if(numpy.isnan(dados[i,j+k_up]) == False and numpy.isinf(dados[i,j+k_up]) == False and dados[i,j+k_up] != 0 and dados[i,j+k_up] < 1 and dados[i,j+k_up] > -1):
                        dado_up =  dados[i,j+k_up]
                        break

                if (dado_up != 0 and dado_down != 0):
                    dados_zeros[i,j] = (dado_up + dado_down)/2
                elif(dado_up + dado_down != 0):
                    dados_zeros[i,j] =  (dado_up + dado_down) 
                else:
                    dados_zeros[i,j] = 0
                    #print('[%d, %d] dado_up: %.2f dado_down: %.2f' %(i, j, dado_up, dado_down))
                    #sys.exit(1)
            else:
                dados_zeros[i,j] = dados[i,j]
            
            if(dados_zeros[i,j] >= 0.9 or dados_zeros[i,j] <= -0.9):
                dados_zeros[i,j] == 0 
    filename = 'data_ds' + str(dateStart) + '_de' + str(dateEnd) + '_' + typeData + '_zeros.csv'
    print('--- Save output file: %s'%(filename))
    numpy.savetxt(outputDir + filename, dados_zeros, fmt='%1.3f', delimiter = ' ', newline='\n')

    # Apply Savtisk Golae filter
    for i in range(len(dados)):
        dados_savgol[i] = savgol_filter(dados_zeros[i], 19, 16, mode ='nearest')
    filename = 'data_ds' + str(dateStart) + '_de' + str(dateEnd) + '_' + typeData + '_savgol.csv'
    print('--- Save output file: %s'%(filename))
    numpy.savetxt(outputDir + filename, dados_savgol, fmt='%1.3f', delimiter = ' ', newline='\n')
               
