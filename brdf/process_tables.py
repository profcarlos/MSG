import sys
sys.path.insert(0, 'C:\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
import numpy
import matplotlib.pyplot as plt
import itertools
import kernels
import os
import gc
import glob
import sys
from utils import * 
import pandas as pd

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
        [43, 177] Natural Savânica   (Terraclass 12) 
""" % argv[0]

# E O HISTOGRAMA DO HORÁRIO DA MELHOR AMOSTRA DO DIA

CULTURE = [[143,48],[178,14],[155,32],[144,38],[146,33],[146,58],[139,57],[144,67],[134,79],[130,63]]
PASTURE = [[35,81],[8,83],[32,118],[34,123],[50,83],[123,81],[130,88],[156,108], [15, 80], [44,188]]
TYPE_FILE = ['BRDF2', '6S', 'SMAC']
use_cols = ['DOY', 'HOUR_1330', 'MAX_DAY', 'MEAN_DAY','MAX_SPARSE_THIN', 'MAX_SPARSE_THICK']
#use_cols = ['DOY', 'MAX_SPARSE_THIN', 'MAX_SPARSE_THICK', 'MAX_DENSE_THIN',  'MAX_DENSE_THICK', 'MAX_ROUJEAN_THIN', 'MAX_ROUJEAN_THICK', 'MAX_GAO', 'MAX_PROUD', 'MAX_SILVEIRA']

def insert_in_cols (type, cols):
    # get end of cols
    end = len(cols)
    if(type == 'BRDF2'): 
        type = 'BRDF'
    # start in third col
    for x in range(end):
        if(cols[x] != 'DOY'):
            cols[x] = str(type) + '_NDVI_' + str(cols[x])
    return cols


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

    for typeFile in TYPE_FILE:
        ALL = []
        for pos_class in range(len(CLASS)):
            MSG_pix_x = CLASS[pos_class][0]
            MSG_pix_y = CLASS[pos_class][1]

            # Create name files
            mod_filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_MOD13.csv'
            msg_filename = 'data_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_typ_'+ str(typeData) + '_NDVI_' + str(typeFile) + '.csv'
            out_filename = 'info_' + str(dateStart) +'_'+ str(dateEnd) + '_pix_' + str(MSG_pix_x) +'_'+ str(MSG_pix_y) + '_typ_'+ str(typeData) + '_NDVI_' + str(typeFile) + '.csv'
            cmp_filename = out_filename.replace('.csv', '_CMP.csv')
            all_filename =  'info_' + str(dateStart) +'_'+ str(dateEnd) + '_class_' + str(classData) + '_typ_'+ str(typeData) + '_NDVI_' + str(typeFile) + '_ALL.csv'
            
            # Read data files
            # MOD13  = [DOY NDVI RED NIR]
            if(os.path.isfile(inputDir + mod_filename)):
                MOD13 = pd.read_csv(inputDir + mod_filename, sep = " ", index_col = False)
            else:
                print("Error! Verify data to read: %s"%(inputDir + mod_filename))
                sys.exit() 

            if(os.path.isfile(inputDir + msg_filename)):
                type_cols = insert_in_cols (typeFile, use_cols)
                MSG_IN = pd.read_csv(inputDir + msg_filename, sep = " ", index_col = False, usecols = type_cols)
            else:
                print("Error! Verify data to read: %s"%(inputDir + msg_filename))
                sys.exit()
        
            # Process data
            print("Process data: %s"%(msg_filename))
            MSG_OUT = MSG_IN[MSG_IN.DOY > 0]
            MSG_OUT = MSG_IN[MSG_IN.DOY < 366]
            MSG_OUT.to_csv(outputDir + out_filename, sep=',', index = False, encoding='utf-8')
            MSG_COMP = MSG_OUT[MSG_OUT.DOY == 0]
            mod_doy_ant = 0
            for mod_doy in MOD13["DOY"]:
                # Verify most near value in msg_doy
                dif = 30
                if(mod_doy_ant > mod_doy):
                    break
                for msg_doy in MSG_OUT["DOY"]:
                    if(mod_doy < msg_doy):
                        break
                    if(mod_doy - msg_doy < dif):
                        dif = mod_doy - msg_doy
                        if (dif == 0):
                            break
                mod_doy_ant = mod_doy
                #print('mod_doy: %d msg_doy: %d'%(mod_doy, msg_doy))
                #print(MSG_SAV[MSG_SAV.DOY == msg_doy])
                MSG_COMP = MSG_COMP.append(MSG_OUT[MSG_OUT.DOY == msg_doy], ignore_index=True)
            #print(MSG_COMP)
            MSG_COMP = MSG_COMP.rename(columns={'DOY': 'DOY_MSG'})
            MOD13 = MOD13.rename(columns={'DOY': 'DOY_MOD'})
            #print('-----------COMP')
            #print(MOD13)
            #print(MSG_COMP)
            MOD13 = pd.concat([MSG_COMP, MOD13], axis = 1)
            MOD13 = MOD13.sort_index(axis=1)
            #print(MOD)
            MOD13.to_csv(outputDir + cmp_filename, sep=',', index = False, encoding='utf-8')
            if(len(ALL) == 0):
                ALL = MOD13
            else:
                ALL = ALL.append(MOD13, ignore_index=True)
            id = ALL.shape[0]
            ALL = ALL[:id-1]
        ALL.to_csv(outputDir + all_filename, sep=',', index = False, encoding='utf-8')
        del ALL

        
 
