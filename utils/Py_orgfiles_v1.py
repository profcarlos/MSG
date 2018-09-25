import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
import glob
import os.path
import shutil
from sys import argv
from utils import mapDict

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to output data
        -fi     folter to input data
        -ty     type input data (grb, tar, hrf, h5)
""" % argv[0]


yearStart  = 2011
yearEnd    = 2016
monthStart = 1
monthEnd   = 13
vet_month = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

class main():

    argDict = mapDict(argv, usage)

    if "-ty" in argDict and "-fo" in argDict and "-fi" in argDict:
        typeData = argDict["-ty"]
        dirData = argDict["-fi"]
        dirOutput = argDict["-fo"]

    else:
        exit(usage)     
    if not(typeData == 'grb' or typeData == 'tar' or typeData == 'hdf' or typeData == 'h5'):
        exit(usage)

    print('.Organize data in folder: ' + dirData)
    for year in range(yearStart, yearEnd):
        print('...........Organize data in year: ' + str(year))
        folder = dirOutput + '\\' + str(year)
        if not(os.path.isdir(folder)):
            os.makedirs(folder)
        for month in range(monthStart, monthEnd):
            print('......Organize data in month: ' + vet_month[month])
            folder = dirOutput + str(year) + '\\'  + vet_month[month]
            print('.Organize files in folder: ' + folder)
            files = glob.glob(dirData + '*-' + str(year) + vet_month[month] + '*.' + typeData)
            rows = len(files)
            if(rows == 0):
                print('. No files to folder: ' + folder)
                continue
            if not(os.path.isdir(folder)):
                os.makedirs(folder)
            id = len(dirData)
            for i in range (0, rows, 1):
                print('.Move file: ' + files[i][id:])      
                shutil.move(files[i], folder + '\\' + files[i][id:]) 
                
