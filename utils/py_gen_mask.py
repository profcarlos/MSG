import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import gdal
from utils import *

def createMask(inputFile, outputFile):
    
    DATA = readFileBand(inputFile, 1)
    if(DATA == []):
        print("Error reading input file: %s"%(inputFile))
        return False
    else:
        ds = gdal.Open(inputFile)
        y_123 = ds.RasterYSize
        x_123 = ds.RasterXSize
        DATA = numpy.zeros((2, y_123, x_123))
        for x in range(x_123):
            for y in range(y_123):
                DATA[0][y][x] = x
                DATA[1][y][x] = y
        saveRASTERfile(outputFile, inputFile, DATA)
    return True

class main():
    inputFile = 't:\eumetsat\\MOD13_250M\\2013\\01\\20130101_MOD13.tif'
    outputFile= 't:\\eumetsat\\analyses\\maskVectorGoias_250m.tif'
    createMask(inputFile, outputFile)
    
