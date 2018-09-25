import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

from sys import argv
from utils import *
import subprocess
import os
import numpy
from scipy.special import factorial

csv_file   = 'C:\\Users\\carlos.silveira\\Dropbox\\Gabriel\\pontos_coord_v2.csv'
input_dir  = 't:\\modis\\'


def convert (x_input, y_input, type):

    '''
    Files: T:\MODIS\MOD15_500m\2000161_MOD15.tif  
    Size is 6096, 2624
    Origin = (-63.850666342091017,-9.999999999104965)
    Pixel Size = (0.003810953930405,-0.003810953930405)

    Files: T:\MODIS\MOD13_1000m\2000049_MOD13.tif
    Size is 3048, 1312
    Origin = (-63.850666342091017,-9.999999999104965)
    Pixel Size = (0.007621907860809,-0.007621907860809)

    Files: T:\MODIS\MOD13_250m\2000049_MOD13.tif
    Size is 12193, 5248
    Origin = (-63.850666342091017,-9.999999999104965)
    Pixel Size = (0.001905476965202,-0.001905476965202)
    '''
    x_orig = -63.850666342091017
    y_orig = -9.999999999104965

    
    if(type == 'MOD15_500m'):
        size = 0.003810953930405
    elif(type == 'MOD13_1000m'):
        size = 0.007621907860809
    elif(type == 'MOD13_250m'):
        size = 0.001905476965202
    else:
        print('...Error in type convert: %s'%(type))
        return (0,0)

    x_sample = int(abs((x_orig - x_input)/size))
    y_sample = int(abs((y_orig - y_input)/size))

    return(y_sample, x_sample)

# MOD15 bands
#['":MOD_Grid_MOD15A2H:Fpar_500m', '":MOD_Grid_MOD15A2H:FparExtra_QC', '":MOD_Grid_MOD15A2H:FparLai_QC', '":MOD_Grid_MOD15A2H:FparStdDev_500m']

# MOD13 bands
#['":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days NDVI"', '":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days composite day of the year"','":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days pixel reliability"',  '":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days VI Quality"']


class main():
  

    COORD = numpy.loadtxt(csv_file, delimiter = ',', skiprows = 1)
   
    DATA_OUT = datafile('', input_dir, '_ndvi_fpar' + '_.csv', "ID LAT LON YEAR DOY NDVI_MOD13_250m PIX_MOD13_250m QUAL_MOD13_250m MODLAND_QA_MOD13_250m VI_USEFUL_MOD13_250m AEROSOL_MOD13_250m BRDF_MOD13_250m CLOUD_MOD13_250m NDVI_MOD13_1000m PIX_MOD13_1000m QUAL_MOD13_1000m MODLAND_QA_MOD13_1000m VI_USEFUL_MOD13_1000m AEROSOL_MOD13_1000m BRDF_MOD13_1000m CLOUD_MOD13_1000m FPAR_500m_1 FPAR_EXTRA_QC_1 AEROSOL_1 FPAR_LAI_QC_1 MODLAND_QC_1 CLOUD_1 FPAR_500m_2 FPAR_EXTRA_QC_2 AEROSOL_1 FPAR_LAI_QC_2 MODLAND_QC_2 CLOUD_2", len(COORD)*(365*(2015-2000)/8), 33)
    for year in range(2000, 2016, 1):
        for doy in range(1, 366, 16):
            # Read MODIS data files
            file_MOD13 = input_dir + '\\MOD13_250m\\' + str(year) + str('%03d'%doy) + '_MOD13.tif'
            if(os.path.isfile(file_MOD13)):
                NDVI_MOD13_250m = readFileBand(file_MOD13, 1)
                PIX_MOD13_250m  = readFileBand(file_MOD13, 3)
                QUAL_MOD13_250m = readFileBand(file_MOD13, 4)
                print('Reading file: %s'%(file_MOD13))
            else:
                print('Error to read file: %s'%(file_MOD13))
                continue

            file_MOD13 = input_dir + '\\MOD13_1000m\\' + str(year) + str('%03d'%doy) + '_MOD13.tif'
            if(os.path.isfile(file_MOD13)):
                NDVI_MOD13_1000m = readFileBand(file_MOD13, 1)
                PIX_MOD13_1000m  = readFileBand(file_MOD13, 3)
                QUAL_MOD13_1000m = readFileBand(file_MOD13, 4)
                print('Reading file: %s'%(file_MOD13))
            else:
                print('Error to read file: %s'%(file_MOD13))
                continue

            file_MOD15 = input_dir + '\\MOD15_500m\\' + str(year) + str('%03d'%doy) + '_MOD15.tif'
            if(os.path.isfile(file_MOD15)):
                FPAR_500m_1     = readFileBand(file_MOD15, 1)
                FPAR_EXTRA_QC_1 = readFileBand(file_MOD15, 2)
                FPAR_LAI_QC_1   = readFileBand(file_MOD15, 3)
                print('Reading file: %s'%(file_MOD15))
            else:
                print('Error to read file: %s'%(file_MOD15))
                #continue
            if(doy != 1):
                file_MOD15 = input_dir + '\\MOD15_500m\\' + str(year) + str('%03d'%(doy-8)) + '_MOD15.tif'
            else:
                file_MOD15 = input_dir + '\\MOD15_500m\\' + str(year-1) + str('%03d'%(361)) + '_MOD15.tif'
            if(os.path.isfile(file_MOD15)):
                FPAR_500m_2     = readFileBand(file_MOD15, 1)
                FPAR_EXTRA_QC_2 = readFileBand(file_MOD15, 2)
                FPAR_LAI_QC_2   = readFileBand(file_MOD15, 3)
                print('Reading file: %s'%(file_MOD15))
            else:
                print('Error to read file: %s'%(file_MOD15))
                #continue

            print('Save data file: %s'%(str(year) + str('%03d'%doy)))
            # Get data for all coord   
            for n in range(len(COORD)):
                data_out = [COORD[n][0], COORD[n][1], COORD[n][2], year, doy]
                # NDVI_MOD13_250m PIX_MOD13_250m QUAL_MOD13_250m MODLAND_QA_MOD13_250m VI_USEFUL_MOD13_250m AEROSOL_MOD13_250m BRDF_MOD13_250m CLOUD_MOD13_250m
                (y, x) = convert(COORD[n][1], COORD[n][2], 'MOD13_250m')
                #FILL -3000
                if(NDVI_MOD13_250m[y][x]>0):
                    data_out.append(NDVI_MOD13_250m[y][x]*0.0001)
                else:
                    data_out.append(0)
                data_out.append(PIX_MOD13_250m[y][x])
                data_out.append(QUAL_MOD13_250m[y][x])
                data_out.append(QUAL_MOD13_250m[y][x]  &(2**0 + 2**1))
                data_out.append((QUAL_MOD13_250m[y][x] &(2**2 + 2**3 + 2**4 + 2**5))>>2) 
                data_out.append((QUAL_MOD13_250m[y][x] &(2**6 + 2**7))>>6)  
                data_out.append((QUAL_MOD13_250m[y][x] &(2**9))>>9)  
                data_out.append(int((QUAL_MOD13_250m[y][x] &(2**8))>>8) + int((QUAL_MOD13_250m[y][x] &(2**10))>>9) + int((QUAL_MOD13_250m[y][x] &(2**15))>>13) )

                (y, x) = convert(COORD[n][1], COORD[n][2], 'MOD13_1000m')
                # FILL -3000
                if(NDVI_MOD13_1000m[y][x]>0):
                    data_out.append(NDVI_MOD13_1000m[y][x]*0.0001)
                else:
                    data_out.append(0)
                data_out.append(PIX_MOD13_1000m[y][x])
                data_out.append(QUAL_MOD13_1000m[y][x])
                data_out.append(QUAL_MOD13_1000m[y][x] &(2**0 + 2**1))
                data_out.append((QUAL_MOD13_1000m[y][x]&(2**2 + 2**3 + 2**4 + 2**5))>>2)
                data_out.append((QUAL_MOD13_1000m[y][x] &(2**6 + 2**7))>>6)
                data_out.append((QUAL_MOD13_1000m[y][x] &(2**9))>>9) 
                data_out.append(int((QUAL_MOD13_1000m[y][x] &(2**8))>>8) + int((QUAL_MOD13_1000m[y][x] &(2**10))>>9) + int((QUAL_MOD13_1000m[y][x] &(2**15))>>13)) 

                (y, x) = convert(COORD[n][1], COORD[n][2], 'MOD15_500m')

                try:
                    # FILL 249-255
                    if(FPAR_500m_1[y][x] < 249):
                        data_out.append(FPAR_500m_1[y][x]*0.01)
                    else:
                        data_out.append(0)
                    data_out.append(FPAR_EXTRA_QC_1[y][x])
                    data_out.append((FPAR_EXTRA_QC_1[y][x]&(2**3 + 2**4))>>3) # Aerosol e Cirrus
                    data_out.append(FPAR_LAI_QC_1[y][x])
                    data_out.append(FPAR_LAI_QC_1[y][x]&(2**0)) # MODLAND_QC
                    data_out.append((FPAR_LAI_QC_1[y][x]&(2**3 + 2**4))>>3) # CloudState
                except:
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                try:
                    # FILL 249-255
                    if(FPAR_500m_2[y][x] < 249):
                        data_out.append(FPAR_500m_2[y][x]*0.01)
                    else:
                        data_out.append(0)                        
                    data_out.append(FPAR_EXTRA_QC_2[y][x])
                    data_out.append((FPAR_EXTRA_QC_2[y][x]&(2**3 + 2**4))>>3) # Aerosol e Cirrus
                    data_out.append(FPAR_LAI_QC_2[y][x])
                    data_out.append(FPAR_LAI_QC_2[y][x]&(2**0)) # MODLAND_QC
                    data_out.append((FPAR_LAI_QC_2[y][x]&(2**3 + 2**4))>>3) # CloudState
                except:
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                    data_out.append(0)
                #print(data_out)
                DATA_OUT.addAtrib(data_out)
                DATA_OUT.addSample()
            # Save data file
            DATA_OUT.save()

    	







