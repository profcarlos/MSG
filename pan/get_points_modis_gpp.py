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

    Files: T:\OUTROS\STAT_GAB\GPP_LAPIG\2014001_gpp.tif
    Size is 3858, 3730
    Origin = (-53.252403461635829,-12.393279057154817)
    Pixel Size = (0.001905476965202,-0.001905476965202)

    Files: T:\OUTROS\STAT_GAB\GPP_clm\2014001_gpp_clm.tif
    Size is 3858, 3730
    Origin = (-53.252403461635794,-12.393279067399060)
    Pixel Size = (0.001905476965202,-0.001905476965202)

    Files: T:\MODIS\MOD17\2014001_GPP.tif
    Size is 1929, 1865
    Origin = (-53.252403461635794,-12.393279067399060)
    Pixel Size = (0.003810953930405,-0.003810953930405)

    '''
    x_orig = -53.252403461635794
    y_orig = -12.393279067399060

    if(type == 'gpp_calc'):
        size = 0.001905476965202
    if(type == 'gpp_mod'):
        size = 0.003810953930405

    x_sample = int(abs((x_orig - x_input)/size))
    y_sample = int(abs((y_orig - y_input)/size))

    return(y_sample, x_sample)

# MOD15 bands
#['":MOD_Grid_MOD15A2H:Fpar_500m', '":MOD_Grid_MOD15A2H:FparExtra_QC', '":MOD_Grid_MOD15A2H:FparLai_QC', '":MOD_Grid_MOD15A2H:FparStdDev_500m']

# MOD13 bands
#['":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days NDVI"', '":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days composite day of the year"','":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days pixel reliability"',  '":MODIS_Grid_16DAY_250m_500m_VI:"250m 16 days VI Quality"']


class main():
  

    COORD = numpy.loadtxt(csv_file, delimiter = ',', skiprows = 1)
   
    DATA_OUT = datafile('', input_dir, 'verif_gpp_comp.csv', "ID LAT LON YEAR DOY GPP_CALC GPP_MOD_1 GPP_CLM_1 GPP_MOD_2 GPP_CLM_2", len(COORD)*int(365/8), 10)
    for year in range(2014, 2015, 1):
        for doy in range(1, 366, 16):
            # Read MODIS data files
            # GPP sem nuvens !!! pasta GPP
            file_GPP_CALC = 'T:\\OUTROS\\STAT_GAB\\GPP\\' + str(year) + str('%03d'%doy) + '_gpp.tif'
            if(os.path.isfile(file_GPP_CALC)):
                GPP_CALC = readFileBand(file_GPP_CALC, 1)
                print('Reading file: %s'%(file_GPP_CALC))
            else:
                print('Error to read file: %s'%(file_GPP))
                continue

            file_GPP_MOD17 = 'T:\\MODIS\\MOD17\\' + str(year) + str('%03d'%doy) + '_gpp.tif'
            if(os.path.isfile(file_GPP_MOD17)):
                GPP_MOD_1 = readFileBand(file_GPP_MOD17, 1)
                GPP_CLM_1 = readFileBand(file_GPP_MOD17, 2)
                print('Reading file: %s'%(file_GPP_MOD17))
            else:
                print('Error to read file: %s'%(file_GPP_MOD17))

            file_GPP_MOD17 = 'T:\\MODIS\\MOD17\\' + str(year) + str('%03d'%(doy+8)) + '_gpp.tif'
            if(os.path.isfile(file_GPP_MOD17)):
                GPP_MOD_2 = readFileBand(file_GPP_MOD17, 1)
                GPP_CLM_2 = readFileBand(file_GPP_MOD17, 2)
                print('Reading file: %s'%(file_GPP_MOD17))
            else:
                print('Error to read file: %s'%(file_GPP_MOD17))

            print('Save data file: %s'%(str(year) + str('%03d'%doy)))
            # Get data for all coord   
            for n in range(len(COORD)):
                data_out = [COORD[n][0], COORD[n][1], COORD[n][2], year, doy]
                # NDVI_MOD13_250m PIX_MOD13_250m QUAL_MOD13_250m MODLAND_QA_MOD13_250m VI_USEFUL_MOD13_250m AEROSOL_MOD13_250m BRDF_MOD13_250m CLOUD_MOD13_250m
                (y, x) = convert(COORD[n][1], COORD[n][2],'gpp_calc')
                data_out.append(GPP_CALC[y][x]*86400/1000000)

                (y, x) = convert(COORD[n][1], COORD[n][2],'gpp_mod')
                data_out.append(GPP_MOD_1[y][x]/8*1000*0.0001)
                data_out.append(GPP_CLM_1[y][x])

                data_out.append(GPP_MOD_2[y][x]/8*1000*0.0001)
                data_out.append(GPP_CLM_2[y][x])
                #print(data_out)
                DATA_OUT.addAtrib(data_out)
                DATA_OUT.addSample()
            # Save data file
            DATA_OUT.save()

    	







