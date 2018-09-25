import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import numpy
import sys
import subprocess
import os
import glob
from utils import *

outputDir = 't:\\OUTROS\\STAT_GAB\\'

fpar_file = 't:\\OUTROS\\FPAR_CALC\\date_fpar.tif'
efic_file = 't:\\OUTROS\\STAT_GAB\\EFIC\\date_efic.tif'
par_file  = 't:\\OUTROS\\STAT_GAB\\PAR\\date_par.tif'
fpar_file = 't:\\OUTROS\\STAT_GAB\\FPAR\\date_fpar.tif'
apar_file = 't:\\OUTROS\\STAT_GAB\\APAR\\date_apar.tif'
gpp_file  = 't:\\OUTROS\\STAT_GAB\\GPP\\date_gpp.tif'
exeDir = 'C:\\osgeo4w\\bin\\'

year = 2014
for doy in range(1,365,16):

	#print('...process to doy: %03d'%doy)
	cmd_array = exeDir + 'gdal_calc.py -A ' + fpar_file + ' -B ' + par_file + ' --outfile=' + apar_file + ' --calc="A*B"'
	cmd_array = cmd_array.replace('date',str(year)+str('%03d'%doy))
	print(cmd_array)
	#subprocess.call(cmd_array)

	cmd_array = exeDir + 'gdal_calc.py -A ' + apar_file + ' -B ' + efic_file + ' --outfile=' + gpp_file + ' --calc="A*B"'
	cmd_array = cmd_array.replace('date',str(year)+str('%03d'%doy))
	print(cmd_array)
	#subprocess.call(cmd_array)
