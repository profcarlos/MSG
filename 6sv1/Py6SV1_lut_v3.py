import sys
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')

from Py6S import *

import os
import numpy
import time
import gc
from sys import argv
from utils import * 
import gc

usage = """\
Usage: %s [OPTIONS]
		-ty     type data (VZA [550 650])
""" % argv[0]

'''
	|Variable | Min value | Max value | Step   | Num Step| Mult |
	| VZA     | 55        | 65        | 2.5    | 4       | 10   |
	| VAZ     | 100       | 110       | 2.5    | 4       | 10   | 100 a 112.5 ?
	| SZA     | 0         | 60        | 10     | 6       | 10   |  0 a 700
	| SAZ     | 0         | 360       | 60     | 6       | 10   |  0 a 420
	| WV      | 0         | 7.5       | 1.5    | 5       | 10   |  0 a  90
	| OZ      | 0         | 0.30      | 0.06   | 5       | 1000 |  0 a 360 
	| AOT     | 0         | 0.4       | 0.1    | 4       | 1000 |  0 a 800
'''

START = [[550, 651, 25, 10], [1000, 1101, 25, 10], [0, 601, 100, 10], [0, 361, 60, 10], [0, 76, 15, 10],[0, 301, 60, 1000], [600, 801, 100, 1000]]


def sixsv1(VZA, VAZ, SZA, SAZ, WV, OZ, AOT):
	
	# Config 6SV1
	try:
		s = SixS(os.getcwd() + "\\sixsV1_1_lab.exe")
	except:
		print( ' Except in sixv1 path.')
		s.produce_debug_report()
	#print('%.4f %.4f %.4f %.4f %.4f %.4f %.4f'%(VZA, VAZ, SZA, SAZ, WV, OZ, AOT))
	s.altitudes.set_sensor_satellite_level()
	s.altitudes.set_target_sea_level()
	s.geometry = Geometry.User()
	s.geometry.solar_z = SZA
	s.geometry.solar_a = SAZ
	s.geometry.view_z = VZA
	s.geometry.view_a = VAZ
	s.geometry.day = 1
	s.geometry.month = 1
	s.atmos_profile = AtmosProfile.UserWaterAndOzone(WV, OZ)
	s.aot550 = AOT

	# Calc to RED band
	s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(100) 
	s.wavelength = Wavelength (0.56, 0.71)
	try:
		s.run()
	except:
		print('. Error in Py6S (RED)')
		time.sleep(5)
		gc.collect()
		s.run()
	xa = s.outputs.values['coef_xa']
	xb = s.outputs.values['coef_xb']
	xc = s.outputs.values['coef_xc']
	RED = [xa, xb, xc]

	# Calc to NIR band    
	s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(100)  
	s.wavelength = Wavelength (0.74, 0.88)
	try:
		s.run()
	except:
		print('. Error in Py6S (NIR)')
		time.sleep(5)
		gc.collect()
		s.run()
	xa = s.outputs.values['coef_xa']
	xb = s.outputs.values['coef_xb']
	xc = s.outputs.values['coef_xc']
	NIR = [xa, xb, xc]
	return (RED, NIR)


def drange(data):
	[start, stop, step] = data
	rng = start
	while rng < stop:
		yield rng
		rng += step
	yield rng
	return

class main():

	argDict = mapDict(argv, usage)
	gc.collect()
	
	if ("-ty" in argDict):
		VZA = int(argDict["-ty"])
	else:
		exit(usage)

	#for VZA in drange(START[0][:3]):
	lut_file = os.getcwd() + '\\LUT_6S_2\\LUT_6S_VZA_' + str(VZA) + '.csv'
	num = 0
	if(os.path.isfile(lut_file)):
		print('... Loading data file: %s' %(lut_file))
		try:
			DATA = numpy.loadtxt(lut_file, delimiter = ' ', skiprows = 1)
		except:
			readFile = open(lut_file)
			lines = readFile.readlines()
			readFile.close()
			w = open(lut_file,'w')
			w.writelines([item for item in lines[:-1]])
			w.close
			DATA = numpy.loadtxt(lut_file, dtype = int, delimiter = ' ', skiprows = 1)
	else:
		print('... Creating data file: %s' %(lut_file))
		DATA = numpy.zeros((1, 13))
	for VAZ in drange(START[1][:3]):
		for SZA in drange(START[2][:3]):
			for SAZ in drange(START[3][:3]):
				for WV in drange(START[4][:3]):
					for OZ in drange(START[5][:3]):
						for AOT in drange(START[6][:3]):
							if(len(DATA) > num):
								#print(DATA[num])
								if(DATA[num][0] == VZA and DATA[num][1] == VAZ and DATA[num][2] == SZA and DATA[num][3] == SAZ and DATA[num][4] == WV and DATA[num][5] == OZ and DATA[num][6] == AOT):
									print('. DATA in : %d %d %d %d %d %d %d'%(VZA, VAZ, SZA, SAZ, WV, OZ, AOT))
									num = num + 1
									continue
							
							RED, NIR = sixsv1(VZA/START[0][3], VAZ/START[1][3], SZA/START[2][3], SAZ/START[3][3], WV/START[4][3], OZ/START[5][3], AOT/START[6][3])
							RED = [int(float(RED[0])*100000), int(float(RED[1])*100000), int(float(RED[2])*100000)]
							NIR = [int(float(NIR[0])*100000), int(float(NIR[1])*100000), int(float(NIR[2])*100000)]
							newrow = [VZA, VAZ, SZA, SAZ, WV, OZ, AOT, RED[0], RED[1], RED[2], NIR[0], NIR[1], NIR[2]]
							print('. 6SV1 at : %d %d %d %d %d %d %d %d %d %d %d %d %d'%(VZA, VAZ, SZA, SAZ, WV, OZ, AOT, RED[0], RED[1], RED[2], NIR[0], NIR[1], NIR[2]))	
							if(len(DATA) == num):
								DATA  = numpy.vstack([DATA, newrow])
							else:
								try:
									DATA[num] = newrow
								except:
									print("...Error to save DATA")
									print('num: %d'%(num))
									print(DATA.size)
									input()
							num = num + 1
					if(not(len(DATA) > num)):
						print('. Save data')
						numpy.savetxt(lut_file, DATA[:num], fmt='%d', comments = '', header = 'VZA VAZ SZA SAZ WV OZ AOT XA_RED XB_RED XC_RED XA_NIR XB_NIR XC_NIR', delimiter = ' ', newline='\n')

