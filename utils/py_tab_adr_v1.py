import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

from sys import argv                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
import subprocess
import os
import numpy
from scipy.special import factorial
import random
from scipy.signal import savgol_filter


DATA = numpy.genfromtxt('c:\\tmp\\tab_data.csv', dtype='int32', delimiter = ';', skiprows = 1)

y_len = len(DATA)
x_len = len(DATA[0])
print(x_len)
print(y_len)
OUT = numpy.zeros((y_len, x_len))

for y in range(y_len):
	#print(DATA[y])
	for x in range(x_len):
		if(x == 0):
			OUT[y,x] = DATA[y,x]
			continue
		if(DATA[y,x] > 0):
			date = str(DATA[y,x])
			
			if(len(date) == 7):
				date = '0' + date
			#print('date: %s'%date)
			dia = int(date[0:2])
			mes = int(date[2:4])
			ano = int(date[4:8])
			#print(dia)
			#print(mes)
			#print(ano)
			out = False
			if(ano == 2015 and mes >= 6):
				if(mes == 6):
					if(dia >= 10):
						out = True
				else:
					out = True
			if(ano == 2016 and mes <= 6):
				if(mes == 6):
					if(dia < 10):
						out = True
				else:
					out = True
			OUT[y,x] = int(out)
	#print(OUT[y])
	#sys.exit()
numpy.savetxt('c:\\tmp\\' + 'tab_format.csv', OUT, header = "", delimiter=",",fmt='%d')