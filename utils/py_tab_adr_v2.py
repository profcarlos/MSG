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



DATA = numpy.genfromtxt('c:\\tmp\\BFAST_corrigido.csv', delimiter = ';', skiprows = 1)
HEADER = ['ID','NBK','PBK_1','PBK_2','PBK_3','PBK_4','PBK_5','PBK_6','PBK_7','PBK_8','PBK_9','PBK_10','PBK_11','PBK_12','PBK_13','PBK_14','PBK_15','PBK_16','PBK_17','PBK_18','PBK_19','PBK_20','PBK_21','PBK_22','PBK_23','PBK_24','PBK_25','PBK_26','PBK_27','PBK_28','PBK_29','PBK_30','DBK_1','DBK_2','DBK_3','DBK_4','DBK_5','DBK_6','DBK_7','DBK_8','DBK_9','DBK_10','DBK_11','DBK_12','DBK_13','DBK_14','DBK_15','DBK_16','DBK_17','DBK_18','DBK_19','DBK_20','DBK_21','DBK_22','DBK_23','DBK_24','DBK_25','DBK_26','DBK_27','DBK_28','DBK_29','DBK_30','MAG_1','MAG_2','MAG_3','MAG_4','MAG_5','MAG_6','MAG_7','MAG_8','MAG_9','MAG_10','MAG_11','MAG_12','MAG_13','MAG_14','MAG_15','MAG_16','MAG_17','MAG_18','MAG_19','MAG_20','MAG_21','MAG_22','MAG_23','MAG_24','MAG_25','MAG_26','MAG_27','MAG_28','MAG_29','MAG_30','INT_1','INT_2','INT_3','INT_4','INT_5','INT_6','INT_7','INT_8','INT_9','INT_10','INT_11','INT_12','INT_13','INT_14','INT_15','INT_16','INT_17','INT_18','INT_19','INT_20','INT_21','INT_22','INT_23','INT_24','INT_25','INT_26','INT_27','INT_28','INT_29','INT_30','SLO_1','SLO_2','SLO_3','SLO_4','SLO_5','SLO_6','SLO_7','SLO_8','SLO_9','SLO_10','SLO_11','SLO_12','SLO_13','SLO_14','SLO_15','SLO_16','SLO_17','SLO_18','SLO_19','SLO_20','SLO_21','SLO_22','SLO_23','SLO_24','SLO_25','SLO_26','SLO_27','SLO_28','SLO_29','SLO_30']

y_len = len(DATA)
x_len = len(DATA[0])
print(x_len)
print(y_len)
OUT = numpy.zeros((y_len, 6))

for y in range(y_len):
	OUT[y,0] = DATA[y,0]
	
	for x in range(HEADER.index('DBK_1'),HEADER.index('DBK_1')+int(DATA[y,1])):
		date = str(DATA[y,x])
		#print(len(date))
		if(len(date) == 9):
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
		if(out == True):
			# Substitui o breakpoint anterior se o atual é positivo e o posterior é negativo ou o atual é negatimo mas maior do que o posterior
			if((OUT[y,1] == 0) or (OUT[y,1] != 0 and ((OUT[y,5] > 0  and DATA[y,x+90] < 0) or (OUT[y,5] < 0  and DATA[y,x+90] < OUT[y,5])))):
				OUT[y,1] = mes
				OUT[y,2] = ano
				OUT[y,3] = DATA[y,x+30]
				OUT[y,4] = DATA[y,x+60]
				OUT[y,5] = DATA[y,x+90]
				print(OUT[y])
	#if(DATA[y,0] > 5):
	#	sys.exit()
numpy.savetxt('c:\\tmp\\' + 'tab_format_bfast_corrigido.csv', OUT, header = 'ID;MONTH;YEAR;MAG;INT;SLOP', delimiter=";", fmt='%.4f')