import numpy
import time
import os
import sys
import itertools

'''
	|Variable | Min value | Max value | Step   | Num Step| Mult |
	| VZA     | 55        | 65        | 2.5    | 4       | 10   |
	| VAZ     | 100       | 110       | 2.5    | 4       | 10   |
	| SZA     | 30        | 60        | 5.0    | 6       | 10   |
	| SAZ     | 120       | 230       | 20     | 6       | 10   |
	| WV      | 1.0       | 7.0       | 1.0    | 7       | 10   |
	| OZ      | 0.20      | 0.30      | 0.025  | 4       | 1000 |
	| AOT     | 0.00      | 0.4       | 0.1    | 4       | 1000 |
'''

class corr_atm():

	def __init__(self):
		print('...__init__')
		self.VARS = ['VZA', 'VAZ', 'SZA', 'SAZ', 'WV', 'OZ', 'AOT', 'XA_RED', 'XB_RED', 'XC_RED', 'XA_NIR', 'XB_NIR', 'XC_NIR']
		self.START = [[550, 651, 25, 10], [1000, 1101, 25, 10], [300, 601, 50, 10], [1200, 2301, 200, 10], [10, 71, 10, 10],[200, 301, 25, 1000], [0, 401, 100, 1000]]
		# Order VZA VAZ SZA SAZ WV OZ AOT XA_RED XB_RED XC_RED XA_NIR XB_NIR XC_NIR
		# Order 550 1050 600 2200 30 250 200 516 10170 9516 710 6016 6533

		tac = time.clock()
		for vza in (self.drange(self.START[0][:3])):
			path = os.getcwd() + '\\LUT_6S_VZA_' + str(vza) + '.csv'
			print(vza)
			if(vza == 550):
				DATA = numpy.loadtxt(path, dtype='int32', delimiter = ' ', skiprows = 1)
			else:
				DATA = numpy.append(DATA, numpy.loadtxt(path, dtype='int32', delimiter = ' ', skiprows = 1), axis = 0)
			print('DATA [%d][%d]' %(len(DATA), len(DATA[0])))
		print('time to get data %.2f'%(time.clock()-tac))
		DATA = DATA.T
		for name, i in zip(self.VARS, range(len(self.VARS))):
			globals()[name] = numpy.vstack(DATA[i])
		DATA = []
		print('VZA [%d][%d]' %(len(VZA), len(VZA[0])))
		print('time to org data %.2f'%(time.clock()-tac))

	def adjust_parameters(self, PAR):
		print('...adjust_parameters')
		print('PAR input:')
		print(PAR) 
		# Put input parameters in range of data values of look-up tables
		for i in range(len(PAR)):
			par = PAR[i]*self.START[i][3]
			new_par = 0
			print('i: %d par: %d'%(i, par))
			for lut in self.drange(self.START[i][:3]):
				#print('i: %d par: %.2f lut: %.2f'%(i, par, lut))
				print('lut: %d (par - lut): %d dif: %d' %(lut, numpy.abs(par - lut), self.START[i][2]))
				if(numpy.abs(par - lut) <= self.START[i][2]):
					new_par = lut
					break
			if(i != 6 and new_par == 0):
				print(". Error! Verify data in look-up tables (variable %d): %(i)")
				print(PAR)
				sys.exit()
			PAR[i] = new_par
		print('PAR processed:')
		print(PAR)
		return(PAR)

	def calc_corr(self, PAR):
		print('...calc_corr')
		tac = time.clock()
		VARS_sample = []
		for i in range(len(self.VARS[:7])):
			VARS_sample.append(self.VARS[i] + '_sample')
		print(PAR)
		for name, i in zip(VARS_sample, range(len(VARS_sample))):
			globals()[name] = PAR[i]
		paser = numpy.logical_and(VZA == VZA_sample, numpy.logical_and(VAZ == VAZ_sample, numpy.logical_and(SZA == SZA_sample, \
				 numpy.logical_and(SAZ == SAZ_sample, numpy.logical_and(WV == WV_sample, numpy.logical_and(OZ == OZ_sample, \
				 AOT  == AOT_sample))))))
		samples = len(paser[paser == True])
		print('samples: %d' %(samples))
		if(samples == 0):
			print('Error! not have sample to this data!')
		else:
			print('sample_data')
			print('RED_DATA: %d %d %d' %(XA_RED[paser], XB_RED[paser], XC_RED[paser]))
			print('NIR_DATA: %d %d %d' %(XA_NIR[paser], XB_NIR[paser], XC_NIR[paser]))
		print('time to search %.2f'%(time.clock()-tac))

	def drange(self, data):
		[start, stop, step] = data
		rng = start
		while rng < stop:
			yield rng
			rng += step
		return rng

class main():
	PAR = [55.3, 107.8, 63.4, 223.5, 3.2, 0.18, 0.09]
	six = corr_atm()
	six.calc_corr(six.adjust_parameters(PAR))