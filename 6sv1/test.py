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

	def __init__():

		self.VARS = ['VZA', 'VAZ', 'SZA', 'SAZ', 'WV', 'OZ', 'AOT', 'XA_RED', 'XB_RED', 'XC_RED', 'XA_NIR', 'XB_NIR', 'XC_NIR']
		self.START = [[550, 625, 25, 10], [1000, 1101, 25, 10], [300, 601, 50, 10], [1200, 2301, 200, 10], [10, 71, 10, 10],[200, 301, 25, 1000], [0, 401, 100, 1000]]
		self.DATA = []
		# Order VZA VAZ SZA SAZ WV OZ AOT XA_RED XB_RED XC_RED XA_NIR XB_NIR XC_NIR
		# Order 550 1050 600 2200 30 250 200 516 10170 9516 710 6016 6533

		tac = time.clock()
		for vza in (drange(self.START[0][:3])):
			path = os.getcwd() + '\\LUT_6S_VZA_' + str(PAR[0]) + '.csv'
			DATA.append(numpy.loadtxt(path, dtype='int32', delimiter = ' ', skiprows = 1))
		print('. DATA[%d][%d]'%(len(DATA), len(DATA[0])))
		DATA = DATA.T
		for name, i in zip(VARS, range(len(VARS))):
			globals()[name] = numpy.vstack(DATA[i])
		print('time to get data %.2f'%(time.clock()-tac))

	def drange(data):
		[start, stop, step] = data
		rng = start
		while rng < stop:
			yield rng
			rng += step
		return rng

	def adjust_parameters(PAR):
		print('PAR input:')
		print(PAR) drange(self.START[i][:3]):
		# Put input parameters in range of data values of look-up tables
		for i in range(len(PAR)):
			par = PAR[i]*self.START[i][3]
			new_par = 0
			for lut in drange(self.START[i][:3]):
				#print('i: %d par: %.2f lut: %.2f'%(i, par, lut))
				if(numpy.abs(par - lut) <= self.START[i][2]):
					new_par = lut
			if(new_par == 0):
				print(". Error! Verify data in look-up tables (variable %d): %(i)")
				print(PAR)
				sys.exit()
			PAR[i] = new_par
		print('PAR processed:')
		print(PAR)

	def calc_corr():

		VARS_sample = []
		for i in range(len(self.VARS[:7])):
			VARS_sample.append(self.VARS[i] + '_sample')
		for name, i in zip(VARS_sample, range(len(VARS_sample))):
			globals()[name] = PAR[i]
		print('time to org data %.2f'%(time.clock()-tac))
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

	class main():
		PAR = [55.3, 107.8, 63.4, 223.5, 3.2, 0.27, 0.29]
		six = corr_atm()
		six.adjust_parameters(PAR)
		six.calc_corr()