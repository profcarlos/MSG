import numpy


class main:
	folder = 'C:\\Users\\carlos\\Dropbox\\ANALYSES\\OTTO_1010\\'
	tab_data = folder + 'data_trmm_msg_ndvi_ndwi_v2.csv'
	NAME = ['CLASS', 'MONTH', 'TRMM', 'NDVI', 'NDWI', 'NDVI_1', 'NDWI_1', 'NDVI_2', 'NDWI_2', 'NDVI_3', 'NDWI_3']
	DATA = numpy.loadtxt(tab_data, delimiter = ' ', skiprows = 1)

	print(DATA[0])
	DATA_T = numpy.transpose(DATA)
	print(DATA_T[0])
	for i, var in zip(range(len(NAME)), NAME):
		#print(i)
		#print(var)
		globals()[var] = DATA_T[i]
	n_class = int(numpy.amax(CLASS))
	OUT_DATA = numpy.zeros((n_class+1,12))

	for basin in numpy.arange(1,n_class+1,1):
		paser = numpy.logical_and(CLASS==basin, 1)
		#print(paser)
		OUT_DATA[basin][0] = basin
		OUT_DATA[basin][3] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI[paser]))
		OUT_DATA[basin][5] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI_1[paser]))
		OUT_DATA[basin][7] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI_2[paser]))
		OUT_DATA[basin][9] = numpy.min(numpy.corrcoef(TRMM[paser], NDVI_3[paser]))

		OUT_DATA[basin][4] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI[paser]))
		OUT_DATA[basin][6] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI_1[paser]))
		OUT_DATA[basin][8] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI_2[paser]))
		OUT_DATA[basin][10] = numpy.min(numpy.corrcoef(TRMM[paser], NDWI_3[paser]))
		OUT_DATA[basin][1] = numpy.argmax([OUT_DATA[basin][3], OUT_DATA[basin][5], OUT_DATA[basin][7], OUT_DATA[basin][9]])
		OUT_DATA[basin][2] = numpy.argmax([OUT_DATA[basin][4], OUT_DATA[basin][6], OUT_DATA[basin][8], OUT_DATA[basin][10]])
		
		print(OUT_DATA[basin])
	numpy.savetxt(folder + 'corrcoef.csv', OUT_DATA, header = "BASIN, NDVI_MAX, NDWI_MAX, NDVI_0, NDWI_0, NDVI_1, NDWI_1, NDVI_2, NDWI_2, NDVI_3, NDWI_3", delimiter=",",fmt='%.4f')