import numpy
import sys
import os
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.metrics import classification_report,confusion_matrix

outputDir    = os.getcwd() + '\\process\\'
inumpyutFile = os.getcwd() + '\\data\\TAB_CLM_perc_1_to_10_date_20120101_20121231.csv' 

# DATA = LAT LON BAND_1 BAND_2 BAND_3 BAND_4 BAND_5, BAND_9, CLM
DATA = numpy.loadtxt(inumpyutFile, delimiter = ' ', skiprows = 1, usecols=range(1,10))
#DATA = DATA[1:100]
cols = 6 
rows = len(DATA)
print(cols)
print(rows)
X = numpy.zeros((rows, cols))
X_fit = numpy.zeros((cols, rows))
y = numpy.zeros((rows, 1))

numpy.set_printoptions(precision=4)
#print('DATA')
#for i in range (0, rows, 1):
#    print('[%d][%.4f %.4f %.4f %3d %3d %3d %1d]'% (i, DATA[i][2], DATA[i][3], DATA[i][4], DATA[i][5], DATA[i][6], DATA[i][7], DATA[i][8]))

for i in range (0, rows, 1):
    X[i] = DATA[i][2:8]
    if(DATA[i][8] == 3):
        y[i] = 1
    else:
        y[i] = 0

print('[X] [y]')
for i in range (0, rows, 1):
    print('[%d][%.4f %.4f %.4f %3d %3d %3d] [%1d]'%(i, X[i][0], X[i][1], X[i][2], X[i][3], X[i][4], X[i][5], y[i]))

# Transpose data
X_t = X.transpose()

print('----------------------------')
for i in range (0, cols, 1):
    #print('X_t[%d]'%i)
    #print(X_t[i])
    X_fit[i] = StandardScaler().fit_transform(X_t[i].reshape(-1,1)).transpose()
    #print('X_fit[%d]'%i)
    #print(X_fit[i])

y_fit = y#StandardScaler().fit_transform(y)

#print('X_fit')
#print(X_fit)
#print('y_fit')
#print(y_fit)

X_train, X_test, y_train, y_test = train_test_split(X_fit.transpose(), y_fit, test_size=.4, random_state=int(0.4*rows))

print('X_train')
print(X_train)
print('y_train')
print(y_train)

print('X_test')
print(X_test)
print('y_test')
print(y_test)


print('----------------------------')
print('meam ans std data')
print('X_train: median: %.2f mean: %.2f std: %.2f'%(numpy.median(X_train), numpy.mean(X_train), numpy.std(X_train)))
print('X_test : median: %.2f mean: %.2f std: %.2f'%(numpy.median(X_test), numpy.mean(X_test), numpy.std(X_test)))
print('y_train: median: %.2f mean: %.2f std: %.2f'%(numpy.median(y_train), numpy.mean(y_train), numpy.std(y_train)))
print('y_test : median: %.2f mean: %.2f std: %.2f'%(numpy.median(y_test), numpy.mean(y_test), numpy.std(y_test)))

print('----------------------------')
print('Start MLP')
max_score = 0
results = []
for layer_1 in range(6, 36, 6):
    for layer_2 in range(12, 64, 12):
        for layer_3 in range(3, 15, 3):
            # Use last parameters to process
            warm = False
            for iter in range(200, 1200, 400):
                id_mlp = 'MLP ' + str(layer_1) + ':' + str(layer_2) + ':' + str(layer_3) + ' iter '+ str(iter)
                MLP = MLPClassifier(hidden_layer_sizes=(layer_1, layer_2, layer_3),max_iter = iter, learning_rate_init = 0.01, warm_start = warm, momentum = 0.9)
                warm = True
                MLP.fit(X_train, y_train)
                score = MLP.score(X_test, y_test)
                predictions = MLP.predict(X_test)
                print('ID: %s SCORE: %.4f'%(id_mlp, score))
                #print(confusion_matrix(y_test,predictions))
                #print(classification_report(y_test,predictions))
                results.append([layer_1, layer_2, layer_3, iter, score])
                if(score > max_score):
                    #print('----- save data')
                    coefs = MLP.coefs_
                    intercepts = MLP.intercepts_
                    max_score = score
                    id_max_score = id_mlp

print('best results: %.4f in: %s'%(max_score, id_max_score))
for i in range(len(results)):
    print(results[i])
numpy.save(outputDir + 'MLP_coefs', coefs)
numpy.save(outputDir + 'MLP_intercepts', intercepts)
numpy.savetxt(outputDir + 'MLP_results', results, fmt='%1.4f', comments='', header = "layer_1 layer_2 layer_3 iter score", delimiter = ' ', newline='\n')
