#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=====================
Classifier comparison
=====================

A comparison of a several classifiers in scikit-learn on synthetic datasets.
The point of this example is to illustrate the nature of decision boundaries
of different classifiers.
This should be taken with a grain of salt, as the intuition conveyed by
these examples does not necessarily carry over to real datasets.

Particularly in high-dimensional spaces, data can more easily be separated
linearly and the simplicity of classifiers such as naive Bayes and linear SVMs
might lead to better generalization than is achieved by other classifiers.

The plots show training points in solid colors and testing points
semi-transparent. The lower right shows the classification accuracy on the test
set.
"""
print(__doc__)


# Code source: Gaël Varoquaux
#              Andreas Müller
# Modified for documentation by Jaques Grobler
# License: BSD 3 clause
import numpy
import sys
import os
import gc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.metrics import classification_report,confusion_matrix
h = .02  # step size in the mesh

names = ["Nearest Neighbors", "Linear SVM", "RBF SVM", 
         "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
         "Naive Bayes", "QDA", "Gaussian Process"]

classifiers = [
    KNeighborsClassifier(3),
    SVC(kernel="linear", C=0.25),
    SVC(gamma=2, C=0.25),
    DecisionTreeClassifier(max_depth=20),
    RandomForestClassifier(max_depth=20, n_estimators=40, max_features=1),
    MLPClassifier(hidden_layer_sizes=(12,36,3), alpha=1),
    AdaBoostClassifier(),
    GaussianNB(),
    QuadraticDiscriminantAnalysis(),
    GaussianProcessClassifier(1.0 * RBF(1.0), warm_start=False)]

outputDir    = os.getcwd() + '\\process\\'
inumpyutFile = os.getcwd() + '\\data\\TAB_CLM_perc_1_to_100_date_20120101_20121231.csv' 

gc.collect()

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

#print('[X] [y]')
#for i in range (0, rows, 1):
#    print('[%d][%.4f %.4f %.4f %3d %3d %3d] [%1d]'%(i, X[i][0], X[i][1], X[i][2], X[i][3], X[i][4], X[i][5], y[i]))

# Transpose data
X_t = X.transpose()

print('----------------------------')
for i in range (0, cols, 1):
    #print('X_t[%d]'%i)
    #print(X_t[i])
    X_fit[i] = StandardScaler().fit_transform(X_t[i].reshape(-1,1)).transpose()
    #print('X_fit[%d]'%i)
    #print(X_fit[i])

y = y.transpose()[0]

#print('X_fit')
#print(X_fit)
#print('y_fit')
#print(y_fit)

X_train, X_test, y_train, y_test = train_test_split(X_fit.transpose(), y, test_size=.3, random_state=int(0.4*rows))

print('----------------------------')
print('meam ans std data')
print('X_train: median: %.2f mean: %.2f std: %.2f num: %d'%(numpy.median(X_train), numpy.mean(X_train), numpy.std(X_train), len(X_train)))
print('X_test : median: %.2f mean: %.2f std: %.2f num: %d'%(numpy.median(X_test), numpy.mean(X_test), numpy.std(X_test), len(X_test)))
print('y_train: median: %.2f mean: %.2f std: %.2f num: %d'%(numpy.median(y_train), numpy.mean(y_train), numpy.std(y_train), len(y_train)))
print('y_test : median: %.2f mean: %.2f std: %.2f num: %d'%(numpy.median(y_test), numpy.mean(y_test), numpy.std(y_test), len(y_test)))

X = None
y = None
X_t = None
DATA = None

# iterate over classifiers
for name, clf in zip(names, classifiers):
    clf.fit(X_train, y_train)
    score = clf.score(X_test, y_test)
    print('----------------------------')
    print("%s: %.4f"%(name, score))
    predictions = clf.predict(X_test)
    print(confusion_matrix(y_test,predictions))
    print(classification_report(y_test,predictions))
