import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.metrics import classification_report,confusion_matrix

outputDir = 'c:\\tmp\\'
            
X, y = make_classification(n_features=7, n_redundant=0, n_informative=2,
                           random_state=0, n_clusters_per_class=1)
rng = np.random.RandomState(2)
X += 2 * rng.uniform(size=X.shape)
linearly_separable = (X, y)

X, y = linearly_separable
X = StandardScaler().fit_transform(X)
X_train, X_test, y_train, y_test = \
        train_test_split(X, y, test_size=.4, random_state=40)
print('pre fit')
print('X_train: median: %.2f mean: %.2f std: %.2f'%(np.median(X_train), np.mean(X_train), np.std(X_train)))
print('X_test : median: %.2f mean: %.2f std: %.2f'%(np.median(X_test), np.mean(X_test), np.std(X_test)))

scaler = StandardScaler()
# Fit only to the training data
scaler.fit(X_train)
X_train = scaler.transform(X_train)
scaler.fit(X_test)
X_test = scaler.transform(X_test)
print('pos fit')
print('X_train: median: %.2f mean: %.2f std: %.2f'%(np.median(X_train), np.mean(X_train), np.std(X_train)))
print('X_test : median: %.2f mean: %.2f std: %.2f'%(np.median(X_test), np.mean(X_test), np.std(X_test)))

max_score = 0
results = []
for iter in range(100, 3000, 100):
    print('--------------- iter: %d'%(iter))
    MLP = MLPClassifier(hidden_layer_sizes=(28, 14, 7),max_iter = iter, learning_rate_init = 0.01, warm_start = True, momentum = 0.9)

    MLP.fit(X_train, y_train)

    score = MLP.score(X_test, y_test)

    print(score)
    predictions = MLP.predict(X_test)
    #print(confusion_matrix(y_test,predictions))
    print(classification_report(y_test,predictions))
    results.append([iter, score])
    if(score > max_score):
        print('----- save data')
        coefs = MLP.coefs_
        intercepts = MLP.intercepts_
        max_score = score
        id_max_score = iter

print('best results: %.4f iter: %d'%(max_score, id_max_score))
print(results)
np.save(outputDir + 'MLP_coefs', coefs)
np.save(outputDir + 'MLP_intercepts', intercepts)
    


