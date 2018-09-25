import numpy as np
import matplotlib.pyplot as plt

N = 10
xTrue = np.linspace(0, 1, N)
yTrue = 3 * xTrue
xData = xTrue + np.random.normal(0, 1, N)
yData = yTrue + np.random.normal(0, 1, N)
xData = np.reshape(xData, (N, 1))
yData = np.reshape(yData, (N, 1))
data = np.hstack((xData, yData))
print('data')
print(data)
mu = data.mean(axis=0)
data = data - mu
# data = (data - mu)/data.std(axis=0)  # Uncommenting this reproduces mlab.PCA results
eigenvectors, eigenvalues, V = np.linalg.svd(data.T, full_matrices=False)
projected_data = np.dot(data, eigenvectors)
sigma = projected_data.std(axis=0).mean()
print('mu')
print(mu)
print('data')
print(data)
print('eigenvectors')
print(eigenvectors)
print('eigenvalues')
print(eigenvalues)
print('V')
print(V)


fig, ax = plt.subplots()
ax.scatter(xData, yData)
for axis in eigenvectors:
    start, end = mu, mu + sigma * axis
    ax.annotate(
        '', xy=end, xycoords='data',
        xytext=start, textcoords='data',
        arrowprops=dict(facecolor='red', width=2.0))
ax.set_aspect('equal')
plt.show()
