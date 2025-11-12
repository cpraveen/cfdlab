import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

filePath = sys.argv[1]
Z = np.loadtxt(filePath, delimiter=' ')

x_len = float(sys.argv[4])
y_len = float(sys.argv[5])

x = np.linspace(0.0, x_len, Z.shape[1])
y = np.linspace(0.0, y_len, Z.shape[0])
X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap=cm.viridis)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel(sys.argv[3])
ax.set_title(sys.argv[2])

ax.view_init(elev=30, azim=230)

plt.savefig(filePath.replace('.dat', '.png'))
