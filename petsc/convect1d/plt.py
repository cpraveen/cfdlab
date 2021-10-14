from matplotlib import rcParams
rcParams['font.size'] = 14
rcParams['font.family'] = 'serif'
rcParams['figure.autolayout'] = True
rcParams['lines.linewidth'] = 2
rcParams['lines.markersize'] = 6
rcParams['lines.markerfacecolor'] = 'none'
rcParams['axes.titlesize'] = 14
rcParams['axes.labelsize'] = 14
rcParams['text.usetex'] = True    # This will use Latex fonts (slow)

import matplotlib.pyplot as plt
import numpy as np

d0 = np.loadtxt("sol0.dat")
d1 = np.loadtxt("sol.dat")

plt.plot(d0[:,0],d0[:,1],'-')
plt.plot(d1[:,0],d1[:,1],'o')
plt.legend(("Initial","Final"))
plt.xlabel("x")
plt.ylabel("u")
plt.show()
