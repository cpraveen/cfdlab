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
