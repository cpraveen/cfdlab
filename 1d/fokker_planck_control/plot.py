import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.autolayout'] = True


data = np.loadtxt("sol.dat")

lw = 2
fs = 18

plt.figure()
plt.plot(data[:,0],data[:,1],'r-',linewidth=lw)
plt.plot(data[:,0],data[:,2],'b--',linewidth=lw)
plt.xlabel("$x$",fontsize=fs)
plt.ylabel("$f_1$, $f_2$",fontsize=fs)
plt.legend(("$f_1$","$f_2"))
plt.title("Solution at $t=10$ using WENO5, 500 cells")
plt.savefig("sol.pdf")
