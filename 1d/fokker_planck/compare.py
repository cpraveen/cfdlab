import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.autolayout'] = True

lw = 2
fs = 18

first = np.loadtxt("sol_first.dat")
tvd = np.loadtxt("sol_tvd.dat")
weno5 = np.loadtxt("sol_weno5.dat")

plt.figure()

plt.plot(first[:,0],first[:,1],'r-.',linewidth=lw)
plt.plot(tvd[:,0],tvd[:,1],'r--',linewidth=lw)
plt.plot(weno5[:,0],weno5[:,1],'r-',linewidth=lw)

plt.plot(first[:,0],first[:,2],'b-.',linewidth=lw)
plt.plot(tvd[:,0],tvd[:,2],'b--',linewidth=lw)
plt.plot(weno5[:,0],weno5[:,2],'b-',linewidth=lw)

plt.xlabel("$x$",fontsize=fs)
plt.ylabel("$f$",fontsize=fs)
plt.legend(("First","TVD","WENO5"),loc="upper right")
plt.title("Fokker-Planck solution at $t=5$, 100 cells")
plt.savefig("fp_compare.pdf")
