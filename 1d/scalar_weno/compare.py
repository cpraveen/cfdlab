import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.autolayout'] = True

lw = 2
fs = 18

init = np.loadtxt("init.dat")
first = np.loadtxt("sol_first.dat")
tvd = np.loadtxt("sol_tvd.dat")
#weno3 = np.loadtxt("sol_weno3.dat")
weno5 = np.loadtxt("sol_weno5.dat")

plt.figure()
plt.plot(first[:,0],first[:,1],'r-',linewidth=lw)
plt.plot(tvd[:,0],tvd[:,1],'b--',linewidth=lw)
#plt.plot(weno3[:,0],weno3[:,1],'c-.',linewidth=lw)
plt.plot(weno5[:,0],weno5[:,1],'k-',linewidth=lw)
plt.plot(init[:,0],init[:,1],'c:',linewidth=lw)
plt.xlabel("$x$",fontsize=fs)
plt.ylabel("$u$",fontsize=fs)
plt.legend(("First","TVD","WENO5","Exact"),loc="upper center")
plt.title("Solution of $u_t+u_x=0$ at $t=5$, 100 cells")
plt.savefig("compare.pdf")
