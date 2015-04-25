import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.autolayout'] = True


data = np.loadtxt("sol.dat")
targ = np.loadtxt("target.dat")

lw = 2
fs = 18

plt.figure()
plt.plot(data[:,0],data[:,1],'r-',linewidth=lw)
plt.plot(data[:,0],data[:,2],'b--',linewidth=lw)
plt.plot(targ[:,0],targ[:,1],'k-',linewidth=lw)
plt.plot(targ[:,0],targ[:,2],'k--',linewidth=lw)
plt.xlabel("$x$",fontsize=fs)
plt.ylabel("$f_1$, $f_2$",fontsize=fs)
plt.legend(("$f_1$","$f_2$","$f_1^T$","$f_2^T$"))
plt.title("Solution at $t=5$ using WENO5, 500 cells")
plt.savefig("f_opt.pdf")

#--------------------------------------------------------------------------
data = np.loadtxt("opt.dat")

plt.figure()
plt.plot(data[:,0],data[:,1],'r-',linewidth=lw)
plt.plot(data[:,0],data[:,2],'b--',linewidth=lw)
plt.xlabel("$t$")
plt.ylabel("$u$")
plt.legend(("$u_1$","$u_2$"))
plt.savefig("control.pdf")
