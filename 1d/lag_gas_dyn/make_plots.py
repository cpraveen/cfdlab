import numpy as np
import matplotlib.pyplot as plt

con = np.loadtxt("sol_con.txt")
noncon1 = np.loadtxt("sol_noncon1.txt")
noncon2 = np.loadtxt("sol_noncon2.txt")

plt.figure()
plt.plot(noncon1[:,0], noncon1[:,1], lw=2, label="Non-cons")
plt.plot(noncon2[:,0], noncon2[:,1], lw=2, label="Non-cons, energy cons")
plt.plot(con[:,0], con[:,1], 'k--', lw=2, label="Conservative")
plt.xlabel("x")
plt.ylabel("Density")
plt.legend()
plt.savefig("lgd1.pdf")

plt.figure()
plt.plot(noncon1[:,0], noncon1[:,1], lw=2, label="Non-cons")
plt.plot(noncon2[:,0], noncon2[:,1], lw=2, label="Non-cons, energy cons")
plt.plot(con[:,0], con[:,1], 'k--', lw=2, label="Conservative")
plt.xlim(0.4,0.8)
plt.xlabel("x")
plt.ylabel("Density")
plt.legend()
plt.savefig("lgd2.pdf")

plt.show()
