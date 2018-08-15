import numpy as np
import matplotlib.pyplot as plt
from params import *

x = np.linspace(0,1,200)
y = np.sin(10*x)
z = np.cos(10*x)

plt.figure(figsize=(8,5)) # (width,height)
plt.plot(x,y,c='red',ls='dashed',marker='o',markevery=10,lw=lw,ms=ms)
plt.plot(x,z,c='blue',ls='dashdot',marker='s',markevery=10,lw=lw,ms=ms)
plt.legend(('f','g'))
plt.xlabel('$\\theta$')
plt.ylabel('Function')
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()
