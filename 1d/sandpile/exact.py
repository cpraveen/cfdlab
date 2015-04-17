import numpy as np

x = np.linspace(0.0,1.0,101)
u = 0*x
v = 0*x
for i in range(len(x)):
    u[i] = np.min([x[i], 1-x[i]])
    v[i] = np.max([0.5-x[i], x[i]-0.5])

f = open('exact.dat','w')
for i in range(len(x)):
    f.write(str(x[i])+"  "+str(u[i])+"  "+str(v[i])+"\n")
