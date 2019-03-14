"""
Second order ADER for u_t + a u_x = 0
First predictor : Taylor expansion and use PDE
Second predictor: Correct time slope by upwind flux
"""
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

uinit = lambda x: np.sin(2*np.pi*x)

nc = 100

xmin, xmax = 0.0, 1.0
dx = (xmax - xmin)/nc

x = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nc)

u0 = np.zeros((nc,2))
u1 = np.zeros((nc,2))
B  = np.zeros(nc)

# Set initial condition
for i in range(nc):
    val = integrate.quad(uinit, x[i]-0.5*dx, x[i]+0.5*dx)
    u0[i,0] = val[0]/dx
    foo = lambda xx: ((xx-x[i])/dx)**2
    val0 = integrate.quad(foo, x[i]-0.5*dx, x[i]+0.5*dx)
    foo = lambda xx: uinit(xx)*(xx-x[i])/dx
    val1 = integrate.quad(foo, x[i]-0.5*dx, x[i]+0.5*dx)
    u0[i,1] = val1[0]/val0[0]

plt.figure(figsize=(10,5))
plt.subplot(121)
plt.plot(x,u0[:,0])
plt.subplot(122)
plt.plot(x,u0[:,1])

Tf  = 1.0
cfl = 0.33
a   = 1.0
dt  = cfl*dx/a

c0 = 1.0/2.0 - 2.0*cfl/3.0

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, u0[:,0], 'r')
line2, = ax.plot(x, u0[:,0], 'b')
ax.set_xlabel('x'); ax.set_ylabel('u')
plt.draw()
    
t = 0.0
while t < Tf:
    # First predictor
    for i in range(nc):
        B[i] = - cfl*u0[i,1]
    # Second predictor
    for i in range(nc):
        B[i] = -cfl*((u0[i,0] + c0*u0[i,1]) - (u0[i-1,0] + c0*u0[i-1,1]))
    # Update conserved variable
    for i in range(nc):
        fl = u0[i-1,0] + 0.5*u0[i-1,1] + 0.5*B[i-1]
        fr = u0[i,0] + 0.5*u0[i,1] + 0.5*B[i]
        u1[i,0] = u0[i,0] - cfl*(fr - fl)
        u1[i,1] = u0[i,1] + 12*cfl*(u0[i,0]+0.5*B[i]) \
                  - 6*cfl*(fl+fr)
    t += dt
    u0[:,:] = u1
    line2.set_ydata(u1[:,0])
    plt.draw(); plt.pause(0.1);
