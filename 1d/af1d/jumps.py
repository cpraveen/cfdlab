'''
Gradient jump indicator
'''
from pylab import *

def indicator(xmin,xmax,fun,nc):
    nv = nc + 1
    n = 2*nc + 1

    dx = (xmax - xmin) / nc
    x = linspace(xmin, xmax, n)
    xv = linspace(xmin, xmax, nv)

    u = fun(x)

    # loop over faces
    s = zeros(nv)
    for f in range(1,nv-1):
        i = 2*f
        uxb = (u[i-2] - 4 * u[i-1] + 3 * u[i]) / dx
        uxf = (-3 * u[i] + 4 * u[i+1] - u[i+2]) / dx
        uxxb = (u[i-2] - 2 * u[i-1] + u[i]) / (0.5*dx)**2
        uxxf = (u[i] - 2 * u[i+1] + u[i+2]) / (0.5*dx)**2
        s[f] = abs(uxb - uxf) + dx * abs(uxxb - uxxf)

    return dx,x,u,xv,s


xmin, xmax = 0.0, 1.0
fun = lambda x: sin(10*pi*x)
nc = 50

dx,x,u,xv,s = indicator(xmin,xmax,fun,nc)

figure()
plot(x,u)

figure()
plot(xv, s)

ncells = array([50,100,200,400,800])
ind = zeros(len(ncells))
h   = zeros(len(ncells))
for i,nc in enumerate(ncells):
    dx,x,u,xv,s = indicator(xmin,xmax,fun,nc)
    ind[i] = s.max()
    h[i] = dx

figure()
loglog(ncells, ind, "-o")
xlabel("ncells")
ylabel("Max jump indicator")

figure()
loglog(ncells, ind/h**2, "-o")
xlabel("ncells")
ylabel("Max jump indicator")

show()
