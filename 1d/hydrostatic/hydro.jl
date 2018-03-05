using ApproxFun
using Plots

gam = 1.4;

x = Fun(identity,0..1);
d = domain(x);

# potential
phi = x;

# Hydrostatic temperature Temperature
T = 1 - (gam-1)/gam*phi;

D = Derivative(d);
L = D - (1/T)*(D*phi);
