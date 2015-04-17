function [ws,epsn,c,us,gs] = stationarysol(x,nu)

epsn = 0.6;
c = atan(1/(1 + epsn)) - pi*(1+epsn)/4;

ws = -nu*pi*(1+epsn)*tan(pi*(1+epsn)*x/4 + c)/2;
us = ws(1);
gs = -nu*pi/2*(1+epsn)*sec(pi*(1+epsn)/4 + c)^2*pi*(1+epsn)/4;
