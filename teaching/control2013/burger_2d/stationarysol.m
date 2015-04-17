function [ws,wsx] = stationarysol(x,y)
global a b nu

epsn = 0.6;
c = atan(1/(1 + epsn)) - pi*(1+epsn)/4;

%epsn=0.323745909529;
%c=-pi/8;

fun = @(x) -0.5*nu*pi*(1+epsn)*tan(0.25*pi*(1+epsn)*x/a + c);

g  = sin(pi*y/b);
ws = fun(x) .* g;

wsx= imag(fun(x+1i*eps))/eps .* g;

%us = ws(1);
%gs = -nu*pi/2*(1+epsn)*sec(pi*(1+epsn)/4 + c)^2*pi*(1+epsn)/4;
