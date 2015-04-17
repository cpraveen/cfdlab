function f = dflux(u)

eps = 1.0e-20;
f = imag(flux(u + i*eps))/eps;
