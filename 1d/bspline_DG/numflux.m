% Local Lax-Friedrich flux
function f = numflux(u, v)

a = max(abs(dflux(u)), abs(dflux(v)));
f = 0.5*( flux(u) + flux(v) ) - 0.5*a*(v - u);
