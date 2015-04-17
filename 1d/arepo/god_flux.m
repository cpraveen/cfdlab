function f = god_flux(ul, ur, w)

u1 = max(w,ul);
g1 = 0.5*u1*u1 - u1*w;

u2 = min(w,ur);
g2 = 0.5*u2*u2 - u2*w;

f = max(g1,g2);
