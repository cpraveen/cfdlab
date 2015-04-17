function f = densityeq(r, x, s1)

global gamma frate H0

if r <= 0.0
   disp('Density is negative')
   disp([r x s1])
   pause
   f = 1.0e20;
   return
end

a = nozarea(x);
p = s1*r^gamma;
c = sqrt(gamma*p/r);
u = frate/(r*a);
f = c^2/(gamma-1) + 0.5*u^2 - H0;
