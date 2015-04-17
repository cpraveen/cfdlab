function mach = machnum(as, x, flow)

global xthroat gamma

a = nozarea(x);

ar= as/a;
pp= (gamma+1)/(gamma-1);
fun = @(m) ( ar - m*( 2*(1 + 0.5*(gamma-1)*m^2)/(gamma+1) )^(-0.5*pp) );

if flow == 'subsonic'
   mach = fzero(fun, [0,1]);
else
   mach = fzero(fun, [1,10]);
end
