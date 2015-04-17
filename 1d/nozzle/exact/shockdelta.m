function d = shockdelta(a1,c1,c2)

global xs gamma s1 r1 p0 as1

xsm = xs - 1.0e-10;
xsp = xs + 1.0e-10;

a = nozarea(xs);

[rl ul pl ml] = solution(xs,s1,r1,xsm);
[rr ur pr mr] = solution(xs,s1,r1,xsp);

EPS = 1.0e-20;
dfdM1 = imag( fM1(ml + i*EPS) )/EPS;

m    = rl*ul;
dMdm = (ml/m)*( (1+0.5*(gamma-1)*ml^2)/(1-ml^2) );

% dM/dx at x=xsm
dadx = imag( nozarea(xs+i*EPS) )/EPS;
fun=@(M)((1/M^2)*(2*(1+0.5*(gamma-1)*M^2)/(gamma+1))^((gamma+1)/(gamma-1)));
dfundM = imag( fun(ml+i*EPS) )/EPS;
dMdx = 2*(a/as1)*(dadx/as1)/dfundM;

d = c2 - c1*fM1(ml) - p0*dfdM1*a1*dMdm/a;
d = d/(p0*dfdM1*dMdx);
