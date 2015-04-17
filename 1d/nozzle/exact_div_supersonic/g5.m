% f = h(x)*dm/dp0, m = rho*u
function f = g5(x)

global xs s1 r1 gamma H0

[r u p m] = solution(xs,s1,r1,x);

for j=1:length(x)

   a = nozarea(x(j));
   sonic = sqrt(gamma*p(j)/r(j));
   drhodp0 = gamma/(gamma-1)/H0*(1 + 0.5*(gamma-1)*m(j)^2)^(-1/(gamma-1));
   dmdp0 = sonic*m(j)*drhodp0;
   f(j) = a*dmdp0;

end
