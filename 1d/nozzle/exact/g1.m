% f = (1/h(x))*dp/dm, m = rho*u
function f = g1(x)

global xs s1 r1

[r u p m] = solution(xs,s1,r1,x);

for j=1:length(x)

   a = nozarea(x(j));
   f(j) = -u(j)/a/(1-m(j)^2);

end
