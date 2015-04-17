% f = dp/dp0
function f = g2(x)

global xs s1 r1 gamma

[r u p m] = solution(xs,s1,r1,x);

for j=1:length(x)

   a = nozarea(x(j));
   f(j) = (1 + 0.5*(gamma-1)*m(j)^2)^(-gamma/(gamma-1));

end
