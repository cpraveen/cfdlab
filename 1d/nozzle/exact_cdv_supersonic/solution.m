function [r u p m] = solution(xs, r1, s1, x)

global gamma frate p0 as1 xthroat

a  = nozarea(xs);

% left state
ml = machnum(as1,xs,'supsonic');
pl = p0*(1 + 0.5*(gamma-1)*ml^2)^(-gamma/(gamma-1));
rl = (pl/s1)^(1/gamma);
cl = sqrt(gamma*pl/rl);
ul = ml*cl;

   % right state
   pr = pl;
   rr = rl;
   ur = ul;
   cr = cl;
   mr = ml;

% post-shock entropy
s2 = pr/rr^gamma;
p02= pr*(1 + 0.5*(gamma-1)*mr^2)^(gamma/(gamma-1)); % total pressure
ff = (1 + 0.5*(gamma-1)*mr^2)*2/(gamma+1);
ff = (ff)^((gamma+1)/(gamma-1));
ff = ff/mr/mr;
as2= a/sqrt(ff);

for j=1:length(x)
   a = nozarea(x(j));
   if x(j) <= xthroat
      m(j) = machnum(as1, x(j), 'supsonic');
      p(j) = p0*(1 + 0.5*(gamma-1)*m(j)^2)^(-gamma/(gamma-1));
      r(j) = (p(j)/s1)^(1/gamma);
   elseif x(j) > xthroat & x(j) < xs
      m(j) = machnum(as1, x(j), 'supsonic');
      p(j) = p0*(1 + 0.5*(gamma-1)*m(j)^2)^(-gamma/(gamma-1));
      r(j) = (p(j)/s1)^(1/gamma);
   else
      m(j) = machnum(as2, x(j), 'supsonic');
      p(j) = p02*(1 + 0.5*(gamma-1)*m(j)^2)^(-gamma/(gamma-1));
      r(j) = (p(j)/s2)^(1/gamma);
   end
   u(j) = frate/(r(j)*a);
end
