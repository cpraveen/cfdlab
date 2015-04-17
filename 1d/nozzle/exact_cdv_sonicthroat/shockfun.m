function f = shockfun(xs, x2, r1, s1, p2)

global gamma frate as1 p0

a  = nozarea(xs);

% left state
ml = machnum(as1,xs,'supsonic');
pl = p0*(1 + 0.5*(gamma-1)*ml^2)^(-gamma/(gamma-1));
rl = (pl/s1)^(1/gamma);
cl = sqrt(gamma*pl/rl);
ul = ml*cl;

if ml>1.0
   % right state
   pr = pl*(1 + 2*gamma*(ml^2-1)/(gamma+1));
   rr = rl*(gamma+1)*ml^2/((gamma-1)*ml^2 + 2);
   ur = ul*rl/rr;
   cr = sqrt(gamma*pr/rr);
   mr = ur/cr;
else
   pr = pl;
   rr = rl;
   ur = ul;
   cr = cl;
   mr = ml;
end

% post-shock entropy
p02= pr*(1 + 0.5*(gamma-1)*mr^2)^(gamma/(gamma-1)); % total pressure
ff = (1 + 0.5*(gamma-1)*mr^2)*2/(gamma+1);
ff = (ff)^((gamma+1)/(gamma-1));
ff = ff/mr/mr;
as2= a/sqrt(ff);

% outflow conditions
m2 = machnum(as2,x2,'subsonic');
pc = p02*(1 + 0.5*(gamma-1)*m2^2)^(-gamma/(gamma-1));

% function
f = abs(pc - p2)/p2;
