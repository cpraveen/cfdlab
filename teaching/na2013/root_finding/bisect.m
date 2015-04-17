% [a,b] = initial nontrivial bracket containing a root
a  = 0.0;
b  = 2.0;
epsln = 1.0e-6; % tolerance on root
delta = 1.0e-6; % tolerance on function value
N = 50;         % maximum number of iterations

f = @(x) x^2 - 2;

fa = f(a); fb = f(b);
sa = sign(fa); sb = sign(fb);
assert(sa ~= sb);

for n=1:N
   e = b-a;
   c = a + 0.5*e;
   fc = f(c);
   fprintf('%5d %20.12e %20.12e %20.12e\n',n,c,fc,c-sqrt(2))
   if abs(e) < epsln*abs(c) || abs(fc) < delta
      break
   end
   sc = sign(fc);
   if sc ~= sb
      a = c; fa = fc; sa = sc;
   else
      b = c; fb = fc; sb = sc;
   end
end

if n==N
   fprintf('Reached limit of iterations\n')
end
