function f = derivative(U, s)
globals;

N = length(s);
p = length(U) - 1;

f = zeros(1,N);
for j=1:N
   f(j) = 0;
   for n=0:p
      f(j) = f(j) + U(n+1) * dbernstein(p,n,s(j));
   end
end

f = f/h;
