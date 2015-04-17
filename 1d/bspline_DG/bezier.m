function f = bezier(U, s)

N = length(s);

f = zeros(1,N);
for j=1:N
   f(j) = decas(U, s(j));
end
