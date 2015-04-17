function f = decas(U, s)

p = length(U) - 1;
s1 = 1 - s;

Ua = U;

for r=1:p
   for i=1:p+1-r
      Ua(i) = s1 * Ua(i) + s * Ua(i+1);
   end
end

f = Ua(1);
