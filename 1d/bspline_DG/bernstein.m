function f = bernstein(N,n,x)

if n>N || n<0
   f = 0;
elseif N==0
   f = 1;
else
   f = nchoosek(N,n) * x.^n .* (1-x).^(N-n);
end
