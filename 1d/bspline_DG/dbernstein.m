function f = dbernstein(N,n,x)

f = N*( bernstein(N-1, n-1, x) - bernstein(N-1, n, x) );
