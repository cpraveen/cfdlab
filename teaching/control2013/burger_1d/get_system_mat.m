function [M,A1,A2,D1,d1,d2,H] = get_system_mat(N)

h = 1/N;

d1 = sparse(N,1);
d1(1) = -1/6;
    
d2 = sparse(N,1);
d2(1) = -1/h;

e = ones(N,1);
M = h*spdiags([e/6 2*e/3 e/6], -1:1, N,N);
M(N,N) = h/3;

A1 = (1/h)*spdiags([-e 2*e -e], -1:1, N, N);
A1(N,N) = (1/h);

A2 = sparse(N,N);
A2(1,1) = -1/6;

% tensor stored as a list
Dij = sparse(N,N);
for k = 1 : N
   Dij = 0*Dij; 
   if k==1 
      at = [k, k+1];
      Dij(at,at) = [0, 1/3; -1/6, 1/6];
   end
   if (1<k && k<N)
       at = [k-1, k, k+1];
       Dij(at,at) = [-1/6, 1/6, 0; -1/3, 0, 1/3; 0, -1/6, 1/6];
   end
   if k==N
       at = [k-1, k];
       Dij(at,at) = [-1/6, 1/6; -1/3, 1/3];
   end
   D1{k} = Dij;
end

H = sparse(1,N);
H(1,N) = 1;