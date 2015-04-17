% Hautus Criterion
function [nu] = hautus(A,B)

n = length(A(:,1));

[V,D] = eig(A');
D = real(diag(D));

% Determining the unstable eigenvalues

count = 0; % Number of eigenvalues with non-negative real parts satisfying the Hautus Criterion
nu = 0; % Number of eigenvalues with positive real part
nz = 0; % Number of eigenvalues with zero real part
tol = 1e-13; % To check Hautus criterion

for j = 1:n
   if(D(j) >= -tol)
       if(D(j) > tol)
           nu = nu + 1;
       else
           nz = nz + 1;
       end
       if( abs(B'*V(:,j)) > tol)
          count = count + 1;
       end
   end   
end

fprintf('Number of positive eigenvalues = %d \n',nu)
fprintf('Number of zero eigenvalues = %d\n\n',nz)

if count == nu + nz
   disp('system is stabilizable')
else
   disp('system is not stabilizable')
end

