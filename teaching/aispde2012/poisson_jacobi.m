% Solve 2-D poisson equation using Jacobi method
n=25;
xmin=0; xmax=1; ymin=xmin; ymax=xmax;
h = (xmax-xmin)/(n-1); TOL=1e-6; res=TOL+1; iter=0;
x=linspace(xmin,xmax,n);
y=linspace(ymin,ymax,n);
[X,Y]=meshgrid(x,y);
f=2*(2*pi)^2*sin(2*pi*X).*sin(2*pi*Y);
u=zeros(n,n);
while res>TOL
   res = 0; iter=iter+1;
   uo  = u;
   for i=2:n-1
      for j=2:n-1
         u(i,j) = 0.25*(uo(i-1,j) + uo(i+1,j) + uo(i,j-1) + uo(i,j+1)) + ...
                  0.25*h^2*f(i,j);
         res = res + (u(i,j) - uo(i,j))^2;
       end
   end
   res = sqrt(res / n^2);
   fprintf(1,'Iter=%d,  residue=%e\n',iter,res);
end
figure(1); contourf(X,Y,u,25); title('Numerical solution'); colorbar;
% Exact solution
ue=sin(2*pi*X).*sin(2*pi*Y);
figure(2); contourf(X,Y,u-ue,25); title('Error'); colorbar;
