% p = degree of bspline
% N = number of cells = no. of knot intervals
% cflmode = 'zhang' or 'praveen'

function [Ul, ubar, mintheta] = project(U)

globals;

[N,p] = size(U);
p = p - 1;

mintheta = 1.0e20;

for j=1:N
   ubar(j) = sum(U(j,:))/(p+1);
end

% limiter
for j=1:N

   % Theta limiter of Zhang-Shu
   ff = bezier(U(j,:), xx);
   mj = min(ff);
   Mj = max(ff);
   a1 = abs(umax(j) - ubar(j))/(abs(Mj - ubar(j)) + 1.0e-14);
   a2 = abs(umin(j) - ubar(j))/(abs(mj - ubar(j)) + 1.0e-14);
   theta = min( [a1, a2, 1.0] );
   mintheta = min([mintheta, theta]);
   Ul(j,:) = theta*(U(j,:) - ubar(j)) + ubar(j);

   uxc = 0; uxb = 0; uxf = 0;
   if j>1 && j<N
      uxb= (ubar(j) - ubar(j-1))/h;
      uxc= (ubar(j+1) - ubar(j-1))/(2*h);
      uxf= (ubar(j+1) - ubar(j))/h;
   else
      if j==1  && periodic==yes
         uxb= (ubar(1) - ubar(N))/h;
         uxc= (ubar(2) - ubar(N))/(2*h);
         uxf= (ubar(2) - ubar(1))/h;
      elseif j==N && periodic==yes
         uxb= (ubar(N) - ubar(N-1))/h;
         uxc= (ubar(1) - ubar(N-1))/(2*h);
         uxf= (ubar(1) - ubar(N))/h;
      end
   end

   ux = minmod( uxc, 2*uxb, 2*uxf );

   if ux ~= 0.0 % then cell average solution is monotone
      ux1 = derivative(Ul(j,:), xx);
      ux  = [ux ux1];
      tolimit = minmodtest(ux);

      if tolimit % if yes, solution is oscillatory, so project
         U1 = zeros(2,1);
         U1(1) = ubar(j) - 0.5*h*ux(1);
         U1(2) = ubar(j) + 0.5*h*ux(1);
         % elevate U1 to degree p
         for n=0:p
            Ul(j,n+1) = U1(1) + (n/p)*(U1(2) - U1(1));
         end
      end

   end

end
