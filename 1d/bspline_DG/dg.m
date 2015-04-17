% p = degree of bspline
% N = number of cells = no. of knot intervals
% cflmode = 'zhang' or 'praveen'

function [ndof, err] = dg(p,N,cflmode)

globals;

assert(p>=1,'Error: p < 1'); assert(N>=1,'Error: N < 1');

% Default initial values
ndof = 0; err = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control variables
U    = zeros(N,p+1);
Uold = zeros(N,p+1);
res  = zeros(N,p+1);

% burgers with initial shock -- step function
if testcase==burger_step
   fluxfun = burger;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 1.0;
   mmconst = 1e20;

   periodic = no
   for j=1:N/2
      U(j,:) = 1.0;
   end
elseif testcase==burger_sine
   fluxfun = burger;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 2.0;
   mmconst = 50;

   periodic = yes;
elseif testcase==burger_rare
   fluxfun = burger;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 0.6;

   periodic = no;
   U(:,:) = 1.0;
   for j=1:N/2
      U(j,:) = -1.0;
   end
elseif testcase==lincon_step
   fluxfun = lincon;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 0.5;

   periodic = no
   for j=1:N/2
      U(j,:) = 1.0;
   end
   mmconst=1e20;
elseif testcase==lincon_sine
   fluxfun = lincon;
   % domain size
   xmin = 0.0
   xmax = 1.0
   tfinal = 2.0;
   periodic = yes;
   mmconst = 50;
elseif testcase==lincon_hat
   fluxfun = lincon;
   % domain size
   xmin = 0.0
   xmax = 1.0
   tfinal = 2.0;
   periodic = yes;
   mmconst = 1e20;
elseif testcase==lincon_zalesak
   fluxfun = lincon;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 2.0;
   periodic = yes;
   mmconst = 1e20;
else
   fprintf(1,'Unknown test case\n');
   return;
end

% Number of gauss quadrature points
if fluxfun==lincon
   Ng = p; % Exact for linear convection
elseif fluxfun==burger
   Ng = ceil(3*p/2); % Exact for quadratic flux
else
   fprintf(1,'Dont know Ng\n');
   return
end

assert(Ng>=1 && Ng<=10, 'We need 1 <= Ng <=10\n');

% Cell size
h = (xmax - xmin)/N;

% weights and quadrature points for limiting - zhang-shu
if strcmp(cflmode,'zhang')
   if p==1
      xx = [-0.5 +0.5];
      cfl= 1.0/4.0;
   elseif p==2 || p==3
      xx = [-0.5 0.0 +0.5];
      cfl= 1.0/6.0;
   elseif p==4 || p==5
      xx = [-0.5 -1/sqrt(20) +1/sqrt(20) +0.5];
      cfl= 1.0/12.0;
   else
      fprintf(1,'Dont know cfl for this value of p\n');
      return
   end
   xx = xx + 0.5; % shift to [0,1]
else
   xx = linspace(0,p,p+1)/p;
   if p==1
      cfl= 1.0/4.0;
   elseif p==2
      cfl= 1.0/6.0;
   elseif p==3
      cfl= 1.0/8.0;
   elseif p==4
      cfl= 1.0/12.85714285;
   elseif p==5
      cfl= 1.0/15.15789473;
   else
      fprintf(1,'Dont know cfl for this value of p\n');
      return
   end
end

% cell centers
for j=1:N
   xc(j) = xmin + (j-1)*h + 0.5*h;
end

% mass matrix
for m=0:p
   for n=0:p
      M(m+1,n+1) = nchoosek(p,m) * nchoosek(p,n) * ...
                   gamma(2*p + 1 - (m+n)) * gamma(m+n+1) / gamma(2+2*p);
   end
end
invM = inv(M);

% initial condition
if testcase ~= burger_step && ...
   testcase ~= burger_rare && ...
   testcase ~= lincon_step
   %for j=1:N
   %   x1 = xmin + (j-1)*h;
   %   U(j,1) = initcond(x1);
   %   U(j,p+1) = initcond(x1+h);
   %end

   if p>=1
      AA = M(2:p,2:p); invAA = inv(AA);
      %b = zeros(p-1,1);
      b = zeros(p+1,1);
      for j=1:N
         x1 = xmin + (j-1)*h;
         for n=1:p+1
            fun = @(x)( initcond(x1+h*x) .* bernstein(p,n-1,x) );
            %b(n) = quadgk(fun,0,1) - U(j,1)*M(1,n+1) - U(j,p+1)*M(p+1,n+1);
            b(n) = quadgk(fun,0,1);
         end
         %U(j,2:p) = invAA * b;
         U(j,:) = invM * b;
      end
   end
end

umin1 = 1.0e20;
umax1 =-1.0e20;

% plot initial condition
figure(1)
hold on
for j=1:N
   x1 = xmin + (j-1)*h;
   s  = linspace(0,1,10);
   up = bezier(U(j,:), s);
   xp = x1 + h*s;
   plot(xp,up)
   umin1 = min([umin1, min(up)]);
   umax1 = max([umax1, max(up)]);
end
title('Initial condition');

umin1
umax1

umin = umin1 * ones(N,1);
umax = umax1 * ones(N,1);

% Parameters for runge-kutta
ark = [0.0 3.0/4.0 1.0/3.0];
brk = 1 - ark;

amax  = 1; % max speed
dt    = cfl*h/amax;
niter = tfinal/dt;


% Project initial condition
[Ul,ubar,mintheta] = project(U);
if mintheta < 1.0
   fprintf(1,'Initial condition is being projected\n');
   fprintf(1,'Minimum theta = %e\n', mintheta);
end
U = Ul;

plot(xc,ubar,'o')

time = 0.0;
iter = 0;

% time integration
while time < tfinal
%for iter=1:niter
   % Exactly match final time
   if time + dt > tfinal
      dt = tfinal - time;
   end

   Uold = U;
   %minmax(ubar);

   % Runge-kutta stages
   for rks=1:3

      res(:,:) = 0.0;

      % Inter-element fluxes
      if periodic==no
         flx = numflux(U(1,1), U(1,1));
      else
         flx = numflux(U(N,p+1), U(1,1));
      end
      res(1,1) = res(1,1) - flx;

      for j=1:N-1
         flx = numflux(U(j,p+1), U(j+1,1));
         res(j,  p+1) = res(j,  p+1) + flx;
         res(j+1,1  ) = res(j+1,1  ) - flx;
      end

      if periodic==no
         flx = numflux(U(N,p+1), U(N,p+1));
      else
         flx = numflux(U(N,p+1), U(1,1));
      end
      res(N,p+1) = res(N,p+1) + flx;

      % Flux integral
      for j=1:N
         for n=0:p
            fun = @(x)( flux( bezier(U(j,:), x) ).*dbernstein(p,n,x) );
            %res(j,n+1) = res(j,n+1) - quadgk(fun,0,1);
            fg = fun(x_values(Ng,:));
            res(j,n+1) = res(j,n+1) - sum( fg.*c_values(Ng,:) );
         end
      end

      % Update control variables
      for j=1:N
         res(j,:) = invM * (res(j,:))';
         U(j,:) = ark(rks)*Uold(j,:) + ...
                  brk(rks)*(U(j,:) - (dt/h)*res(j,:));
      end

      % limiter
      [Ul,ubar,mintheta] = project(U);
      U = Ul;

   end % end of RK

   % compute L2 norm of solution
   energy = 0.0;
   for j=1:N
      energy = energy + U(j,:) * M * (U(j,:))';
   end
   energy = h * energy;

   time = time + dt;
   iter = iter + 1;
   fprintf(1,'%d %f %f %f %f\n', iter, dt, time, mintheta, energy);

   figure(2)
   xp = []; up = [];
   for j=1:N
      x1 = xmin + (j-1)*h;
      s  = linspace(0,1,10);
      up1= bezier(U(j,:), s);
      xp1= x1 + h*s;
      xp = [xp xp1]; up = [up up1];
      %plot(x1+h*linspace(0,p,p+1)/p,U(j,:),'-*g')
      %hold on
   end
   plot(xp,up,'LineWidth',1.5); hold on
   plot(xc,ubar,'o')
   if testcase==lincon_sine || testcase==lincon_hat || ...
      testcase==lincon_zalesak
      x = linspace(xmin,xmax,200);
      plot(x,initcond(x-time),'--r')
   elseif testcase==lincon_step
      plot([xmin,time,time,xmax],[1,1,0,0],'--r')
   elseif testcase==burger_step
      plot([xmin,0.5*time,0.5*time,xmax],[1,1,0,0],'--r')
   end
   axis tight
   hold off

end

% Plot legend
if testcase==lincon_sine || testcase==lincon_hat || ...
   testcase==lincon_step || testcase==lincon_zalesak || ...
   testcase==burger_step
   legend('DGFEM', 'Cell avg', 'Exact');
end

% Compute error for linear convection equation
err = 0.0;
if testcase==lincon_sine
   for j=1:N
      x1  = xmin + (j-1)*h;
      fun = @(x)( (bezier(U(j,:), x) - initcond(x1+h*x-time)).^2 );
      err1= quadgk(fun,0,1,'AbsTol',1e-10,'RelTol',0);
      err = err + err1;
   end
   err = sqrt(h*err);
else
   fprintf(1,'Error computation not implemented\n');
end

fprintf(1,'L2 error = %e\n', err);

ndof = N*(p+1);

% save cell averages to file
fp = fopen('solavg.dat','w');
fprintf(fp,'%e %e\n', [xc; ubar]);
fclose(fp);
