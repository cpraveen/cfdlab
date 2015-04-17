%------------------------------------------------------------------------------
% Solve linear convection equation with periodic BC
% function lin_hyp_1d_periodic(N, cfl, scheme)
% N      = number of grid points
% cfl    = cfl number
% scheme = BD, FD, LW
%------------------------------------------------------------------------------
function lin_hyp_1d_periodic(N, cfl, scheme)

xmin = 0;
xmax = 1;
a    = 1;
Tf   = 10;

h  = (xmax - xmin)/(N-1);
dt = cfl * h / abs(a);
nu = a * dt / h;

fprintf(1,'N      = %d\n', N);
fprintf(1,'cfl    = %f\n', cfl);
fprintf(1,'scheme = %s\n', scheme);
fprintf(1,'h      = %f\n', h);
fprintf(1,'dt     = %f\n', dt);
fprintf(1,'nu     = %f\n', nu);

% Make grid
x  = linspace(xmin, xmax, N);

% Initial condition
f = @(x) sin(2*pi*x);

% Set initial condition
u = f(x);

t = 0;
while t < Tf
   if strcmp(scheme,'bd')
      u = update_bd(nu, u);
   elseif strcmp(scheme, 'fd')
      u = update_fd(nu, u);
   elseif strcmp(scheme, 'lw')
      u = update_lw(nu, u);
   else
      fprintf(1,'Unknown scheme = %s\n', scheme);
      return
   end
   t = t + dt;
   fe = f(x - a * t);
   plot(x, u, 'r-', x, fe, 'b--', 'LineWidth', 2)
   legend('Numerical', 'Exact')
   pause(0.1);
end

%------------------------------------------------------------------------------
% Backward difference in space
%------------------------------------------------------------------------------
function u = update_bd(nu, u)

uold = u;
N = length(u);

for j=2:N
   u(j) = (1-nu)*uold(j) + nu*uold(j-1);
end

u(1) = u(N);

%------------------------------------------------------------------------------
% Forward difference in space
%------------------------------------------------------------------------------
function u = update_fd(nu, u)

uold = u;
N = length(u);

for j=1:N-1
   u(j) = (1+nu)*uold(j) - nu*uold(j+1);
end

u(N) = u(1);

%------------------------------------------------------------------------------
% Lax-Wendroff
%------------------------------------------------------------------------------
function u = update_lw(nu, u)

uold = u;
N = length(u);

j = 1;
u(j) = uold(j) - 0.5*nu*(uold(j+1) - uold(N-1)) ...
     + 0.5*nu^2*(uold(N-1) - 2*uold(j) + uold(j+1));

for j=2:N-1
   u(j) = uold(j) - 0.5*nu*(uold(j+1) - uold(j-1)) ...
        + 0.5*nu^2*(uold(j-1) - 2*uold(j) + uold(j+1));
end

u(N) = u(1);
