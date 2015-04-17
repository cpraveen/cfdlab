M      = 51;  % number of grid points
xmin   = 0;   % left  boundary
xmax   = 1;   % right boundary
mu     = 1;   % coefficient
Tf     = 1;   % final time
lambda = 0.6;

h  = (xmax - xmin)/(M-1);
dt = lambda*h^2/mu;
x  = linspace(xmin, xmax, M);

% Set initial condition (note transpose on x)
U = sin(pi*x');

% Perform time integration
t = 0;
while t < Tf
   Uold = U;
   for j=2:M-1
      U(j) = lambda*Uold(j-1) + (1-2*lambda)*Uold(j) + lambda*Uold(j+1);
   end
   t = t + dt;
   Ue = exp(-mu*pi^2*t)*sin(pi*x);
   plot(x, U, 'o', x, Ue, '-')
   axis([xmin, xmax, 0, 1])
   legend('Numerical', 'Exact')
   pause(0.1)
end
