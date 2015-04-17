M      = 20;  % number of grid points
xmin   = 0;   % left  boundary
xmax   = 1;   % right boundary
mu     = 1;   % coefficient
Tf     = 1;   % final time
lambda = 2.0;

h  = (xmax - xmin)/(M-1);
dt = lambda*h^2/mu;
x  = linspace(xmin, xmax, M);

dm1 = -0.5*lambda*ones(M-2,1);
dp1 = -0.5*lambda*ones(M-2,1);
d0  = 1+lambda*ones(M-2,1);
B  = spdiags([dm1, d0, dp1], -1:1, M-2, M-2);

dm1 = 0.5*lambda*ones(M-2,1);
dp1 = 0.5*lambda*ones(M-2,1);
d0  = 1-lambda*ones(M-2,1);
A  = spdiags([dm1, d0, dp1], -1:1, M-2, M-2);

% Set initial condition (note transpose on x)
U = sin(pi*x');

% Perform time integration
t = 0;
while t < Tf
   Uold = U;
   U(2:M-1) = B \ (A*Uold(2:M-1));
   t = t + dt;
   Ue = exp(-mu*pi^2*t)*sin(pi*x);
   plot(x, U, 'o', x, Ue, '-')
   axis([xmin, xmax, 0, 1])
   legend('Numerical', 'Exact')
   pause(0.1)
end
