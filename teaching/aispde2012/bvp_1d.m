ua=0.0;
ub=0.0;
a = 0.0;
b = 2*pi;
N = 10;
h = (b - a)/(N+1);
x = linspace(a, b, N+2);

% Construct matrix
A = diag(-ones(N-1,1),-1) + diag(2*ones(N,1),0) + diag(-ones(N-1,1),1);

% Construct rhs vector
B = h*h*sin(x(2:N+1));  B = B';
B(1) = B(1) + ua;
B(N) = B(N) + ub;

% Solve
U = A \ B;

% Plot solution
xx = linspace(a,b,100);
plot(x(2:N+1),U,'o',xx,sin(xx),'LineWidth',2);
legend('Numerical','Exact')
xlabel('x');
ylabel('u(x)')
axis tight

% Compute error norm
u = sin(x(2:N+1))';
e = max(abs(U-u));
fprintf(1,'Max norm of error = %e\n', e)

% Save figure to pdf file
%print -dpdf ode1.pdf
