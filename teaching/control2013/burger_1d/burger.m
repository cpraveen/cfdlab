% on-linear Burgers' without control

clear all
close all

nu = 1;
ni = 101; % number of nodes
N = ni -1;
h = 1/N;

% Space and time discretization
x = (0:h:1);
nT=1000; dt=0.01; tspan=0:dt:nT*dt;

% Obtaining stationary solution and other stationary conditions
[ws,epsn,c,us,gs] = stationarysol(x,nu);

% Get system matrices
[M,A1,A2,D1] = get_system_mat(N);
[L,U] = lu(M);

% Set initial condition
delta = 0;
z0(1:N) = delta*sin(pi*x(2:N+1)/2);

% Solving the system
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,z] = ode15s(@rhs_burger,tspan,z0,[],L,U,A1,A2,D1,N,nu,us,ws);

for i = 1:nT
    if mod(i,10)==0
        zz = [0, z(i,:)];
        figure(2)
        plot(x,zz+ws,'-','LineWidth',2)
        xlabel('x') 
        ylabel('z')
        titletext = sprintf('Solution at t=%f',tspan(i));
        title(titletext)
        legend('Numerical Soln')
        hold off
    end
end
