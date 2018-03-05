clear all
close all

% Equation of state: rho = rho(p,T)
rho = @(p,T) p./T;

dom = [0,1];
r = chebfun(@(r) r, dom);

% Gravitational potential
phi = r;

% Adiabatic index
gam = 1.4;

% Hydrostatic temperature
T = 1 - (gam-1)*phi/gam;

% Define the problem: dp/dr + rho*dphi/dr = 0
N = chebop(@(p) diff(p,1) + rho(p,T)*diff(phi,1), dom);
N.lbc = [1.0]; % pressure at r=0

% Solve
p = N\0;

figure(1)
plot(p)
xlabel('r')
ylabel('p')

figure(2)
plot(rho(p,T))
xlabel('r')
ylabel('\rho')

figure(3)
plot(T)
xlabel('r')
ylabel('T')
