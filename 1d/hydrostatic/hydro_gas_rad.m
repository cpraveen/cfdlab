% Pressure            : p = rho*R*T + (a/3)*T^4
% Energy per unit mass: e = R*T/(gamma-1) + a*T^4/rho
% Entropy change      : T*ds = de - (p/rho^2)*drho
clear all
close all

a = 0.1;

% Equation of state: rho = rho(p,T)
rho = @(p,T) (p - (a/3)*T.^4)./T;

dom = [0,1];
r = chebfun(@(r) r, dom);

% Gravitational potential
phi = r;

% Adiabatic index
gam = 1.4;

% Hydrostatic temperature
T = 1 - 0.1*r;

% Define the problem: dp/dr + rho*dphi/dr = 0
N = chebop(@(p) diff(p,1) + rho(p,T)*diff(phi,1), dom);
N.lbc = [2.0]; % pressure at r=0

% Solve
p = N\0;

figure(1)
plot(p)
xlabel('r')
ylabel('p')
print -dpdf 'pre.pdf'

figure(2)
plot(rho(p,T))
xlabel('r')
ylabel('\rho')
print -dpdf 'rho.pdf'

figure(3)
plot(T)
xlabel('r')
ylabel('T')
print -dpdf 'T.pdf'


den  = rho(p,T);
dsdr = (1/(gam-1) + 4*a*T^3/den)/T * diff(T,1) - (p + a*T^4)/(T * den^2) * diff(den,1);
figure(4)
plot(dsdr)
xlabel('r')
ylabel('T ds/dr')
print -dpdf 'dsdr.pdf'

% Save solution on mesh of 2^n points
n = 10;
N = 2^n;
x = linspace(0,1,N);
p1 = p(x);
r1 = den(x);
data = [x', p1', r1'];
save('hydro.txt','-ascii','data')
