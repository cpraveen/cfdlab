clear all
close all

global gamma H0 frate L xthroat athroat as1 p0 noztyp
global xs s1 r1

gamma=1.4;

noztyp = 3;

% x location of nozzle inlet and outlet
x1 = 0.0; % inlet
a1 = nozarea(x1); % inlet area
x2 = L; % L is set when nozarea is called above

as1 = athroat;

% Inflow conditions
m1=machnum(as1, 0, 'subsonic'); % mach number
r1=1.0; % density
p1=1.0;
p0=p1*(1 + 0.5*(gamma-1)*m1^2)^(gamma/(gamma-1));

% pressure ratio: Poutlet/Pinlet
prat=0.80;

c1=sqrt(gamma*p1/r1);
u1=m1*c1;
s1=p1/r1^gamma;

% Enthalpy: this is constant
H0 = c1^2/(gamma-1.0) + 0.5*u1^2;

% outflow pressure
p2=p1*prat;

% flow rate: this is constant
frate = r1*u1*a1;

% Find shock location xs
fun = @(x) shockfun(x, x2, r1, s1, p2);
xs = fminbnd(fun, athroat, x2, optimset('TolX',1e-12))
fun(xs)

% integral of pressure
fun = @(x) pressure(xs, r1, s1, x);
integral = quadgk(fun, x1, x2, 'Waypoints', [xs], 'AbsTol', 1e-10, 'RelTol', 0);
fprintf(1,'Pressure Integral = %24.15e\n', integral);

% Load numerical solution
%load flow.dat

% plots
x = linspace(x1,x2,1000);
[r u p m] = solution(xs, s1, r1, x);
for j=1:length(x)
   a(j) = nozarea(x(j));
end

figure(1)
%plot(x,r,flow(:,1),flow(:,2),'o');
plot(x,r)
%legend('Exact', 'Numerical')
title('Density')

figure(2)
%plot(x,u,flow(:,1),flow(:,3),'o');
plot(x,u)
%legend('Exact', 'Numerical')
title('Velocity')

figure(3)
%plot(x,p,flow(:,1),flow(:,4),'o');
plot(x,p)
%legend('Exact', 'Numerical')
title('Pressure')

figure(4)
%plot(x,m,flow(:,1),flow(:,5),'o');
plot(x,m)
%legend('Exact', 'Numerical')
title('Mach number')

figure(5)
plot(x,a)
title('Area variation')
