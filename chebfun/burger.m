% Find solution using cole-hopf tranformation
clear all

nu = 0.05;

x = chebfun('x',[0,2*pi]);
u0 = sin(x);
plot(u0); hold on

phi0 = exp(0.5*cos(x)/nu - 0.5/nu);
%plot(phi0)

L = chebop(0,2*pi);
L.op = @(x,u) nu*diff(u,2);
L.lbc = @(u) diff(u);
L.rbc = @(u) diff(u);

dt = 1.0;
N  = 5;
for t = linspace(dt,N*dt,N)
   phi = expm(L, t, phi0);
   %plot(phi); hold on
   u = -2*nu*diff(log(phi));
   plot(u); hold on
end

hold off
