% 1d burgers code on moving grids
% written by Praveen. C
np = 100;
xmin = 0;
xmax = 1;
tf = 0.5;
cfl = 1.0;

% amount of velocity correction
% alpha = 0 is pure lagrangian case
alpha = 0.1;

% riemann data
ul = 1;
ur = 0.5;
s = 0.5*(ul+ur);

% mesh generating points
xp = linspace(xmin,xmax,np);
xp = xp';

% no. of faces and their locations
nf = np + 1;
xf = zeros(nf,1);
xf(2:nf-1) = 0.5*(xp(1:np-1) + xp(2:np));
xf(1) = xp(1) - (xf(2) - xp(1));
xf(nf) = xp(np) + (xp(np) - xf(nf-1));
xc = 0.5*(xf(1:nf-1) + xf(2:nf));
h = xf(2:nf) - xf(1:nf-1);
u = zeros(size(h));

% initial condition
for j=1:np
   if xp(j) < 0.5
      u(j) = ul;
   else
      u(j) = ur;
   end
end

t = 0;
while(t < tf)
   dt = cfl * min(h ./ abs(u));
   w = (1-alpha)*u + alpha*(xc - xp)/dt;
   res = residual(u,w);
   hu = h.*u - dt*res;
   % update points
   xp = xp + w*dt;
   % sort points
   [y,is] = sort(xp,'ascend');
   xp = xp(is);
   hu = hu(is);
   % compute xf, xc, h
   xf(2:nf-1) = 0.5*(xp(1:np-1) + xp(2:np));
   xf(1) = xp(1) - (xf(2) - xp(1));
   xf(nf) = xp(np) + (xp(np) - xf(nf-1));
   xc = 0.5*(xf(1:nf-1) + xf(2:nf));
   h = xf(2:nf) - xf(1:nf-1);
   % compute u
   u = hu ./ h;
   t = t + dt;
   [dt, t, min(h)]
   plot(xp,u,'o-',...
       [xp(1), 0.5+s*t, 0.5+s*t, xp(end)],[ul,ul,ur,ur],'r-',...
       'LineWidth',1.5);
   legend('AREPO','Exact')
   axis([0.4+0.75*t, 0.6+0.75*t, 0.4, 1.1])
   pause(1)
end
