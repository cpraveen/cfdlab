% Solve heat equation using FEM and BDF
% Try with 
%      alpha = 0             | solution is stable
%      alpha = 0.4 + pi^2*mu | solution is unstable, one unstable eigenvalue
clear all

a     = 0; 
b     = 1; 
ni    = 100; 
mu    = 1/60; 
alpha = 0.4 + pi^2*mu;
%alpha = 0

n = ni - 1;
h = (b-a)/ni;

% Generate the system matrices
[M,A,B,C,D,N] = matrix_fem(ni,mu,alpha);

% eigenvalues
eo=eig(full(A),full(M));

figure(1)
plot(real(eo),imag(eo),'o')
title('Eigenvalues')

% Function to compute energy
compute_energy = @(z) z'*M*z;

nT=800; dt=0.1; t=0:dt:nT*dt;

z = zeros(n, nT+1);
energy = zeros(nT+1,1);

% Space mesh: only interior points
x = linspace(a+h,b-h,n);

% Set initial condition
z0(1:n) = (x.^2).*(1-x).^3;
i=1;
z(1:n,i) = z0; % Initial flux is zero
energy(i) = compute_energy(z(1:n,i));
figure(2)
plot(x,z(1:n,i),'o-')
xlabel('x')
ylabel('z')
title('Initial condition')

% First time step: use BDF1 (Backward Euler)
beta = 1; a1   = -1;
Am = (1/(beta*dt))*M - A;
[L1,U1,P1,Q1] = lu(Am);

i = 2;
rhs = -(a1/(beta*dt))*M*z(:,i-1);
z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
energy(i) = compute_energy(z(1:n,i));

% Remaining time steps: use BDF2
beta = 2/3; a1 = -4/3; a2 = 1/3;
Am = (1/(beta*dt))*M - A;
[L1,U1,P1,Q1] = lu(Am);

for i=3:nT+1
   rhs = -(1/(beta*dt))*M*(a1*z(:,i-1) + a2*z(:,i-2));
   z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
   energy(i) = compute_energy(z(1:n,i));
   if mod(i,10)==0
      zz = [0, z(1:n,i)', 0];
      xx = [a, x, b];
      figure(3)
      subplot(1,2,1), plot(xx,zz,'-','LineWidth',2)
      xlabel('x'); ylabel('z')
      title('Solution')
      subplot(1,2,2), semilogy(1:i,energy(1:i),'-','LineWidth',2)
      xlabel('Time step');
      title('Energy')
   end
end
