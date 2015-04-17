% Heat equation using FEM and BDF
% Feedback control

clear all
close all

a     = 0; 
b     = 1; 
ni    = 100; 
mu    = 1; 
alpha = 0.4 + pi^2*mu;
n = ni - 1;
h = (b-a)/ni;

% Generate the system matrices
[M,A,B,Q,R] = matrix_fem(ni,mu,alpha);

% uncontrolled eigenvalues
eo=eigs(A,M,10,'LR');

K = feedback_matrix(M,A,B,Q,R);
A=A-B*sparse(K); 

% eigenvalues with feedback control
ec=eigs(A,M,10,'LR');

figure(2)
plot(real(eo),imag(eo),'o',real(ec),imag(ec),'*')
title('Eigenvalues')
legend('Original','With feedback')
grid on

% Function to compute energy
compute_energy = @(z) z'*M*z;

nT=1000; dt=0.01; t=0:dt:nT*dt;

z = zeros(n, nT+1);
energy = zeros(nT+1,1);
u = zeros(nT+1,1);

% Space mesh: only interior points
x = linspace(a+h,b-h,n);

% Set initial condition
z0(1:n) = (x.^2).*(1-x).^3;

i=1;
z(1:n,i) = z0; % Initial flux is zero
energy(i) = compute_energy(z(1:n,i));
u(i) = -K*z(1:n,i);
figure(3)
plot(x,z(1:n,i),'o-')
title('Initial condition')
xlabel('x')

% First time step: use BDF1 (Backward Euler)
beta = 1; a1   = -1;
Am = (1/(beta*dt))*M - A;
[L1,U1,P1,Q1] = lu(Am);

i = 2;
rhs = -(a1/(beta*dt))*M*z(:,i-1);
z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
energy(i) = compute_energy(z(1:n,i));
u(i) = -K*z(1:n,i);

% Remaining time steps: use BDF2
beta = 2/3; a1 = -4/3; a2 = 1/3;
Am = (1/(beta*dt))*M - A;
[L1,U1,P1,Q1] = lu(Am);

for i=3:nT+1
   rhs = -(1/(beta*dt))*M*(a1*z(:,i-1) + a2*z(:,i-2));
   z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
   energy(i) = compute_energy(z(1:n,i));
   u(i) = -K*z(1:n,i);
   if mod(i,10)==0
      figure(4)
      xx=[a,x,b]; zz=[0,z(1:n,i)',u(i)];
      subplot(1,3,1), plot(xx,zz,'-','LineWidth',2)
      title('Solution')
      xlabel('x'); ylabel('z')
      subplot(1,3,2), semilogy(1:i,energy(1:i),'-','LineWidth',2)
      title('Energy')
      xlabel('Time step')
      subplot(1,3,3), plot(1:i, u(1:i), 'LineWidth', 2)
      title('Control')
      xlabel('Time step')
   end
end
