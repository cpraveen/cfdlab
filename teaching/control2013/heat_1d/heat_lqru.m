% Heat equation using FEM and BDF
% Feedback control
% Stabilizing only the unstable components

clear all
close all

a     = 0; 
b     = 1; 
ni    = 100; 
mu    = 1; 
alpha = 0.4 + pi^2;

n = ni - 1;
h = (b-a)/ni;

% Generate the system matrices
[M,A,B,Q] = matrix_fem(ni,mu,alpha);

% uncontrolled eigenvalues
[V,D]=eigs(A,M,10,'LR');
eo=diag(D);

% Find number of unstable eigenvalues
nu = 1;
V = V(:,1:nu);
D = D(1:nu, 1:nu);

% Bu is the matrix acting on the control corresponding to the unstable 
% portion of the diagonalized system, 
Bu = V'*B;

% Solving the ARE for the unstable part of the diagonalized system
Qu = zeros(nu);
Ru = 1;
[Pu,L,G] = care(D,Bu,Qu,Ru);

disp('The initial unstable eigenvalues were')
eig(D)
disp('The unstable eigenvalues are modified to')
eig(D - Bu*Bu'*Pu)

% Matrix P for the initial system.
P = V*Pu*V'; 

% Feedback matrix for the original system
K = sparse(B'*P*M);
A = A - B*K;

% eigenvalues with feedback control
ec=eigs(A,M,10,'LR');

figure(1)
plot(real(eo),imag(eo),'o',real(ec),imag(ec),'*')
title('Eigenvalues')
legend('Original','With feedback')
grid on

% Function to compute energy
compute_energy = @(z) z'*M*z;

nT=2000; dt=0.01; t=0:dt:nT*dt;

z = zeros(n, nT+1);
energy = zeros(nT+1,1);
u = zeros(nT+1,1);

% Space mesh: only interior points
x = linspace(a+h,b-h,n);

% Set initial condition
z0(1:n) = sqrt(2)*sin(pi*x);
i=1;
z(1:n,i) = z0; % Initial flux is zero
energy(i) = compute_energy(z(1:n,i));
u(i) = -K*z(1:n,i);

figure(2)
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
      hold off
      subplot(1,3,2), semilogy(1:i,energy(1:i),'-','LineWidth',2)
      title('Energy')
      xlabel('Time step')
      subplot(1,3,3), plot(1:i, u(1:i), 'b-','LineWidth', 2)
      title('Control')
      xlabel('Time step')
      hold off
   end
end
