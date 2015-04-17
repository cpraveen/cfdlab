% Heat equation using FEM and BDF
% Partial observation and noise

clear all
%close all

a     = 0; 
b     = 1; 
ni    = 100; 
mu    = 1; 
alpha = 10;

n = ni - 1;
h = (b-a)/ni;

% Generate the system matrices
[M,A,B,Q,R,H] = matrix_fem(ni,mu,alpha);

% uncontrolled eigenvalues
eo=eigs(A,M,10,'LA');

% Obtaining control feedback matrix K
K = feedback_matrix(M,A,B,Q,R);
ec=eigs(A-B*K,M,10,'LR');

figure(2)
plot(real(eo),imag(eo),'o',real(ec),imag(ec),'*')
legend('Eigenvalue of (A,M)','Eigenvalue of (A-BK,M)')
grid on

% Parameters for time integration. nT = number of time steps
nT=6000; dt=0.01; t=0:dt:nT*dt;

% Noise in state
w = (5e-2)^2 * randn(nT+1,n);
Rw = diag(std(w).^2);
Rw = sparse(Rw);

% Noise in observation
v = (5e-2)^2 * randn(nT+1,1);
Rv = diag(std(v).^2);
Rv = sparse(Rv);

% Obtaining gain matrix L
[S,T,L] = care(full(A'), full(H'), full(Rw), full(Rv),[],full(M));
L = real(L');
L = sparse(L);

% Check that (A-LH,M) is stable
el=eigs(A-L*H, M, 10, 'LR');
figure(3)
plot(real(el),imag(el),'o')
grid on
legend('Eigenvalue of (A-LH,M)')

% Function to compute energy
compute_energy = @(z) z'*M*z;

% Mass matrix  
Me = [M ,          sparse(n,n); ...
      sparse(n,n), M];

% Coupled system
Ae = [A,   -B*K; ...
      L*H,  A-L*H-B*K];

% Compute eigenvalue of extended system
ee = eigs(Ae,Me,10,'LR');
figure(4)
plot(real(ee),imag(ee),'o')
grid on
legend('Eigenvalue of (Ae,Me)')

% Noise matrix
Be = [speye(n),   sparse(n,1); ...
      sparse(n,n), L];

% Observation matrix
He = [H,          sparse(1,n); ...
      sparse(1,n), H];

% Space mesh: only interior points
x = linspace(a+h,b-h,n);

% N : number of state variables, no : number of observation
N = n; no=1;
z0 = sin(pi*x'); % Initial profile
%z0(1:N) = sin(pi*x).*(1 + 0.01*randn(1,N)); % Initial profile

z = zeros(2*N,nT+1);
y = zeros(2*no,nT+1);
energy = zeros(nT+1,1);
u = zeros(nT+1,1);
zerr = zeros(nT+1,1);

% Set initial condition
i=1;
z(:,i) = [z0+0.1*sin(2*pi*x'); z0]; 
energy(i) = compute_energy(z(1:n,i));
u(i) = -K*z(N+1:2*N,i);
zerr(i) = norm(z(1:N,i) - z(N+1:2*N,i));

% Plot initial condition
figure(5)
plot(x,z(1:N,1),'o-',x,z(N+1:end,1),'*--')
legend('True system','Estimator')
title('Initial condition')
xlabel('x'); ylabel('z')

% First time step: use BDF1 (Backward Euler)
beta = 1; a1   = -1;
Am = (1/(beta*dt))*Me - Ae;
[L1,U1,P1,Q1] = lu(Am);

i = 2;
rhs = -(a1/(beta*dt))*Me*z(:,i-1) + Be*[w(i,:)';v(i,:)'];
z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
energy(i) = compute_energy(z(1:N,i));
u(i) = -K*z(N+1:2*N,i);
y(:,i) = He*z(:,i) + [v(i,:)';0];

% Remaining time steps: use BDF2
beta = 2/3; a1 = -4/3; a2 = 1/3;
Am = (1/(beta*dt))*Me - Ae;
[L1,U1,P1,Q1] = lu(Am);

for i=3:nT+1
   rhs = -(1/(beta*dt))*Me*(a1*z(:,i-1) + a2*z(:,i-2)) + Be*[w(i,:)';v(i,:)'];
   z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
   energy(i) = compute_energy(z(1:N,i));
   u(i) = -K*z(N+1:2*N,i);
   zerr(i) = norm(z(1:N,i) - z(N+1:2*N,i));
   if mod(i,50)==0
      figure(6)
      xx=[a,x,b]; zz=[0,z(1:N,i)',u(i)];
      zze = [0,z(N+1:2*N,i)',u(i)];
      subplot(1,3,1), plot(xx,zz,'-','LineWidth',1)
      hold all
      plot(xx,zze,'o','LineWidth',1),
      title('Solution')
      xlabel('x'); ylabel('z')
      legend('z','ze')
      hold off
      subplot(1,3,2), semilogy(1:i,energy(1:i),'-','LineWidth',1)
      title('Energy')
      xlabel('Time step')
      subplot(1,3,3)
      [AX,H1,H2] = plotyy(1:i, u(1:i), 1:i, log10(zerr(1:i)));
      set(get(AX(1),'Ylabel'),'String','Control u') 
      set(get(AX(2),'Ylabel'),'String','log||z-z_e||')
      title('Control and log||z-z_e||')
      xlabel('Time step')
   end
end
