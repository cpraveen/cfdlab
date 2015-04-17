% Partial observation with noise

clear all
close all

a = 0;
b = 1;
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
[M,A1,A2,D1,d1,d2,H] = get_system_mat(N);

% Assembling final system matrices
A =  zeros(N,N);
for i = 1:N
    A(i,:) = -nu*A1(i,:) - us*A2(i,:) - ws(2:N+1)*(D1{i} + D1{i}');
end
A = sparse(A);
B = -2*us*d1 - A2*ws(2:N+1)' - nu*d2;
B= sparse(B);

% uncontrolled eigenvalues
eo=eig(full(A),full(M));

% Finding feedback matrix
Q = sparse(N,N);
R = 1;
K = feedback_matrix(M,A,B,Q,R);

ec=eig(full(A-B*K),full(M));

figure(1)
plot(real(eo),imag(eo),'o',real(ec),imag(ec),'*')
grid on
legend('A','A-BK')
title('Eigenvalues')

% Noise in state
eta = (5e-2)^2 * randn(nT+1,N);
Reta = diag(std(eta).^2);
Reta = sparse(Reta);

% Noise in observation
mu = (5e-2)^2 * randn(nT+1,1);
Rmu = diag(std(mu).^2);
Rmu = sparse(Rmu);

% Obtaining gain matrix L
[S,T,L] = care(full(A)', full(H)', full(Reta), full(Rmu),[],full(M));
L = real(L');
L = sparse(L);

el=eig(full(A-L*H),full(M));
figure(2)
plot(real(el),imag(el),'o')
grid on
title('Eigenvalues of (A-LH,M)')

% Coupled system
Ae = [A,   -B*K; ...
      L*H,  A-L*H-B*K];

% Noise matrix
Be = [speye(N),   sparse(N,1); ...
      sparse(N,N), L];

% Observation matrix
He = [H,          sparse(1,N); ...
      sparse(1,N), H];
  
Me = [M , sparse(N,N); ...
      sparse(N,N), M];

ee=eig(full(Ae),full(Me));
figure(3)
plot(real(ee),imag(ee),'o')
grid on
title('Eigenvalues of coupled system')

% Set initial condition
delta = 0.1; % This controls size of initial perturbation
z0(1:N) = delta*sin(pi*x(2:N+1)/2);

z = zeros(2*N,nT+1);
y = zeros(2,nT+1);
u = zeros(nT+1,1);

% Set initial condition
i=1;
z(:,1) = [z0';zeros(N,1)]; 
u(i) = -K*z(N+1:2*N,i);

figure(4)
plot(x(2:N+1),z(1:N,1),'o-',x(2:N+1),z(N+1:end,1),'*--')
legend('True model','Estimator')
title('Initial condition')
xlabel('x')

% First time step: use BDF1 (Backward Euler)
beta = 1; a1 = -1;
Am = (1/(beta*dt))*Me - Ae;
[L1,U1,P1,Q1] = lu(Am);

i = 2;
rhs = -(a1/(beta*dt))*Me*z(:,i-1) + Be*[eta(i,:)';mu(i,:)'];
z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
u(i) = -K*z(N+1:2*N,i);
y(:,i) = He*z(:,i) + [mu(i,:)';0];

% Remaining time steps: use BDF2
beta = 2/3; a1 = -4/3; a2 = 1/3;
Am = (1/(beta*dt))*Me - Ae;
[L1,U1,P1,Q1] = lu(Am);

for i=3:nT+1
   rhs = -(1/(beta*dt))*Me*(a1*z(:,i-1) + a2*z(:,i-2)) + Be*[eta(i,:)';mu(i,:)'];
   z(:,i) = Q1 * (U1 \ (L1 \ (P1 * rhs)));
   u(i) = -K*z(N+1:2*N,i);
   y(:,i) = He*z(:,i) + [mu(i,:)';0];
   if mod(i,5)==0
      figure(5)
      zz=[u(i),z(1:N,i)'];
      zze = [u(i),z(N+1:2*N,i)'];
      subplot(2,2,1), plot(x,zz,'-','LineWidth',1)
      hold all
      plot(x,zze,'o','LineWidth',1),
      grid on
      title('Solution')
      xlabel('x'); ylabel('z')
      legend('z','ze')
      hold off
      subplot(2,2,2), plot(1:i, u(1:i), 'LineWidth', 1)
      grid on
      title('Control')
      xlabel('Time step')
      subplot(2,2,3)
      plot(1:i,y(1,1:i),'r-',1:i,y(2,1:i),'b--')
      grid on
      legend('True observation','Estimated observation')
      xlabel('Time')
   end
end
