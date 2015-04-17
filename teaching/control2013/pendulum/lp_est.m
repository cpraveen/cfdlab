% Linear model for the inverted pendulum, partial observation with noise
% Solving the coupled estimation and control problem
% Time integration using BDF scheme

clear all
%close all
parameters;
[A,B] = get_system_matrices();

Q = eye(4);
Ru = 1/3^2;

[K,X] = lqr(A, B, Q, Ru);
disp('Eigenvalues of A-B*K')
eig(A-B*K)

% Observation operator
H = [1, 0, 0, 0;
     0, 0, 1, 0];
%H = eye(4);

% Number of time steps
Nt = 2500; N = 4; No = length(H(:,1));

% Initial condition
z0 = [0; 0; 20*pi/180; 0];

% Noise in state
w = (5e-2)^2 * randn(Nt,N);
Rw = diag(std(w).^2);

% Noise in observation
v = (5e-2)^2 * randn(Nt,No);
Rv = diag(std(v).^2);

[L,S] = lqr(A', H', Rw, Rv);
L = real(L');

disp('Eigenvalues of A-L*H')
eig(A-L*H)

% Coupled system
Ae = [A,   -B*K; ...
      L*H,  A-L*H-B*K];

disp('Eigenvalues of coupled system')
coupled_ev = eig(full(Ae))

% Noise matrix
Be = [eye(N,N),   zeros(N,No); ...
      zeros(N,N), L];

% Observation matrix
He = [H,          zeros(No,N); ...
      zeros(No,N), H];
  
% Function to compute energy
compute_energy = @(x) 0.5*(M*x(:,2).^2 + m*(x(:,2) ...
                     + l*x(:,4).*cos(x(:,3))).^2 + ...
                      (I+m*l^2).*x(:,4).^2 ) + m*g*l*cos(x(:,3));

% N : number of state variables, no : number of observation
dt = 0.01;
t = 0:dt:dt*(Nt-1);

zc = zeros(2*N,Nt);
u = zeros(1,Nt);
yc = zeros(2*No,Nt);
energy = zeros(1,Nt);

% Set initial condition
i=1;
zc(:,i) = [z0; zeros(N,1)] ; % Initial condition
u(i) = -K*zc(N+1:2*N,i);
energy(i) = compute_energy(zc(N+1:2*N,i)');

% First time step: use BDF1 (Backward Euler)
beta = 1; a1   = -1;
Am = (1/(beta*dt))*eye(2*N) - Ae;

i = 2;
rhs = -(a1/(beta*dt))*zc(:,i-1) + Be*[w(i,:)';v(i,:)'];
zc(:,i) = Am\rhs;
u(i) = -K*zc(N+1:2*N,i);
yc(:,i) = He*zc(:,i) + [v(i,:)';zeros(No,1)];
energy(i) = compute_energy(zc(N+1:2*N,i)');


% Remaining time steps: use BDF2
beta = 2/3; a1 = -4/3; a2 = 1/3;
Am = (1/(beta*dt))*eye(2*N) - Ae;

for i=3:Nt
   rhs = -(1/(beta*dt))*(a1*zc(:,i-1) + a2*zc(:,i-2)) + Be*[w(i,:)';v(i,:)'];
   zc(:,i) = Am\rhs;
   u(i) = -K*zc(N+1:2*N,i);
   yc(:,i) = He*zc(:,i) + [v(i,:)';zeros(No,1)];    
   energy(i) = compute_energy(zc(N+1:2*N,i)');
end

% Plotting all state variables
figure(1),
title('Noisy state')
subplot(2,2,1)
hold all
plot(t,zc(1,:),t,zc(5,:))
legend('z','z_e')
title('Position of cart')
grid on
subplot(2,2,2)
hold all
plot(t,zc(2,:),t,zc(6,:))
legend('z','z_e')
title('Speed of cart')
grid on
subplot(2,2,3)
hold all
plot(t,zc(3,:),t,zc(7,:))
legend('z','z_e')
title('Angle of Pendulum')
grid on
subplot(2,2,4)
hold all
plot(t,zc(4,:),t,zc(8,:))
legend('z','z_e')
title('Anglular speed of Pendulum')
grid on

% Plotting control
figure(2)
hold all
plot(t,u)
legend('Control from estimation problem')

% Plotting energy
figure(3)
hold all
plot(t,energy)
legend('Evolution of energy')
