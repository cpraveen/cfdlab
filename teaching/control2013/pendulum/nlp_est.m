% Solving coupled estimation and control problem for the nonlinear system

clear all
close all

parameters;
[A,B] = get_system_matrices();

C = [1 0 0 0; 
     0 0 1 0];
Q = C'*C; 
Ru = 1/3^2;

[K,X] = lqr(A, B, Q, Ru);
disp('Eigenvalues of A-B*K')
eig(A-B*K)

% Observation operator
H = [1, 0, 0, 0;
     0, 0, 1, 0];

% Initial condition
x0 = [0; 0; 0; 0];

% Covariance matrix for noise in state
Rw = (5e-2)^2 * diag([1, 1, 1, 1]);

% Covariance matrix for noise in observation
Rv = (5e-2)^2 * diag([1, 1]);

% Adding noise to the initial conditions
x0c = [x0; 0*x0];

[L,S] = lqr(A', H', Rw, Rv);
L = real(L');

disp('Eigenvalues of A-L*H')
eig(A-L*H)

% Time dicretizition
tspan = [0:0.01:10];

% Nonlinear system without estimation and control with noise full
% information
%options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,zc] = ode15s(@rhs_nlpest,tspan,x0c,[],A,B,H,L,M,l,m,g,k,c,I,K,Alpha,Beta);

figure(1),
subplot(2,2,1)
hold all
plot(t,zc(:,1),t,zc(:,5))
legend('z','z_e')
title('Position of cart')
grid on
subplot(2,2,2)
hold all
plot(t,zc(:,2),t,zc(:,6))
legend('z','z_e')
title('Speed of cart')
grid on
subplot(2,2,3)
hold all
plot(t,zc(:,3),t,zc(:,7))
legend('z','z_e')
title('Angle of Pendulum')
grid on
subplot(2,2,4)
hold all
plot(t,zc(:,4),t,zc(:,8))
legend('z','z_e')
title('Anglular speed of Pendulum')
grid on

u = -K*zc(:,5:8)';

% Plotting control
figure(2)
plot(t,u)
legend('Control from estimation problem')
