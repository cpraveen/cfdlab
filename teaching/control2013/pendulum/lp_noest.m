% Linear model for inverted pendulum partial observation with noise
% We use the noisy observation to generate the control and solve without
% estimation
% We use implicit Euler integration


clear all
close all

parameters;
[A,B] = get_system_matrices();

% Finding matrix K
Q = eye(4);
Ru = 1/3^2;
[K,X] = lqr(A, B, Q, Ru);
disp('Eigenvalues of A-B*K')
eig(A-B*K)

% Observation operator; for full observation taking H = identity matrix
H = [1, 0, 0, 0;
     0, 0, 1, 0];
%H = eye(4);

% Function to compute energy
compute_energy = @(x) 0.5*(M*x(:,2).^2 + m*(x(:,2) ...
                     + l*x(:,4).*cos(x(:,3))).^2 + ...
                    (I+m*l^2).*x(:,4).^2 ) + m*g*l*cos(x(:,3));

% Number of time steps; N = number of states; No = number of observations
Nt = 25000; N = 4; No = length(H(:,1));

% Initial condition
x0 = [0; 0; 20*pi/180; 0];

% Noise in state
w = (5e-2)^2 * randn(Nt+1,N);
Rw = diag(std(w).^2);

% Noise in observation
v = (5e-2)^2 * randn(Nt+1,No);
Rv = diag(std(v).^2);

% Time discretization
dt = 0.01;
t = 0:dt:dt*(Nt-1);

% State variables, control variable and observation variables
z = zeros(N,Nt); 
u = zeros(1,Nt);
y = zeros(No,Nt);
energy = zeros(1,Nt);

% Initial conditions 
z(:,1) = x0;

Am = eye(4) - dt*(A);

nt = 1;

y(:,nt) = H*z(:,nt) + v(nt,:)';   
u(nt) = 0; %initial control is zero
energy(nt) = compute_energy(z(:,nt)');

% Implicit Euler for time integration
for nt=1:Nt-1    
    
    % calculate state with noise without estimation
    z(:,nt+1) = Am\(z(:,nt) + dt*w(nt+1,:)' + dt*B*u(nt));
    
    % observation without estimation
    y(:,nt+1) = H*z(:,nt+1) + v(nt+1,:)';   
    
    % calculate control without estimation
    u(nt+1) = -K*[y(1,nt+1); ...
                  (y(1,nt+1)-y(1,nt))/dt; ...
                  y(2,nt+1); ...
                  (y(2,nt+1)-y(2,nt))/dt];
    
    % calculating energy
    energy(nt) = compute_energy(z(:,nt)');
    
end;

% Plotting z
figure('Name','Evolution of state variables with control for nonlinear system')
subplot(2,2,1)
plot(t,z(1,:),'LineWidth',1)
legend('z without estimation')
title('Position of cart')
hold all
grid on
subplot(2,2,2)
plot(t,z(2,:),'LineWidth',1)
legend('z without estimation')
title('Speed of cart')
hold all
grid on
subplot(2,2,3)
plot(t,z(3,:),'LineWidth',1)
legend('z without estimation')
title('Angle of pendulum')
hold all
grid on
subplot(2,2,4)
plot(t,z(2,:),'LineWidth',1)
legend('z without estimation')
title('Angular speed of pendulum')
hold all
grid on

% Plotting control
figure(2)
hold all
plot(t,u,'LineWidth',1)
legend('Evolution of control without estimation')
grid on

% Plotting energy
figure(3)
hold all
plot(t,energy)
legend('Evolution of energy')
