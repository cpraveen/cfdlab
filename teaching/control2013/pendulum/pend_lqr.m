% Solving the linear and nonlinear inverted pendulum models with feedback
% control from the LQR approach

close all

% File declaring the values of the various parameters to be used
parameters;

% Evaluating the various matrices for the linearized system
[A,B] = get_system_matrices();

% Initial condition
x0 = [0; 0; 0; 0];

% Time dicretizition
tspan = [0:0.01:3];

% Obtaining the feedback matrix L
% Measurement matrix
C = [1 0 0 0; 
     0 0 1 0];
%C = eye(4);
Q = C'*C;
R = 0.01;
[K,S,E] = lqr(A, B, Q, R);

% Solving nonlinear system with control
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,x] = ode15s(@rhs_nlpc,tspan,x0,options,M,m,l,g,k,c,I,K,Alpha,Beta) ;

% Plotting the evolution of control
figure(1)
plot(t,-K*x','LineWidth',1)
title('Evolution of control for nonlinear system')
hold all

% Function to compute energy
compute_energy = @(x) 0.5*(M*x(:,2).^2 + m*(x(:,2) + ...
                           l*x(:,4).*cos(x(:,3))).^2 + ...
                           (I+m*l^2).*x(:,4).^2 ) + m*g*l*cos(x(:,3));

% Plotting the evolution of energy     
energy = compute_energy(x);
figure(2)
plot(t,energy','LineWidth',1)
title('Evolution of energy for nonlinear system')
hold all


% Plotting the solution at final time
%figure('Name','Evolution of state variables with control for nonlinear system')
figure(3)
subplot(2,2,1), plot(t,x(:,1), 'LineWidth',1), title('Position of cart'),
hold all;
subplot(2,2,2) ; plot(t,x(:,2), 'LineWidth',1) ; title('Speed of cart') ;
hold all;
subplot(2,2,3) ; plot(t,x(:,3), 'LineWidth',1) ; title('Angle of pendulum') ;
hold all;
subplot(2,2,4), plot(t,x(:,4),'LineWidth',1), title('Angular speed of pendulum') ;
hold all;

% Solving linear model with control
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,x] = ode15s(@rhs_lpc,tspan,x0,options,A,B,K);

% Plotting the solutions
%figure('Name','Evolution of state variables with control for linear system')
figure(4)
subplot(2,2,1), plot(t,x(:,1), 'LineWidth',1), title('Position of cart (linear)'),
hold all;
subplot(2,2,2) ; plot(t,x(:,2), 'LineWidth',1) ; title('Speed of cart (linear)') ;
hold all;
subplot(2,2,3) ; plot(t,x(:,3), 'LineWidth',1) ; title('Angle of pendulum (linear)') ;
hold all;
subplot(2,2,4), plot(t,x(:,4),'LineWidth',1), title('Angular speed of pendulum (linear)') ;
hold all;
