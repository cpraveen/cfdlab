% Nonlinear model for the inverted pendulum, without any feedback control

close all

% File declaring the values of the various parameters to be used
parameters;

% Initial condition
x0 = [0; 0; 0; 0];

% Time dicretizition
tspan = [0:0.01:100];

% Nonlinear system without control
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x] = ode15s(@rhs_nlp,tspan,x0,options,M,m,l,g,k,c,I) ;

% Plotting the solution at final time
figure('Name','Evolution of state variables without control for nonlinear system')
subplot(2,2,1), plot(t,x(:,1), 'LineWidth',1), title('Position of cart');
subplot(2,2,2) ; plot(t,x(:,2), 'LineWidth',1) ; title('Speed of cart') ;
subplot(2,2,3) ; plot(t,x(:,3), 'LineWidth',1) ; title('Angle of pendulum') ;
subplot(2,2,4), plot(t,x(:,4), 'LineWidth',1), title('Angular speed of pendulum') ;


