clear all
close all

parameters;

% Initial condition
x0 = [0.5; 0.2; 0.4; 1];

% Time interval and times at which solution desired
tspan = [0:0.01:20];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% without control
[t,x] = ode15s(@fbo,tspan,x0,options,M,m,l,g,k,c,I) ;

figure(1),  title('Evolution of state variables without control'),
         subplot(2,2,1), plot(t,x(:,1)), title('Position of cart'),
         subplot(2,2,2) ; plot(t,x(:,2)) ; title('Speed of cart') ;
         subplot(2,2,3) ; plot(t,x(:,3)) ; title('Angle of pendulum') ;
         subplot(2,2,4), plot(t,x(:,4)), title('Angular speed of pendulum') ;
