function dxdt = rhs_nlpest(t,x,A,B,H,L,M,l,m,g,k,c,I,K,Alpha,Beta)
             
dxdt = zeros(8,1);

% Compute feedback using estimator
u = -K*x(5:8);
F = Alpha*u - Beta*x(2);

% nonlinear model
D = (M + m)*(I + m*l^2) - m^2*l^2*cos(x(3))^2; 

dxdt = [ x(2);...
         (m*l*cos(x(3))*(c*x(4) - m*g*l*sin(x(3))) + (I + m*l^2)*(F - k*x(2) + m*l*x(4)^2*sin(x(3))))/D;...
         x(4);...
        ((M+m)*(-c*x(4) + m*g*l*sin(x(3))) - m*l*cos(x(3))*(F-k*x(2) + m*l*x(4)^2*sin(x(3))))/D] ;

% linear estimator
dxdt(5:8,1) = L*H*x(1:4) + (A-L*H-B*K)*x(5:8);
