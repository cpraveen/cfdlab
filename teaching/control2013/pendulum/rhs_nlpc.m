% RHS for nonlinear model with feedback
% no noise
function dxdt = rhs_nlpc(t,x,M,m,l,g,k,c,I,K,Alpha,Beta)

dxdt = zeros(4,1);

u = -K*x;
F = Alpha*u - Beta*x(2);


D = (M + m)*(I + m*l^2) - m^2*l^2*cos(x(3))^2; 

dxdt = [ x(2);...
         (m*l*cos(x(3))*(c*x(4) - m*g*l*sin(x(3))) + (I + m*l^2)*(F - k*x(2) + m*l*x(4)^2*sin(x(3))))/D;...
         x(4);...
        ((M+m)*(-c*x(4) + m*g*l*sin(x(3))) - m*l*cos(x(3))*(F-k*x(2) + m*l*x(4)^2*sin(x(3))))/D] ;

