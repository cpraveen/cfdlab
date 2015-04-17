% Nonlinear model
% no control, no noise
function dxdt = fbo(t,x,M,m,l,g,k,c,I)                  

dxdt = zeros(4,1);

D = (M + m)*(I + m*l^2) - m^2*l^2*cos(x(3))^2; 

dxdt = [ x(2);...
         (m*l*cos(x(3))*(c*x(4) - m*g*l*sin(x(3))) + (I + m*l^2)*( - k*x(2) + m*l*x(4)^2*sin(x(3))))/D;...
         x(4);...
        ((M+m)*(-c*x(4) + m*g*l*sin(x(3))) - m*l*cos(x(3))*(-k*x(2) + m*l*x(4)^2*sin(x(3))))/D] ;
