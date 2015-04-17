% Stabilizing the inverted pendulum. Equations of the linearized model are
%
%        x_t = A x + B u + w
%          y = H x + v
%
% u : state
% w : noise in the model
% v : noise in the observation
%--------------------------------------------------------------------------
% INPUT : m : mass of the pendulum (kg)
%         M : mass of the cart (kg)
%         k : coefficient of friction between the cart and the floor
%         d : coefficient of friction due to air on the rotating pendulum
%         I : moment of inertie kg/m^2
%         l : half the length of the pendulum axis 
%         g : acceleration due to gravity
%         (alpha, beta) used to express the force applied as a function of
%         voltage V(t) and dx / dt by the relation: 
%         F(t) = alpha*V(t) - beta*x'(t)
%
% OUTPUT : Matrices for the model describing the linearized system
% 
%         x1 : position of cart
%         x2 : speed of cart
%         x3 : angle of the pendulum wrt tht vertical
%         x4 : angular speed of the pendulum
%
% The outputs are measured to estimate y1 = x1 and y2 = x3.
%
% The values of the input paramters have been declared in parameter.m
%--------------------------------------------------------------------------


function [A,B] = get_system_matrices()

parameters;

% auxiliary variables of the model
v1 = (M+m)/(I*(M+m)+l^2*m*M);
v2 = (I+l^2*m)/(I*(M+m)+l^2*m*M);

% Matrix A
A = [0   1                      0                       0; ...
     0 -k*v2           -(l*m)^2*g*v2/(I+l^2*m)  l*m*c*v2/(I+l^2*m);...
     0   0                      0                       1; ...
     0 l*m*k*v1/(M+m)        l*m*g*v1                 -c*v1];

% Matrix B
B = [0; v2 ; 0 ;  -l*m*v1/(M+m)];

% CHECK: B contains alpha which should not be added to A
A = A + B*[0 -Beta 0 0];

B = B*Alpha;