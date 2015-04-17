function [M,A,B,Q,R,H]=matrix_fem(n,nu,alpha)
% INPUT : n : number of elements
%         nu : viscosity coefficient
%
% OUTPUT : Matrices for P1 FEM, Dirichlet BC at x=0 and x=1
%          Dirichlet control at x=1
%
%   dz/dt = A*z + B*u
%          M: Mass matrix
%          A: state matrix
%          B: control matrix
%
%  J = (1/2) int(z' * Q * z) dt + (1/2) int( u' * R * u) dt
%
% Mesh point spacing
h = 1/n;

e = ones(n-1,1);

% mass matrix
M = h*spdiags([e/6 2*e/3 e/6], -1:1, n-1, n-1);

A = (nu/h)*spdiags([e -2*e e], -1:1, n-1, n-1) + alpha*M;

% matrix in the LQR problem
Q = sparse(n-1,n-1);
R = 1;

B = sparse(n-1,1);
B(n-1,1) = nu/h;

% Observation: Heat flux at x=0 : nu*dz/dx(0,t)
H = sparse(1,n-1);

% CHOOSE ONE OF THE FOLLOWING

% First order finite difference z_x = (z1 - z0)/h
%H(1,1) = nu/h;

% Second order finite difference
% z_x(x0) = (1/(2h))*( -3*z0 + 4*z1 - z2 )
H(1,1) =  2*(nu/h);
H(1,2) = -0.5*(nu/h);
