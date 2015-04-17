% RHS of linear model
% no control, no noise
function dxdt = rhs_lp(t,x,A)

dxdt = A*x;
