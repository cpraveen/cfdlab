% RHS for linear model with feedback
% no noise
function dxdt = rhs_lpc(t,x,A,B,K)

dxdt = (A-B*K)*x;
