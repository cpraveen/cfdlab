function f = fM1(M1)

global gamma

p2byp1 = 1 + (2*gamma/(gamma+1))*(M1^2 - 1);
M22 = (1 + 0.5*(gamma-1)*M1^2)/(gamma*M1^2 - 0.5*(gamma-1));
f = p2byp1*((1+0.5*(gamma-1)*M22)/(1+0.5*(gamma-1)*M1^2))^(gamma/(gamma-1));
