% RHS of nonlinear model
% no control, no noise
function dxdt = rhs_burger(t,x,L,U,A1,A2,D1,N,nu,us,ws)                  

dxdt = zeros(N,1);

for i = 1:N
    dxdt(i,1) = -x'*D1{i}*x ...
                - (nu*A1(i,:) + us*A2(i,:)  + ws(2:N+1)*(D1{i} + D1{i}'))*x;
end

dxdt = U\(L\dxdt);
