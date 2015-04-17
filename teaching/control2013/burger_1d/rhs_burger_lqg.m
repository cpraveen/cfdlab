% RHS of nonlinear model
% no control, no noise
function dxdt = rhs_burger_lqg(t,x,Lc,Uc,A1,A2,d1,d2,D1,H,A,B,K,L,N,nu,us,ws)                  

dxdt = zeros(2*N,1);

u = - K*x(N+1:2*N);

for i = 1:N
   dxdt(i,1) = -x(1:N)'*D1{i}*x(1:N) ...
             - (nu*A1(i,:) + us*A2(i,:)  + ws(2:N+1)*(D1{i} + D1{i}'))*x(1:N);
end

dxdt(1:N,1) = dxdt(1:N,1) - u*(A2*(x(1:N) + ws(2:N+1)') + 2*us*d1 + nu*d2) ...
              - u^2*d1; 

dxdt(N+1:2*N,1) = L*H*x(1:N) + (A-L*H-B*K)*x(N+1:2*N);

dxdt(1:N,1)     = Uc \ (Lc \ dxdt(1:N,1));
dxdt(N+1:2*N,1) = Uc \ (Lc \ dxdt(N+1:2*N,1));
