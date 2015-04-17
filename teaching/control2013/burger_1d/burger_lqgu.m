clear all
close all

nu = 1;
ni = 101; % number of nodes
N = ni -1;
h = 1/N;

% Space and time discretization
x = (0:h:1);
nT=500; dt=0.01; tspan=0:dt:nT*dt;

% Obtaining stationary solution and other stationary conditions
[ws,epsn,c,us,gs] = stationarysol(x,nu);

% Get system matrices
[M,A1,A2,D1,d1,d2,H] = get_system_mat(N);

% Assembling final system matrices
A =  zeros(N,N);
for i = 1:N
    A(i,:) = -nu*A1(i,:) - us*A2(i,:) - ws(2:N+1)*(D1{i} + D1{i}');
end
A = sparse(A);
B = -2*us*d1 - A2*ws(2:N+1)' - nu*d2;
B= sparse(B);

Ne = 5; % no of eigenvalues to find

% uncontrolled eigenvalues
[V,D]=eigs(A,M,Ne,'LR');
[W,D]=eigs(A',M,Ne,'LR');
eo=diag(D)
Nu = length(find(real(eo)>0))
assert(Ne > Nu) % did we get all unstable eigenvalues
V = V(:,1:Nu);
W = W(:,1:Nu);
for j=1:Nu
   a = W(:,j)'*M*V(:,j);
   V(:,j) = V(:,j)/a;
end
Au= D(1:Nu,1:Nu);
Bu= W'*B;

% Finding feedback matrix
Qu = zeros(Nu,Nu);
Ru= 1;
[Pu,E,K1] = care(Au,Bu,Qu,Ru)
P = W*Pu*V';
K = B'*P*M;
eigs(A-B*K,M,5,'LR')

% Noise in state
eta = (5e-2)^2 * randn(nT+1,N);
Reta = diag(std(eta).^2);
Reta = sparse(Reta);

% Noise in observation
mu = (5e-2)^2 * randn(nT+1,1);
Rmu = diag(std(mu).^2);
Rmu = sparse(Rmu);

Hu = H*V;
Qeta = V'*Reta*V;
[Pe,E,K2] = care(Au',Hu',Qeta,Rmu)
Pe = W*Pe*V';
L = (M*Pe*H') / Rmu;
eigs(A-L*H,M,5,'LR')

% Obtaining gain matrix L
%[S,T,L] = care(full(A)', full(H)', full(Reta), full(Rmu),[],full(M));
%L = real(L');

% Set initial condition
delta = 0.1; % This controls size of initial perturbation
z0(1:N) = delta*sin(pi*x(2:N+1)/2);
z0c = [z0(1:N)';zeros(N,1)];

Me = [M , sparse(N,N); ...
      sparse(N,N), M];

[Lc,Uc] = lu(M);

% Solving the system
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,zc] = ode15s(@rhs_burger_lqg,tspan,z0c,[],Lc,Uc,A1,A2,d1,d2,D1,H,A,B,K,L,N,nu,us,ws);

u = -K*zc(:,N+1:2*N)';

for i = 1:nT
    if mod(i,5)==0
        zz = [u(i), zc(i,1:N)];
        zze = [u(i), zc(i,N+1:2*N)];
        figure(2)
        subplot(1,2,1)
        plot(x,zz,'-','LineWidth',1)
        hold all
        plot(x,zze,'o','LineWidth',1)
        xlabel('x') 
        ylabel('z')
        titletext = sprintf('Solution at t=%f',tspan(i));
        title(titletext)
        legend('z','z_e')
        hold off
        
        subplot(1,2,2)
        plot(1:i, u(1:i), 'LineWidth', 1)
        xlabel('time') 
        ylabel('control')
    end
end
