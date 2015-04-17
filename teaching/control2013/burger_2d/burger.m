% Program to solve 2d Burger equation
%
% ws = stationary solution
%      -nu*(ws_xx + ws_yy) + ws*ws_x = fs
%
% z  = w - ws
%
% PDE for perturbation z
%         z_t + ws*z_x + z*ws_x + z*z_x - nu*(z_xx + z_yy) = 0
%
% FEM approximation
%        M z' = A z + N(z)
%
globals;
parameters;

read_data();

[M,A,B,Q,H] = get_matrices();

za = sparse ( nNodes, 1 );        % Complete solution
z  = sparse ( nFreeNodes, 1 );    % Interior solution
u  = sparse ( nControlNodes, 1 ); % Control vector

dt = 0.01;
Nt = 1000;
A1 = M - dt*A;
[L1,U1,P1,Q1] = lu(A1);
za(:)= -sin(0.5*pi*coordinates(:,1)/a) .* sin(pi*coordinates(:,2)/b);
z(:) = za(FreeNodes);
u(:) = sin(pi*coordinates(ControlNodes,1)/a); % sin(pi*x/a)
show ( elements3, coordinates, full(za) );
pause(2)

% Time loop
t = 0;
energy = zeros(Nt,1);
time = zeros(Nt,1);
for it = 1:Nt
   rhs = M*z;% + dt*burgN(za);;
   z = Q1 * (U1 \ (L1 \ (P1 * rhs) ) );
   t = t + dt;
   za(FreeNodes) = z;
   v = 0;
   za(ControlNodes) = v*u;
   time(it) = t;
   energy(it) = z'*M*z;
   y = H*z;
   fprintf(1,'Time = %e, Energy = %e, Obs = %e %e %e\n', ...
           t, energy(it), y(1), y(2), y(3))
   if mod(it,5) == 0
      show ( elements3, coordinates, full(za) );
      pause(1)
   end
end
show ( elements3, coordinates, full(za) );

% Plot energy
figure(10)
semilogy(time, energy)
xlabel('Time')
ylabel('Energy')
