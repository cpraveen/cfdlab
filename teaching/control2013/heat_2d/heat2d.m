% Program to solve 2d heat equation without or with feedback control
% using full state information
%
% To run without feedback control
% >> heat2d(0)
% To run with feedback control
% >> heat2d(1)
%
function heat2d(to_control)
   globals;
   parameters;

   read_data();

   [M,A,B,Q,H] = get_matrices();
   
   if to_control==1 
%--------------------------------------------------------------------------
% Feedback matrix based on only the unstable components
%--------------------------------------------------------------------------
      % Unstable eigenvalues and eigenvectors
      Nu = 1;
      [V,D] = eigs(A,M,Nu,'lr');
      % Bu is the matrix acting on the control corresponding to the unstable 
      % portion of the diagonalized system, 
      Bu = V'*B;

      % Solving the ARE for the unstable part of the diagonalized system
      Qu = zeros(Nu);
      Ru = 1;
      [Pu,L,G] = care(D,Bu,Qu,Ru);

      disp('The initial unstable eigenvalues were')
      eig(D)
      disp('The unstable eigenvalues are modified to')
      eig(D - Bu*Bu'*Pu)

      % Matrix P for the initial system.
      % P = V*Pu*V'
      % Feedback matrix for the original system
      % K = B'*P*M;   

      % We directly compute K to avoid computing P which will be large matrix
      Pu= sparse(Pu);
      V = sparse(V);
      K = ((B'*V)*Pu)*(V'*M);

   else

      K = sparse(1, nFreeNodes);

   end
%--------------------------------------------------------------------------   

  za = sparse ( nNodes, 1 );        % Complete solution
  z  = sparse ( nFreeNodes, 1 );    % Interior solution
  u  = sparse ( nControlNodes, 1 ); % Control vector

  dt = 0.01;
  Nt = 10000;
  A1 = M - dt*(A-B*K);
  [L1,U1,P1,Q1] = lu(A1);
  za(:)= -cos(0.5*pi*coordinates(:,1)) .* sin(pi*coordinates(:,2));
  z(:) = za(FreeNodes);
  u(:) = sin(pi*coordinates(ControlNodes,2)); % sin(pi*y)
  show ( elements3, coordinates, full(za) );
  pause(2)

% Time loop
  t = 0;
  energy = zeros(Nt,1);
  time = zeros(Nt,1);
  for it = 1:Nt
     z = Q1 * (U1 \ (L1 \ (P1 * (M*z) ) ) );
     t = t + dt;
     za(FreeNodes) = z;
     v = -K*z;
     za(ControlNodes) = v*u;
     time(it) = t;
     energy(it) = z'*M*z;
     y = H*z;
     fprintf(1,'Time = %e, Energy = %e, Obs = %e %e %e\n', ...
             t, energy(it), y(1), y(2), y(3))
     if mod(it,500) == 0
      show ( elements3, coordinates, full(za) );
      pause(1)
     end
     v = -K*z;
  end
  show ( elements3, coordinates, full(za) );

  % Plot energy
  figure(10)
  semilogy(time, energy)
  xlabel('Time')
  ylabel('Energy')

end
