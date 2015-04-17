function heat2d_est ( )
   
   globals;
   parameters;

   read_data();

   [M,A,B,Q,H] = get_matrices();
   
   dt = 0.01;
   Nt = 1000;
 
%==========================================================================
% EVALUATION OF THE FEEDBACK MATRIX
%==========================================================================
   
   %Unstable eigenvalues and eigenvectors
   %For omega=0.4 there is only one unstable eigenvalue
   nu = 1;
   [V,D] = eigs(A,M,nu,'lr');
   
   % Bu is the matrix acting on the control corresponding to the unstable 
   % portion of the diagonalized system, 
   Bu = V'*B;
   
   % Solving the ARE for the unstable part of the diagonalized system
   Qu = zeros(nu);
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
   K = ((B'*V)*Pu)*(V'*M);

%==========================================================================
   
   % Noise in state
   w = (5e-3)^2 * randn(Nt+1,nFreeNodes);
   Rw = diag(std(w).^2);
   Rw = sparse(Rw);
 
   % Noise in observation
   v = (5e-3)^2 * randn(Nt+1,3);
   Rv = diag(std(v).^2);
   Rv = sparse(Rv);  
%==========================================================================
% EVALUATION OF THE FILTERING GAIN MATRIX
%==========================================================================
   Hu = H*V;
   Rwu = V'*Rw*V;
   Peu = care(D',full(Hu)',full(Rwu),full(Rv));
   Pe = sparse(V*Peu*V');
   L = sparse((M*Pe*H')/Rv);

  zac = sparse ( 2*nNodes, 1 );        % Complete solution
  zc  = sparse ( 2*nFreeNodes, 1 );    % Interior solution
  u  = sparse ( nControlNodes, 1 ); % Control vector

  % Coupled system
  Ae = [A,   -B*K; ...
        L*H,  A-L*H-B*K];

  % Noise matrix  
  Be = [speye(nFreeNodes),   sparse(nFreeNodes,3); ...
        sparse(nFreeNodes,nFreeNodes), L];

  % Observation matrix
  He = [H,          sparse(3,nFreeNodes); ...
        sparse(3,nFreeNodes), H];
  
  Me = [M , sparse(nFreeNodes,nFreeNodes);  
        sparse(nFreeNodes,nFreeNodes), M];
 
  A1 = Me - dt*Ae;
  [L1,U1,P1,Q1] = lu(A1);
  zac(1:nNodes)= -cos(0.5*pi*coordinates(:,1)) .* sin(pi*coordinates(:,2));
  zac(nNodes+1:2*nNodes) = 0*zac(1:nNodes);
  zc(1:nFreeNodes) = zac(FreeNodes);
  zc(nFreeNodes +1:2*nFreeNodes) = 0*zac(nNodes + FreeNodes);
  u(:) = sin(pi*coordinates(ControlNodes,2));
  show_c ( elements3, coordinates, full(zac),nNodes);
  pause(2)

% Time loop
  t = 0;
  energy = zeros(Nt,1);
  time = zeros(Nt,1);
  for it = 1:Nt
     zc = Q1*(U1\(L1\(P1*Me*zc)));
     t = t + dt;
     zac(FreeNodes) = zc(1:nFreeNodes);
     zac(nNodes + FreeNodes) = zc(nFreeNodes+1:2*nFreeNodes);
     v = -K*zc(nFreeNodes+1:2*nFreeNodes);
     zac(ControlNodes) = v*u;
     zac(nNodes + ControlNodes) = v*u;
     time(it) = t;
     energy(it) = zc(1:nFreeNodes)'*M*zc(1:nFreeNodes);
     y = He*zc;
     fprintf(1,'Time = %e, Energy = %e, Obs = %e %e %e\n', ...
             t, energy(it), full(y(1)), full(y(2)), full(y(3)))
     if mod(it,10) == 0
      show_c ( elements3, coordinates, full(zac),nNodes);
     end
  end
  show_c ( elements3, coordinates, full(zac),nNodes);
  % Plot energy
  figure(10)
  semilogy(time, energy)
  xlabel('Time')
  ylabel('Energy')

end
