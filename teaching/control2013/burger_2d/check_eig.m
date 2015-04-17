   globals;
   parameters;

   read_data();

   [M,A,B,Q,H] = get_matrices();

%  number of eigenvalues to find
   ne = 5;
   
% Compute eigenvalues and eigenfunctions
  opts.p = 50;
  [V,D] = eigs(A, M, ne, 'LR', opts);
  e = diag(D)
  figure(11)
  plot(real(e),imag(e),'o')
  grid on

% Same for adjoint
  [W,D] = eigs(A', M, ne, 'LR', opts);
  E = diag(D)

% Number of unstable eigenvalues
  Nu = length(find(real(E) > 0))
  assert(Nu < ne)

% Get only unstable eigenvectors
  V = V(:,1:Nu);
  W = W(:,1:Nu);

% stabilizability
  disp('Stabilizability')
  B'*W

% Detectability
  disp('Detectability')
  H*V

% Orthonormalize
  for j=1:Nu
     a = V(:,j)'*M*W(:,j);
     V(:,j) = sign(a)*V(:,j)/sqrt(abs(a));
     W(:,j) =         W(:,j)/sqrt(abs(a));
  end

% Check: This should be identity matrix
  V'*M*W
