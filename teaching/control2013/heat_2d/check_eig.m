   globals;
   parameters;

   read_data();

   [M,A,B,Q,H] = get_matrices();
   
% Compute eigenvalues and eigenfunctions
  [V,D] = eigs(A, M, 5, 'LA');
  e = diag(D);
  figure(11)
  plot(real(e),imag(e),'o')
  grid on
  fprintf(1,'Most unstable eigenvalue (numerical) = %f\n', e(1))
  fprintf(1,'Most unstable eigenvalue (exact)     = %f\n', -pi^2/40 + omega)
  ef = zeros(nNodes,1);
  ef(FreeNodes) = V(:,1);
  show ( elements3, coordinates, full ( ef ) );
