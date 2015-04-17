function show_c ( elements3, coordinates, u, N )
%
%  Parameters:
%
%    Input, integer ELEMENTS3(N3,3), the nodes that make up each triangle.
%
%
%    Input, real COORDINATES(N,1:2), the coordinates of each node.
%
%    Input, real U(N), the finite element coefficients which represent the solution.
%    There is one coefficient associated with each node.
%
%
%  Display the information associated with triangular elements.
%
  figure(1)
  hold off
  subplot(2,1,1)
  trisurf ( elements3, coordinates(:,1), coordinates(:,2), u(1:N)' );
  view ( -67.5, 30 );
  title('z')
  subplot(2,1,2)
  trisurf ( elements3, coordinates(:,1), coordinates(:,2), u(N+1:2*N)');
  view ( -67.5, 30 );
  title ( 'z_e' )

%   figure(2)
%   hold off
%   subplot(2,1,1)
%   contour(coordinates(:,1), coordinates(:,2), u(1:N)');
%   title('z')
%   subplot(2,1,2)
%   contour(coordinates(:,1), coordinates(:,2), u(N+1:2*N)');
%   title ( 'z_e' )