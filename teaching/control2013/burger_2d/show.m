function show ( elements3, coordinates, u )
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
  trisurf ( elements3, coordinates(:,1), coordinates(:,2), u' );
  xlabel('x')
  ylabel('y')
  zlabel('z')
%
%  Define the initial viewing angle.
%
  view ( -67.5, 30 );

  title ( 'Solution to the Problem' )

  %figure(2)
  %hold off
  %pdecont(coordinates', elements3', u');
