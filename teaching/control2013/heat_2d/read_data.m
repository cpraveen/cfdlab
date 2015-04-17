function read_data ( )
   globals;
%
%  Read the nodal coordinate data file.
%
  load coordinates.dat;
%
%  Read the triangular element data file.
%
  eval ( 'load elements3.dat;', 'elements3=[];' );
%
%  Read the Neumann boundary condition data file.
%
  eval ( 'load neumann.dat;', 'neumann=[];' );
%
%  Read the Dirichlet boundary condition data file.
%
  eval ( 'load dirichlet.dat;', 'dirichlet=[];' );

end
