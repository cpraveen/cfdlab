function create_boundary_edges()

global problem n_eb eb X Y e_dir e_neu
global xmin xmax ymin ymax

if problem==1 | problem==2 | problem==4
   e_dir = [1:n_eb];
   e_neu = [];
elseif problem==3
   xe = 0.5*(X(eb(:,1)) + X(eb(:,2)));
   ye = 0.5*(Y(eb(:,1)) + Y(eb(:,2)));
   e_dir = find( xe-xmin<eps | xe-xmax>-eps);
   e_neu = find( ye-ymin<eps | ye-ymax>-eps);
else
   fprintf(1,'Unknown problem\n');
   pause
end

fprintf(1,'Number of dirichlet edges = %d\n', length(e_dir));
fprintf(1,'Number of neumann   edges = %d\n', length(e_neu));
assert(length(e_dir)+length(e_neu) == length(eb));
