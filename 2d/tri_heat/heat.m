% prob = 1: sine initial condition
% prob = 2: (1/pi)*atan(y/x)
function [error_inf, error_l2] = heat(prob,nx,ny,randf)

% define global variables
global X Y u ei eb n_ei n_eb nbr_ei nbr_eb v_ei v_eb tri gradu e_type
global norm_ei norm_eb np nt tarea carea rhs mu tri_normal
global problem cpen e_dir e_neu
global xmin xmax ymin ymax

problem = prob;

% BEGIN PARAMS================================================================
% 1 = explicit, 2 = implicit
timescheme=2

% Nitsche penalty
cpen = 10;

% Whether to use strong BC
strong = false;
% END PARAMS================================================================

if strong == true
   cpen = 0;
end

fprintf(1,'----------------------------------------------\n');
% domain size
assert(nx==ny);
assert(randf>=0 && randf<0.5);
if problem==1
   xmin=0; xmax=1;
   ymin=0; ymax=1;
   tf=0.05;
elseif problem==2
   xmin=-0.5; xmax=0.5;
   ymin=0; ymax=1;
   tf=1e20;
elseif problem==3 || problem==4
   xmin=0; xmax=1;
   ymin=0; ymax=1;
   tf=1e20;
   assert(mod(nx,2)==1);
end
h = (xmax-xmin)/(nx-1);

% make grid
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
c=1;
X=zeros(nx*ny,1);
Y=zeros(nx*ny,1);
for j=1:nx
   for k=1:ny
      X(c) = x(j);
      Y(c) = y(k);
      if j~=1 && j~=nx && k~=1 && k~=ny
         X(c) = X(c) + (2*rand()-1)*h*randf;
         Y(c) = Y(c) + (2*rand()-1)*h*randf;
      end
      c    = c + 1;
   end
end
tri=DelaunayTri([X,Y]);
triplot(tri)
e=tri.edges;           % all edges
eb=tri.freeBoundary;   % boundary edges

np   = length(X);   % no of vertices
nt   = size(tri,1); % no of triangles
fprintf(1,'Number of points         = %d\n', np);
fprintf(1,'Number of triangles      = %d\n', nt);

n_e  = size(e ,1); % no of all edges
n_eb = size(eb,1); % no of boundary edges

hmin = 1e20;

% find neighbouring cell of boundary edges
nbr_eb = zeros(n_eb,1);
norm_eb= zeros(n_eb,2);
for j=1:n_eb
   nbr = tri.edgeAttachments(eb(j,:));
   nbr = nbr{:};
   nbr_eb(j) = nbr(1);
   dx = X(eb(j,2)) - X(eb(j,1));
   dy = Y(eb(j,2)) - Y(eb(j,1));
   % order vertices 
   v  = tri.Triangulation(nbr(1),:);
   dx1= sum(X(v))/3 - X(eb(j,1));
   dy1= sum(Y(v))/3 - Y(eb(j,1));
   if(dx1*dy - dx*dy1 > 0)
      tmp = eb(j,:);
      eb(j,1) = tmp(2);
      eb(j,2) = tmp(1);
   end
   % compute normal
   dx = X(eb(j,2)) - X(eb(j,1));
   dy = Y(eb(j,2)) - Y(eb(j,1));
   norm_eb(j,:) = [dy, -dx];
   hmin = min(hmin, sqrt(dx^2 + dy^2));
end

% order all the edges
for j=1:n_e
   nbr = tri.edgeAttachments(e(j,:));
   nbr = nbr{:};
   dx = X(e(j,2)) - X(e(j,1));
   dy = Y(e(j,2)) - Y(e(j,1));
   % order vertices 
   v  = tri.Triangulation(nbr(1),:);
   dx1= sum(X(v))/3 - X(e(j,1));
   dy1= sum(Y(v))/3 - Y(e(j,1));
   if(dx1*dy - dx*dy1 > 0)
      tmp = e(j,:);
      e(j,1) = tmp(2);
      e(j,2) = tmp(1);
   end
end

ei=setdiff(e,eb,'rows'); % interior edges
n_ei = size(ei,1); % no of boundary edges
fprintf(1,'Number of all      edges = %d\n', n_e);
fprintf(1,'Number of boundary edges = %d\n', n_eb);
fprintf(1,'Number of interior edges = %d\n', n_ei);

assert(n_eb + n_ei == n_e);

% find neighbouring cell of interior edges
nbr_ei = zeros(n_ei,2);
norm_ei= zeros(n_ei,2);
for j=1:n_ei
   nbr = tri.edgeAttachments(ei(j,:));
   nbr = nbr{:};

   dx = X(ei(j,2)) - X(ei(j,1));
   dy = Y(ei(j,2)) - Y(ei(j,1));
   % order vertices 
   v  = tri.Triangulation(nbr(1),:);
   dx1= sum(X(v))/3 - X(ei(j,1));
   dy1= sum(Y(v))/3 - Y(ei(j,1));
   if(dx1*dy - dx*dy1 > 0)
      tmp = ei(j,:);
      ei(j,1) = tmp(2);
      ei(j,2) = tmp(1);
   end

   nbr_ei(j,1) = nbr(1);
   nbr_ei(j,2) = nbr(2);
   dx = X(ei(j,2)) - X(ei(j,1));
   dy = Y(ei(j,2)) - Y(ei(j,1));
   norm_ei(j,:) = [dy, -dx];
   hmin = min(hmin, sqrt(dx^2 + dy^2));
end

%figure(10)
%for j=1:n_eb
%   plot(X(eb(j,:)), Y(eb(j,:)), '-r','Linewidth',2);
%   hold on
%end
%for j=1:n_ei
%   plot(X(ei(j,:)), Y(ei(j,:)), '-b');
%   hold on
%end

% find vertex oppsite to interior edge
v_ei = zeros(n_ei,2);
for j=1:n_ei
   ve = ei(j,:);
   % left vertex
   vt = tri(nbr_ei(j,1),:);
   v_ei(j,1) = setdiff(vt,ve);
   % right vertex
   vt = tri(nbr_ei(j,2),:);
   v_ei(j,2) = setdiff(vt,ve);
end

% find vertex oppsite to boundary edge
v_eb = zeros(n_eb,1);
for j=1:n_eb
   ve = eb(j,:);
   % left vertex
   vt = tri(nbr_eb(j),:);
   v_eb(j) = setdiff(vt,ve);
end

% compute triangle area and dual cell area
tarea = zeros(nt,1);
carea = zeros(np,1);
tri_normal = zeros(nt,3,2);
for j=1:nt
   v = tri.Triangulation(j,:);
   dx1 = X(v(2)) - X(v(1));
   dy1 = Y(v(2)) - Y(v(1));
   dx2 = X(v(3)) - X(v(1));
   dy2 = Y(v(3)) - Y(v(1));
   tarea(j) = (dx1 * dy2 - dx2 * dy1) / 2.0;
   assert(tarea(j) > 0);
   carea(v) = carea(v) + tarea(j) / 3.0;
   tri_normal(j,3,:) = [+dy1, -dx1];
   tri_normal(j,2,:) = [-dy2, +dx2];
   tri_normal(j,1,:) = -(tri_normal(j,2,:) + tri_normal(j,3,:));
end
fprintf(1,'Min/max triangle  area = %e  %e\n', min(tarea), max(tarea));
fprintf(1,'Min/max dual cell area = %e  %e\n', min(carea), max(carea));
fprintf(1,'Mininimum edge length  = %e\n', hmin);

% Identity dirichlet and neumann edges
create_boundary_edges();

% Coefficients for RK3 CHECK
ark = [0, 1/2];
brk = 1 - ark;

u     = zeros(np,1);
rhs   = zeros(np,1);
gradu = zeros(nt,2);

% set viscosity
mu = zeros(nt,1);
compute_mu()

% Compute jacobian matrix
if timescheme==20
   fprintf(1,'Find matrices ...\n');
   eps=1e-20;
   A=zeros(np,np);
   for j=1:np
      u(:) = 0;
      u(j) = complex(u(j), eps);
      compute_rhs(0);
      Jcol = imag(rhs)/eps;
      A(:,j) = Jcol;
      fprintf(1,'%d ', j);
   end
   A=sparse(A);
   u(:) = 0;
   fprintf(1,'Find source vector ...\n');
   compute_rhs(0);
   b = rhs;
end

if timescheme==2 && strong==true
   [A,b] = compute_matrix_strong();
elseif timescheme==2 && strong==false
   [A,b] = compute_matrix_weak();
end
%[A1,b1] = compute_matrix();
%full(A)
%full(A1)
%AA = A - A1
%A-A'
%A1-A1'
%A+A1
%pause

% set initial condition
set_initial_condition();
figure(2)
%trisurf(tri,X,Y,u);
pdecont([X';Y'],tri.Triangulation',u,30)
title('Initial condition')

% set bc type for edges
%e_type = zeros(n_eb, 1);
%set_edge_type()

cfl = 0.1;

t = 0;
it=0;
if timescheme ==1
   dt= cfl * hmin^2 / max(mu);
else
   dt = 0.01*hmin;
end
res=1e20; res_tol=1e-8;

if timescheme==2
   A1 = sparse(diag(carea)) - 0.5*dt*A;
   A2 = sparse(diag(carea)) + 0.5*dt*A;
end

while t<tf && res>res_tol
   if(t+dt > tf) % adjust dt for last step to reach t=tf
      dt = tf - t;
   end
   if timescheme==1
      compute_rhs(t);
      u = u + dt*rhs./carea;
      res = norm(rhs);
   else
      u = A1 \ (A2*u + dt*b);
      res = norm(A*u+b);
   end
   t = t + dt;
   it= it+ 1;
   if strong == true
      ue = exact_solution(t,X,Y);
      u(eb(:,1))=ue(eb(:,1)); u(eb(:,2))=ue(eb(:,2));
   end
   energy(it,1) = t; energy(it,2) = 0.5*sum(u.^2 .* carea);
   fprintf(1,'%e %e %e %e %e\n', t, min(u), max(u), energy(it,2), res);
end
fprintf(1,'%e %e %e %e %e\n', t, min(u), max(u), energy(it,2), res);

% Get exact solution
ue = exact_solution(t,X,Y);
%ue = exact_cell_average(t);

% Uncomment this to not include boundary errors in norm calculation
%ue(eb(:,1))=u(eb(:,1)); ue(eb(:,2))=u(eb(:,2));

% compute errors
error_inf = max( abs(u-ue) );
error_l2  = sqrt(sum((u-ue).^2 .* carea));
fprintf(1,'Linf error = %e\n', error_inf);
fprintf(1,'L2   error = %e\n', error_l2);

figure(3)
pdecont([X';Y'],tri.Triangulation',u,30)
title('Numerical Solution','FontSize',16)
axis square

figure(4)
pdecont([X';Y'],tri.Triangulation',ue,30)
title('Exact Solution')

figure(5)
%pdecont([X';Y'],tri.Triangulation',u-ue,30)
pdeplot([X';Y'],[],tri.Triangulation','xydata',u-ue,'levels',30,'colorbar','on')
title('Error')

figure(6)
plot(energy(:,1), energy(:,2), 'LineWidth', 2);
title('Energy')
