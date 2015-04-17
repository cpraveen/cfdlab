function square(n)

parameters;

xmin=0;xmax=a;
ymin=0;ymax=b;
nx=n;
ny=n;
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
c=1;
X=zeros(nx*ny,1);
Y=zeros(nx*ny,1);
for j=1:nx
   for k=1:ny
      X(c) = x(j);
      Y(c) = y(k);
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

% find neighbouring cell of boundary edges
for j=1:n_eb
   nbr = tri.edgeAttachments(eb(j,:));
   nbr = nbr{:};
   dx = X(eb(j,2)) - X(eb(j,1));
   dy = Y(eb(j,2)) - Y(eb(j,1));
   % order vertices so that domain is to left side
   v  = tri.Triangulation(nbr(1),:);
   dx1= sum(X(v))/3 - X(eb(j,1));
   dy1= sum(Y(v))/3 - Y(eb(j,1));
   if(dx1*dy - dx*dy1 > 0)
      tmp = eb(j,:);
      eb(j,1) = tmp(2);
      eb(j,2) = tmp(1);
   end
end

fid=fopen('coordinates.dat','w');
fprintf(fid,'%24.14e %24.14e\n', tri.X');
fclose(fid);

fid=fopen('elements3.dat','w');
fprintf(fid,'%8d %8d %8d\n', tri.Triangulation');
fclose(fid);

% We want bottom/top/right boundaries as dirichlet
xmid = 0.5*(X(eb(:,1)) + X(eb(:,2)));
id = find(xmid > eps);

fid=fopen('dirichlet.dat','w');
fprintf(fid,'%8d %8d\n', eb(id,:));
fclose(fid);

% Left boundary is neumann
in = find(xmid < eps);
fid=fopen('neumann.dat','w');
fprintf(fid,'%8d %8d\n', eb(in,:));
fclose(fid);
