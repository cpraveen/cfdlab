function [M,A,B,Q,H] = get_matrices ( )
  globals;

  nNodes = size(coordinates,1)
  nTri = size(elements3,1)

  Af = sparse ( nNodes, nNodes );
  Mf = sparse ( nNodes, nNodes );
  area = zeros(nTri,1);
%
%  Assembly.
%
  for j = 1 : nTri
%    area of triangle
     area(j) = 0.5*det([1,1,1;coordinates(elements3(j,:),:)']);
%    Mass matrix
     Mloc = 2*area(j)*[2,1,1;1,2,1;1,1,2]/24;
     Mf(elements3(j,:),elements3(j,:)) = Mf(elements3(j,:), elements3(j,:)) ...
        + Mloc;
%    centroid of triangle
     coor = coordinates(elements3(j,:),:);
     cmid = sum(coor,1)/3.0;
     [ws,wsx] = stationarysol(cmid(1),cmid(2));
%    Dloc = int( phi_i * d(phi_j)/dx )
     Dloc = coordinates(elements3(j,1),2) * [ 0 -1  1; 0 -1  1; 0 -1  1] ...
          + coordinates(elements3(j,2),2) * [ 1  0 -1; 1  0 -1; 1  0 -1] ...
          + coordinates(elements3(j,3),2) * [-1  1  0;-1  1  0;-1  1  0];
%   Stiffness matrix
    Af(elements3(j,:),elements3(j,:)) = Af(elements3(j,:),elements3(j,:)) ...
      + nu*stima3(coordinates(elements3(j,:),:)) ...
      + ws*Dloc/6.0 ...
      + wsx*Mloc;
  end
%
%  Determine which nodes are associated with Dirichlet conditions.
%  and which are free
  BoundNodes = unique ( dirichlet );
  FreeNodes = setdiff ( 1:nNodes, BoundNodes );

% Find nodes for which y = b
  ControlNodes = find ( coordinates(:,2) - b > -eps );

% Get matrix for free nodes only
  M = Mf(FreeNodes, FreeNodes);
  A = -Af(FreeNodes, FreeNodes) + omega * M;
  B = -Af(FreeNodes, ControlNodes) + omega * Mf(FreeNodes, ControlNodes);

% One dimensional control u(x,t) = v(t) sin(pi*x)
  B = sparse(B * sin(pi*coordinates(ControlNodes,1)));

  nBoundNodes = size(BoundNodes,1)
  nFreeNodes = size(FreeNodes,2)
  nControlNodes = size(ControlNodes,1)

% Observations on left boundary
  nObs = 3;
  H = sparse ( nObs, nNodes );
  ymid = 0.5*( coordinates(neumann(:,1),2) + coordinates(neumann(:,2),2) );
  ymin = [0.20, 0.50, 0.80];
  ymax = [0.25, 0.55, 0.85];
  for j=1:nObs
    ie = find(ymid > ymin(j) & ymid < ymax(j));
    ds = abs(coordinates(neumann(ie,2),2) - coordinates(neumann(ie,1),2));
    H(j,neumann(ie,1)) = H(j,neumann(ie,1)) + 0.5*ds';
    H(j,neumann(ie,2)) = H(j,neumann(ie,2)) + 0.5*ds';
    ds = sum(ds);
    H(j,:) = H(j,:)/ds;
  end
  H = H(:, FreeNodes);

% LQR problem
  Q = sparse(nFreeNodes,nFreeNodes);

end
