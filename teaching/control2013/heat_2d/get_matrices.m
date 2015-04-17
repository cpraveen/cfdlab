function [M,A,B,Q,H] = get_matrices ( )
  globals;

  nNodes = size(coordinates,1)
  nTri = size(elements3,1)


  Af = sparse ( nNodes, nNodes );
  Mf = sparse ( nNodes, nNodes );
%
%  Assembly.
%
  for j = 1 : nTri
%   Stiffness matrix
    Af(elements3(j,:),elements3(j,:)) = Af(elements3(j,:),elements3(j,:)) ...
      + nu*stima3(coordinates(elements3(j,:),:));
%   Mass matrix
    Mf(elements3(j,:),elements3(j,:)) = Mf(elements3(j,:), elements3(j,:)) ...
       + det([1,1,1;coordinates(elements3(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
  end
%
%  Determine which nodes are associated with Dirichlet conditions.
%  and which are free
  BoundNodes = unique ( dirichlet );
  FreeNodes = setdiff ( 1:nNodes, BoundNodes );

% Find nodes for which x = a
  ControlNodes = find ( coordinates(:,1) - a > -eps );

% Get matrix for free nodes only
  M = Mf(FreeNodes, FreeNodes);
  A = -Af(FreeNodes, FreeNodes) + omega * M;
  B = -Af(FreeNodes, ControlNodes) + omega * Mf(FreeNodes, ControlNodes);

% One dimensional control u(y,t) = v(t) sin(pi*y)
  B = sparse(B * sin(pi*coordinates(ControlNodes,2)));

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
