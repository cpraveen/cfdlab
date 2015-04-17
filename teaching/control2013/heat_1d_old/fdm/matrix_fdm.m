function [M,A,B,C,D,N]=matrice_chaleur_diri_DF(n,nu,alpha)
% INPUT : n : nombre d'intervalles -> (n+1) points de discrétisation
%         nu : viscosité
%         alpha : shift
% OUTPUT : Matrices de l'eq de la chaleur sur [0,1] avec CL de Dirichlet,
%          contrôle frontière en x=0.
%
% Schéma DF à 3 points.
%          M: matrice de masse (identité)
%          A: matrice de rigidité
%          B: matrice de contrôle
%          C: matrice d'observation
%          D: matrice de feedforward dans l'expression de l'observation à contrôler
%          N: matrice de coiplage entre l'état et le contrôle dans la fct
%          cout
%          R: matrice devant le terme quadratique du contrôle dans la fct
%          cout


% pas de la subdivision
h = 1/n;

% matrice de masse
M = speye(n-1);

e = ones(n-1,1);
A = nu*spdiags([e -2*e e], -1:1, n-1, n-1)/h^2 + alpha*M;

% matrice de contrôle
B = sparse(n-1,1);
B(1,1) = nu/h^2;

% matrice de sortie
C = h*ones(1,n-1);

% matrice provenant du couplage
D = h/2; 
D = sparse(D);

% matrice de couplage entre état et contrôle.
N = C'*D;
N = sparse(N);

% matrice de poids devant le terme quadratique du contrôle
%R = 1;
