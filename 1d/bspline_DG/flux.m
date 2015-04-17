% Flux function
function f = flux(u)

globals;

if fluxfun==lincon
   f = u;        % Linear flux
elseif fluxfun==burger
   f = 0.5*u.^2; % Burgers flux
else
   fprintf(1,'Unknown flux function\n');
   exit
end
