function [] = set_initial_condition()

global X Y u problem

fprintf(1,'Setting initial condition\n');

if problem == 1
   u = sin(2*pi*X) .* sin(2*pi*Y) + sin(4*pi*X) .* sin(4*pi*Y);
   u = 50 * u;
elseif problem == 2
   %u = exact_cell_average(0);
   u = exact_solution(0,X,Y);
   %u(size(X)') = 0;
elseif problem == 3 || problem == 4
   u = exact_solution(0,X,Y);
else
   fprintf(1,'Unknown problem %d\n', problem);
   pause
end

fprintf(1,'   min = %e\n', min(u));
fprintf(1,'   max = %e\n', max(u));
