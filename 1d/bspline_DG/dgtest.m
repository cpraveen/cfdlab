clear all
close all

globals;
init();

% Test case
testcase = burger_step;
%testcase = burger_sine;
%testcase = burger_rare;
%testcase = lincon_sine;
%testcase = lincon_step;
%testcase = lincon_hat;
%testcase = lincon_zalesak;

% CFL condition
%cflmode = 'zhang';
cflmode = 'praveen';

% Number of elements
N = [20 40 80 160];
N = [40];

% Degree of bspline
p = 2;

ndof = [];
err  = [];
for j=1:length(N)
   [a1, a2] = dg(p,N(j),cflmode);
   ndof = [ndof a1];
   err  = [err a2];
   if length(N) > 1
      figure(10)
      loglog(ndof,err,'o-');
   end
end

% Compute convergence order
for j=2:length(ndof)
   pp =-(log10(err(j)) - log10(err(j-1)))/(log10(ndof(j)) - log10(ndof(j-1)))
end
