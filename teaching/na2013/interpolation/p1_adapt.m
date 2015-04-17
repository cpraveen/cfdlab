clear all

xmin = 0;
xmax = 1;
fun = @(x) exp(-100*(x-0.5).^2);
fun = @(x) exp(-100*(x-0.5).^2) .* sin(4*pi*x);

N = 10;
x = linspace(xmin,xmax,N);
f = fun(x);

ne = 100;
xe = linspace(xmin,xmax,ne);
fe = fun(xe);
plot(xe,fe,'-',x,f,'or--','LineWidth',2);
title('Initial approximation')
pause

err = 1e-2;
nadapt = 50;

for n=1:nadapt
   h = x(2:end) - x(1:end-1);
   H = zeros(1,N-1);
   for j=2:N-2
      P = polyfit(x(j-1:j+2), f(j-1:j+2), 3);
      H(j) = max(abs(6*P(1)*x(j:j+1) + 2*P(2)));
   end
   elem_err = (1.0/8.0) * h.^2 .* abs(H);
   [current_err,i] = max(elem_err);
   if current_err < err
      fprintf('Satisfied error tolerance\n')
      break
   end
   % Following two lines for adaptation
   %x = [x(1:i), 0.5*(x(i)+x(i+1)), x(i+1:end)];
   %f = [f(1:i), fun(x(i+1)), f(i+1:end)];
   % Following two lines do uniform refinement
   x = linspace(xmin,xmax,2*N);
   f = fun(x);
   plot(xe,fe,'-',x,f,'or--','LineWidth',2);
   N = length(x);
   pause
end
fprintf('Number of points = %d\n',N)
