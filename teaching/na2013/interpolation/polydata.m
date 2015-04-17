% effect of uniform and chebychev spacing
% plots the function appearing in the interpolation error formula

N=16;
xu=linspace(-1,+1,N+1);
xc=cos( (2*(N - (0:N)) + 1)/(2*N + 2) * pi );

np=200;
x = linspace(-1,+1,np);
for i=1:np
   dx    = x(i) - xu;
   wu(i) = prod(dx);
   dx    = x(i) - xc;
   wc(i) = prod(dx);
end

plot(x, wu, '--', x, wc, 'LineWidth', 2);
set(gca,'FontSize',16)
legend('Uniform', 'Chebyshev')
title('N=16')
xlabel('x')
ylabel('\omega_N(x)')
%print -dpdf polydata.pdf
