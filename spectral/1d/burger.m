% Solution for viscous burgers using spectral method
% u_t + (u^2/2)_x = nu u_xx in (0,1)
% IC: u(x,0) = sin(2*pi*x)
% BC: periodic
%
% The formulation is based on
%     Trefethen: Spectral Methods in Matlab
%
% Semi-implicit scheme: explicit for convection, implicit for diffusion
%
% Arguments: nu, N, cfl, nT
function burger(nu, N, cfl, nT)

  if nargin == 0
    nu = 0.1;
    N  = 2^5;
    nT = 500;
    cfl= 0.5;
  end

  fprintf(1,"nu = %e\n", nu);
  fprintf(1,"N  = %d\n", N);
  fprintf(1,"cfl= %f\n", cfl);

  xf = linspace(0,1,N+1);
  x  = xf(2:end);
  k  = [0:N/2-1 0 -N/2+1:-1];
  h  = 1/N;

  u  = sin(2*pi*x);
  uh = fft(u);
  plot(uh,'o');

  t  = 0;
  for it=1:nT
     dt = cfl*h/max(abs(u));
     w  = u.^2/2;
     wh = fft(w);
     uh = (uh - dt * 1i * k .* wh) ./ (1 + dt * nu * k.^2 );
     u  = real(ifft(uh));
     t  = t + dt; ts = strcat('Time = ', num2str(t));
     if(mod(it,10) == 0)
        subplot(2,1,1)
        plot(xf,[u(end) u],'o-');
        xlabel('x'); ylabel('u')
        title(ts); grid on
        subplot(2,1,2)
        plot(k(1:N/2+1), abs(uh(1:N/2+1)), 'o')
        xlabel('k'); ylabel('|u_k|')
        pause(0.5)
     end
  end

end
