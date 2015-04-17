%
% Compare trignometric and polynomial interpolation of a periodic function
%
fun = @(x) max(0,1-abs(x-pi)/2);       % cont, diff a.e.
fun = @(x) exp(sin(x));                % inf diff, periodic
fun = @(x) ( abs(x-pi) < 0.5*pi ) + 0; % Hat function

ng = 500;
xg = linspace(0,2*pi,ng);
fe = fun(xg);

N = [6 8 10 12 14 16 18 20 22 24]; 
data = [];
for n=N
   h = 2*pi/n; 
   x = h*(0:n-1)';
   v = fun(x);

   v_hat = fft(v);
   k = [0:n/2, -n/2+1:-1];

   fg = exp(1i*xg'*k) * v_hat / n;
   fg = real(fg);
   err_trig = max(abs(fe - fg'));

   % polynomial interpolation at chebyshev points
   xp = cos( (2*(n - (0:n)) + 1)/(2*n + 2) * pi );
   xp = pi*(1+xp);
   vp = fun(xp);
   P = polyfit(xp,vp,n-1);
   fp = polyval(P,xg);
   err_poly = max(abs(fe - fp));

   figure(1)
   plot(x,v,'o',xg,fg,'-',xg,fp,'--',xg,fe,':')
   legend('Data trig','Trigonometric','Polynomial','Exact')
   data = [data; n, err_trig, err_poly];
   pause
end

figure(2)
semilogy(data(:,1), data(:,2), 'o-', data(:,1), data(:,3), '*--')
legend('Trigonometric','Polynomial')
