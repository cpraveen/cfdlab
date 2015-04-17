% Try two functions. For each try uniform data and chebyshev data
xmin=-1; xmax=+1;
fun = @(x) 1.0 ./ (1+16*x.^2);
%fun = @(x) exp(-5*x.^2);

xx = linspace(xmin,xmax,100);
ye = fun(xx);

for i=1:6
   N=2*i;
   subplot(3,2,i)
   %x=linspace(xmin,xmax,N+1);
   x=cos( (2*(N - (0:N)) + 1)/(2*N + 2) * pi );
   y=fun(x);
   P=polyfit(x,y,N);
   yy=polyval(P,xx);
   plot(x,y,'o',xx,ye,'--',xx,yy,'LineWidth',2);
   axis([xmin xmax -1 +1])
end

%print -dpdf polintbad.pdf
