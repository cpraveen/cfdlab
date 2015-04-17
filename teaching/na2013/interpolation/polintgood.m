xmin = 0; xmax = 2*pi;
fun = @(x) cos(x);

xx = linspace(xmin,xmax,100);
ye = fun(xx);

for i=1:6
   N=2*i;
   subplot(3,2,i)
   x=linspace(xmin,xmax,N+1);
   y=fun(x);
   P=polyfit(x,y,N);
   yy=polyval(P,xx);
   plot(x,y,'o',xx,ye,'--',xx,yy,'LineWidth',2);
   axis([xmin xmax -1.1 +1.1])
end

%print -dpdf polintgood.pdf
