% interpolate f(x) = |x|, -1 <= x <= +1 using polynomials
% Try uniform and chebyshev data
% Press Enter key to try with increasing number of data points

for N=2:20
   x = linspace(-1,+1,N+1);
   %t = linspace(0,pi,N+1); x = -cos(t);
   y = abs(x);
   P = polyfit(x,y,N);

   xs= linspace(-1,+1,1000);
   ys= polyval(P, xs);
   plot(x,y,'o',xs,ys,xs,abs(xs),'LineWidth',1.5);
   set(gca,'FontSize',16)
   axis([-1, 1, -0.1, 1.1])
   legend('Data','Interpolation','Exact','Location','North')
   pause
end
