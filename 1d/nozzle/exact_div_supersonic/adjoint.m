
q1dnoz

global xlower xupper

xsm = xs - 1.0e-10;
[rl ul pl ml] = solution(xs,s1,r1,xsm);
xsp = xs + 1.0e-10;
[rr ur pr mr] = solution(xs,s1,r1,xsp);

xi = linspace(x1,x2,100);
xi = [xi xsm xsp];
xi = sort(xi);

I1 = zeros(size(xi));
I2 = zeros(size(xi));
I3 = zeros(size(xi));

% limits for integral
%xlower = x1;
%xupper = x2;
xlower = 2.0;
xupper = 4.0;

gg1 = @(x)( charfun(x).*g1(x) );
gg2 = @(x)( charfun(x).*g2(x) );

% Computation of I1

for j=1:length(xi)

   if xi(j) <= xs

      a1 = 1; c1 = 0;
      sol = inv([1,  g5(xsp); g1(x2),  g2(x2)])*[1; 0];
      a2=sol(1); c2=sol(2);

      if xs > xlower & xs < xupper
         delta = shockdelta(a1,c1,c2);
      else
         delta = 0.0;
      end

      I1(j) = quadgk(gg1,xi(j),xsm) + a2*quadgk(gg1,xsp,x2) + ...
              c2*quadgk(gg2,xsp,x2) - (pr-pl)*delta;

   elseif xi(j) > xs

      a1 = 0; b1 = 0; c1 = 0;
      b2 = 0; 
      b3 = 0;
      AA = [-1, 1, 0; 1, 0, g5(xsp); 0, g1(x2), g2(x2)];
      sol=inv(AA)*[1; 0; 0];
      a2=sol(1); a3=sol(2); c2=sol(3); c3=sol(3);

      if xs > xlower & xs < xupper
         delta = shockdelta(a1,c1,c2);
      else
         delta = 0.0;
      end

      I1(j) = a2*quadgk(gg1,xsp,xi(j)) + a3*quadgk(gg1,xi(j),x2) + ...
              c2*quadgk(gg2,xsp,x2) - (pr-pl)*delta;
   end

end

% Computation of I2
% This is zero

% Computation of I3

for j=1:length(xi)

   if xi(j) <= xs

      a1 = 0; b1 = 0; c1 = 1;
      b2 = 0;
      sol = inv([1,  g5(xsp); g1(x2),  g2(x2)])*[g5(xsm); 0];
      a2=sol(1); c2=sol(2);

      if xs > xlower & xs < xupper
         delta = shockdelta(a1,c1,c2);
      else
         delta = 0.0;
      end


      I3(j) = quadgk(gg2,xi(j),xsm) + a2*quadgk(gg1,xsp,x2) + ...
              c2*quadgk(gg2,xsp,x2) - (pr-pl)*delta;

   elseif xi(j) > xs

      a1 = 0; b1 = 0; c1 = 0;
      b2 = 0; 
      b3 = 0;
      AA = [0, -1, 1; 1, g5(xsp), 0; g1(x2), 0, g2(x2)];
      sol=inv(AA)*[1; 0; 0];
      a2=sol(1); c2=sol(2); c3=sol(3);
      a3=a2;

      if xs > xlower & xs < xupper
         delta = shockdelta(a1,c1,c2);
      else
         delta = 0.0;
      end


      I3(j) = a2*quadgk(gg1,xsp,x2) + c2*quadgk(gg2,xsp,xi(j)) + ...
              c3*quadgk(gg2,xi(j),x2) - (pr-pl)*delta;
   end

end

figure(6)
plot(xi,I1,'r',xi,I2,'g',xi,I3,'b','LineWidth',2)
legend('I1','I2','I3')

% Compute adjoint
v = zeros(length(xi),3);

for j=1:length(xi)

   a = nozarea(xi(j));
   [r u p m] = solution(xs,s1,r1,xi(j));
   totp = p*(1 + 0.5*(gamma-1)*m^2)^(gamma/(gamma-1));
   f1 = [1; u; H0];
   f2 = (a/2/H0)*[-r*u; 0; r*u*H0];
   f3 = (a/totp)*[r*u; p+r*u*u; r*u*H0];

   vv = [I1(j) I2(j) I3(j)]*inv([f1 f2 f3]);
   v(j,1) = vv(1);
   v(j,2) = vv(2);
   v(j,3) = vv(3);
end

figure(7)
subplot(3,1,1)
plot(xi,v(:,1),'LineWidth',2)
axis square
subplot(3,1,2)
plot(xi,v(:,2),'LineWidth',2)
axis square
subplot(3,1,3)
plot(xi,v(:,3),'LineWidth',2)
axis square

