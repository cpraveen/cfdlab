q1dnoz

xsm = xs - 1.0e-10;
[rl ul pl ml] = solution(xs,s1,r1,xsm);
xsp = xs + 1.0e-10;
[rr ur pr mr] = solution(xs,s1,r1,xsp);

xi = linspace(x1,x2,100);
xi = [xi xthroat-1e-5 xthroat+1e-5 xsm xsp];
xi = sort(xi);

I1 = zeros(size(xi));
I2 = zeros(size(xi));
I3 = zeros(size(xi));

% Computation of I1

for j=1:length(xi)

      I1(j) = quadgk(@g1,xi(j),x2);

end

% Computation of I2
% This is zero

% Computation of I3

for j=1:length(xi)

   I3(j) = quadgk(@g2,xi(j),x2);

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
plot(xi,v(:,1),'r',xi,v(:,2),'g',xi,v(:,3),'b','LineWidth',2)
legend('v_1','v_2','v_3')

