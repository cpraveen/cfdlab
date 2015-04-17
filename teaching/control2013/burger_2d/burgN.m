function res = burgN(za)
globals;

res = zeros(nNodes,1);

MM = [2,1,1;1,2,1;1,1,2]/24;

for j=1:nTri
   v  = elements3(j,:);
   zx = za(v(1))*(coordinates(v(2),2) - coordinates(v(3),2)) ...
      + za(v(2))*(coordinates(v(3),2) - coordinates(v(1),2)) ...
      + za(v(3))*(coordinates(v(1),2) - coordinates(v(2),2));
   %zx = 0.5*zx/area(j);
   %Mloc = 2*area(j)*[2,1,1;1,2,1;1,1,2]/24;
   %res(v) = res(v) - zx*Mloc*za(v);
   res(v) = res(v) - zx*MM*za(v);
end
res = res(FreeNodes);
