function u = boundary_value(time, iedge, ipoint)

global X Y eb problem

xp = X(eb(iedge,ipoint));
yp = Y(eb(iedge,ipoint));

if problem == 1
   u = 0;
elseif problem == 2
   if yp < eps
      xm = 0.5*(X(eb(iedge,1)) + X(eb(iedge,2)));
      ym = 0.5*(Y(eb(iedge,1)) + Y(eb(iedge,2)));
      if xm < 0
         u = 1.0;
      elseif xm > 0
         u = 0.0;
      else
         u = 0.5;
      end
   else
      u = (1.0/pi)*atan2(yp, xp);
   end
elseif problem == 3
   k1=1; k2=4;
   if xp <= 0.5
      u = (k2*xp + 2*k1*k2) / (0.5*(k1+k2) + 4*k1*k2);
   else
      u = (k1*xp + 2*k1*k2 + 0.5*(k2-k1)) / (0.5*(k1+k2) + 4*k1*k2);
   end
elseif problem == 4
   k1=1; k2=4;
   a=1; b=1; c=1;
   if xp <= 0.5
      u = a + b*xp + c*yp;
   else
      u = a - b*(k1-k2)/(2*k2) + b*(k1/k2)*xp + c*yp;
   end
else
   fprintf(1,'Unknown problem = %d\n', problem);
end
