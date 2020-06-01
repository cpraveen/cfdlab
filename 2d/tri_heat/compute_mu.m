function compute_mu()

global problem mu X Y tri nt

if problem==1 || problem==2
   mu(:) = 1.0;
elseif problem==3 || problem==4 % test case from Subramanian and Perot (2006)
   k1=1; k2=4;
   for j=1:nt
      % get three vertices of tri(j)
      v = tri(j,:);
      % compute centroid coordinates
      xc = sum(X(v))/3;
      yc = sum(Y(v))/3;
      if xc <= 0.5
         mu(j) = k1;
      else
         mu(j) = k2;
      end
   end
else
   fprintf(1,'compute_mu: unknown problem\n');
   pause
end
