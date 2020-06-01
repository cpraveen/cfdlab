function ue = exact_solution(t,X,Y)

global problem
%global X Y problem

if problem==1
   ue = exp(-2*4*pi^2*t)  * sin(2*pi*X) .* sin(2*pi*Y) ...
      + exp(-2*16*pi^2*t) * sin(4*pi*X) .* sin(4*pi*Y);
   ue = 50 * ue;
elseif problem==2
   ue = (1/pi)*atan2(Y,X);
   ii = find(abs(X)+abs(Y) < 1.0e-12);
   assert(size(ii,1)<=1);
   if size(ii,1) == 1
      ue(ii(1)) = 0.5;
   end
elseif problem==3
   k1=1; k2=4;
   ii = find(X <= 0.5);
   ue(ii) = (k2*X(ii) + 2*k1*k2) / (0.5*(k1+k2) + 4*k1*k2);
   ii = find(X > 0.5);
   ue(ii) = (k1*X(ii) + 2*k1*k2 + 0.5*(k2-k1)) / (0.5*(k1+k2) + 4*k1*k2);
   ue=ue';
elseif problem==4
   k1=1; k2=4;
   a = 1; b = 1; c = 1;
   ii = find(X <= 0.5);
   ue(ii) = a + b*X(ii) + c*Y(ii);
   ii = find(X > 0.5);
   ue(ii) = a - b*(k1-k2)/(2*k2) + b*(k1/k2)*X(ii) + c*Y(ii);
   ue=ue';
end
