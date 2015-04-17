% initial condition
function f = initcond(x)

globals;

n = length(x);

if testcase==lincon_sine
   for j=1:n
      f(j) = sin(2*pi*x(j));
   end
elseif testcase==lincon_hat
   for j=1:n
      if x(j) < 0
         nn = ceil( -x(j) ); x(j) = x(j) + nn;
      end
      if x(j) > 1
         nn = floor(1-x(j)); x(j) = x(j) + nn;
      end
      if x(j) < 0.25 || x(j) > 0.75
         f(j) = 0.0;
      elseif x(j) >= 0.25 && x(j) < 0.5
         f(j) = 4*(x(j) - 0.25);
      else
         f(j) = 1.0 - 4*(x(j) - 0.5);
      end
   end
elseif testcase==lincon_zalesak
   a=0.5; z=-0.7; delta=0.005; alpha=10; beta=log(2)/(36*delta^2);
   G = inline('exp(-beta*(x-z)^2)', 'x', 'beta', 'z');
   F = inline('sqrt(max([1-alpha^2*(x-a)^2, 0]))', 'x', 'alpha', 'a');
   for j=1:n
      if x(j) < -1
         nn = ceil( -(1+x(j))/2 ); x(j) = x(j) + 2*nn;
      end
      if x(j) > 1
         nn = floor((1-x(j))/2); x(j) = x(j) + 2*nn;
      end
      if x(j) <=-0.6 && x(j) >= -0.8
         f(j) = ( G(x(j),beta,z-delta) + G(x(j),beta,z+delta) + ...
                  4*G(x(j),beta,z) )/6;
      elseif x(j) <= -0.2 && x(j) >= -0.4
         f(j) = 1;
      elseif x(j) <=0.2 && x(j) >= 0
         f(j) = 1 - abs( 10*x(j) - 1 );
      elseif x(j) >= 0.4 && x(j) <= 0.6
         f(j) = (F(x(j),alpha,a-delta) + F(x(j),alpha,a+delta) + ...
                 4*F(x(j),alpha,a))/6;
      else
         f(j) = 0;
      end
   end

elseif testcase==burger_sine
   for j=1:n
      f(j) = 0.25 + 0.5 * sin(pi*x(j));
   end
else
   fprintf(1,'Unknown initial condition\n');
   pause
end
