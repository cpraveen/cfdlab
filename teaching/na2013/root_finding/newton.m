r     = 2.0;    % initial guess for root
delta = 1.0e-6; % tolerance on function value
N     = 50;     % maximum number of iterations

f = @(x) x^2 - 2;
df= @(x) 2*x;
re= sqrt(2);

for n=1:N
   fr  = f(r);
   dfr = df(r);
   dr  = - fr/dfr;
   r   = r + dr;
   fprintf('%5d %20.12e %20.12e %20.12e\n',n,r,fr,r-re)
   if abs(fr) < delta
      break
   end
end

if n==N
   fprintf('Reached limit of iterations\n')
end
