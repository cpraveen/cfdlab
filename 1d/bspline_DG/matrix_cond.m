% Condition number of mass-matrix
%N = degree of bernstein
for N=1:20
   A=zeros(N+1,N+1);
   for m=0:N
      for n=0:N
         A(m+1,n+1) = nchoosek(N,m) * nchoosek(N,n) * ...
                      gamma(2*N + 1 - (m+n)) * gamma(m+n+1) / gamma(2+2*N);
      end
   end
   %inv(A)
   %pause

   x(N)=N; y(N)=cond(A);

   %weigths
   w=[];
   for n=0:N
      w1 = nchoosek(N,n)*gamma(N+1-n)*gamma(n+1)/gamma(2+N);
      w = [w w1];
   end
   fprintf(1,'%d %e %e\n', x(N), y(N), min(w));
end

plot(x,log10(y),'-o')
xlabel('Degree')
ylabel('log10(condition)')

%weights

for N=1:10

   %weigths for control points
   w=[];
   for n=0:N
      w1 = nchoosek(N,n)*gamma(N+1-n)*gamma(n+1)/gamma(2+N);
      w = [w w1];
   end

   x=linspace(0,N,N+1)/N;
   A=zeros(N+1,N+1);
   for m=0:N
      for n=0:N
         A(m+1,n+1) = bernstein(N,n,x(m+1));
      end
   end

   wt = w * inv(A);
   [N 1/min(wt) max(wt) sum(wt)]

end
