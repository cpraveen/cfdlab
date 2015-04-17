for N=2:20
   x = linspace(-1,+1,N+1); x = x';
   V=[];
   for j=0:N
      V = [V x.^j];
   end
   fprintf(1,'%d %e\n', N, cond(V));
end
