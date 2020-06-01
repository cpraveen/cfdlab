function set_edge_type()

global e_type n_eb eb X Y problem

if problem==1 || problem==2
   e_type(:) = 0;
elseif problem==3 % bottom and top are neumann bc
   e_type(:)=0;
   for j=1:n_eb
      v = eb(j,:);
      ym = sum(Y(v))/2;
      if ym < 1e-10 || ym-1 > -1e-10
         e_type(j)=1;
      end
   end
end
