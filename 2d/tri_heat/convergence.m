n = [21 41 81 161];

error_inf = zeros(length(n),1);
error_l2  = zeros(length(n),1);
ord_inf   = zeros(length(n),1);
ord_l2    = zeros(length(n),1);
h         = zeros(length(n),1);
for j=1:length(n)
   [a,b] = heat(2,n(j),n(j),0.0);
   error_inf(j) = a;
   error_l2(j) = b;
   h(j) = 1/(n(j)-1);
   if j>=2
      ord_inf(j) = log(error_inf(j-1)/error_inf(j))/log(2);
      ord_l2(j)  = log(error_l2(j-1)/error_l2(j))/log(2);
      [h(j)  error_inf(j)  ord_inf(j)  error_l2(j)  ord_l2(j)]
   end
end

fprintf(1,'----------------------------------------------------------------\n')
[h  error_inf  ord_inf  error_l2  ord_l2]
fprintf(1,'----------------------------------------------------------------\n')

figure(50)
loglog(h, error_inf, 'o-', h, error_l2 , '*-', 'LineWidth', 2)
legend('L^\infty error','L^2 error','Location','NorthWest')
