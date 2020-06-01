uni=load('uni.dat');
figure(1)
h = uni(:,1);
error_inf = uni(:,2);
error_l2  = uni(:,4);
loglog(h, error_inf, 'o-', h, error_l2 , '*-', 'LineWidth', 2)
leg=legend('L^\infty error','L^2 error','Location','NorthWest');
set(leg,'FontSize',16)
xlabel('h','FontSize',16)
ylabel('Error','FontSize',16)
set(gca,'FontSize',16)
j=1;
res = [h(j), error_inf(j), 0, error_l2(j), 0];
for j=2:length(h)
   rate_inf = log(error_inf(j-1)/error_inf(j))/log(2);
   rate_l2  = log(error_l2(j-1)/error_l2(j))/log(2);
   res = [res; h(j), error_inf(j), rate_inf, error_l2(j), rate_l2];
end
res
print -dpdf prob1_uni.pdf


ran=load('rand.dat');
figure(2)
h = ran(:,1);
error_inf = ran(:,2);
error_l2  = ran(:,4);
loglog(h, error_inf, 'o-', h, error_l2 , '*-', 'LineWidth', 2)
leg=legend('L^\infty error','L^2 error','Location','NorthWest');
set(leg,'FontSize',16)
xlabel('h','FontSize',16)
ylabel('Error','FontSize',16)
set(gca,'FontSize',16)
j=1;
res = [h(j), error_inf(j), 0, error_l2(j), 0];
for j=2:length(h)
   rate_inf = log(error_inf(j-1)/error_inf(j))/log(2);
   rate_l2  = log(error_l2(j-1)/error_l2(j))/log(2);
   res = [res; h(j), error_inf(j), rate_inf, error_l2(j), rate_l2];
end
res
print -dpdf prob1_rand.pdf
