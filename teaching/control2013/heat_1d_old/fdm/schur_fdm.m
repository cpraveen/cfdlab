function [K,X] = schur_fdm(A,B,C,R,D,N)

npts = size(B,1);
n = npts+1;
h=1/n;

% avec lqr
[K,X] = lqr(A,B,C'*C,R+D^2,N);

% noyau sur [0,1]
k_DF = [0 K 0];
x = 0:h:1;

figure(1),
plot(x,k_DF,'.-g'),hold on
xlabel('x','fontsize',24),ylabel('k(x)','fontsize',24),
legend('Evolution du noyau k du gain statique en fct de x')

end
