function K = feedback_matrix(M,A,B,Q)

npts = size(B,1);
n = npts+1;
h = 1/n;

[K,X] = lqr(full(M\A), full(M\B), full(Q), 1);
K = real(K)/M;

% plot control operator
k_EF = [0 K 0];
x = 0:h:1;
figure(1),
plot(x,k_EF,'o-r'),hold on
xlabel('x','fontsize',24),ylabel('k(x)','fontsize',24),
legend('K')

end
