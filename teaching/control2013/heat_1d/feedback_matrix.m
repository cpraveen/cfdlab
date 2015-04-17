function K = feedback_matrix(M,A,B,Q,R)

npts = size(B,1);
n = npts+1;
h = 1/n;

%[K,X] = lqr(full(M\A), full(M\B), full(Q), R);
%K = real(K);
[X,L,K] = care(full(A),full(B), full(Q),R,[],full(M));
K = sparse(K);

% plot control operator
k_EF = [0 K 0];
x = 0:h:1;
figure(1)
plot(x,k_EF,'o-')
hold all
xlabel('x','fontsize',24)
ylabel('k(x)','fontsize',24)
legend('K')


end
