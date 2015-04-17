function K = feedback_matrix(M,A,B,Q,R)

N = size(B,1);
h = 1/N;

%[K,X] = lqr(full(M\A), full(M\B), full(Q), R);
%K = real(K);
[X,L,K] = care(full(A),full(B), full(Q),R,[],full(M));

% plot control operator
x = h:h:1;
figure(10)
plot(x,K,'o-')
hold all
xlabel('x','fontsize',24)
ylabel('k(x)','fontsize',24)
legend('K')

K = sparse(K);

end
