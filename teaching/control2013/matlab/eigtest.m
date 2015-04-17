n = 100;
xmin = 0;
xmax = 1;
h = (xmax-xmin)/n;
x = linspace(xmin+h, xmax-h, n-1);

e = ones(n-1,1);
A = -(1/h^2)*spdiags([e,-2*e,e], -1:1, n-1, n-1);

% Compute all eigenvalues
[V,D] = eig(full(A));
l = diag(D);

% Exact eigenvalues
le = pi^2 * ((1:(n-1)).^2)';

% Compute error in eigenvalues
err = abs(l - le);
figure(200)
plot(err,'o')
grid on
xlabel('n')
ylabel('Error in eigenvalue \lambda_n')

% Plot exact and numerical eigenfunctions
% Exact eigenfunctions are sin(j*pi*x), j=1,2,...
figure(1)
ve = sin(pi*x); ve = ve/norm(ve);
plot(x,ve,'--',x,V(:,1)','o','LineWidth',2)
xlabel('x'), ylabel('u_1')
title('First eigenfunction')
legend('Exact','Finite difference')

figure(2)
ve = -sin(2*pi*x); ve = ve/norm(ve);
plot(x,ve,'--',x,V(:,2)','o','LineWidth',2)
xlabel('x'), ylabel('u_2')
title('Second eigenfunction')
legend('Exact','Finite difference')

figure(3)
ve = -sin(3*pi*x); ve = ve/norm(ve);
plot(x,ve,'--',x,V(:,3)','o','LineWidth',2)
xlabel('x'), ylabel('u_3')
title('Third eigenfunction')
legend('Exact','Finite difference')

figure(4)
ve = sin(4*pi*x); ve = ve/norm(ve);
plot(x,ve,'--',x,V(:,4)','o','LineWidth',2)
xlabel('x'), ylabel('u_4')
title('Fourth eigenfunction')
legend('Exact','Finite difference')
