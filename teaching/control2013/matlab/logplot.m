x = linspace(0,1,100);
y = exp(-10*x);

subplot(1,2,1)
plot(x,y,'LineWidth',2)
title('Linear scale for x and y')
xlabel('x')
ylabel('y')

subplot(1,2,2)
semilogy(x,y,'LineWidth',2)
title('Linear scale for x, log scale for y')
xlabel('x')
ylabel('y')
