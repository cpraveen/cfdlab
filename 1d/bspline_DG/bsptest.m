clear all
close all

set(0,'DefaultAxesColorOrder',[1 0 1; 0 1 1; 0 1 0])

k=2; % degree
n=5; % n+1 = no of control points
m=k + n + 1; % m+1 = no of knots

% knot vector
knot=[0 0 0 0.5 0.5 0.5 1 1 1];

u=linspace(0,1,200);
for j=0:n

% control points
cont=zeros(1,n+1); cont(1,j+1)=1.0;

p=bspeval(k,cont,knot,u);

plot(u,p)
hold on

end
