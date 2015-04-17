clear all
close all

k=2; % degree
nk=4; % no of distinct knots

m = (k+1) + (nk-2)*(k+1) + (k+1) - 1;

n = m - k - 1; % n+1 = no of control points

% knot vector
knot=[zeros(1,k+1)];
for j=2:nk-1
   knot=[knot (j-1)*ones(1,k+1)];
end
knot=[knot (nk-1)*ones(1,k+1)];

u=linspace(0,nk-1,200);
for j=0:n

% control points
cont=zeros(1,n+1); cont(1,j+1)=1.0;

p=bspeval(k,cont,knot,u);

disp(j)
plot(u,p)
hold on
pause(2)

end

figure(2)
cont=zeros(1,n+1); 
cont(1,1:3)=1.0;
cont(1,4:6)=2.0;
cont(1,7:9)=3.0;
p=bspeval(k,cont,knot,u);
plot(u,p)
