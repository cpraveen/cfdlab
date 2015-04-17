% NOTE: Interior knots are repeated k+2 times. This is not necessary
% Enough to repeat k+1 times, just like the end knots
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

u=linspace(0,nk-1,1000);
for j=0:n

% control points
cont=zeros(1,n+1); cont(1,j+1)=1.0;

p=bspeval(k,cont,knot,u);

disp(j)
if j>=0 && j<=2
   subplot(3,1,1)
   if j==0
      plot(u,p,'-r','LineWidth',2)
   elseif j==1
      plot(u,p,'-b','LineWidth',2)
   else
      plot(u,p,'-m','LineWidth',2)
   end
   hold on
elseif j>=3 && j<=5
   subplot(3,1,2)
   if j==3
      plot(u,p,'-r','LineWidth',2)
   elseif j==4
      plot(u,p,'-b','LineWidth',2)
   else
      plot(u,p,'-m','LineWidth',2)
   end
   hold on
elseif j>=6 && j<=8
   subplot(3,1,3)
   if j==6
      plot(u,p,'-r','LineWidth',2)
   elseif j==7
      plot(u,p,'-b','LineWidth',2)
   else
      plot(u,p,'-m','LineWidth',2)
   end
   hold on
end
pause(1)

end

print -dpdf basis.pdf

figure(2)
cont=zeros(1,n+1); 
cont(1,1:3)=1.0;
cont(1,4:6)=2.0;
cont(1,7:9)=3.0;
p=bspeval(k,cont,knot,u);
plot(u,p,'LineWidth',2)

print -dpdf discfun.pdf
