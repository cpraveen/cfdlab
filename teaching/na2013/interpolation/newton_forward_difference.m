% f(x) = sqrt(x) rounded to seven decimal places
% Only six decimals are accurate
x=2.0:0.1:2.7;
f=[1.414214,1.449138,1.483240,1.516575,1.549193,1.581139,1.612452,1.643168];
% Following is accurate to about 14 digits
%f=sqrt(x);

for d=1:length(f)-1
   fprintf('Divided difference of order %d\n',d)
   for i=1:length(f)-d
      f(i) = f(i+1) - f(i);
      fprintf('%5d %24.14f\n',i,f(i))
   end
end
