function f = charfun(x)

global xlower xupper

if x < xlower | x > xupper
   f = zeros(size(x));
else
   f = ones(size(x));
end
