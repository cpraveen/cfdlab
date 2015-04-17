% the minmod function
function f = minmod(a,b,c)

globals;

if a*b > 0 && b*c > 0
   f = sign(a) * min( abs([a b c]) );
else
   f = 0.0;
end
