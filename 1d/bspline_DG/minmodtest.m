function tolimit = minmodtest(ux)

globals;

tolimit = false;

if abs(ux(1)) < mmconst*h*h
   return
end

sgn = sign(ux(1));

for j=2:length(ux)
   if sign(ux(j)) ~= sgn
      tolimit = true;
      return
   end
end
