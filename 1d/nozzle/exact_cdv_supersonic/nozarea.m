function [a] = nozarea(x)

global L xthroat athroat noztyp

if noztyp == 1

   % straight diverging nozzle
   ain = 1.0512;
   aout= 1.75;
   L   = 10.0;
   a = ain + x*(aout-ain)/L;

elseif noztyp == 2

   % another diverging nozzle, S-shaped
   ain = 1.0512;
   aout= 1.75;
   L   = 10.0;
   dc  = 0.8;
   dd  = 4.0;
   db  = (aout - ain)/(tanh(10.0*dc - dd) - tanh(-dd));
   da  = ain - db*tanh(-dd);
   a = da + db*tanh(dc*x - dd);

elseif noztyp == 3

   % Converging-diverging nozzle
   L=2.0;
   xthroat=L/2;
   aratio=2.0; % = ain/athroat
   athroat=1.0;
   ain = aratio*athroat;
   if x<0.25*L || x>0.75*L
      a = ain;
   else
      a = athroat + (ain-athroat)*(sin(2*pi*(x-L/2)/L))^2;
   end

elseif noztyp == 4

   % Converging-diverging nozzle
   L=2.0;
   if x<0.5 || x>1.5
      a = 2.0;
   else
      a = 1.0 + (sin(pi*(x-1.0)))^2;
   end

elseif noztyp == 5

   % Polynomial converging diverging nozzle
   L=4.0;
   ain=2.0;
   athroat=1.0;
   xx=x-2.0;
   if xx<-1 || xx>+1
      a = ain;
   else
      a = athroat + (ain-athroat)*(3*xx^2 - 3*xx^4 + xx^6);
   end

else

   fprintf(1,'Unknown nozzle type given\n');
   return

end
