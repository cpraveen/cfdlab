function init()

globals;

lincon = 1;
burger = 2;
no=0; yes=1;

% Test cases
burger_step = 11;
burger_sine = 12;
burger_rare = 13;
lincon_sine = 21;
lincon_step = 22;
lincon_hat  = 23;
lincon_zalesak  = 24;

% Constant for minmod
% Dont change this; set it in dg.m in the testcase selection step
mmconst = 1.0e20;

% Setup Gaussian Coefficients and Evaluation Points
x_values = zeros(10,10) ;
c_values = zeros(10,10) ;

i=1 ;
x_values(i,1) = 0.0 ;
c_values(i,1) = 2.0 ;

i=2 ;
x_values(i,1) = -0.5773502691896260 ; 
x_values(i,2) = -x_values(i,1) ;
c_values(i,1) = 1.0 ;
c_values(i,2) = 1.0 ;

i=3 ;
x_values(i,1) = -0.7745966692414830 ;
x_values(i,2) = 0.0 ;
x_values(i,3) = -x_values(i,1) ;
c_values(i,1) = 0.5555555555555560 ;
c_values(i,2) = 0.8888888888888890 ;
c_values(i,3) =c_values(i,1) ;

i=4 ;
x_values(i,1) = -0.8611363115940530 ;
x_values(i,2) = -0.3399810435848560 ;
x_values(i,3) = -x_values(i,2) ;
x_values(i,4) = -x_values(i,1) ;
c_values(i,1) = 0.3478548451374540 ;
c_values(i,2) = 0.6521451548625460 ;
c_values(i,3) =c_values(i,2) ;
c_values(i,4) =c_values(i,1) ;

i=5 ;
x_values(i,1) = -0.9061798459386640 ;
x_values(i,2) = -0.5384693101056830 ;
x_values(i,3) = 0.0 ;
x_values(i,4) = -x_values(i,2) ;
x_values(i,5) = -x_values(i,1) ;
c_values(i,1) = 0.2369368850561890 ;
c_values(i,2) = 0.4786386704993660 ;
c_values(i,3) = 0.5688888888888890 ;
c_values(i,4) =c_values(i,2) ;
c_values(i,5) =c_values(i,1) ;

i=6 ;
x_values(i,1) = -.9324695142032520 ;
x_values(i,2) = -.6612093864662650 ;
x_values(i,3) = -.2386191860831970 ;
x_values(i,4) = -x_values(i,3) ;
x_values(i,5) = -x_values(i,2) ;
x_values(i,6) = -x_values(i,1) ;
c_values(i,1) = 0.1713244923791700 ;
c_values(i,2) = 0.3607615730481390 ;
c_values(i,3) = 0.4679139345726910 ;
c_values(i,4) =c_values(i,3) ;
c_values(i,5) =c_values(i,2) ;
c_values(i,6) =c_values(i,1) ;

i=7 ;
x_values(i,1) = -0.9491079123427590 ;
x_values(i,2) = -0.7415311855993940 ;
x_values(i,3) = -0.4058451513773970 ;
x_values(i,4) = 0.0 ;
x_values(i,5) = -x_values(i,3) ;
x_values(i,6) = -x_values(i,2) ;
x_values(i,7) = -x_values(i,1) ;
c_values(i,1) = 0.1294849661688700 ;
c_values(i,2) = 0.2797053914892770 ;
c_values(i,3) = 0.3818300505051190 ;
c_values(i,4) = 0.4179591836734690 ;
c_values(i,5) =c_values(i,3) ;
c_values(i,6) =c_values(i,2) ;
c_values(i,7) =c_values(i,1) ;

i=8 ;
x_values(i,1) = -0.9602898564975360 ;
x_values(i,2) = -0.7966664774136270 ;
x_values(i,3) = -0.5255324099163290 ;
x_values(i,4) = -0.1834346424956500 ;
x_values(i,5) = -x_values(i,4) ;
x_values(i,6) = -x_values(i,3) ;
x_values(i,7) = -x_values(i,2) ;
x_values(i,8) = -x_values(i,1) ;
c_values(i,1) = 0.1012285362903760 ;
c_values(i,2) = 0.2223810344533740 ;
c_values(i,3) = 0.3137066458778870 ;
c_values(i,4) = 0.3626837833783620 ;
c_values(i,5) =c_values(i,4) ;
c_values(i,6) =c_values(i,3) ;
c_values(i,7) =c_values(i,2) ;
c_values(i,8) =c_values(i,1) ;

i=9 ;
x_values(i,1) = -0.9681602395076260 ;
x_values(i,2) = -0.8360311073266360 ;
x_values(i,4) = -0.6133714327005900 ;
x_values(i,4) = -0.3242534234038090 ;
x_values(i,5) = 0.0 ;
x_values(i,6) = -x_values(i,4) ;
x_values(i,7) = -x_values(i,3) ;
x_values(i,8) = -x_values(i,2) ;
x_values(i,9) = -x_values(i,1) ;
c_values(i,1) = 0.0812743883615740 ;
c_values(i,2) = 0.1806481606948570 ;
c_values(i,3) = 0.2606106964029350 ;
c_values(i,4) = 0.3123470770400030 ;
c_values(i,5) = 0.3302393550012600 ;
c_values(i,6) =c_values(i,4) ;
c_values(i,7) =c_values(i,3) ;
c_values(i,8) =c_values(i,2) ;
c_values(i,9) =c_values(i,1) ;

i=10 ;
x_values(i,1)  = -0.9739065285171720 ;
x_values(i,2)  = -0.8650633666889850 ;
x_values(i,3)  = -0.6794095682990240 ;
x_values(i,4)  = -0.4333953941292470 ;
x_values(i,5)  = -0.1488743389816310 ;
x_values(i,6)  = -x_values(i,5) ;
x_values(i,7)  = -x_values(i,4) ;
x_values(i,8)  = -x_values(i,3) ;
x_values(i,9)  = -x_values(i,2) ;
x_values(i,10) = -x_values(i,1) ;
c_values(i,1)  = 0.0666713443086880 ;
c_values(i,2)  = 0.1494513491505810 ;
c_values(i,3)  = 0.2190863625159820 ;
c_values(i,4)  = 0.2692667193099960 ;
c_values(i,5)  = 0.2955242247147530 ;
c_values(i,6)  = c_values(i,5) ;
c_values(i,7)  = c_values(i,4) ;
c_values(i,8)  = c_values(i,3) ;
c_values(i,9)  = c_values(i,2) ;
c_values(i,10) = c_values(i,1) ;

% Above values for interval [-1,+1]. Hence normalize to [0,1]
x_values = (x_values + 1.0)/2.0;
c_values = c_values/2.0;
