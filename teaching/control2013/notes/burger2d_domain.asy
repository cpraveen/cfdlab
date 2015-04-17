size(200,200);

real b = 0.5; 

pair p0 = (0,0);
pair p1 = (1,0);
pair p2 = (1,b);
pair p3 = (0,b);

draw(p0--p1--p2--p3--cycle,linewidth(1.0));
draw((1,0)--(1,b), blue+linewidth(1.0));

label("$\Gamma_{d_1}$", (0.5*(p0+p3)), W);
label("$\Gamma_c$", (0.5*(p2+p3)), N);
label("$\Gamma_{d_2}$", (0.5*(p0+p1)), S);
label("$\Gamma_n$", (0.5*(p1+p2)), E);
label("$\Omega$", 0.25*(p0+p1+p2+p3));
