size(200,200);

real b = 0.5; 

pair p0 = (0,0);
pair p1 = (1,0);
pair p2 = (1,b);
pair p3 = (0,b);

draw(p0--p1--p2--p3--cycle,linewidth(1.0));

real y0=b/2.0;
real y1=b/6.0;

draw((1,y0)--(1,y1), red+linewidth(2));

label("$\Gamma_o$", (1,0.5*(y0+y1)), E);
label("$\Omega$", 0.25*(p0+p1+p2+p3));
