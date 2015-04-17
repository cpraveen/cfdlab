size(200,200);

pair p0 = (0,0);
pair p1 = (1,0);
pair p2 = (1,1);
pair p3 = (0,1);

draw(p0--p1--p2--p3--cycle,linewidth(1.0));

real y0=0.2;
real y1=0.25;
real y2=0.5;
real y3=0.55;
real y4=0.8;
real y5=0.85;

draw((0,y0)--(0,y1), red+linewidth(2));
draw((0,y2)--(0,y3), red+linewidth(2));
draw((0,y4)--(0,y5), red+linewidth(2));

label("$y_1$", (0, 0.5*(y0+y1)), W);
label("$y_2$", (0, 0.5*(y2+y3)), W);
label("$y_3$", (0, 0.5*(y4+y5)), W);
label("$\Omega$", 0.25*(p0+p1+p2+p3));
