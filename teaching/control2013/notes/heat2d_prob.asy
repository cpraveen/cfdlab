size(200,200);

pair p0 = (0,0);
pair p1 = (1,0);
pair p2 = (1,1);
pair p3 = (0,1);

draw(p0--p1--p2--p3--cycle,linewidth(2));

label("$z=0$", 0.5*(p0+p1), S);
label("$z=0$", 0.5*(p2+p3), N);
label("$z=u$", 0.5*(p1+p2), E);
label("$z_x=0$", 0.5*(p0+p3), W);
label("$\Omega$", 0.25*(p0+p1+p2+p3));
