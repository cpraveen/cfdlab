size(150,150);

pair p0 = (0,0);
pair p1 = (1, 1);
pair p2 = (-0.5,1.5);

draw(p0--p1--p2--cycle, linewidth(2));

label("$1$", p0, S);
label("$2$", p1, E);
label("$3$", p2, NW);
label("$K$", (p0+p1+p2)/3);
