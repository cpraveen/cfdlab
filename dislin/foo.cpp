#include <iostream>
#include <cmath>
#include "discpp.h"

int main ()
{ int n = 100, i, ic;
  double fpi = 3.1415926 / 180.0, step, x;
  double xray[100], y1ray[100], y2ray[100];
  Dislin g;

  step = 360. / (n - 1);

  for (i = 0; i < n; i++)
  { xray[i] = i * step;
    x = xray[i] * fpi;
    y1ray[i] = sin (x);
    y2ray[i] = cos (x);
  }

  g.metafl ("xwin");
  g.scrmod ("revers");
  g.disini ();
  g.pagera ();
  g.complx ();
  g.axspos (450, 1800);
  g.axslen (2200, 1200);

  g.name   ("X-axis", "x");
  g.name   ("Y-axis", "y");

  g.labdig (-1, "x");
  g.ticks  (9, "x");
  g.ticks  (10, "y");

  g.titlin ("Demonstration of CURVE", 1);
  g.titlin ("SIN(X), COS(X)", 3);

  ic=g.intrgb (0.95,0.95,0.95);
  g.axsbgd (ic);

  g.graf   (0.0, 360.0, 0.0, 90.0, -1.0, 1.0, -1.0, 0.5);
  g.setrgb (0.7, 0.7, 0.7);
  g.grid   (1, 1);

  g.color  ("fore");
  g.height (50);
  g.title  ();

  g.color  ("red");
  g.curve  (xray, y1ray, n);
  g.color  ("green");
  g.curve  (xray, y2ray, n);
  g.disfin ();
  return 0;
}
