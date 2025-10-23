use VTK;
use Math;

proc fun1(x, y) {return sin(2*pi*x) * sin(2*pi*y);}
proc fun2(x, y) {return (sin(2*pi*x) * sin(2*pi*y), cos(2*pi*x) * cos(2*pi*y));}
proc fun3(x, y) {return [sin(2*pi*x) * sin(2*pi*y), cos(2*pi*x) * cos(2*pi*y)];}

proc main()
{
   param nx = 50, ny = 50;
   const dx = 1.0 / (nx-1);
   const dy = 1.0 / (ny-1);
   var x : [1..nx] real;
   var y : [1..ny] real;
   for i in x.domain do
      x[i] = (i-1)*dx;
   for j in y.domain do
      y[j] = (j-1)*dy;

   const D = {1..nx,1..ny};
   var u1 : [D] real;        // array of real
   var u2 : [D] 2*real;      // array of tuples
   var u3 : [D] [0..1] real; // array of array

   forall (i,j) in D
   {
      u1[i,j] = fun1(x[i],y[j]);
      u2[i,j] = fun2(x[i],y[j]);
      u3[i,j] = fun3(x[i],y[j]);
   }

   const names1 = ["sol"];
   write_vtk(x, y, 0.0, 0, names1, u1, "test1.vtk");

   const names2 = ["sol1", "sol2"];
   write_vtk(x, y, 0.0, 0, names2, u2, "test2.vtk");
   write_vtk(x, y, 0.0, 0, names2, u3, "test3.vtk");

   // Write in two steps
   // diff test1.vtk and test4.vtk
   const filename = "test4.vtk";
   write_vtk(x, y, 0.0, 0, filename);
   write_vtk(names1, u1, filename);
}
