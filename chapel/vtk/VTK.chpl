module VTK
{

use IO;

// This function can be removed once it is added to Chapel, once this
// https://github.com/chapel-lang/chapel/pull/27936 is released.
pragma "no copy return"
pragma "return not owned"
proc type _array._dom do return chpl__domainFromArrayRuntimeType(this);

// Write 2d rectilinear grid format
// Write only grid
proc write_vtk(x : [?Dx] real, 
               y : [?Dy] real,
               time : real,
               cycle : int,
               filename : string) where Dx.rank == 1 &&
                                        Dy.rank == 1
{
   const nx = x.size;
   const ny = y.size;

   var fw = try! open(filename, ioMode.cw).writer();

   try! fw.writeln("# vtk DataFile Version 3.0");
   try! fw.writeln("Cartesian grid");
   try! fw.writeln("ASCII");
   try! fw.writeln("DATASET RECTILINEAR_GRID");
   try! fw.writeln("FIELD FieldData 2");
   try! fw.writeln("TIME 1 1 double");
   try! fw.writeln(time);
   try! fw.writeln("CYCLE 1 1 int");
   try! fw.writeln(cycle);
   try! fw.writef ("DIMENSIONS %i %i 1\n", nx, ny);

   try! fw.writef ("X_COORDINATES  %i  float\n", nx);
   for v in x do
      try! fw.writeln(v);

   try! fw.writef ("Y_COORDINATES  %i  float\n", ny);
   for v in y do
      try! fw.writeln(v);

   try! fw.writef ("Z_COORDINATES  1  float\n");
   try! fw.writeln(0.0);

   try! fw.writef("POINT_DATA  %i\n", nx * ny);
   try! fw.close();
}

// Write 2d rectilinear grid format
// Write only solution, assumes grid is already written to same filename
proc write_vtk(names : [] string,
               const ref u : [?D] ?T,
               filename : string) where D.rank == 2
{
   // Append to this file
   var fw = try! open(filename, ioMode.a).writer();

   if isReal(T) // u is array of reals
   {
      try! fw.writef("SCALARS %s float\n", names[0]);
      try! fw.writef("LOOKUP_TABLE default\n");
      for j in D.dim(1) do
         for i in D.dim(0) do
            try! fw.writeln(u[i,j]);
   }
   else if isTuple(T) // u is array of tuples
   {
      const nc = T.size; //u[1,1].size;
      assert(nc == names.size, "Inconsistent size of names");
      for (c,name) in zip(0..nc-1,names)
      {
         try! fw.writef("SCALARS %s float\n", name);
         try! fw.writef("LOOKUP_TABLE default\n");
         for j in D.dim(1) do
            for i in D.dim(0) do
               try! fw.writeln(u[i,j][c]);
      }
   }
   else // u is array of arrays
   {
      const dom1d = T.domain;
      assert(dom1d.size == names.size, "Inconsistent size of names");
      for (c,name) in zip(dom1d,names)
      {
         try! fw.writef("SCALARS %s float\n", name);
         try! fw.writef("LOOKUP_TABLE default\n");
         for j in D.dim(1) do
            for i in D.dim(0) do
               try! fw.writeln(u[i,j][c]);
      }
   }

   try! fw.close();
}

// Write 2d rectilinear grid format
proc write_vtk(x : [?Dx] real, 
               y : [?Dy] real,
               time : real,
               cycle : int,
               names : [] string,
               const ref u : [?D] ?T,
               filename : string) where D.rank == 2 && 
                                        Dx.rank == 1 &&
                                        Dy.rank == 1
{
   const nx = D.dim(0).size;
   const ny = D.dim(1).size;
   assert(nx == x.size, "Size of x is wrong");
   assert(ny == y.size, "Size of y is wrong");

   write_vtk(x, y, time, cycle, filename);
   write_vtk(names, u, filename);
}

}
