# Examples of writing vtk files

## Cartesian and structured grids
```
c++ vtk_struct.cc
./a.out
```
This outputs a Cartesian grid and a structured grid in vtk format.

## Gmsh unstructured grid
```
gmsh -2 -format msh2 cylinder.geo
c++ gmsh.cc
./a.out
```
This outputs a vtk file in unstructured format.
