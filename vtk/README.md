# Examples of writing vtk files

## Cartesian and structured grids
```
c++ vtk_struct.cc
./a.out
```
This outputs a Cartesian grid and a structured grid in vtk format. A file in Tecplot format is also written. Open them in Visit.

## Time dependent solution
```
c++ vtk_anim.cc
./a.out
```
generate a sequence of files with time dependent solution. View them in Visit
```
visit -o sol*.vtk
```
and plot pseudocolor or contour plots.

## Gmsh unstructured grid
```
gmsh -2 -format msh2 cylinder.geo
c++ gmsh.cc
./a.out
```
This outputs a vtk file in unstructured format. Open it in Visit.

## Reading vtk files

The two notebooks show how to read vtk files using vtk Python API and meshio. These were taken from https://github.com/pnavaro/plot_vtk_with_matplotlib
