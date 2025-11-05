# 2d vorticity-streamfunction solver

The method is described in my notes `vte2d.tm`. There is a Julia version [here](https://github.com/cpraveen/numpde/tree/master/vte2d) but it may not be as good as the Chapel version.

## Test Poisson solver

```shell
./test_poisson  -h    # shows available command line args
./test_poisson        # run with default options
visit -o poisson.vtk  # see the solution
```

## Test VTE solver

See available options

```shell
./vte -h
```

Run the code

```shell
./vte --n 128 --Re 1000 > log.txt &
```

See the solution

```shell
python stream.py           # plot streamlines, needs pyvista
python vel.py --Re 1000    # compare velocity with Ghia
```

[Ghia et al.](https://doi.org/10.1016/0021-9991(82)90058-4) results were taken from [here](https://github.com/CliMA/Oceananigans.jl/blob/main/validation/lid_driven_cavity/plot_lid_driven_cavity.py). The figures show the streamlines for Re=1000 and comparison of velocity profile along center line of the domain.

<p align="center">
<img src="output/stream_Re1000.svg" width=512>
<img src="output/vel_Re1000.svg">
</p>

## Looping over red-black

An alternate way to do the loops over red and black points in the Poisson solver is,

```chapel
forall i in inner.dim(0) do
forall j in inner.dim(1) by 2 align i
{
    // apply SOR
}

forall i in inner.dim(0) do
forall j in inner.dim(1) by 2 align i+1
{
    // apply SOR
}
```

## Running on multiple locales

I have not tested this but the code should work on multiple locales, except for the solution output functions in the VTK module which are serial.

If you want to run with a single locale, then there is no need to use `stencilDist` as we dont use any ghost points in this code. You can declare a normal domain instead

```chapel
const D = {1..n, 1..n};
```

and remove calls to `updateFluff()` function.
