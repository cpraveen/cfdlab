# 2d vorticity-streamfunction solver

The method is described in my notes `vte2d.tm`.

## Test Poisson solver

```shell
./test_poisson  -h  # shows available command line args
./test_poisson      # run with default options
```

## Test VTE solver

```shell
./vte -h                   # shows available command line args
./vte --n 128 --Re 1000
python plot.py             # needs pyvista
python line.py -Re 1000
```

Ghia et al. results taken from [here](https://github.com/CliMA/Oceananigans.jl/blob/main/validation/lid_driven_cavity/plot_lid_driven_cavity.py).

## Looping over red-black

An alternate way to do the loops is to use a conditional

```chapel
forall (i,j) in inner do
if (i+j)%2 == 0 // red points
{
    // apply SOR
}

forall (i,j) in inner do
if (i+j)%2 == 1 // black points
{
    // apply SOR
}
```

But my timing tests seem to show this is bit slower.
