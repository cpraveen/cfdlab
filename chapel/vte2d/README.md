# 2d vorticity-streamfunction solver

The method is described in my notes `vte2d.tm`.

## Test Poisson solver

```shell
./test_poisson  # shows available command line args
./test_poisson
```

## Test VTE solver

```shell
./vte -h  # shows available command line args
./vte
python plot.py
```

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
