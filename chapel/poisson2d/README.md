# Solve 2d Poisson equation

The code has red-black SOR and V-cycle multigrid. It requires the VTK writer, see the `makefile` to set the path.

Compile the code

```shell
make
```

Run the code with multigrid (default)

```shell
./main --nx 128 --ny 128 --levels 7
python ./plot.py
```

Increase grid size

```shell
./main --nx 256 --ny 256 --levels 8
python ./plot.py
```

For best multigrid convergence rate, choose `nx = ny = 2^levels`.

Try with SOR

```shell
./main --nx 256 --ny 256 --method sor
python ./plot.py
```
