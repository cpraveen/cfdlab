# 1d active flux

Compile the code

```shell
make
```

## Burgers

Vanilla jacobian splitting, compare to Duan, Fig., 3 (top left)

```shell
./burgers --n 200 --ic 2  --Tf 0.5 --diff 1 --jac 0
gnuplot plot.gnu
open sol.pdf
```

Modified jacobian splitting

```shell
./burgers --n 200 --ic 2  --Tf 0.5 --diff 1 --jac 1
gnuplot plot.gnu
open sol.pdf
```

Vanilla flux splitting

```shell
./burgers --n 200 --ic 2  --Tf 0.5 --diff 3
gnuplot plot.gnu
open sol.pdf
```
