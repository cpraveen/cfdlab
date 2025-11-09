# 1d Poisson equation

Solves the problem

```text
    -u''(x) = f(x)
```

with Dirichlet bc. There are two methods available SOR and v-cycle multigrid.

```shell
./main --N 128 --method mg  --niter 50
./main --N 128 --method sor --niter 500
```
