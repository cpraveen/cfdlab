# Solve 2d Euler equations

* `isentropic_forall`: We compute all the fluxes and then update solution. This is needed to avoid race condition.
* `isentropic_coforall`: We dont store flux, but directly compute residual. To avoid race condition, we use coforall.

## forall version

If we compute `res` in forall version, it would look like this

```c
forall (i,j) in Dx
{
   ...
   res[i-1,j] += flux;
   res[i  ,j] -= flux;
}
```

There is danger that different threads may simultaneously modify the same element of `res`. Hence we compute and store all the fluxes, instead of computing `res`.

## Timing

Here are timing obtained on M1 macmini.

```shell
time ./isentropic
./isentropic  56.51s user 0.14s system 378% cpu 14.953 total

time ./isentropic_coforall
./isentropic_coforall  70.43s user 0.13s system 396% cpu 17.797 total
```

The coforall version is slower. The forall version requires more storage to store all the fluxes.
