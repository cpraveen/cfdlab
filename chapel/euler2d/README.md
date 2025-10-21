# Solve 2d Euler equations

> This will work on shared memory computers, it may not work correctly in distributed systems.

* `euler_flux`: We compute and store all the fluxes and then update the solution. This is avoids race condition as each loop is parallel.
* `euler_res`: We dont store flux, but directly compute and store residual. To avoid race condition, we use forall only in one direction.

## euler_res version

If we compute `res` using a single forall, it would look like this

```c
forall (i,j) in Dx
{
   ...
   res[i-1,j] += flux;
   res[i  ,j] -= flux;
}
```

These loops are not independent; for each `(i,j)` we modify two values of `res`, at `(i-1,j)` and `(i,j)`. There is danger that different threads may simultaneously modify the same element of `res`. Hence we use forall in only one direction

```c
forall j in 1..ny
{
   for i in 1..nx+1
   {
      ...
      res[i-1,j] += flux;
      res[i  ,j] -= flux;
   }
}
```

In `i` direction we use a serial for loop so there is no race condition.

## Timing

Here are timing obtained on M1 macmini.

```shell
time ./euler_flux
./euler_flux  55.90s user 0.13s system 378% cpu 14.791 total

time ./euler_res
./euler_res  58.43s user 0.13s system 373% cpu 15.672 total
```

The times are somewhat similar, though `euler_flux` seems 1-2 seconds faster. The `euler_flux` version requires more storage to store all the fluxes and this storage increases in 3-D.
