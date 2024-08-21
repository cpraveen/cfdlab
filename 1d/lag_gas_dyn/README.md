# 1d gas dynamics in lagrangian coordinates

Solves using HLL, first order in space and SSPRK3 in time

* `lag_gas_con.f`: conservation form
* `lag_gas_noncon1.f`: non-conservation form, does not conserve total energy
* `lag_gas_noncon2.f`: non-conservation form, conserves total energy

See 

Abgrall and Karni, JCP, 2010

Run as 

```
make
./con
./noncon1
./noncon2
```

and plot in gnuplot using `plot.gnu` file.
