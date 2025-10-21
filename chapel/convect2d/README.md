# Linear advection equation in 2-D

Finite volume method with WENO5 reconstruction, upwind flux and SSPK3 scheme

* `convect2d.chpl`: constant speed case, there is an MPI/PETSc version of this
* `convect2d_variable.chpl`: variable speed case

We declare `res` as an array of atomic reals in to order avoid race condition. For other strategies, see the `euler2d` example.
