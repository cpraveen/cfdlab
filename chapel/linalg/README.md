# Linear Algebra

The makefile assumes that `BLAS_DIR` and `LAPACK_DIR` are set. I am also using gfortran from homebrew.

```shell
make test1
DYLD_LIBRARY_PATH=$LAPACK_DIR/lib ./test1
```
