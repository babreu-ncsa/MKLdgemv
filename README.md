# MKLdgemv

Compilation on Bridges-2:

`ml purge`

`ml openmpi/3.1.5-nvhpc21.7`

`ml mkl/2020.4.304`

`ml cuda/11.1.1`

CUDA is not really necessary, but we have more ambitious goals here...

Regardless, a typical MKL call is:

```fortran
call dgemv(m, n, alpha, a, lda, x, incx, beta, y, incy)
```

and the point is that the array `a` can be much larger than the matrix A that we use to perform `y -> alpha*A*x + beta*y`. When that's the case, you can specify the initial element `a(start_row, start_column)` from which `m`, `n` are going to be counted. The code in this repo does that for a simple 3x3 array `a` that multiplies 2-entries vectors using all possible blocks from `a`. Results can be verified quite easily.
