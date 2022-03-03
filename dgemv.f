      program dgemv_analysis
        !use mkl_service
        implicit none
        include "mkl_lapack.fi"

        !! order of arrays
        integer, parameter :: ord1 = 2
        integer, parameter :: ord2 = 2

        !! regular arrays
        real*8, dimension(:), allocatable :: x, y      
        real*8, dimension(:,:), allocatable :: a

        !! scalars
        real*8 :: alpha, beta
        
        !! integer that handle dimensions
        integer :: m, n, lda
        integer :: incx, incy

        !! helpers
        integer :: i, j

        !! build vectors
        allocate(x(ord2))
        allocate(y(ord1))
        y = 0.0d0
        do i = 1, ord2
                x(i) = 1.0d0
        enddo
        write(*,*) 'x:', x
        write(*,*) 'y:', y

        !! build matrix (0 outside diagonal)
        allocate(a(ord1,ord2))
        a = 0.d0
        !! this will be weird if ord1 /= ord2, but whatever...
        do i = 1, ord1
                a(i,i) = 1.0d0
        enddo
        write(*,*) 'A:', a


        !! time to dgmev these
        m = ord1
        n = ord2
        lda = ord1
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0

        call dgemv('N', m, n, alpha, a, lda, x, incx, beta, y, incy)

        write(*,*) 'y = alpha*A*x + beta*y:', y

        deallocate(x,y,a)
      end program dgemv_analysis
