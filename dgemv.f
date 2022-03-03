      program dgemv_analysis
        !use mkl_service
        implicit none
        include "mkl_lapack.fi"

        !! order of arrays
        integer, parameter :: ord1 = 2
        integer, parameter :: ord2 = 2
        integer, parameter :: ord1_large = 3
        integer, parameter :: ord2_large = 3

        !! regular arrays
        real*8, dimension(:), allocatable :: x, y      
        real*8, dimension(:,:), allocatable :: a
        real*8, dimension(:,:), allocatable :: a_large

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

        !! build extended matrix
        allocate(a_large(ord1_large, ord2_large))
        a_large = 0.0d0
        do j = ord2+1, ord2_large
                do i = 1, ord1_large
                        a_large(i,j) = 1.0d0
                enddo
        enddo
        do j = 1, ord2_large
                do i = ord1+1, ord1_large
                        a_large(i,j) = 1.0d0
                enddo
        enddo
        do i = 1, ord1
                a_large(i,i) = 1.0d0
        enddo
        write(*,*) 'A_large:'
        write(*,*) a_large


        !! standard dgemv call
        m = ord1
        n = ord2
        lda = ord1
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0
        write(*,*)
        write(*,*) 'Standard call to MKL'
        call dgemv('N', m, n, alpha, a, lda, x, incx, beta, y, incy)
        write(*,*) 'y = alpha*A*x + beta*y:'
        write(*,*) y

        

        !! call with A extended
        y = 0.0d0
        m = ord1
        n = ord2
        lda = ord2_large
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0
        write(*,*)
        write(*,*) 'Call with A extended, no spec'
        call dgemv('N', m, n, alpha, a_large, lda, x, incx, beta, 
     &                  y, incy)
        write(*,*) 'y = alpha*A_large*x + beta*y:'
        write(*,*) y


        !! I can also specify where to depart from...
        y = 0.0d0
        m = ord1
        n = ord2
        lda = ord2_large
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0
        write(*,*)
        write(*,*) 'Call with A extended, top left 2x2 block'
        call dgemv('N', m, n, alpha, a_large(1,1), lda, x, incx, beta, 
     &                  y, incy)
        write(*,*) 'y = alpha*A_large*x + beta*y:'
        write(*,*) y
        
        !! bottom left 2x2 block
        y = 0.0d0
        m = ord1
        n = ord2
        lda = ord2_large
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0
        write(*,*)
        write(*,*) 'Call with A extended, bottom left 2x2 block'
        call dgemv('N', m, n, alpha, a_large(2,1), lda, x, incx, beta, 
     &                  y, incy)
        write(*,*) 'y = alpha*A_large*x + beta*y:'
        write(*,*) y

        !! top right 2x2 block of A_large
        y = 0.0d0
        m = ord1
        n = ord2
        lda = ord2_large
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0
        write(*,*)
        write(*,*) 'Call with A extended, top right 2x2 block'
        call dgemv('N', m, n, alpha, a_large(1,2), lda, x, incx, beta, 
     &                  y, incy)
        write(*,*) 'y = alpha*A_large*x + beta*y:'
        write(*,*) y



        y = 0.0d0
        m = ord1
        n = ord2
        lda = ord2_large
        incx = 1
        incy = 1
        alpha = 1.0d0
        beta = 0.0d0
        write(*,*)
        write(*,*) 'Call with A extended, bottom right 2x2 block'
        call dgemv('N', m, n, alpha, a_large(2,2), lda, x, incx, beta, 
     &                  y, incy)
        write(*,*) 'y = alpha*A_large*x + beta*y:'
        write(*,*) y

        deallocate(x,y,a)
        deallocate(a_large)
      end program dgemv_analysis
