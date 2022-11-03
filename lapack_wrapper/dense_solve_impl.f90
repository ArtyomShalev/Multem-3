!=======================================================================
! Split to call from f77 codes compiled with f90
subroutine zgetrs_wrap (mat_array, bx, permutation_indices)
    ! call Lapack zgetrs
    !
    ! ZGETRS solves a system of linear equations
    !    A * X = B,  A**T * X = B,  or  A**H * X = B
    ! with a general N-by-N matrix A using the LU factorization computed
    ! by ZGETRF.
    integer, parameter:: dp=kind(0.d0)
    complex(dp) ::mat_array(:,:), bx(:)
    integer :: permutation_indices(:)
    integer :: info, lda, ldb, m, n!, array_rank
    integer :: array_shape(2)

    !allocate(array_shape(1:rank(mat_array)))
    array_shape = shape(mat_array)
    m = array_shape(1)
    n = array_shape(2)
    lda = size(permutation_indices)
    ldb = size(bx)
    info = 0
    if( m < 0 ) then
        info = -1
    else if( n < 0 ) then
        info = -2
    else if( m /= n ) then
        info = -3
    else if( lda < max( 1, m ) ) then
        info = -4
    else if( lda /= ldb ) then
        info = -5

    end if
    if( info /= 0 ) then
        write(*,*) "Erorr in zgetrf_wrap", -info
        stop 1
        return
    end if

    ! Call Lapack worker
    call zgetrs( 'no transpose', n, 1, mat_array, &
            & lda, permutation_indices, bx, ldb, info )

    if( info /= 0 ) then
        write(*,*) "Erorr in zgetrs_wrap", -info
        stop 1
        return
    end if
end subroutine
!=======================================================================
subroutine zgetrf_wrap ( mat_array, permutation_indices )
    ! call Lapack zgetrf
    !
    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    !
    ! The factorization has the form
    !    A = P * L * U
    ! where P is a permutation matrix, L is lower triangular with unit
    ! diagonal elements (lower trapezoidal if m > n), and U is upper
    ! triangular (upper trapezoidal if m < n).
    !
    ! This is the right-looking Level 3 BLAS version of the algorithm.
    integer, parameter:: dp=kind(0.d0)
    complex(dp) ::mat_array(:,:)
    integer :: permutation_indices(:)
    integer :: info, lda, m, n!, array_rank
    integer :: array_shape(2)

    !allocate(array_shape(1:rank(mat_array)))
    array_shape = shape(mat_array)
    m = array_shape(1)
    n = array_shape(2)
    lda = size(permutation_indices)
    info = 0
    if( m < 0 ) then
        info = -1
    else if( n < 0 ) then
        info = -2
    else if( m /= n ) then
        info = -3
    else if( lda < max( 1, m ) ) then
        info = -4
    end if
    if( info /= 0 ) then
        write(*,*) "Erorr in zgetrf_wrap", -info
        stop 1
        return
    end if

    ! Call Lapack worker
    call zgetrf( m, n, mat_array, lda, permutation_indices, info )

    if( info /= 0 ) then
        write(*,*) "Erorr in zgetrf_wrap", -info
        stop 1
    end if
    !stop
    return
end subroutine zgetrf_wrap



