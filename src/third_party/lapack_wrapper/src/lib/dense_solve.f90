module dense_solve
    implicit none
    integer, parameter:: dp=kind(0.d0)
    private
    public zgetrf_wrap, zgetrs_wrap, zgeev_wrap, zgeevx_wrap, zgebal_wrap, zgehrd_wrap
    interface
        ! TODO add intent(in) and intent(out) for wrapped subroutines
        subroutine zgetrs_wrap (mat_array, bx, permutation_indices)
            integer, parameter:: dp=kind(0.d0)
            complex(dp) ::mat_array(:,:), bx(:)
            integer :: permutation_indices(:)
        end subroutine

        subroutine zgetrf_wrap ( mat_array, permutation_indices )
            integer, parameter:: dp=kind(0.d0)
            complex(dp) ::mat_array(:,:)
            integer :: permutation_indices(:)
        end subroutine zgetrf_wrap

    end interface

contains
    !=======================================================================
    subroutine zgehrd_wrap(a, w, vr, scale, ilo, ihi)
        CHARACTER          job, compz, eigsrc, initv
        INTEGER            ihi, ilo, info, lda, ldz, n, lwork, shapes(2), ldvl, ifaill(1)
        integer, allocatable :: ifailr(:)
        REAL(dp)               scale(:)
        real(dp), allocatable:: rwork(:)
        COMPLEX(dp)            a(:,: ), w(:),z(1,1), vl(1,1), vr(:,:)
        logical, allocatable :: selector(:)
        complex(dp), allocatable :: tau(:), work(:)
        !-------------------------------------------------------------------
        shapes = shape(a)
        n = shapes(1)
        lda = shapes(2)
        job = 'B'

        allocate(tau(1:n-1))
        !workspace query
        lwork = -1
        allocate(work(1:2))
        call zgehd2(n, ilo, ihi, a, lda, tau, work, info)
!        call zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
!        lwork = work(1)
!        deallocate(work)
!        allocate(work(1:lwork))
!        call zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
        if (info /= 0) stop 'error in zeghrd_wrap'

    end subroutine
    !=======================================================================
    !=======================================================================
    subroutine zgebal_wrap(a, w, vr, scale, ilo, ihi)
        CHARACTER          job, compz, eigsrc, initv
        INTEGER            ihi, ilo, info, lda, ldz, n, lwork, shapes(2), ldvl, ifaill(1)
        integer, allocatable :: ifailr(:)
        REAL(dp)               scale(:)
        real(dp), allocatable:: rwork(:)
        COMPLEX(dp)            a(:,: ), w(:),z(1,1), vl(1,1), vr(:,:)
        logical, allocatable :: selector(:)
        complex(dp), allocatable :: tau(:), work(:)
        !-------------------------------------------------------------------
        shapes = shape(a)
        n = shapes(1)
        lda = shapes(2)
        job = 'B'
        call zgebal(job, n, a, lda, ilo, ihi, scale, info)
        if (info /= 0) stop 'error in zegbal_wrap'

!        allocate(tau(1:n-1))
!        !workspace query
!        lwork = -1
!        allocate(work(1:2))
!        call zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
!        lwork = work(1)
!        deallocate(work)
!        allocate(work(1:lwork))
!        call zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
!        if (info /= 0) stop 'error in zegbal_wrap'
!
!        job='E'
!        compz='N'
!        ldz = 1
!        deallocate(work)
!        allocate(work(1:n))
!        call zhseqr( job, compz, n, ilo, ihi, a, lda, w, z, ldz, work, lwork, info )
!        if (info /= 0) stop 'error in zegbal_wrap'
!
!        job='R'
!        eigsrc = 'Q'
!        initv = 'N'
!        allocate(selector(1:n))
!        selector(:) = .true.
!        allocate(rwork(1:n))
!        ldvl = 1
!        allocate(ifailr(1:n))
!        call zhsein( job, eigsrc, initv, selector, n, a, lda, w, vl, ldvl, vr, &
!                     lda, lda, lda, work, rwork, ifaill, ifailr, info )
!        if (info /= 0) stop 'error in zegbal_wrap'
!        deallocate(tau)
!        deallocate(work)
!        deallocate(rwork)
!        deallocate(select)
!        deallocate(ifailr)

    end subroutine
    !=======================================================================
    subroutine zgeev_wrap (a, w, vr)
        ! call Lapack zgeevx
        ! A * vr(:,j) = w(j) * vr(:,j)
        !
        ! ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
        ! eigenvalues and, optionally, the left and/or right eigenvectors.
        !
        ! The right eigenvector v(j) of A satisfies
        !                  A * v(j) = lambda(j) * v(j)
        ! where lambda(j) is its eigenvalue.
        ! The left eigenvector u(j) of A satisfies
        !               u(j)**H * A = lambda(j) * u(j)**H
        ! where u(j)**H denotes the conjugate transpose of u(j).
        !
        ! The computed eigenvectors are normalized to have Euclidean norm
        ! equal to 1 and largest component real.
        character          jobvl, jobvr
        integer            info, lda, ldvl, ldvr, lwork, n
        real(dp), allocatable :: rwork(:)
        complex(dp)        a( :, : ), vl( 1, 1 ), vr( :, : ), &
                w( : )
        complex(dp), allocatable :: work(:)
        !-------------------------------------------------------------------
        lwork = -1
        jobvl = 'N'
        jobvr = 'V'
        n = size(w)
        lda = n
        ldvr = n
        ldvl = 1
        allocate(rwork(1:2*n))
        !workspace query
        lwork = -1
        allocate(work(1:2))
        call zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        if (info /= 0) stop 'workspace query problem in zgeev_wrap'
        lwork = work(1)
        deallocate(work)
        allocate(work(1:lwork))
        call zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        if (info /= 0) stop 'eigen solver problem in zgeev_wrap'
    end subroutine
    !=======================================================================
    subroutine zgeevx_wrap (a, w, vr)
        ! call Lapack zgeev
        ! A * vr(:,j) = w(j) * vr(:,j)
        !
        ! ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
        ! eigenvalues and, optionally, the left and/or right eigenvectors.
        !
        ! The right eigenvector v(j) of A satisfies
        !                  A * v(j) = lambda(j) * v(j)
        ! where lambda(j) is its eigenvalue.
        ! The left eigenvector u(j) of A satisfies
        !               u(j)**H * A = lambda(j) * u(j)**H
        ! where u(j)**H denotes the conjugate transpose of u(j).
        !
        ! The computed eigenvectors are normalized to have Euclidean norm
        ! equal to 1 and largest component real.
        character          jobvl, jobvr, balanc, sense
        integer            info, lda, ldvl, ldvr, lwork, n, ilo, ihi
        real(dp) :: abnrm
        real(dp), allocatable :: rwork(:), scale(:), rconde(:), rcondv(:)
        complex(dp)        a( :, : ), vr( :, : ), &
                           w( : )
        complex(dp), allocatable :: work(:), vl( :, : )
        !-------------------------------------------------------------------
        lwork = -1
        balanc = 'B'
        sense = 'B'
        jobvl = 'V'
        jobvr = 'V'
        n = size(w)
        lda = n
        ldvr = n
        ldvl = n
        allocate(rwork(1:2*n))
        allocate(scale(1:n))
        allocate(rconde(1:n))
        allocate(rcondv(1:n))
        allocate(vl(1:n, 1:n))
        !workspace query
        lwork = -1
        allocate(work(1:2))
        call zgeevx( balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, &
            ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info)
        if (info /= 0) stop 'workspace query problem in zgeevx_wrap'
        lwork = work(1)
        deallocate(work)
        allocate(work(1:lwork))
        call zgeevx( balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, &
                ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info)
        if (info /= 0) stop 'eigen solver problem in zgeevx_wrap'
    end subroutine
    !=======================================================================
end module



