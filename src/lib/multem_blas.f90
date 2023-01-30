module multem_blas
    use dense_solve
    implicit none
    private
    integer, parameter, private:: dp=kind(0.0D0)
    complex(dp), parameter, private :: ci    = (0.0_dp, 1.0_dp)
    public comlr2, comhes, cbabk2, cnaa, zsu, zge

contains
    !=======================================================================
    complex(dp) function cmplx_dp(re, im)
        real(dp), intent(in) :: re, im
        cmplx_dp = cmplx(re,im, kind=dp)
    end function cmplx_dp
    !=======================================================================
    subroutine zge(a, int, n, nc, emach)

        !     ------------------------------------------------------------------
        !     zge is a standard subroutine to perform gaussian elimination on
        !     a nc*nc matrix 'a' prior  to inversion, details stored in 'int'
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer n, nc
        real(dp) emach
        !
        ! ..  array arguments  ..
        !
        integer    int(nc)
        complex(dp) a(nc, nc)
        !
        ! ..  local scalars  ..
        !
        integer    i, ii, in, j, k
        complex(dp) yr, dum
        !     ------------------------------------------------------------------
        !
        do ii = 2, n
            i = ii - 1
            yr = a(i, i)
            in = i
            do j = ii, n
                if(abs(yr) - abs(a(j, i)) >= 0) cycle
                yr = a(j, i)
                in = j
            end do
            int(i) = in
            if(in - i /= 0) then
                do j = i, n
                    dum = a(i, j)
                    a(i, j) = a(in, j)
                    a(in, j) = dum
                end do
            end if
            if(abs(yr) - emach > 0) then
                do j = ii, n
                    if(abs(a(j, i)) - emach<=0) cycle
                    a(j, i) = a(j, i) / yr
                    do k = ii, n
                        a(j, k) = a(j, k) - a(i, k) * a(j, i)
                    end do
                end do
            end if
        end do
        return
    end subroutine
    !=======================================================================
    subroutine zsu(a, int, x, n, nc, emach)

        !     ------------------------------------------------------------------
        !     zsu  is  a standard back-substitution  subroutine  using the
        !     output of zge to calculate  a-inverse times x, returned in x
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer n, nc
        real(dp) emach
        !
        ! ..  array arguments  ..
        !
        integer    int(nc)
        complex(dp) a(nc, nc), x(nc)
        !
        ! ..  local scalars  ..
        !
        integer    i, ii, in, j, ij
        complex(dp) dum
        !     ------------------------------------------------------------------
        !
        do ii = 2, n
            i = ii - 1
            if(int(i) - i /= 0) then
                in = int(i)
                dum = x(in)
                x(in) = x(i)
                x(i) = dum
            end if
            do j = ii, n
                if(abs(a(j, i)) - emach >0) x(j) = x(j) - a(j, i) * x(i)
            end do
        end do
        do ii = 1, n
            i = n - ii + 1
            ij = i + 1
            if(i - n /= 0) then
                do j = ij, n
                    x(i) = x(i) - a(i, j) * x(j)
                end do
            end if
            if(abs(a(i, i)) - emach * 1.0d-7 < 0) then
                a(i, i) = emach * 1.0d-7 * (1.0_dp, 1.0_dp)
            else
                x(i) = x(i) / a(i, i)
            endif

        end do
        return
    end subroutine
    !=======================================================================
    subroutine cnaa(ndim, n, ar, ai, evr, evi, vecr, veci, ierr)

        !     ------------------------------------------------------------------
        !     'eispack'  is a  collection  of codes for  solving  the algebraic
        !     eigenvalue  problem.  the original  algol  codes were  written by
        !     j. h. wilkinson, et.al., and subsequently  translated to  fortran
        !     and tested at argonne national laboratory.
        !
        !     this   subroutine  computes  all  eigenvalues  and  corresponding
        !     eigenvectors  of  an  arbitrary   complex  matrix.  the matrix is
        !     balanced by exact norm  reducing  similarity  transformations and
        !     then  is  reduced  to  complex  hessenberg   form  by  stabilized
        !     elementary similarity transformations. a modified lr algorithm is
        !     used to compute the eigenvalues of the hessenberg matrix.
        !
        !       on input--->
        !          ndim     must be the row dimension of the arrays ar,ai,vecr,
        !                   and veci in the calling program dimension statement
        !          n        is the order of the matrix. n must not exceed ndim.
        !                   n*ndim  must not exceed 22500=150*150=53744(octal).
        !                   n must not exceed 150.  n may be 1.
        !          ar,ai    arrays with  exactly  ndim  rows  and  at  least  n
        !                   columns.  the leading n by n subarrays must contain
        !                   the real and  imaginary  parts  respectively of the
        !                   arbitrary complex matrix whose eigensystem is to be
        !                   computed.
        !
        !        on output--->
        !          evr,evi    contain the real and imaginary parts respectively
        !                     of the computed eigenvalues.  the eigenvalues are
        !                     not ordered in any way.
        !          vecr,veci  contain in the leading n by n  subarrays the real
        !                     and imaginary parts respectively  of the computed
        !                     eigenvectors.  the j-th columns  of vecr and veci
        !                     contain the  eigenvector  associated  with evr(j)
        !                     and  evi(j).  the eigenvectors are not normalized
        !                     in any way.
        !          ierr       is a status code.
        !                   --normal code.
        !                     0 means the lr iterations converged.
        !                   --abnormal codes.
        !                     j means the j-th eigenvalue has not been found in
        !                     30 iterations. the first j-1 elements of evr  and
        !                     evi contain those eigenvalues  already  found. no
        !                     eigenvectors are computed.
        !                    -1 means the input values of n, ndim are too large
        !                     or inconsistent.
        !          ar,ai      are destroyed.
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer ierr, n, ndim
        !
        ! ..  array arguments  ..
        !
        real(dp)ar(ndim, n), ai(ndim, n), evr(n), evi(n)
        complex(dp) :: a(ndim, n), w(n), vr(ndim, n)
        real(dp)vecr(ndim, n), veci(ndim, n)
        !
        ! ..  local scalars  ..
        !
        integer   i, ierrpi, igh, low, nmierr
        !
        ! ..  local arrays  ..
        !
        integer int(270)
        real(dp)scale(270)
        !     ------------------------------------------------------------------
        !
        if(ndim<n .or. n<1) go to 10
        if(n * ndim > 72900) go to 10
        a = ar + ci * ai
        call zgebal_wrap(a, w, vr, scale, low, igh)
        ar = dble(a)
        ai = aimag(a)
        !       call zgehrd_wrap(a, w, vr, scale, low, igh)
        !     ierr = 0
        !     ar = dble(a)
        !     ai = aimag(a)
        !     evr = dble(w)
        !     evi = aimag(w)
        !     vecr = dble(vr)
        !     veci = dble(vr)

        call comhes(ndim, n, low, igh, ar, ai, int)
        call comlr2(ndim, n, low, igh, int, ar, ai, evr, evi, vecr, veci, ierr)
        if(ierr==0) go to 2
        !     call errchk(54,54hin cnaa  , some eigenvalue not found in 30 itera
        !    1tions.)
        if(ierr==n) go to 20
        nmierr = n - ierr
        do i = 1, nmierr
            ierrpi = ierr + i
            evr(i) = evr(ierrpi)
            evi(i) = evi(ierrpi)
        end do
        go to 20
        2    call cbabk2(ndim, n, low, igh, scale, n, vecr, veci)
        go to 20
        10    write(*, *)"in cnaa input dim in error or matrix is too big."

        ierr = -1
        20    if(ierr > 0) ierr = n - ierr + 1
        return
    end subroutine
    !=======================================================================
    subroutine cbabk2(nm, n, low, igh, scale, m, zr, zi)

        !     ------------------------------------------------------------------
        !     this subroutine forms the eigenvectors of a complex general
        !     matrix by back transforming those of the corresponding
        !     balanced matrix determined by  cbal.
        !
        !     on input--->
        !        nm must be set to the row dimension of two-dimensional
        !          array parameters as declared in the calling program
        !          dimension statement,
        !
        !        n is the order of the matrix,
        !
        !        low and igh are integers determined by  cbal,
        !
        !        scale contains information determining the permutations
        !          and scaling factors used by  cbal,
        !
        !        m is the number of eigenvectors to be back transformed,
        !
        !        zr and zi contain the real and imaginary parts,
        !          respectively, of the eigenvectors to be
        !          back transformed in their first m columns.
        !
        !     on output--->
        !        zr and zi contain the real and imaginary parts,
        !          respectively, of the transformed eigenvectors
        !          in their first m columns.
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer nm, n, low, igh, m
        !
        ! ..  array arguments  ..
        !
        real(dp) scale(n), zr(nm, m), zi(nm, m)
        !
        ! ..  local scalars  ..
        !
        integer i, j, k, ii
        real(dp)s
        !     ------------------------------------------------------------------
        !
        if (m==0) return
        if (igh/=low) then
            do i = low, igh
                s = scale(i)
                !     ********** left hand eigenvectors are back transformed
                !                if the foregoing statement is replaced by
                !                s=1.0/scale(i). **********
                do j = 1, m
                    zr(i, j) = zr(i, j) * s
                    zi(i, j) = zi(i, j) * s
                end do
            end do
        end if
        !     ********** for i=low-1 step -1 until 1,
        !                igh+1 step 1 until n do -- **********
        do ii = 1, n
            i = ii
            if (i >= low .and. i <= igh) cycle
            if (i < low) i = low - ii
            k = scale(i)
            if (k == i) cycle
            !
            do j = 1, m
                s = zr(i, j)
                zr(i, j) = zr(k, j)
                zr(k, j) = s
                s = zi(i, j)
                zi(i, j) = zi(k, j)
                zi(k, j) = s
            end do
        end do
        !
        return
    end subroutine
    !=======================================================================
    subroutine comhes(nm, n, low, igh, ar, ai, int)

        !     ------------------------------------------------------------------
        !     given a  complex  general  matrix, this  subroutine  reduces  a
        !     submatrix situated in rows and columns low through igh to upper
        !     hessenberg form by stabilized elementary similarity transforms.
        !
        !     on input--->
        !        nm       must be set to the row dimension of two-dimensional
        !                 array parameters as declared in the calling program
        !                 dimension statement
        !        n        is the order of the matrix
        !        low,igh  are integers determined by the balancing subroutine
        !                 cbal. if  cbal  has not been used, set low=1, igh=n
        !        ar,ai    contain the real and imaginary parts, respectively,
        !                 of the complex input matrix.
        !
        !     on output--->
        !        ar,ai    contain the real and imaginary parts, respectively,
        !                 of the hessenberg matrix.the multipliers which were
        !                 used in the  reduction  are stored in the remaining
        !                 triangles under the hessenberg matrix,
        !        int      contains information on the rows and columns inter-
        !                 changed in the reduction. only elements low through
        !                 igh are used.
        !
        !     arithmetic is real except for the replacement of the algol
        !     procedure cdiv by complex division using subroutine cmplx.
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer nm, n, low, igh
        !
        ! ..  array arguments  ..
        !
        integer int(igh)
        real(dp)ar(nm, n), ai(nm, n)
        !
        ! ..  local scalars  ..
        !
        integer    i, j, m, la, kp1, mm1, mp1
        real(dp)   xr, xi, yr, yi
        complex(dp) z3
        !     ------------------------------------------------------------------
        !
        la = igh - 1
        kp1 = low + 1
        if (la < kp1) go to 200
        !
        do m = kp1, la
            mm1 = m - 1
            xr = 0.0_dp
            xi = 0.0_dp
            i = m
            !
            do j = m, igh
                if (abs(ar(j, mm1)) + abs(ai(j, mm1))&
                        <= abs(xr) + abs(xi)) go to 100
                xr = ar(j, mm1)
                xi = ai(j, mm1)
                i = j
                100    continue
            end do
            !
            int(m) = i
            if (i == m) go to 130
            !     ********** interchange rows and columns of ar and ai **********
            do j = mm1, n
                yr = ar(i, j)
                ar(i, j) = ar(m, j)
                ar(m, j) = yr
                yi = ai(i, j)
                ai(i, j) = ai(m, j)
                ai(m, j) = yi
            end do
            !
            do j = 1, igh
                yr = ar(j, i)
                ar(j, i) = ar(j, m)
                ar(j, m) = yr
                yi = ai(j, i)
                ai(j, i) = ai(j, m)
                ai(j, m) = yi
            end do
            !     ********** end interchange **********
            130      if (xr == 0.0_dp .and. xi == 0.0_dp) go to 180
            mp1 = m + 1
            !
            do i = mp1, igh
                yr = ar(i, mm1)
                yi = ai(i, mm1)
                if (yr == 0.0_dp .and. yi == 0.0_dp) go to 160
                z3 = cmplx_dp(yr, yi) / cmplx_dp(xr, xi)
                yr = dble(z3)
                yi = aimag (z3)
                ar(i, mm1) = yr
                ai(i, mm1) = yi
                !
                do j = m, n
                    ar(i, j) = ar(i, j) - yr * ar(m, j) + yi * ai(m, j)
                    ai(i, j) = ai(i, j) - yr * ai(m, j) - yi * ar(m, j)
                end do
                !
                do j = 1, igh
                    ar(j, m) = ar(j, m) + yr * ar(j, i) - yi * ai(j, i)
                    ai(j, m) = ai(j, m) + yr * ai(j, i) + yi * ar(j, i)
                end do
                !
                160    continue
            end do
            !
            180 continue
        end do
        !
        200 return
    end subroutine
    !=======================================================================
    subroutine comlr2(nm, n, low, igh, int, hr, hi, wr, wi, zr, zi, ierr)

        !     ------------------------------------------------------------------
        !     this subroutine finds  the  eigenvalues and  eigenvectors  of  a
        !     complex upper hessenberg  matrix by the modified  lr method. the
        !     eigenvectors  of a complex  general matrix  can also be found if
        !     comhes has been used to reduce this general matrix to hessenberg
        !     form.
        !
        !     on input--->
        !        nm      must  be set to the row dimension  of two-dimensional
        !                array  parameters as  declared in the calling program
        !                dimension statement
        !        n       is the order of the matrix
        !        low,igh are integers determined by the  balancing  subroutine
        !                cbal.  if  cbal  has not been used,  set low=1, igh=n
        !        int     contains information on the rows and  columns  inter-
        !                changed in the reduction by comhes,if performed. only
        !                elements low through igh are used.if the eigenvectors
        !                of the hessenberg matrix are desired,set int(j)=j for
        !                these elements
        !        hr,hi   contain the real and imaginary parts, respectively,of
        !                the complex upper hessenberg matrix. their lower tri-
        !                angles  below the subdiagonal contain the multipliers
        !                which   were  used  in  the  reduction by  comhes, if
        !                performed.  if  the  eigenvectors  of  the hessenberg
        !                matrix are desired,these elements must be set to zero
        !
        !      on output--->
        !                the   upper hessenberg portions of hr and hi have been
        !                destroyed, but  the location hr(1,1) contains the norm
        !                of the triangularized matrix,
        !        wr,wi   contain the real and imaginary parts, respectively, of
        !                the   eigenvalues.  if  an  error  exit  is  made, the
        !                eigenvalues should be correct for indices ierr+1,...,n
        !        zr,zi   contain the real and imaginary parts, respectively, of
        !                the eigenvectors.the eigenvectors are unnormalized. if
        !                an error exit is  made, none of the  eigenvectors  has
        !                been found
        !        ierr    is set to  zero for normal return,
        !          j     if the j-th  eigenvalue has not been  determined after
        !                30 iterations.
        !
        !     arithmetic  is  real  except  for the  replacement  of  the algol
        !     procedure cdiv by  complex division and  use of  the  subroutines
        !     csqrt and cmplx in computing complex square roots.
        !     ------------------------------------------------------------------
        !
        ! ..  scalar arguments  ..
        !
        integer nm, n, low, igh, ierr
        !
        ! ..  array arguments  ..
        !
        integer int(igh)
        real(dp) hr(nm, n), hi(nm, n), wr(n), wi(n), zr(nm, n), zi(nm, n)
        !
        ! ..  local scalars  ..
        !
        integer    i, j, k, l, m, en, ii, jj, ll, mm, nn, im1, ip1, its, mp1, enm1, iend
        real(dp)   si, sr, ti, tr, xi, xr, yi, yr, zzi, zzr, norm, machep
        complex(dp) z3
        !     ------------------------------------------------------------------
        !
        !     ********** machep is a machine dependent parameter specifying
        !                the relative precision of floating point arithmetic.
        !
        machep = 2.0_dp**(-47)  !TODO: convergence constant in comlr2()
        !
        ierr = 0
        !     ********** initialize eigenvector matrix **********
        do i = 1, n
            !
            do j = 1, n
                zr(i, j) = 0.0_dp
                zi(i, j) = 0.0_dp
                if (i == j) zr(i, j) = 1.0_dp
            end do
        end do
        !     ********** form the matrix of accumulated transformations
        !                from the information left by comhes **********
        iend = igh - low - 1
        if (iend <= 0) go to 180
        !     ********** for i=igh-1 step -1 until low+1 do -- **********
        do ii = 1, iend
            i = igh - ii
            ip1 = i + 1
            !
            do k = ip1, igh
                zr(k, i) = hr(k, i - 1)
                zi(k, i) = hi(k, i - 1)
            end do
            !
            j = int(i)
            if (i == j) go to 160
            !
            do k = i, igh
                zr(i, k) = zr(j, k)
                zi(i, k) = zi(j, k)
                zr(j, k) = 0.0_dp
                zi(j, k) = 0.0_dp
            end do
            !
            zr(j, i) = 1.0_dp
            160 continue
        end do
        !     ********** store roots isolated by cbal **********
        180 do i = 1, n
            if (i >= low .and. i <= igh) go to 200
            wr(i) = hr(i, i)
            wi(i) = hi(i, i)
            200 continue
        end do
        !
        en = igh
        tr = 0.0_dp
        ti = 0.0_dp
        !     ********** search for next eigenvalue **********
        220 if (en < low) go to 680
        its = 0
        enm1 = en - 1
        !     ********** look for single small sub-diagonal element
        !                for l=en step -1 until low do -- **********
        240 do ll = low, en
            l = en + low - ll
            if (l == low) go to 300
            if (abs(hr(l, l - 1)) + abs(hi(l, l - 1)) <=&
                    machep * (abs(hr(l - 1, l - 1)) + abs(hi(l - 1, l - 1))&
                            + abs(hr(l, l)) + abs(hi(l, l)))) go to 300
        end do
        !     ********** form shift **********
        300 if (l == en) go to 660
        if (its == 30) go to 1000
        if (its == 10 .or. its == 20) go to 320
        sr = hr(en, en)
        si = hi(en, en)
        xr = hr(enm1, en) * hr(en, enm1) - hi(enm1, en) * hi(en, enm1)
        xi = hr(enm1, en) * hi(en, enm1) + hi(enm1, en) * hr(en, enm1)
        if (xr == 0.0_dp .and. xi == 0.0_dp) go to 340
        yr = (hr(enm1, enm1) - sr) / 2.0_dp
        yi = (hi(enm1, enm1) - si) / 2.0_dp
        z3 = sqrt(cmplx_dp(yr**2 - yi**2 + xr, 2.0_dp * yr * yi + xi))
        zzr = dble(z3)
        zzi = aimag(z3)
        if (yr * zzr + yi * zzi >= 0.0_dp) go to 310
        zzr = -zzr
        zzi = -zzi
        310 z3 = cmplx_dp(xr, xi) / cmplx_dp(yr + zzr, yi + zzi)
        sr = sr - dble(z3)
        si = si - aimag(z3)
        go to 340
        !     ********** form exceptional shift **********
        320 sr = abs(hr(en, enm1)) + abs(hr(enm1, en - 2))
        si = abs(hi(en, enm1)) + abs(hi(enm1, en - 2))
        !
        340 do i = low, en
            hr(i, i) = hr(i, i) - sr
            hi(i, i) = hi(i, i) - si
        end do
        !
        tr = tr + sr
        ti = ti + si
        its = its + 1
        !     ********** look for two consecutive small
        !                sub-diagonal elements **********
        xr = abs(hr(enm1, enm1)) + abs(hi(enm1, enm1))
        yr = abs(hr(en, enm1)) + abs(hi(en, enm1))
        zzr = abs(hr(en, en)) + abs(hi(en, en))
        !     ********** for m=en-1 step -1 until l do -- **********
        do mm = l, enm1
            m = enm1 + l - mm
            if (m == l) go to 420
            yi = yr
            yr = abs(hr(m, m - 1)) + abs(hi(m, m - 1))
            xi = zzr
            zzr = xr
            xr = abs(hr(m - 1, m - 1)) + abs(hi(m - 1, m - 1))
            if (yr <= machep * zzr / yi * (zzr + xr + xi)) go to 420
        end do
        !     ********** triangular decomposition h=l*r **********
        420 mp1 = m + 1
        !
        do i = mp1, en
            im1 = i - 1
            xr = hr(im1, im1)
            xi = hi(im1, im1)
            yr = hr(i, im1)
            yi = hi(i, im1)
            if (abs(xr) + abs(xi) >= abs(yr) + abs(yi)) go to 460
            !     ********** interchange rows of hr and hi **********
            do j = im1, n
                zzr = hr(im1, j)
                hr(im1, j) = hr(i, j)
                hr(i, j) = zzr
                zzi = hi(im1, j)
                hi(im1, j) = hi(i, j)
                hi(i, j) = zzi
            end do
            !
            z3 = cmplx_dp(xr, xi) / cmplx_dp(yr, yi)
            wr(i) = 1.0_dp
            go to 480
            460      z3 = cmplx_dp(yr, yi) / cmplx_dp(xr, xi)
            wr(i) = -1.0_dp
            480      zzr = dble(z3)
            zzi = aimag(z3)
            hr(i, im1) = zzr
            hi(i, im1) = zzi
            !
            do j = i, n
                hr(i, j) = hr(i, j) - zzr * hr(im1, j) + zzi * hi(im1, j)
                hi(i, j) = hi(i, j) - zzr * hi(im1, j) - zzi * hr(im1, j)
            end do
            !
        end do
        !     ********** composition r*l=h **********
        do j = mp1, en
            xr = hr(j, j - 1)
            xi = hi(j, j - 1)
            hr(j, j - 1) = 0.0_dp
            hi(j, j - 1) = 0.0_dp
            !     ********** interchange columns of hr, hi, zr, and zi,
            !                if necessary **********
            if (wr(j) <= 0.0_dp) go to 580
            !
            do i = 1, j
                zzr = hr(i, j - 1)
                hr(i, j - 1) = hr(i, j)
                hr(i, j) = zzr
                zzi = hi(i, j - 1)
                hi(i, j - 1) = hi(i, j)
                hi(i, j) = zzi
            end do
            !
            do i = low, igh
                zzr = zr(i, j - 1)
                zr(i, j - 1) = zr(i, j)
                zr(i, j) = zzr
                zzi = zi(i, j - 1)
                zi(i, j - 1) = zi(i, j)
                zi(i, j) = zzi
            end do
            !
            580   do i = 1, j
                hr(i, j - 1) = hr(i, j - 1) + xr * hr(i, j) - xi * hi(i, j)
                hi(i, j - 1) = hi(i, j - 1) + xr * hi(i, j) + xi * hr(i, j)
            end do
            !     ********** accumulate transformations **********
            do i = low, igh
                zr(i, j - 1) = zr(i, j - 1) + xr * zr(i, j) - xi * zi(i, j)
                zi(i, j - 1) = zi(i, j - 1) + xr * zi(i, j) + xi * zr(i, j)
            end do
            !
        end do
        !
        go to 240
        !     ********** a root found **********
        660 hr(en, en) = hr(en, en) + tr
        wr(en) = hr(en, en)
        hi(en, en) = hi(en, en) + ti
        wi(en) = hi(en, en)
        en = enm1
        go to 220
        !     ********** all roots found.  backsubstitute to find
        !                vectors of upper triangular form **********
        680 norm = 0.0_dp
        !
        do i = 1, n
            !
            do j = i, n
                norm = norm + abs(hr(i, j)) + abs(hi(i, j))
            end do
        end do
        !
        hr(1, 1) = norm
        if (n == 1 .or. norm == 0.0_dp) go to 1001
        !     ********** for en=n step -1 until 2 do -- **********
        do nn = 2, n
            en = n + 2 - nn
            xr = wr(en)
            xi = wi(en)
            enm1 = en - 1
            !     ********** for i=en-1 step -1 until 1 do -- **********
            do ii = 1, enm1
                i = en - ii
                zzr = hr(i, en)
                zzi = hi(i, en)
                if (i == enm1) go to 760
                ip1 = i + 1
                !
                do j = ip1, enm1
                    zzr = zzr + hr(i, j) * hr(j, en) - hi(i, j) * hi(j, en)
                    zzi = zzi + hr(i, j) * hi(j, en) + hi(i, j) * hr(j, en)
                end do
                !
                760           yr = xr - wr(i)
                yi = xi - wi(i)
                if (yr == 0.0_dp .and. yi == 0.0_dp) yr = machep * norm
                z3 = cmplx_dp(zzr, zzi) / cmplx_dp(yr, yi)
                hr(i, en) = dble(z3)
                hi(i, en) = aimag(z3)
            end do
            !
        end do
        !     ********** end backsubstitution **********
        enm1 = n - 1
        !     ********** vectors of isolated roots **********
        do i = 1, enm1
            if (i >= low .and. i <= igh) go to 840
            ip1 = i + 1
            !
            do j = ip1, n
                zr(i, j) = hr(i, j)
                zi(i, j) = hi(i, j)
            end do
            !
            840 continue
        end do
        !     ********** multiply by transformation matrix to give
        !                vectors of original full matrix.
        !                for j=n step -1 until low+1 do -- **********
        do jj = low, enm1
            j = n + low - jj
            m = min0(j - 1, igh)
            !
            do i = low, igh
                zzr = zr(i, j)
                zzi = zi(i, j)
                !
                do k = low, m
                    zzr = zzr + zr(i, k) * hr(k, j) - zi(i, k) * hi(k, j)
                    zzi = zzi + zr(i, k) * hi(k, j) + zi(i, k) * hr(k, j)
                end do
                !
                zr(i, j) = zzr
                zi(i, j) = zzi
            end do
        end do
        !
        go to 1001
        !     ********** set error -- no convergence to an
        !                eigenvalue after 30 iterations **********
        1000 ierr = en
        1001 return
    end subroutine
end module