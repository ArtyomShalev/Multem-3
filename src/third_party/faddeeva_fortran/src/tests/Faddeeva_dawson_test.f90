program Faddeeva_dawson_test
    use faddeeva_fortran_interface, only: Faddeeva_dawson, find_relerr, dp
    ! -----------------------------------------------------------------------
    implicit none
    complex(dp), allocatable :: res(:)
    real(dp), allocatable :: re_err(:), im_err(:)
    real(dp) :: relerr, errmax
    integer :: i, num_args
    include '/references/reference_data.f90'
    ! -----------------------------------------------------------------------
    num_args = size(Faddeeva_dawson_args)
    allocate(res(1:num_args))
    allocate(re_err(1:num_args))
    allocate(im_err(1:num_args))
    relerr = 0.0_dp
    errmax = 0.0_dp
    do i = 1, num_args
        res(i) = Faddeeva_dawson(Faddeeva_dawson_args(i), relerr)
        re_err(i) = find_relerr(real(res(i)), real(Faddeeva_dawson_wolfram_refs(i)))
        im_err(i) = find_relerr(aimag(res(i)), aimag(Faddeeva_dawson_wolfram_refs(i)))
        print *, "dawson(", real(Faddeeva_dawson_args(i)),"+i*", aimag(Faddeeva_dawson_args(i)),") = ", real(res(i)),"+i*", &
                aimag(res(i)), " (vs.", real(Faddeeva_dawson_wolfram_refs(i)),"+i*", aimag(Faddeeva_dawson_wolfram_refs(i)),"), &
                re/im rel. err. = ", re_err(i),"/", im_err(i)
        if (re_err(i) > errmax) errmax = re_err(i)
        if (im_err(i) > errmax) errmax = im_err(i)
        if (errmax > 1e-13) then
            print *, "FAILURE -- relative error", errmax, "is too large!"
            stop 1
        end if
    end do
    print *, "SUCCESS (max relative error = ", errmax
end program Faddeeva_dawson_test