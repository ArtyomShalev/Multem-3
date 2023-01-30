module faddeeva_fortran_interface
    integer, parameter, public :: dp=kind(0.d0)
    interface
        !> @brief the error function
        !!
        !! computes erf(z), where erf - the error function
        !!
        !! @param[in]   z       a complex value
        !! @param[in]   relerr  demanding relative error
        !! @return      a complex value of the mentioned function
        complex (c_double_complex) function Faddeeva_erf(z, relerr) bind(c, name = 'Faddeeva_erf')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function Faddeeva_erf

        !> @brief the imaginary error function
        !!
        !! computes erfi(z) = -i*erf(i*z), where erf is the error function
        !!
        !! @param[in]   z       a complex value
        !! @param[in]   relerr  demanding relative error
        !! @return      a complex value of the mentioned function
        complex (c_double_complex) function Faddeeva_erfi(z, relerr) bind(c, name = 'Faddeeva_erfi')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function Faddeeva_erfi

        !> @brief the complementary error function
        !!
        !! computes erfc(z) = 1 - erf(z), where erf is the error function
        !!
        !! @param[in]   z       a complex value
        !! @param[in]   relerr  demanding relative error
        !! @return      a complex value of the mentioned function
        complex (c_double_complex) function Faddeeva_erfc(z, relerr) bind(c, name = 'Faddeeva_erfc')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function Faddeeva_erfc

        !> @brief the scaled complementary error function
        !!
        !! computes erfcx(z) = exp(z*z)*erfc(z), where erfc is the complementary error function
        !!
        !! @param[in]   z       a complex value
        !! @param[in]   relerr  demanding relative error
        !! @return      a complex value of the mentioned function
        complex (c_double_complex) function Faddeeva_erfcx(z, relerr) bind(c, name = 'Faddeeva_erfcx')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function Faddeeva_erfcx

        !> @brief the Faddeeva function of complex argument
        !!
        !! computes w(z) = exp(-z*z)*erfc(-i*z), where erfc - the complementary error function
        !!
        !! @param[in]   z       argument
        !! @param[in]   relerr  demanding relative error
        !! @return      a complex value of the mentioned function
        complex (c_double_complex) function Faddeeva_w(z, relerr) bind(c, name = 'Faddeeva_w')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function Faddeeva_w

        !> @brief the dawson function
        !!
        !! computes Dawson(z) = sqrt(pi)/2*exp(-z*z)*erfi(z), the erfi is the imaginary error function
        !!
        !! @param[in]   z       a complex value
        !! @param[in]   relerr  demanding relative error
        !! @return      a complex value of the mentioned function
        complex (c_double_complex) function Faddeeva_dawson(z, relerr) bind(c, name = 'Faddeeva_Dawson')
            use iso_c_binding
            complex(c_double_complex), value :: z
            real(c_double), value :: relerr
        end function Faddeeva_dawson
    end interface

    contains
        !> @brief provide an example of usage
        !!
        !! provide a terminal output of function name, input arguments,
        !! demanding relative error and computational results
        !!
        !! @param[in]   func            a pointer to the function
        !! @param[in]   func_name       function name as a string
        !! @param[in]   z               a complex value
        !! @param[in]   relerr          demanding relative error
        subroutine example(func, func_name, z, relerr)
            character(len=15), intent(in) :: func_name
            complex(dp) :: z
            real(dp) :: relerr

            interface AFunc
                function func(z, relerr)
                    use iso_c_binding
                    complex(c_double_complex), value :: z
                    real(c_double), value :: relerr
                    complex(c_double_complex) :: func
                end function
            end interface

            if (relerr < 2.22e-16_dp) relerr = 2.22e-16_dp !according to https://en.wikipedia.org/wiki/Machine_epsilon
            print *, func_name, ' of argument z =', z, 'with relative error' , relerr, 'equals ', func(z, relerr)

        end subroutine example

        !> @brief define a relative error between two real input arguments
        !!
        !! @param[in]   a       real argument
        !! @param[in]   b       real argument
        !! @return      a relative error
        real(dp) function find_relerr(a, b)
            real(dp), intent(in) :: a, b

            find_relerr = abs((b-a)/a)

        end function find_relerr
end module faddeeva_fortran_interface
