program get_multipole_combination_test
    use libmultem2b, only: get_multipole_combination
    ! ----------------------------------------------
    implicit none
    integer, allocatable :: lmax(:), multipole_type(:, :), multipole_order(:, :), m_projection(:, :), &
                            is_multipole_type_selected(:), is_multipole_order_selected(:),            &
                            is_m_projection_selected(:)
    include '/references/get_multipole_combination_reference.f90'
    ! ----------------------------------------------
    print *, lmax(0)

end program get_multipole_combination_test