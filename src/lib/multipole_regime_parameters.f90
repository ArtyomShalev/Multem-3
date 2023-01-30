module multipole_regime_parameters
    use flap, only : command_line_interface
    use finer, only : file_ini

    implicit none

    type, private :: multipole_parameters
        integer ::  is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected
        integer, allocatable :: multipole_type(:), multipole_order(:), m_projection(:)
    end type multipole_parameters

    type(multipole_parameters), public :: mrp

    character(999), private :: ini_config_file_name  !< Name of INI file.
    type(file_ini), private :: fini       !< INI file handler.
    integer, allocatable :: array(:)

    public cli_parse, ini_parse

    contains

    subroutine ini_parse()

        integer :: error, num
        write(6, *) 'reading config from file: ', ini_config_file_name
        call fini%load(filename = ini_config_file_name)

        call fini%get(section_name = 'selectors', option_name = 'is_multipole_type_selected', val = num, error = error)
        mrp%is_multipole_type_selected = 0
        if ( error==0 .and. (num == 0 .or. num == 1) ) mrp%is_multipole_type_selected = num

        call fini%get(section_name = 'selectors', option_name = 'is_multipole_order_selected', val = num, error = error)
        mrp%is_multipole_order_selected = 0
        if ( error==0 .and. (num == 0 .or. num == 1) ) mrp%is_multipole_order_selected = num

        call fini%get(section_name = 'selectors', option_name = 'is_m_projection_selected', val = num, error = error)
        mrp%is_m_projection_selected = 0
        if ( error==0 .and. (num == 0 .or. num == 1) ) mrp%is_m_projection_selected = num

        allocate(array(1:fini%count_values(section_name='regime', option_name='multipole_type')))
        call fini%get(section_name = 'regime', option_name = 'multipole_type', val = array, error = error)
        if (error==0) mrp%multipole_type = array
        deallocate(array)

        allocate(array(1:fini%count_values(section_name='regime', option_name='multipole_order')))
        call fini%get(section_name = 'regime', option_name = 'multipole_order', val = array, error = error)
        if (error==0) mrp%multipole_order = array
        deallocate(array)

        allocate(array(1:fini%count_values(section_name='regime', option_name='m_projection')))
        call fini%get(section_name = 'regime', option_name = 'm_projection', val = array, error = error)
        if (error==0) mrp%m_projection = array
        deallocate(array)

    end subroutine ini_parse

    subroutine cli_parse()
        !       !< build and parse test cli.
        type(command_line_interface) :: cli  !< command line interface.
        integer :: error !< error trapping flag.

        call cli%init(progname = 'multem2', &
                authors = '', &
                help = 'usage: ', &
                examples = ["multem2 -i config.ini"], &
                epilog = new_line('a') // "all done")

        call cli%add(switch = '--ini', &
                switch_ab = '-i', &
                help = 'name of ini file', &
                required = .false., &
                def = 'multipole_regime_parameters.ini', &
                act = 'store')

        call cli%parse(error = error) ; if (error/=0) stop

        call cli%get(switch = '--ini', val = ini_config_file_name)
    endsubroutine cli_parse

end module multipole_regime_parameters