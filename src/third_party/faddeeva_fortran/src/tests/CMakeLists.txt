include(GNUInstallDirs)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests/${CMAKE_INSTALL_BINDIR}")

function(add_faddeeva_fortran_test test_src libraries)
    get_filename_component(test_name ${test_src} NAME_WE)
    get_filename_component(test_dir ${test_src} DIRECTORY)
    set(result_src "${test_dir}/${test_name}.result")

    add_executable(${test_name} ${test_src})
    foreach(lib IN LISTS libraries)
        target_link_libraries(${test_name} ${lib})
    endforeach()
    add_test(NAME ${test_name} COMMAND ${test_name})
    #
    #    file(STRINGS ${result_src} result_raw)
    #    string(REGEX REPLACE "([+.])([^+.]*)" "\\\\\\1\\2" result ${result_raw})
    #
    #    set_tests_properties(${test_name} PROPERTIES
    #        PASS_REGULAR_EXPRESSION "${result}"
    #    )
endfunction()

file(GLOB faddeeva_fortran_tests *.f90)
foreach(faddeeva_fortran_test IN LISTS faddeeva_fortran_tests)
    add_Faddeeva_fortran_test(${faddeeva_fortran_test} faddeeva_fortran)
endforeach()