# set type specific output defaults
include(GNUInstallDirs)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")


#include_directories(${CMAKE_BINARY_DIR}/amos
#        ${CMAKE_BINARY_DIR}/lapack_wrapper
#        )

#find_package(lapack::lapack)

add_library(multem_blas multem_blas.f90)
#target_link_libraries(multem_blas FiNeR::FiNeR)
target_link_libraries(multem_blas lapack::lapack)
add_library(multipole_regime_parameters multipole_regime_parameters.f90)
target_link_libraries(multipole_regime_parameters FiNeR::FiNeR)
#


add_library(libmultem2b libmultem2b.f90)
#target_link_libraries(libmultem2b amos dense_solve multem_blas errfun::errfun multipole_regime_parameters)
#
target_link_libraries(libmultem2b
        multipole_regime_parameters
        multem_blas
        amos::amos
        faddeeva_fortran::faddeeva_fortran)

add_executable(multem3 multem_mod.f90)

target_link_libraries(multem3
        libmultem2b)
#        dense_solve libmultem2b multem_blas
#        errfun
#        )



#add_custom_command(
#        TARGET multem2 POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        #        ${CMAKE_SOURCE_DIR}/multem2mod/input2_band.txt
#        ${CMAKE_SOURCE_DIR}/multem2mod/input_inoue.txt
#        ${CMAKE_CURRENT_BINARY_DIR}/../../bin/fort.10)
#
#add_custom_command(
#        TARGET multem2 POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_SOURCE_DIR}/multem2mod/multipole_regime_parameters_ref.txt
#        ${CMAKE_CURRENT_BINARY_DIR}/../../bin/multipole_regime_parameters.ini)