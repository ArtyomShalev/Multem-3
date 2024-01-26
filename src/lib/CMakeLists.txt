# set type specific output defaults
include(GNUInstallDirs)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")

# adding libs and setting dependencies
add_library(multem_blas multem_blas.f90)
target_link_libraries(multem_blas lapack::lapack)

add_library(multipole_regime_parameters multipole_regime_parameters.f90)
target_link_libraries(multipole_regime_parameters FiNeR::FiNeR)

add_library(libmultem2b libmultem2b.f90)

# linking all libs together
target_link_libraries(libmultem2b
        multipole_regime_parameters
        multem_blas
        amos::amos
        faddeeva_fortran::faddeeva_fortran)

add_executable(multem3 multem_mod.f90)

target_link_libraries(multem3
        libmultem2b)

add_custom_command(
       TARGET multem3 POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy
       ${CMAKE_SOURCE_DIR}/bin/multem3
       ${CMAKE_SOURCE_DIR}/simulation/multem3)