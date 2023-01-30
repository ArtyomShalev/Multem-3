cmake_minimum_required(VERSION 3.22.1)

SET (PROJECTNAME "Multem-3")
PROJECT (${PROJECTNAME} Fortran)


set(CMAKE_BUILD_TYPE Debug
        CACHE STRING "Build type (Debug, Release)" FORCE)
set(warn "-Wall -Wunused-variable -Wimplicit-interface -Wunused-dummy-argument -Wextra -Wargument-mismatch  -fPIC -fmax-errors=1")
set(CMAKE_Fortran_FLAGS_DEBUG   "${common} ${warn} -g -fbounds-check -fcheck-array-temporaries -fbacktrace")

#set(CMAKE_BUILD_TYPE Release
#            CACHE STRING "Build type (Debug, Release)" FORCE)
##set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops -mtune=native -fopenmp" )
#set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops -mtune=native" )

add_subdirectory(src/lib)
add_subdirectory(src/third_party/FiNeR)
add_subdirectory(src/third_party/lapack_wrapper)
add_subdirectory(src/third_party/amos)
add_subdirectory(src/third_party/faddeeva_fortran)