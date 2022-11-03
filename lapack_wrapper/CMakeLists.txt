project(dense_solve Fortran C)
FIND_PACKAGE(LAPACK REQUIRED)

set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_solve.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_solve_impl.f90

)

ADD_LIBRARY(dense_solve ${SRC})
target_link_libraries(dense_solve ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES} #${HDF5_Fortran_LIBRARIES}
        )
