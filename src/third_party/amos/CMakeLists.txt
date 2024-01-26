cmake_minimum_required(VERSION 3.22.1)
project(amos VERSION 2.0.3 LANGUAGES Fortran)

# seach path for additional cmake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules)

# set export variables needed for building
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
set(NAMESPACE "${PROJECT_NAME}::")

# generate the library and install instructions
add_subdirectory(src/lib)

# generate package config files
include(GNUInstallDirs)
set(project_config "${PROJECT_NAME}Config.cmake")
set(cmake_files_dir "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles")
set(config_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
set(config_build_dir "${CMAKE_CURRENT_BINARY_DIR}/${config_install_dir}")

# export targets for install
install(EXPORT ${TARGETS_EXPORT_NAME}
        NAMESPACE
        ${NAMESPACE}
        DESTINATION
        ${config_install_dir}
        COMPONENT Development
        )

# export targets into build
export(EXPORT ${TARGETS_EXPORT_NAME}
        NAMESPACE
        ${NAMESPACE}
        FILE
        ${config_build_dir}/${TARGETS_EXPORT_NAME}.cmake
        )
