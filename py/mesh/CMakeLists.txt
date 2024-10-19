add_library(PyMesh MODULE pymesh.cc)

find_package(Mesh REQUIRED)
find_package(Math REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)

target_compile_features(PyMesh PRIVATE cxx_std_20)
target_include_directories(PyMesh
        PRIVATE
            ${EIGEN3_INCLUDE_DIR}
            ${Mesh_INCLUDE_DIR}
            ${Math_INCLUDE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}
        )

target_link_libraries(PyMesh PRIVATE pybind11::module Mesh::Mesh Mesh::Math)
set_target_properties(PyMesh PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
set_target_properties(PyMesh PROPERTIES OUTPUT_NAME "pymesh$<$<CONFIG:Debug>:-d>")

set(PyMesh_Version 1.0)

include(GNUInstallDirs)
install(TARGETS PyMesh
        EXPORT PyMeshTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT PyMeshTargets
        FILE PyMeshTargets.cmake
        NAMESPACE Mesh::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("PyMeshConfigVersion.cmake"
        VERSION ${PyMesh_Version}
        COMPATIBILITY SameMajorVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/PyMeshConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/PyMeshConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/PyMeshConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/PyMeshConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

export(TARGETS PyMesh FILE ${CMAKE_CURRENT_BINARY_DIR}/PyMeshTargets.cmake NAMESPACE Mesh::)
export(PACKAGE PyMesh)
