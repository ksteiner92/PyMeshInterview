include(CMakeFindDependencyMacro)

get_filename_component(PYMESH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
list(APPEND CMAKE_MODULE_PATH ${PYMESH_CMAKE_DIR})

find_package(Mesh REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)

if(NOT TARGET Mesh::PyMesh)
    include("${CORE_CMAKE_DIR}/PyMeshTargets.cmake")
endif()
