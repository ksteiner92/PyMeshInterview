include(CMakeFindDependencyMacro)

get_filename_component(MESH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
list(APPEND CMAKE_MODULE_PATH ${MESH_CMAKE_DIR})

find_package(Eigen3 REQUIRED)
find_package(Math REQUIRED)
find_package(Triangle REQUIRED)
find_package(IO REQUIRED)

if(NOT TARGET Mesh::Mesh)
    include("${MESH_CMAKE_DIR}/MeshTargets.cmake")
endif()
