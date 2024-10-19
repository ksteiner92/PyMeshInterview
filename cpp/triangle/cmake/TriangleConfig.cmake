get_filename_component(TRIANGLE_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
list(APPEND CMAKE_MODULE_PATH ${TRIANGLE_CMAKE_DIR})

if(NOT TARGET Mesh::Triangle)
    include("${TRIANGLE_CMAKE_DIR}/TriangleTargets.cmake")
endif()
