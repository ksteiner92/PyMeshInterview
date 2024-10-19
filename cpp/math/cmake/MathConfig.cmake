get_filename_component(MATH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
list(APPEND CMAKE_MODULE_PATH ${MATH_CMAKE_DIR})

if(NOT TARGET Mesh::Math)
    include("${TRIANGLE_CMAKE_DIR}/MathTargets.cmake")
endif()
