get_filename_component(IO_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
list(APPEND CMAKE_MODULE_PATH ${IO_CMAKE_DIR})

if(NOT TARGET Mesh::IO)
    include("${TRIANGLE_CMAKE_DIR}/IOTargets.cmake")
endif()
