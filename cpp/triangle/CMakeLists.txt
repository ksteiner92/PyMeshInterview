add_library(Triangle STATIC triangle.c)
add_library(Mesh::Triangle ALIAS Triangle)

target_compile_definitions(Triangle PUBLIC TRILIBRARY ANSI_DECLARATORS)
if(WIN32)
   target_compile_definitions(Triangle PUBLIC -DNO_TIMER)
else()
   target_compile_definitions(Triangle PUBLIC -DLINUX)
endif()
target_include_directories(Triangle
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	    PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
)
set_target_properties(Triangle PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(Triangle PROPERTIES OUTPUT_NAME "triangle$<$<CONFIG:Debug>:-d>")
set_source_files_properties(triangle.c PROPERTIES LANGUAGE CXX)

set(Triangle_Version 1.0)

include(GNUInstallDirs)
install(TARGETS Triangle
        EXPORT TriangleTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT TriangleTargets
        FILE TriangleTargets.cmake
        NAMESPACE Mesh::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("TriangleConfigVersion.cmake"
        VERSION ${Triangle_Version}
        COMPATIBILITY SameMajorVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/TriangleConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/TriangleConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/TriangleConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/TriangleConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS Triangle FILE ${CMAKE_CURRENT_BINARY_DIR}/TriangleTargets.cmake NAMESPACE Mesh::)
export(PACKAGE Triangle)
