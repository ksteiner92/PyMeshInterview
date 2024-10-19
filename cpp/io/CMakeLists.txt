add_library(IO SHARED logger.cc)
add_library(Mesh::IO ALIAS IO)

target_compile_features(IO PRIVATE cxx_std_20)
target_include_directories(IO
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        )
set_target_properties(IO PROPERTIES OUTPUT_NAME "io$<$<CONFIG:Debug>:-d>")

set(IO_Version 1.0)

include(GNUInstallDirs)
install(TARGETS IO
        EXPORT IOTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT IOTargets
        FILE IOTargets.cmake
        NAMESPACE Mesh::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("IOConfigVersion.cmake"
        VERSION ${IO_Version}
        COMPATIBILITY SameMajorVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/IOConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/IOConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/IOConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/IOConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS IO FILE ${CMAKE_CURRENT_BINARY_DIR}/IOTargets.cmake NAMESPACE Mesh::)
export(PACKAGE IO)
