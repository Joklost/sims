cmake_minimum_required(VERSION 3.7.1)
set(CMAKE_CXX_STANDARD 14)
project(sims VERSION 1.0)

include(GNUInstallDirs)

find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if (GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if (NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif ()
    endif ()
endif ()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/libs/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif ()

add_subdirectory(libs)

add_library(sims
        src/math.cpp src/radiomodel.cpp src/linkmodel.cpp
        src/cholesky.cpp src/clustering.cpp
        src/datagen.cpp src/node.cpp src/link.cpp src/linalg.cpp)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(sims PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        PRIVATE src test benchmark)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS sims EXPORT SimsLibraryConfig
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from install directory.
install(EXPORT SimsLibraryConfig DESTINATION cmake)

# This makes the project importable from the build directory
export(TARGETS sims FILE SimsLibraryConfig.cmake)


target_link_libraries(sims eigen common geo)
#add_subdirectory(benchmark)
enable_testing()

add_subdirectory(test)
