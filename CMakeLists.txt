cmake_minimum_required(VERSION 3.7.1)
set(CMAKE_CXX_STANDARD 14)
project(sims VERSION 1.0)

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

target_link_libraries(sims eigen common geo)
#add_subdirectory(benchmark)
enable_testing()

add_subdirectory(test)