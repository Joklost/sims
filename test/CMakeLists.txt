add_subdirectory(libs)
add_executable(test_sims
        $<TARGET_OBJECTS:sims>
        $<TARGET_OBJECTS:geo>
        $<TARGET_OBJECTS:common>
        main.cpp test.cpp)

target_link_libraries(test_sims PUBLIC sims Catch2)
target_include_directories(test_sims PUBLIC ${PROJECT_SOURCE_DIR}/test)

add_test(NAME TestSims COMMAND test_sims)