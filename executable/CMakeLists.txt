add_executable(MainStokes main.cpp)
add_executable(Test1d 1dtest.cpp)
add_executable(TestMatRedSolveDirect TestMatRedSolveDirect.cpp)

target_link_libraries(MainStokes libsources)
target_link_libraries(Test1d libsources)
target_link_libraries(TestMatRedSolveDirect libsources)
