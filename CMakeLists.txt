project(StreamComparison)
cmake_minimum_required(VERSION 3.1)

add_subdirectory(ThreadPool)
include_directories(ThreadPool/include)

set(WARNING_FLAGS -Wall -Wextra -pedantic)

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

find_package(OpenMP REQUIRED)

add_executable(STREAM STREAM.cc)
target_compile_options(STREAM PUBLIC ${WARNING_FLAGS} ${OpenMP_CXX_FLAGS})
set_target_properties(STREAM PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED True)
target_link_libraries(STREAM ThreadPool ${OpenMP_CXX_FLAGS})

add_executable(STREAM_capture_test STREAM_capture_test.cc)
target_compile_options(STREAM_capture_test PUBLIC ${WARNING_FLAGS} ${OpenMP_CXX_FLAGS})
set_target_properties(STREAM_capture_test PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED True)
target_link_libraries(STREAM_capture_test ThreadPool ${OpenMP_CXX_FLAGS})

add_executable(mmul mmul.cc dummy.cc)
target_compile_options(mmul PUBLIC ${WARNING_FLAGS} ${OpenMP_CXX_FLAGS})
set_target_properties(mmul PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED True)
target_link_libraries(mmul ThreadPool ${OpenMP_CXX_FLAGS})
