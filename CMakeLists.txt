cmake_minimum_required(VERSION 3.14)
project(PiFromScratch CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

include_directories(include)

# Main library
add_library(pi_math src/BigInt.cpp src/NTT.cpp)

# Main executable
add_executable(pi_calc src/PiCalculator.cpp)
target_link_libraries(pi_calc PRIVATE pi_math)

# Tests
enable_testing()
add_executable(test_bigint tests/test_bigint.cpp)
target_link_libraries(test_bigint PRIVATE pi_math)
add_test(NAME test_bigint COMMAND test_bigint)
