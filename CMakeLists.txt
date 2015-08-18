set(BOOST_ROOT ~/Installs/boost)

# CMake settings
cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

# Find boost
find_package(Boost REQUIRED COMPONENTS unit_test_framework test_exec_monitor prg_exec_monitor)

# Include paths
include_directories(BEFORE SYSTEM ${BOOST_INCLUDE_DIRS})

add_definitions(-DED_ENABLE_ASSERTS_NORMALIZATION=1 -DED_ENABLE_ASSERTS_OVERFLOW=1)

# Tester
add_executable(tests main.cpp tests.cpp extended_double.cpp)
target_link_libraries(tests ${Boost_LIBRARIES})
