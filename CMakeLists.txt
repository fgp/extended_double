set(BOOST_ROOT $ENV{HOME}/Installs/boost)
set(BOOST_FRAMEWORK $ENV{HOME}/Installs/Frameworks)

# CMake settings
cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

# Tester
add_executable(tests main.cpp tests.cpp extended_double.cpp code.cpp)

# Find boost
if (EXISTS ${BOOST_FRAMEWORK})
	message("Searching for boost framework in ${BOOST_FRAMEWORK}")
	find_library(Boost NAMES boost PATHS ${BOOST_FRAMEWORK})
	target_link_libraries(tests ${Boost})
elseif (EXISTS ${BOOST_ROOT})
	message("Searching for boost in ${BOOST_ROOT}")
	find_package(Boost REQUIRED COMPONENTS unit_test_framework test_exec_monitor prg_exec_monitor)
	include_directories(BEFORE SYSTEM ${Boost_INCLUDE_DIR})
	target_link_libraries(tests ${Boost_LIBRARIES})
elseif (1)
	message(SEND_ERROR "Neither ${BOOST_FRAMEWORK} nor ${BOOST_ROOT} exists")
endif (EXISTS ${BOOST_FRAMEWORK})

# Enable SSE4.2
add_definitions(-msse4.2)

# Defines
add_definitions(-DED_ENABLE_ASSERTS_NORMALIZATION=1)
add_definitions(-DED_ENABLE_ASSERTS_OVERFLOW=1)
add_definitions(-DED_ENABLE_ASSERTS_STATIC=1)
