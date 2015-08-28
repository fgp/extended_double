set(BOOST_ROOT $ENV{HOME}/Installs/boost)
set(BOOST_FRAMEWORK $ENV{HOME}/Installs/Frameworks)

# CMake settings
cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)
enable_language(CXX)

# Enable SSE4.2
add_compile_options(-msse4.1 -fvisibility=hidden -fvisibility-inlines-hidden)

# Tester
add_executable(tests main.cpp tests.cpp extended_double.cpp code.cpp)

# Defines
target_compile_definitions(tests PUBLIC ED_ENABLE_ASSERTS_STATIC=1)
target_compile_definitions(tests PUBLIC $<$<CONFIG:Debug>:ED_ENABLE_ASSERTS_NORMALIZATION=1>)

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
