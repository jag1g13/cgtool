cmake_minimum_required(VERSION 2.8.4)
project(cgtool)

#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-sign-compare -Wno-unused-variable")

MESSAGE("${CMAKE_CXX_COMPILER}")
MESSAGE("${CMAKE_CXX_COMPILER_ID}")
MESSAGE("${CMAKE_SYSTEM_NAME}")



if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    # Apple OSX
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    # doesn't work
else()
    # Not Apple OSX
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -openmp")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -g -fopenmp")
    else()
        MESSAGE("Not using GNU or Clang")
    endif()
endif()

message("Using " ${CMAKE_CXX_FLAGS})

#find_package(gromacs REQUIRED)
include_directories ("/usr/local/gromacs/include")
link_directories("/usr/local/lib")

set (EXTRA_LIBS ${EXTRA_LIBS} gromacs)

set(SOURCE_FILES
    src/arrays.cpp
    src/field_map.cpp
    src/frame.cpp
    src/bondset.cpp
    src/cg_map.cpp
    src/parser.cpp
    src/cmd.cpp)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)

add_executable(cgtool ${SOURCE_FILES} src/main.cpp)
target_link_libraries (cgtool ${EXTRA_LIBS})

# gtest
#if(GTEST_FOUND)
#message("Found GTest")
set(GTEST_SOURCE_FILES
    ${SOURCE_FILES}
    src/tests/parser_test.cpp)
    #src/parser.cpp)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src)
add_executable(cgtool_gtest ${GTEST_SOURCE_FILES})
target_link_libraries (cgtool_gtest ${GTEST_BOTH_LIBRARIES} gromacs)
add_test(GTestParserAll cgtool_gtest)
#endif()

# manual test - does it run
add_test(runs_1_arg cgtool ../test_data/ALLA)

# names of the tests
set(cpptests
#    csv
    frame)

# how many parts does each test have
set(csv_parts
    1 2 3)
set(frame_parts 1 2 3)

# if we didn't say how many parts a test has, assume 1
foreach(cpptest ${cpptests})
    if(NOT DEFINED "${cpptest}_parts")
        set(${cpptest}_parts "1")
    endif()
endforeach()

foreach(cpptest ${cpptests})
    set(cpptestsrc ${cpptestsrc} src/tests/${cpptest}_test.cpp)
endforeach()

create_test_sourcelist(srclist test_runner.cpp ${cpptestsrc})
add_executable(test_runner ${srclist})
target_link_libraries(test_runner ${libs})

foreach(cpptest ${cpptests})
  foreach(part ${${cpptest}_parts})
    add_test(test_${cpptest}_${part}
             ${TEST_PATH}/test_runner src/tests/${cpptest}_test ${part})
    set_tests_properties(test_${cpptest}_${part} PROPERTIES
      FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
  endforeach()
endforeach()

enable_testing()
