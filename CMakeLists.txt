cmake_minimum_required(VERSION 2.8.4)
project(traj_process2)

#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Og -g -Wall -fopenmp -Wno-sign-compare -Wno-unused-variable")

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
    src/main.cpp)

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

add_executable(tp2 ${SOURCE_FILES})
target_link_libraries (tp2 ${EXTRA_LIBS})

# manual test - does it run
add_test(tp_runs_1_arg tp2 ../test_data/ALLA)

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
