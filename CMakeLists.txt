cmake_minimum_required(VERSION 2.8.4)
project(traj_process2)

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall")

include_directories ("/usr/local/gromacs/include"
                     "/usr/local/gromacs/lib/x86_64-linux-gnu")

set (EXTRA_LIBS ${EXTRA_LIBS} gromacs)

set(SOURCE_FILES
    src/arrays.cpp
    src/field_map.cpp
    src/frame.cpp
    src/bondset.cpp
    src/cg_map.cpp
    src/parser.cpp
    src/main.cpp)

add_executable(tp2 ${SOURCE_FILES})
target_link_libraries (tp2 ${EXTRA_LIBS})


enable_testing()
add_test(tp_runs_1_arg tp2 ../test_data)
