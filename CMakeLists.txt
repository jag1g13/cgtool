cmake_minimum_required(VERSION 2.8.4)
project(traj_process2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall")

include_directories ("/usr/local/gromacs/include"
                     "/usr/local/gromacs/lib/x86_64-linux-gnu")

set (EXTRA_LIBS ${EXTRA_LIBS} gromacs)

set(SOURCE_FILES
    frame.cpp
    bondset.cpp
    cg_map.cpp
    general.cpp
    parser.cpp
    main.cpp)

add_executable(traj_process2 ${SOURCE_FILES})
target_link_libraries (traj_process2 ${EXTRA_LIBS})

enable_testing()
add_test(tp_runs_1_arg traj_process2 ../test_data)
