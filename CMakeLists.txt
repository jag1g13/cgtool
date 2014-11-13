project(traj_process2)

include_directories(/usr/local/gromacs/include)
add_executable(traj_process2 main.cpp)

install(TARGETS traj_process2 RUNTIME DESTINATION bin)
