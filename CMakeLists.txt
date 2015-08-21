cmake_minimum_required(VERSION 2.8.4)

project(CGTOOL)
set(CGTOOL_VERSION_MAJOR 0)
set(CGTOOL_VERSION_MINOR 3)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/")

MESSAGE("CMAKE ${CMAKE_VERSION}")
MESSAGE("CXX   ${CMAKE_SYSTEM_NAME} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE("C     ${CMAKE_SYSTEM_NAME} ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
MESSAGE("MOD   ${CMAKE_MODULE_PATH}")

# TODO add C++11 checks
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -march=native")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Defaults to Release build" FORCE)
endif()

if(NOT CMAKE_FLAGS_Release)
    set(CMAKE_FLAGS_Release "-O3")
endif()
if(NOT CMAKE_FLAGS_Debug)
    set(CMAKE_FLAGS_Release "-g")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_Release}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_Release}")
endif()

# Should we ignore warnings - default ignore a few - only on gcc
option(INGORE_SOME_WARNINGS "Don't ignore any warnings" ON)
set(IGNORED_WARNINGS " -Wno-sign-compare -Wno-unused-variable -Wno-comment -Wno-unused-but-set-variable")
set(IGNORED_WARNINGS_ALWAYS " -Wno-sign-compare")

# Platform specifics
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    # Apple OSX - don't use OPENMP
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # GCC and Clang use different OPENMP flags
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -pthread")

    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        if(IGNORE_SOME_WARNINGS)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${IGNORED_WARNINGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IGNORED_WARNINGS}")
        else()
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${IGNORED_WARNINGS_ALWAYS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IGNORED_WARNINGS_ALWAYS}")
        endif()

    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")

    else()
        MESSAGE("Unrecognised compiler - will attempt to compile")
    endif()
endif()

message("Using " ${CMAKE_CXX_FLAGS})

# Find libxdrfile - compile it if not present - replaces GROMACS library
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
find_library(XDRFILE_LIBRARY NAMES xdrfile
             PATHS /usr/lib /usr/local/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib)
if(XDRFILE_LIBRARY_FOUND)
    MESSAGE("Found libxdrfile: ${XDRFILE_LIBRARY}")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${XDRFILE_LIBRARY})
else()
    MESSAGE("Cannot find libxdrfile, including in compile")
    set(LIBRARY_OUTPUT_DIRCTORY lib)
    add_library(xdrfile STATIC lib/xdrfile/xdrfile.c lib/xdrfile/xdrfile_xtc.c)
    set (EXTRA_LIBS ${EXTRA_LIBS} xdrfile)
endif()

# Find FLENS - required
INCLUDE(ExternalProject)
ExternalProject_Add(flens
                    GIT_REPOSITORY https://github.com/michael-lehn/FLENS.git
                    GIT_TAG v1.1.0
                    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib/flens
                    UPDATE_COMMAND ""
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND ""
                    INSTALL_COMMAND "")
                include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/flens/src/flens)

# Find Boost library - required
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED)
find_package(Boost COMPONENTS REQUIRED program_options)
include_directories(${Boost_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})

# Find Python to create compilation info
find_package(PythonInterp 2.6 REQUIRED)
message("-- Generating compilation info...")
execute_process(COMMAND python
                ${CMAKE_CURRENT_SOURCE_DIR}/src/compile_info.py
                ${CMAKE_CURRENT_SOURCE_DIR}/include/compile_info.inc
                ${CMAKE_CURRENT_SOURCE_DIR}/include/revision_number.inc)

# Add target for documentation if we have Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

# Find source files
set(CGTOOLCORE_FILES
    "src/main/common.cpp"
    "src/frame.cpp"
    "src/cg_map.cpp"
    "src/cmd.cpp"
    "src/parser.cpp"
    "src/residue.cpp"
    "src/small_functions.cpp"
    "src/XTCOutput.cpp")

set(CGTOOL_FILES
    "src/array.cpp"
    "src/boltzmann_inverter.cpp"
    "src/bond_struct.cpp"
    "src/bondset.cpp"
    "src/field_map.cpp"
    "src/histogram.cpp"
    "src/itp_writer.cpp"
    "src/light_array.cpp"
    "src/non_bondset.cpp"
    "src/rdf.cpp")

# Add CGTOOLCORE library
add_library(cgtoolcore ${CGTOOLCORE_FILES})
target_link_libraries(cgtoolcore ${EXTRA_LIBS})

# Add CGTOOL target
add_executable(cgtool "src/main/cgtool.cpp" ${CGTOOL_FILES})
add_dependencies(cgtool flens)
target_link_libraries (cgtool cgtoolcore)

# Add RAMSi target
add_executable(ramsi "src/main/ramsi.cpp" "src/membrane.cpp" "src/array.cpp" ${SOURCE_FILES})
target_link_libraries(ramsi cgtoolcore)

# Testing
add_subdirectory(lib/gtest-1.7.0)
include_directories(lib/gtest-1.7.0/include lib/gtest-1.7.0/include)

# Test parser
add_executable(gtest_parser src/tests/parser_test.cpp)
target_link_libraries (gtest_parser gtest gtest_main cgtoolcore)
add_test(GTestParserAll gtest_parser)
# Test bondset
add_executable(gtest_bondset src/tests/bondset_test.cpp
    src/bondset.cpp src/boltzmann_inverter.cpp src/array.cpp src/histogram.cpp src/bond_struct.cpp)
target_link_libraries (gtest_bondset gtest gtest_main cgtoolcore)
add_test(GTestBondSetAll gtest_bondset)
# Test array
add_executable(gtest_array src/tests/array_test.cpp src/array.cpp)
target_link_libraries (gtest_array gtest gtest_main cgtoolcore)
add_test(GTestArrayAll gtest_array)

# Integration test - does it run
add_test(IntegrationRUNCGTOOL cgtool --dir ../test_data/ALLA -c cg.cfg -x md.xtc -g md.gro -i topol.top)
# Integration test - check output files
add_test(IntegrationITP ${CMAKE_COMMAND} -E compare_files
   ${CMAKE_SOURCE_DIR}/build/ALLA.itp ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.itp)
add_test(IntegrationGRO ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_SOURCE_DIR}/build/ALLA.gro ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.gro)
add_test(IntegrationTOP ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_SOURCE_DIR}/build/ALLA.top ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.top)

add_test(IntegrationRUNRAMSi ramsi -c ../test_data/staph/mem.cfg -x ../test_data/staph/md.xtc -g ../test_data/staph/md.gro)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS gtest_parser gtest_bondset gtest_array cgtool ramsi)
