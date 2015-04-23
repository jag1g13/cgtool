cmake_minimum_required(VERSION 2.8.4)
project(CGTOOL)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Automated coarse-grain mapping for molecular dynamics")
set(CPACK_PACKAGE_CONTACT "James Graham <J.A.Graham@soton.co.uk>")
SET(CPACK_PACKAGE_VENDOR "James Graham")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "x")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "cgtool")
IF(WIN32 AND NOT UNIX)
    # There is a bug in NSI that does not handle full unix paths properly. Make
    # sure there is at least one set of four (4) backlasshes.
    SET(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}\\\\icon.bmp")
    SET(CPACK_NSIS_INSTALLED_ICON_NAME "build\\\\cgtool.exe")
    SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} CGTOOL")
    SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\bitbucket.org/jag1g13/cgtool")
    SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\bitbucket.org/jag1g13/cgtool")
    SET(CPACK_NSIS_CONTACT "J.A.Graham@soton.ac.uk")
    SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(WIN32 AND NOT UNIX)
    SET(CPACK_STRIP_FILES "build/cgtool")
    SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(WIN32 AND NOT UNIX)
SET(CPACK_PACKAGE_EXECUTABLES "cgtool")
INCLUDE(CPack)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

MESSAGE("CMAKE ${CMAKE_VERSION}")
MESSAGE("CXX   ${CMAKE_SYSTEM_NAME} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE("C     ${CMAKE_SYSTEM_NAME} ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
MESSAGE("${CMAKE_MODULE_PATH}")

# TODO add C++11 checks
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -msse2 -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

# Should we ignore warnings - default ignore a few - only on gcc
option(INGORE_SOME_WARNINGS "Don't ignore any warnings" ON)
set(IGNORED_WARNINGS " -Wno-sign-compare -Wno-unused-variable -Wno-comment -Wno-unused-but-set-variable")
set(IGNORED_WARNINGS_ALWAYS " -Wno-sign-compare -Wno-comment")

# Handle conditional compilation of features
option(UPDATE_PROGRESS "Print progress update every few frames" ON)
option(NO_CMD_PARSER "Disable the default command line parser if not using Boost program_options" OFF)

# Platform specifics
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    # Apple OSX - don't use OPENMP
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
else()
    # GCC and CLANG use different OPENMP flags
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
    else()
        MESSAGE("Not using GNU or Clang")
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
    #set(XDR_FILES lib/xdrfile.c lib/xdrfile_xtc.c)
    #set(SOURCE_FILES ${SOURCE_FILES} ${XDR_FILES})
endif()

#include_directories(/usr/local/include/FLENS)
# Find FLENS - required
#find_package(flens REQUIRED)
#include_directories($flens_INCLUDE_DIR)
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
find_package(Boost COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
    MESSAGE("Found Boost Program_Options")
else(Boost_PROGRAM_OPTIONS_FOUND)
    MESSAGE("Boost Program_Options not found - using simple command line")
    set(NO_CMD_PARSER ON)
endif(Boost_PROGRAM_OPTIONS_FOUND)
include_directories(${Boost_INCLUDE_DIRS})
set (EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})

# Add target for documentation if we have Doxygen
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


if(UPDATE_PROGRESS)
    add_definitions(-DUPDATE_PROGRESS)
endif()
if(NO_CMD_PARSER)
    add_definitions(-DNO_CMD_PARSER OR NOT Boost_PROGRAM_OPTIONS_FOUND)
endif()


# Add the main target
FILE(GLOB SOURCE_FILES src/*.cpp)
add_executable(cgtool ${SOURCE_FILES})
add_dependencies(cgtool flens)
target_link_libraries (cgtool ${EXTRA_LIBS})

# Testing
add_subdirectory(lib/gtest-1.7.0)
include_directories(lib/gtest-1.7.0/include lib/gtest-1.7.0/include)

# test parser
add_executable(gtest_parser src/parser.cpp src/tests/parser_test.cpp ${XDR_FILES})
target_link_libraries (gtest_parser gtest gtest_main ${EXTRA_LIBS})
add_test(GTestParserAll gtest_parser)
# test bondset
add_executable(gtest_bondset src/bondset.cpp src/boltzmann_inverter.cpp
               src/bond_struct.cpp src/array.cpp  src/frame.cpp src/parser.cpp
               src/tests/bondset_test.cpp src/small_functions.cpp ${XDR_FILES})
target_link_libraries (gtest_bondset gtest gtest_main ${EXTRA_LIBS})
add_test(GTestBondSetAll gtest_bondset)
# test array
add_executable(gtest_array src/array.cpp src/tests/array_test.cpp ${XDR_FILES})
target_link_libraries (gtest_array gtest gtest_main ${EXTRA_LIBS})
add_test(GTestArrayAll gtest_array)

# integration test - does it run
add_test(IntegrationRUN cgtool --dir ../test_data/ALLA)

# Integration test - check output files
add_test(IntegrationITP ${CMAKE_COMMAND} -E compare_files
   ${CMAKE_SOURCE_DIR}/build/ALLA.itp ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.itp)
add_test(IntegrationGRO ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_SOURCE_DIR}/build/ALLA.gro ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.gro)
add_test(IntegrationTOP ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_SOURCE_DIR}/build/ALLA.top ${CMAKE_SOURCE_DIR}/test_data/ALLA/ALLACG.top)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS gtest_parser gtest_bondset gtest_array cgtool)
