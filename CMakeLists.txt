cmake_minimum_required(VERSION 2.8.4)
project(CGTOOL)

MESSAGE("${CMAKE_SYSTEM_NAME} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-Og -g -Wall")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -g -Wall")

# Should we ignore warnings - default ignore a few
option(INGORE_SOME_WARNINGS "Don't ignore any warnings" ON)


# Platform specifics
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    # Apple OSX - don't use OPENMP
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
else()
    # GCC and CLANG use different OPENMP flags
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -pthread")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Og")
        if(IGNORE_SOME_WARNINGS)
            set(IGNORED_WARNINGS " -Wno-sign-compare -Wno-unused-variable -Wno-comment")
            set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${IGNORED_WARNINGS})
            set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${IGNORED_WARNINGS})
        else()
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-comment -Wno-sign-compare")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-comment -Wno-sign-compare")
        endif()
    else()
        MESSAGE("Not using GNU or Clang")
    endif()
endif()

message("Using " ${CMAKE_CXX_FLAGS})

link_directories("/usr/local/lib")


# Find libxdrfile - compile it if not present - replaces GROMACS library
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
find_library(XDRFILE_LIBRARY NAMES xdrfile
             PATHS /usr/lib /usr/local/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib ${PROJECT_BINARY_DIR})
if(XDRFILE_LIBRARY)
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


# Find Boost library - required
find_package(Boost COMPONENTS program_options REQUIRED)
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


# handle conditional compilation of testing features
option(UPDATE_PROGRESS "Print progress update every few frames" ON)
option(BOLTZMANN_INVERSION "Do a Boltzmann Inversion to calculate force constants" OFF)
option(ELECTRIC_FIELD "Calculate dipoles and the electric field" OFF)

if(UPDATE_PROGRESS)
    add_definitions(-DUPDATE_PROGRESS)
endif()
if(BOLTZMANN_INVERSION)
    add_definitions(-DBOLTZMANN_INVERSION)
endif()
if(ELECTRIC_FIELD)
    add_definitions(-DELECTRIC_FIELD)
endif()


# Add the main target
FILE(GLOB SOURCE_FILES src/*.cpp)
add_executable(cgtool ${SOURCE_FILES})
target_link_libraries (cgtool ${EXTRA_LIBS})

# Testing - do we have GTEST
find_package(GTest)
if(GTEST_FOUND)
    set(GTEST_SOURCE_FILES ${SOURCE_FILES})
    include_directories(${GTEST_INCLUDE_DIRS})

    # test parser
    add_executable(gtest_parser src/parser.cpp src/tests/parser_test.cpp ${XDR_FILES})
    target_link_libraries (gtest_parser ${GTEST_BOTH_LIBRARIES} ${EXTRA_LIBS})
    add_test(GTestParserAll gtest_parser)
    # test bondset
    add_executable(gtest_bondset src/bondset.cpp src/boltzmann_inverter.cpp
                   src/bond_struct.cpp src/array.cpp  src/frame.cpp src/parser.cpp
                   src/tests/bondset_test.cpp ${XDR_FILES})
    target_link_libraries (gtest_bondset ${GTEST_BOTH_LIBRARIES} ${EXTRA_LIBS})
    add_test(GTestBondSetAll gtest_bondset)
    # test array
    add_executable(gtest_array src/array.cpp src/tests/array_test.cpp ${XDR_FILES})
    target_link_libraries (gtest_array ${GTEST_BOTH_LIBRARIES} ${EXTRA_LIBS})
    add_test(GTestArrayAll gtest_array)

    # integration test - does it run
    add_test(Integration1 cgtool ../test_data/ALLA)

    enable_testing()
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                      DEPENDS gtest_parser gtest_bondset cgtool)
endif()
