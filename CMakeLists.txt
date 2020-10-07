cmake_minimum_required(VERSION 3.1)

project(hybridDecoder)
set(CMAKE_CXX_STANDARD 17)

project(matchingBP)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Turn on with 'cmake -Dtest=ON'
option(test "Build all tests." OFF)
# Turn on with 'cmake -Dprofile=ON'
option(profile "Profile using grpof")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(test ON)
    # set(profile ON)
endif()

set(SOURCE_FILES main.cpp)
add_executable(hybridDecoder ${SOURCE_FILES})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(include/blossom5)

include_directories(src)
set(LIB_FILES src/djikstra.h src/djikstra.cpp src/hexLattice.h src/hexLattice.cpp src/decoder.h src/decoder.cpp)
# set(LIB_FILES src/djikstra.h src/djikstra.cpp src/hexLattice.h src/hexLattice.cpp)
add_library(lattice ${LIB_FILES})
target_link_libraries(lattice blossom5)

target_link_libraries(hybridDecoder lattice)
target_link_libraries(hybridDecoder blossom5)

if (test)
    enable_testing()
    add_definitions(-DSCTEST)
    include(CTest)
    
    # Add my tests
    add_subdirectory(test)

    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                    ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                    EXCLUDE_FROM_ALL)

    # The gtest/gtest_main targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
    endif()

    # Add test executables
    add_executable(testHex test/test_hexLattice.cpp)
    add_executable(testDjikstra test/test_djikstra.cpp)
    add_executable(testDecoder test/test_decoder.cpp)
    
    # Standard googletest linking
    target_link_libraries(testHex gtest gtest_main)
    target_link_libraries(testDjikstra gtest gtest_main)
    target_link_libraries(testDecoder gtest gtest_main)

    # Link to libraries
    target_link_libraries(testHex lattice)
    target_link_libraries(testDjikstra lattice)
    target_link_libraries(testDecoder lattice)

    # Enable running tests with 'make test'
    add_test(NAME testHex COMMAND testHex)
    add_test(NAME testDjikstra COMMAND testDjikstra)
    add_test(NAME testDecoder COMMAND testDecoder)
endif()

if (profile)
    # Set gprof flags
    SET(GCC_PROFILE_COMPILE_FLAGS "-pg")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_PROFILE_COMPILE_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()
