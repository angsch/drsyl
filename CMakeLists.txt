cmake_minimum_required(VERSION 3.16)

project(drsyl
    DESCRIPTION "Robust Task-Parallel Solution of the Triangular Sylvester Equation"
    LANGUAGES C)

# Set C standard to gnu99
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

set(CMAKE_C_FLAGS "-Wall -O3")

set(CMAKE_VERBOSE_MAKEFILE ON)

if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -fprefetch-loop-arrays -LNO:prefetch")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops")
else()
    message(WARNING "Untested compiler.")
endif()

add_compile_definitions(ALIGNMENT=64)

include(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -unroll-aggressive")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    endif()
endif()

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(drsyl main.c majorant.c norm.c partition.c random.c reference.c robust.c
    syl-task.c sylvester.c timing.c update-task.c utils.c validation.c)
target_link_libraries(drsyl ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_C_LIBRARIES} m)
