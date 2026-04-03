cmake_minimum_required(VERSION 3.10)
project(k2tree LANGUAGES C)

# Use C11 standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# -----------------------------------------------------------------------------
# Compiler Flags and Options
# -----------------------------------------------------------------------------
# Add common warning and safety flags
add_compile_options(-Wall -Wconversion -Wno-sign-conversion -Wtype-limits)

# Handle Build Types (Debug vs Release)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Optional Sanitizers (e.g. -fsanitize=undefined)
# These are often not available or require extra setup on Windows/MinGW
option(USE_SANITIZERS "Enable sanitizers like UBsan" OFF)

if(USE_SANITIZERS)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=undefined")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=undefined")
    message(STATUS "Sanitizers enabled (Debug mode)")
endif()

# Release flags
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG")

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
# Math library for most executables
find_library(MATH_LIB m)

# OpenMP for bbmmult
find_package(OpenMP)

# Threads (used by k2pagerank.c)
find_package(Threads REQUIRED)

# MinGW specific: Enable ANSI stdio for %zu and other modern format specifiers
if(MINGW)
    add_definitions(-D__USE_MINGW_ANSI_STDIO=1)
endif()

# libsais (already has its own CMakeLists.txt)
add_subdirectory(libsais)

# -----------------------------------------------------------------------------
# Core Source Files
# -----------------------------------------------------------------------------
# K2 Core sources (shared across k2 executables)
# Note: k2ops.c includes other .c files (k2text.c, k2aux.c, k2io.c, minimats.c)
set(K2_CORE_SOURCES 
    k2ops.c 
    bbm.c 
    vu64.c 
    pointers.c 
    rank_0000.c
)

# B128 Core sources
set(B128_CORE_SOURCES 
    b128ops.c 
    bbm.c
)

# -----------------------------------------------------------------------------
# Executables - K2 Suite (Compiled with -DK2MAT)
# -----------------------------------------------------------------------------
set(K2_EXECUTABLES 
    k2sparse 
    k2mult 
    k2showinfo 
    k2subtinfo 
    k2cpdf 
    k2unary   # wait, checked original makefile: k2unary.x is in K2EXECS
    k2sum 
    k2pagerank 
    k2tclosure
)

foreach(target ${K2_EXECUTABLES})
    # Create executable with .x extension as requested
    add_executable(${target}.x ${target}.c ${K2_CORE_SOURCES})
    
    # Define K2MAT for these targets
    target_compile_definitions(${target}.x PRIVATE K2MAT)
    
    # Link with libsais and math library
    target_link_libraries(${target}.x PRIVATE libsais m)

    # Specific dependencies
    if(${target} STREQUAL "k2pagerank")
        target_link_libraries(${target}.x PRIVATE Threads::Threads)
    endif()
endforeach()

# -----------------------------------------------------------------------------
# Executables - B128 Suite (Compiled from k2 sources without -DK2MAT)
# -----------------------------------------------------------------------------
set(B128_EXECUTABLES 
    b128sparse 
    b128mult 
    b128showinfo 
    # b128bbm
    b128unary 
    b128tclosure
)

foreach(target ${B128_EXECUTABLES})
    # Original makefile: b128%.o: k2%.c bbm.h b128.h
    # So we strip 'b128' prefix to find the corresponding k2 source file
    string(REGEX REPLACE "^b128" "k2" src_name ${target})
    
    # Create executable with .x extension
    # It uses the k2 source file but links with B128 core and no K2MAT macro
    add_executable(${target}.x ${src_name}.c ${B128_CORE_SOURCES})
    
    # Link with math library
    target_link_libraries(${target}.x PRIVATE m)
endforeach()

# -----------------------------------------------------------------------------
# Special Executables
# -----------------------------------------------------------------------------
# bbmmult.x: uncompressed matrix multiplication using openmp 
if(OpenMP_C_FOUND)
    add_executable(bbmmult.x bbmmult.c bbm.c)
    target_link_libraries(bbmmult.x PRIVATE OpenMP::OpenMP_C)
endif()

# test_matrixcmp.x: simple textual comparison tool
add_executable(test_matrixcmp.x tests/test_matrixcmp.c)
target_link_libraries(test_matrixcmp.x PRIVATE m)

# -----------------------------------------------------------------------------
# Unit Testing (CTest)
# -----------------------------------------------------------------------------
enable_testing()

add_executable(test_unit.x tests/test_unit.c ${K2_CORE_SOURCES})
target_include_directories(test_unit.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(test_unit.x PRIVATE K2MAT)
target_link_libraries(test_unit.x PRIVATE libsais m Threads::Threads)

add_test(NAME BasicUnitTests COMMAND test_unit.x)

# Consistency (Cross-check) test
if(WIN32)
    add_test(NAME ConsistencyCheck 
        COMMAND powershell.exe -ExecutionPolicy Bypass -File "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_consistency.ps1" -BuildDir "${CMAKE_CURRENT_BINARY_DIR}" -RootDir "${CMAKE_CURRENT_SOURCE_DIR}"
    )
else()
    add_test(NAME ConsistencyCheck 
        COMMAND /bin/bash "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_consistency.sh" "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}"
    )
endif()
