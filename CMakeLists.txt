cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(cassini-parameterization
        VERSION 0.1.0
        LANGUAGES Fortran CXX
        DESCRIPTION "Pashkevich Cassini-oval nuclear shape parameterization"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Fetch GCC compiler options package ---
FetchContent_Declare(
        gcc-compiler-options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/gcc-compiler-options.git
        GIT_TAG 2.0.0
)
FetchContent_MakeAvailable(gcc-compiler-options)

# --- Fetch foundational Fortran modules ---
FetchContent_Declare(
        fortran-foundations
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/fortran-foundations.git
        GIT_TAG 3.0.0
)
FetchContent_MakeAvailable(fortran-foundations)

# --- Preprocessor & build type ---
set(CMAKE_Fortran_PREPROCESS ON)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "RelWithDebInfo" "Release")
endif ()

# --- Compiler flags via gcc-compiler-options ---
include(GCCCompilerOptions/FortranCompilerOptions)
create_fortran_library_interface(TARGET cassini_parameterization_flags)

# --- Source files ---
set(CASSINI_PARAM_SOURCES
        src/cassini_parameterization_workers_mod.f08
        src/cassini_parameterization_mod.f08
        src/cassini_parameterization_c_api_mod.f08
)

# --- No LTO for the library's own objects ---
# One-shot == cached and repeat == first call are BITWISE claims. The
# arithmetic lives in the workers module and the tiers in the
# front module; without LTO no kernel can be inlined into a tier routine, so
# every tier path executes the same machine code. With -flto=auto -ffast-math
# a kernel is inlined into each caller and optimized per call site, and the
# same call then returns different last bits from different call sites. A
# source-file property lands after the gcc-opts interface flags, so it wins.
# Consumers keep LTO for their own code.
set_source_files_properties(${CASSINI_PARAM_SOURCES} PROPERTIES COMPILE_OPTIONS "-fno-lto")

# --- Static library target ---
add_library(cassini_parameterization STATIC ${CASSINI_PARAM_SOURCES})
target_link_libraries(cassini_parameterization
        PUBLIC fortran_foundations
        PRIVATE cassini_parameterization_flags
)
target_compile_definitions(cassini_parameterization PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(cassini_parameterization PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
target_include_directories(cassini_parameterization
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(CassiniParameterization::cassini_parameterization ALIAS cassini_parameterization)

# --- Shared library target (for FFI consumers: Python ctypes, etc.) ---
add_library(cassini_parameterization_shared SHARED ${CASSINI_PARAM_SOURCES})
target_link_libraries(cassini_parameterization_shared
        PUBLIC fortran_foundations
        PRIVATE cassini_parameterization_flags
)
target_compile_definitions(cassini_parameterization_shared PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(cassini_parameterization_shared PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared
        OUTPUT_NAME cassini_parameterization
)
target_include_directories(cassini_parameterization_shared
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(CassiniParameterization::cassini_parameterization_shared
        ALIAS cassini_parameterization_shared)

# --- C++ interface target (consumes the shared library via the C header) ---
add_library(cassini_parameterization_cxx INTERFACE)
target_include_directories(cassini_parameterization_cxx INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(cassini_parameterization_cxx INTERFACE
        cassini_parameterization_shared)
target_compile_features(cassini_parameterization_cxx INTERFACE cxx_std_20)
add_library(CassiniParameterization::cassini_parameterization_cxx
        ALIAS cassini_parameterization_cxx)

# --- Install the shared library into the Python package (wheel builds) ---
install(TARGETS cassini_parameterization_shared
        LIBRARY DESTINATION cassini_parameterization)

# --- Tests (skipped for wheel builds via -DCASSINI_PARAM_BUILD_TESTS=OFF) ---
option(CASSINI_PARAM_BUILD_TESTS "Build the Fortran/C++ test suite" ON)
if (CASSINI_PARAM_BUILD_TESTS)
    enable_testing()

add_library(cassini_param_test_utils STATIC
        tests/test_utils_mod.f08
        tests/cassini_test_fixtures_mod.f08
        tests/cassini_test_reference_mod.f08
)
target_link_libraries(cassini_param_test_utils
        PUBLIC cassini_parameterization fortran_foundations
        PRIVATE cassini_parameterization_flags)
set_target_properties(cassini_param_test_utils PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests)
target_include_directories(cassini_param_test_utils PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests>)

set(CASSINI_PARAM_TEST_SUITES status epsilon map)
foreach (suite IN LISTS CASSINI_PARAM_TEST_SUITES)
    add_executable(cassini_param_${suite}_test tests/cassini_param_${suite}_test.f08)
    target_link_libraries(cassini_param_${suite}_test PRIVATE
            cassini_parameterization cassini_param_test_utils cassini_parameterization_flags)
    set_target_properties(cassini_param_${suite}_test PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests/${suite})
    add_test(NAME ${suite} COMMAND cassini_param_${suite}_test)
endforeach ()
endif ()
