include(CheckLanguage)

set(CMAKE_VERBOSE_MAKEFILE on)
cmake_minimum_required(VERSION 3.18)
if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -v .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()


project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_STANDARD 99)


# Extremely annoying but MSVC requires an additional flag to be set.
if (MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2")
    add_compile_options("/Zc:__cplusplus" "/permissive-")
endif()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#########
#Generate a compile commands file for use by linters.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

##########
#Add non-nanobind submodules.
add_subdirectory(extern/abseil-cpp)


nanobind_add_module(
  antpack_cpp_ext

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  # STABLE_ABI  #Removing this since it leads to compatibility issues.

  NB_STATIC

    extern/antpack_private/src/ant_ext.cpp
    extern/antpack_private/src/annotator_classes/vj_aligner.cpp
    extern/antpack_private/src/annotator_classes/ig_aligner.cpp
    extern/antpack_private/src/annotator_classes/annotator_base_class.cpp
    extern/antpack_private/src/annotator_classes/single_chain_annotator.cpp
    extern/antpack_private/src/annotator_classes/paired_chain_annotator.cpp

    extern/antpack_private/src/vj_assignment/vj_match_counter.cpp

    extern/antpack_private/src/developability_calcs/liability_search_tool.cpp

    extern/antpack_private/src/utilities/utilities.cpp
    extern/antpack_private/src/utilities/cdr_assignment_utilities.cpp
    extern/antpack_private/src/utilities/consensus_file_utilities.cpp
    extern/antpack_private/src/utilities/dna_sequence_utilities.cpp

    extern/antpack_private/src/annotator_classes/prefiltering_tool.cpp

    extern/antpack_private/src/categorical_mixtures/base_categorical_mixture_class.cpp
    extern/antpack_private/src/categorical_mixtures/em_categorical_mixture.cpp
    extern/antpack_private/src/categorical_mixtures/em_update_calcs.cpp
    extern/antpack_private/src/categorical_mixtures/weight_responsibility_calcs.cpp
    extern/antpack_private/src/categorical_mixtures/mixture_utilities.cpp
    extern/antpack_private/src/categorical_mixtures/test_wrappers.cpp

    extern/antpack_private/src/data_management/template_aligner.cpp

    extern/antpack_private/src/license_key_handling/key_check_routines.cpp
    extern/antpack_private/src/license_key_handling/pub_private_utilities.cpp

    extern/antpack_private/src/database_tools/local_database_construction_tool.cpp
    extern/antpack_private/src/database_tools/local_database_management_tool.cpp
    extern/antpack_private/src/database_tools/core_search_utilities.cpp
    extern/antpack_private/src/database_tools/core_sqlite_search_functions.cpp
    extern/antpack_private/src/database_tools/database_utilities.cpp
    extern/antpack_private/src/database_tools/sqlite3.c

    extern/antpack_private/src/license_key_handling/crypto_support/add_scalar.c
    extern/antpack_private/src/license_key_handling/crypto_support/fe.c
    extern/antpack_private/src/license_key_handling/crypto_support/ge.c
    extern/antpack_private/src/license_key_handling/crypto_support/key_exchange.c
    extern/antpack_private/src/license_key_handling/crypto_support/keypair.c
    extern/antpack_private/src/license_key_handling/crypto_support/sc.c
    extern/antpack_private/src/license_key_handling/crypto_support/seed.c
    extern/antpack_private/src/license_key_handling/crypto_support/sha512.c
    extern/antpack_private/src/license_key_handling/crypto_support/sign.c
    extern/antpack_private/src/license_key_handling/crypto_support/verify.c


)


# Ordinarily these should be disabled. Uncomment these if we
# want profiling.
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(antpack_cpp_ext PUBLIC OpenMP::OpenMP_CXX)
    if (MSVC)
	    list(APPEND OpenMP_CXX_FLAGS "-openmp:experimental")
    endif()
    target_compile_options(antpack_cpp_ext PRIVATE ${OpenMP_CXX_FLAGS})
else()
    message( FATAL_ERROR "OpenMP is required and was not found.")
endif()


# Add non-nanobind submodules
target_link_libraries(antpack_cpp_ext PRIVATE absl::base absl::flat_hash_set)

# Target compile flags specific for sqlite.
target_compile_definitions(antpack_cpp_ext PRIVATE
    "SQLITE_DEFAULT_MEMSTATUS=0"
    "SQLITE_OMIT_DECLTYPE"
    "SQLITE_OMIT_PROGRESS_CALLBACK"
    "SQLITE_OMIT_SHARED_CACHE"
)


install(TARGETS antpack_cpp_ext
    LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})
