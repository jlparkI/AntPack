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
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()


project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
#set(QT_DEBUG_FIND_PACKAGE ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_AUTOUIC ON)
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTORCC ON)
# Extremely annoying but for some reason MSVC requires an additional flag
# to be set.
if (MSVC)
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


# Now add requirements for Qt.
#find_package(Qt6 COMPONENTS Widgets PrintSupport Core REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
INCLUDE_DIRECTORIES(src/antpack/antpack_gui)
#INCLUDE_DIRECTORIES(${Qt6Widgets_INCLUDE_DIRS})



nanobind_add_module(
  antpack_cpp_ext

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  # STABLE_ABI  #Removing this since it leads to compatibilit issues.

  NB_STATIC

    src/antpack/cpp_src/ant_ext.cpp
    src/antpack/cpp_src/annotator_classes/vj_aligner.cpp
    src/antpack/cpp_src/annotator_classes/ig_aligner.cpp
    src/antpack/cpp_src/annotator_classes/annotator_base_class.cpp
    src/antpack/cpp_src/annotator_classes/single_chain_annotator.cpp
    src/antpack/cpp_src/annotator_classes/paired_chain_annotator.cpp

    src/antpack/cpp_src/vj_assignment/vj_match_counter.cpp

    src/antpack/cpp_src/developability_calcs/liability_search_tool.cpp

    src/antpack/cpp_src/dna_read_handling/dna_sequence_handling.cpp

    src/antpack/cpp_src/utilities/utilities.cpp
    src/antpack/cpp_src/utilities/cdr_assignment_utilities.cpp
    src/antpack/cpp_src/utilities/consensus_file_utilities.cpp
    src/antpack/cpp_src/utilities/dna_sequence_utilities.cpp

    src/antpack/cpp_src/humanness_calcs/responsibility_calcs.cpp
    src/antpack/cpp_src/annotator_classes/prefiltering_tool.cpp

)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(antpack_cpp_ext PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS antpack_cpp_ext
    LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})
