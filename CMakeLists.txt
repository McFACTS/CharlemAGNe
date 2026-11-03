# CharlemAGNe CMake top file
# General options, not target specifics

cmake_minimum_required(VERSION 3.16)
project(kharma LANGUAGES C CXX)

option(CHARLES_BUILD_TESTS "Build unit tests" OFF)

# We follow Parthenon in requiring C++17 going forward
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX20_STANDARD_COMPILE_OPTION "-std=c++20")

# Set the path to include cmake/ dir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# For including the current git revision in the exe
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe_working_tree(GIT_VERSION --tags)

# Singularity-opac options
set(SINGULARITY_USE_HDF5 ON CACHE BOOL "" FORCE)
# Build Singularity-opac
add_subdirectory(external/singularity-opac)

# cminpack: C port of MINPACK (for lmdif/lmder) rootfinders 
# (consistent with SciPy/pAGN).
add_subdirectory(external/cminpack)

# Build CharlemAGNe
add_subdirectory(src)

if (CHARLES_BUILD_TESTS)
   enable_testing()
   add_subdirectory(tst/unit)
endif()
