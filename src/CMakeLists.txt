
# CharlemAGNe executable

set(EXE_NAME "charlemagne")

AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR} EXE_NAME_SRC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# generate version.cpp file from current git commit in a way that
# 1. Re-generates if hash changes
# 2. Does not re-build everything when hash changes
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/version.cpp" @ONLY)
list(APPEND EXE_NAME_SRC "${CMAKE_CURRENT_BINARY_DIR}/version.cpp" version.hpp)

add_executable(${EXE_NAME} ${EXE_NAME_SRC})

target_link_libraries(${EXE_NAME} 
PUBLIC
    singularity-opac::flags
    cminpack
    hdf5::hdf5
    HDF5::HDF5
)
