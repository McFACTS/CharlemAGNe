project('charlemagne',
    'cpp',
    version: '0.1.0',
    meson_version: '>=1.3.0',
    default_options: ['cpp_std=c++20', 'warning_level=2', 'optimization=3'],
)

py = import('python').find_installation(pure: false)
numpy_dep = dependency('numpy')
hdf5_dep = dependency('hdf5', static: true)
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_19_API_VERSION', language: 'cpp')

# Include directories
inc_dirs = [
    'src',
    'external/singularity-opac',
    'external/singularity-opac/utils',
    'external/singularity-opac/utils/spiner',
    'external/singularity-opac/utils/ports-of-call',
    'external/cminpack/include',
]

incdir = include_directories(inc_dirs)

# Handle cminpack dependency
add_project_arguments('-DSPINER_USE_HDF', '-DCMINPACK_NO_DLL', language: 'cpp')
cminpack_inc = include_directories('external/cminpack/include')
cminpack_lib = static_library('cminpack',
  files(
  ),
  include_directories: cminpack_inc,
  c_args: ['-DCMINPACK_NO_DLL'],
)
cminpack_dep = declare_dependency(link_with: cminpack_lib, include_directories: cminpack_inc)

# Install Python package
subdir('python/charlemagne')
