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
]

incdir = include_directories(inc_dirs)

# Handle cminpack dependency
cmake = import('cmake') # OK we are using cmake from meson I guess
# Seems to initialize an options dict
cminpack_opts = cmake.subproject_options()
cminpack_opts.add_cmake_defines({'BUILD_EXAMPLES': 'OFF', 'CMINPACK_PRECISION': 'd'})
cminpack_opts.set_install(false) # Don't install cminpack into libdir_static
# Set up the targets for meson
cminpack_proj = cmake.subproject('cminpack', options: cminpack_opts)
# We depend on this!
cminpack_dep = cminpack_proj.dependency('cminpack')

# Install Python package
subdir('python/charlemagne')
