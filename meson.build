project('charlemagne',
    'cpp',
    version: '0.1.0',
    cpp_std: 'c++20',
    meson_version: '>=1.3.0',
    default_options: ['cpp_std=c++20', 'warning_level=2', 'optimization=3'],
)

py = import('python').find_installation(pure: false)
numpy_dep = dependency('numpy')
hdf5_dep = dependency('hdf5', static: true)
add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_19_API_VERSION', language: 'cpp')

# Include directories
inc_dirs = [
    'src/charlemagne',
]

incdir = include_directories(inc_dirs)

# Install Python package
subdir('python')
