project('charlemagne',
    'c',
    'cpp',
    version: '0.1.1',
    meson_version: '>=1.3.0',
    default_options: ['cpp_std=c++20', 'warning_level=2', 'optimization=3', 'c_std=c11'],
)

py = import('python').find_installation(pure: false)
cc = meson.get_compiler('c')
numpy_dep = dependency('numpy')
#hdf5_dep = dependency('hdf5', static: true)
#hdf5_hl_dep = dependency('hdf5_hl', static: true)

add_project_arguments('-DNPY_TARGET_VERSION=NPY_1_19_API_VERSION', language: 'cpp')

if host_machine.system() != 'windows'
res1 = run_command('make', 'init')
res2 = run_command('make', 'hdf5')
endif

# Include directories
inc_dirs = [
    'src',
    'external/singularity-opac',
    'external/singularity-opac/utils',
    'external/singularity-opac/utils/spiner',
    'external/singularity-opac/utils/ports-of-call',
    'external/cminpack/include',
    'external/hdf5/include',
]

incdir = include_directories(inc_dirs)

# HDF5
hdf5_lib_path = meson.current_source_dir() / 'external/hdf5/lib'
hdf5_inc_path = 'external/hdf5/include'
hdf5_inc = include_directories(hdf5_inc_path)
hdf5_dep = cc.find_library('hdf5',
  dirs: [hdf5_lib_path],
  required: true,
)
hdf5_hl_dep = cc.find_library('hdf5_hl',
  dirs: [hdf5_lib_path],
  required: true,
)

# Handle cminpack dependency
add_project_arguments('-DSPINER_USE_HDF', '-DCMINPACK_NO_DLL', language: 'cpp')
cminpack_inc = include_directories('external/cminpack/include')
cminpack_lib = static_library('cminpack',
  run_command(py, '-c', 'import glob; print(*glob.glob("external/cminpack/src/*.c"), sep="\\n")'),
  include_directories: cminpack_inc,
  c_args: ['-DCMINPACK_NO_DLL'],
)
cminpack_dep = declare_dependency(link_with: cminpack_lib, include_directories: cminpack_inc)

# Install Python package
subdir('python/charlemagne')
