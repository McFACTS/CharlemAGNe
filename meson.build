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

if host_machine.system() == 'windows'
    add_project_arguments('-D_USE_MATH_DEFINES', language: ['c', 'cpp'])
endif

if host_machine.system() != 'windows'
res1 = run_command('make', 'init', check:false)
res2 = run_command('make', 'hdf5', check:true)
endif

# Include directories
inc_dirs = [
    'src',
    'external/singularity-opac',
    'external/singularity-opac/utils',
    'external/singularity-opac/utils/spiner',
    'external/singularity-opac/utils/ports-of-call',
    'external/cminpack/include',
]
if host_machine.system() != 'windows'
    inc_dirs += 'external/hdf5/include'
endif

incdir = include_directories(inc_dirs)

# HDF5
if host_machine.system() != 'windows'
hdf5_lib_path = meson.current_source_dir() / 'external/hdf5/lib'
hdf5_inc_path = 'external/hdf5/include'
#hdf5_inc = include_directories(hdf5_inc_path)
hdf5_dep = cc.find_library('hdf5',
  dirs: [hdf5_lib_path],
  required: true,
)
hdf5_hl_dep = cc.find_library('hdf5_hl',
  dirs: [hdf5_lib_path],
  required: true,
)
else
# vcpkg: runner sets VCPKG_ROOT; CI must also export VCPKG_DEFAULT_TRIPLET
# matching the `vcpkg install` triplet (x64-windows-static-md / arm64-...)
hdf5_dep = declare_dependency(dependencies: [
    dependency('hdf5', static: true),
    cc.find_library('shlwapi', required: true),
])
hdf5_hl_dep = dependency('hdf5_hl', static: true)
#vcpkg_root = run_command(py, '-c', 'import os; print(os.environ["VCPKG_ROOT"])', check: true).stdout().strip()
#vcpkg_triplet = run_command(py, '-c', 'import os; print(os.environ["VCPKG_DEFAULT_TRIPLET"])', check: true).stdout().strip()
#hdf5_hl_dep = cc.find_library('hdf5_hl', dirs: [vcpkg_root / 'installed' / vcpkg_triplet / 'lib'], required: true)
endif

# Handle cminpack dependency
add_project_arguments('-DSPINER_USE_HDF', '-DCMINPACK_NO_DLL', language: 'cpp')
cminpack_inc = include_directories('external/cminpack/include')
cminpack_lib = static_library('cminpack',
  files(
    'external/cminpack/src/chkder.c',
    'external/cminpack/src/covar1.c',
    'external/cminpack/src/covar.c',
    'external/cminpack/src/dogleg.c',
    'external/cminpack/src/dpmpar.c',
    'external/cminpack/src/enorm.c',
    'external/cminpack/src/fdjac1.c',
    'external/cminpack/src/fdjac2.c',
    'external/cminpack/src/hybrd1.c',
    'external/cminpack/src/hybrd.c',
    'external/cminpack/src/hybrj1.c',
    'external/cminpack/src/hybrj.c',
    'external/cminpack/src/lmder1.c',
    'external/cminpack/src/lmder.c',
    'external/cminpack/src/lmdif1.c',
    'external/cminpack/src/lmdif.c',
    'external/cminpack/src/lmpar.c',
    'external/cminpack/src/lmstr1.c',
    'external/cminpack/src/lmstr.c',
    'external/cminpack/src/qform.c',
    'external/cminpack/src/qrfac.c',
    'external/cminpack/src/qrsolv.c',
    'external/cminpack/src/r1mpyq.c',
    'external/cminpack/src/r1updt.c',
    'external/cminpack/src/rwupdt.c',
  ),
  include_directories: cminpack_inc,
  c_args: ['-DCMINPACK_NO_DLL'],
)
cminpack_dep = declare_dependency(link_with: cminpack_lib, include_directories: cminpack_inc)

# Install Python package
subdir('python/charlemagne')
