py.extension_module(
  '_abi',
  sources: ['_abi.cpp'],
  include_directories: incdir,
  dependencies: [numpy_dep, hdf5_dep, hdf5_hl_dep, cminpack_dep],
  limited_api: '3.9',
  subdir: 'charlemagne',
  install: true
)
