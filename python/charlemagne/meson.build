# Main package meson.build

# Install Python files
py.install_sources([
    '__init__.py',
    '__version__.py',
    'disk.py',
], subdir: 'charlemagne')

# Install data
install_data('data/pagn_combined.sp5', install_dir: py.get_install_dir() / 'charlemagne' / 'data')

# Process subdirectories
subdir('abi')
