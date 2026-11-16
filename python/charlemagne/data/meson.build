# Main package meson.build

# Install Python files
py.install_sources([
    '__init__.py',
], subdir: 'charlemagne/data')

# Install data
install_data('pagn_combined.sp5', install_dir: py.get_install_dir() / 'charlemagne' / 'data')
install_data('pAGN_LICENSE', install_dir: py.get_install_dir() / 'charlemagne' / 'data')
