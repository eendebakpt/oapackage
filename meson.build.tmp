# meson compile -C build

# https://github.com/vislearn/libct/blob/main/meson.build
# https://fuchsia.googlesource.com/third_party/dtc/+/refs/heads/upstream/coverity_scan/meson.build
 

project('oapackage', ['c', 'cpp'], version: '0.1.0', default_options: ['cpp_std=c++11'])

pymod = import('python')
python3 = pymod.find_installation('python')
#message(python3)
#gurobi = dependency('gurobi_c++', fallback: ['gurobi-finder', 'gurobi'])
swig = find_program('swig', required: true)
#message(swig)

cc = meson.get_compiler('c')

swig_include_dir = meson.source_root() / 'src'
include_dir = include_directories(['src', 'src/nauty'])

nautyheaders = ['src\\nauty\\gtools.h',
 'src\\nauty\\naugroup.h',
 'src\\nauty\\nautinv.h',
 'src\\nauty\\naurng.h',
 'src\\nauty\\naugraph.h',
 'src\\nauty\\nausparse.h',
 'src\\nauty\\nautil.h',
 'src\\nauty\\nauty.h',
 'src\\nauty\\schreier.h']
 
sources = ['src/arraytools.cpp',
 'src/arrayproperties.cpp',
 'src/pareto.cpp',
 'src/nonroot.cpp',
 'src/mathtools.cpp',
 'src/oaoptions.cpp',
 'src/tools.cpp',
 'src/md5.cpp',
 'src/strength.cpp',
 'src/graphtools.cpp',
 'src/conference.cpp',
 'src/unittests.cpp',
 'src/Deff.cpp',
 'src/evenodd.cpp',
 'src/lmc.cpp',
 'src/extend.cpp',
 'src/bitarray/bit_array.cpp',
 'src/nauty/nauty.c',
 'src/nauty/nautinv.c',
 'src/nauty/nautil.c',
 'src/nauty/naurng.c',
 'src/nauty/naugraph.c',
 'src/nauty/schreier.c',
 'src/nauty/naugroup.c']
 
oapackage_static = static_library(
  'oapackage', sources,
  include_directories: include_dir,
#  dependencies: [gurobi],
  install: false)

#libct_shared = shared_library(
#  'ct', 'lib/ct.cpp',
#  include_directories: include_dir,
#  version: meson.project_version(),
#  soversion: '0',
#  dependencies: [gurobi],
#  install: true)


    
oapackage_py = custom_target('oapackage_py',
  input: ['oalib.i'],
  output: ['oalib_wrap.py'],
  command: [swig, '-python', '-c++', '-doxygen', '-w503,401,362,302,389,446,509,305', '-Isrc/',
               '-I@0@'.format(swig_include_dir), '-o', '@OUTPUT@', '@INPUT0@','-DSWIGCODE', '-DFULLPACKAGE'])

python3.extension_module('oapackage',
  sources: [oapackage_py],
  include_directories: include_dir,
  link_with: [oapackage_static],
  dependencies: [python3.dependency()],
  install: true,
  install_dir: python3.get_install_dir(pure: false, subdir: '.'))
  
#python3.install_sources([
#    'python/ct/__init__.py',
#    'python/ct/gurobi.py',
#    'python/ct/model.py',
#    'python/ct/primals.py',
#    'python/ct/rounding.py',
#    'python/ct/tracker.py',
#    'python/ct/txt.py',
#    'python/ct/utils.py'
#  ], pure: true, subdir: 'ct')

#install_data(['bin/ct'],
#  install_dir: get_option('bindir'),
#  install_mode: 'rwxr-xr-x')

#install_headers('include/ct.h')

# vim: set ts=8 sts=2 sw=2 et: