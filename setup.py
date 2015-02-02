#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

#%% Load packages
from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from setuptools import Extension 
from codecs import open  # To use a consistent encoding
from os import path
import os
import numpy as np
import platform

here = path.abspath(path.dirname(__file__))


#%% Hack to remove option for c++ code
# see http://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
from setuptools.py31compat import get_path, get_config_vars

(opt,) = get_config_vars('OPT')

#print('OPT %s' % opt)

if not opt is None:
    opt = " ".join(    flag for flag in opt.split() if flag != '-Wstrict-prototypes' )
    os.environ['OPT'] = opt

#print('OPT %s' % opt)

#%% Define sources of the package
oadev=0
srcs=[ 'arrayproperties.cpp', 'pareto.cpp', 'nonroot.cpp','mathtools.cpp', 'oaoptions.cpp', 'tools.cpp', 'arraytools.cpp', 'md5.cpp','strength.cpp']

srcs=srcs+[ 'lmc.cpp', 'extend.cpp']	# code used for extension
if os.path.exists('src/oadevelop.cpp') and 1:
  oadev=1
  print('Building development code')
  srcs=[ 'oadevelop.cpp']+srcs

srcs=[ 'src/' + ff for ff in srcs]

if oadev:
  sources = ['oalib_wrap_dev.cxx'] + srcs + ['bitarray/bit_array.cpp']
  swig_opts=[]
else:
  sources = srcs + ['bitarray/bit_array.cpp']

  if 0:
    sources += ['oalib_wrap.cxx'] 
    swig_opts=[]
  else:
    sources += ['oalib.i']
    swig_opts=['-modern', '-c++', '-w503,401,362' , '-Isrc/']
    if platform.system()=='Windows':
      swig_opts+=['-DWIN32', '-D_WIN32']
    
  
if 'VSC_SCRATCH' in os.environ.keys():
  # we are running on the VSC cluster
  zlibdir=os.environ.get('EBROOTZLIB') # SOFTROOTZLIB
  
  libraries=['z']
  library_dirs=[zlibdir + '/lib']
  include_dirs=['.', 'src', np.get_include(), zlibdir + '/include' ]
  oalib_module = Extension('_oalib',
                           sources=sources,
                           include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries
                           )
else:
  libraries=[]
  library_dirs=[]
  include_dirs=['.', 'src', np.get_include() ]


  oalib_module = Extension('_oalib', sources=sources,
                           include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries, swig_opts=swig_opts
                           )


  progs=['oainfo', 'oasplit', 'oacat']
  progs=[]
  pm=[]
  for ii, p in enumerate(progs):
      prog_module = Extension(p, sources=sources+['utils/%s.cpp' % p],
                           include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries )
      pm.append(prog_module)
    

oalib_module.extra_compile_args = ['-DNEWINTERFACE','-DNOOMP', '-DSWIGCODE'] # '-DHAVE_BOOST'

if platform.system()=='Windows':
	oalib_module.extra_compile_args.append('-DWIN32')
	
if os.name=='nt':
  oalib_module.extra_compile_args += [];  

    # for cygwin/mingw ?
  #oalib_module.extra_compile_args += ['-fpermissive', '-std=gnu++11' ];  
else:
  oalib_module.extra_compile_args += ['-O2', '-Wno-unknown-pragmas', '-Wno-sign-compare', '-Wno-return-type' , '-Wno-unused-variable','-Wno-unused-result','-fPIC'];

if oadev:
  oalib_module.extra_compile_args.append('-DOADEV')

print('find_packages: %s' % find_packages() )

data_files=[]
#data_files+=[('', ['oalib/oahelper.py'])]
#data_files+=[('oalib', ['oalib/test.oa'])]
#data_files+=[ ('', ['scripts/example_python_testing.py'])]

scripts=['scripts/example_python_testing.py']

packages=['oapackage']

setup (name = 'OApackage',
       version = '1.9.89',
       author      = "Pieter Eendebak",
       author_email='pieter.eendebak@gmail.com',
	license="BSD",
       url='http://www.pietereendebak.nl/oapage/software.html',
       description = """Python interface to Orthogonal Array library""",
	keywords=[ "orthogonal arrays, design of experiments"],
       ext_modules = [oalib_module],
      # packages=find_packages(exclude=['oahelper']),
       packages=packages,
       #packages=['oalib'],
       #packages=find_packages(exclude=['oahelper']),
        data_files = data_files,
    scripts=scripts,
       py_modules = ['oalib'],	
	requires=['numpy', 'matplotlib'],
	classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Science/Research', 
	      'Programming Language :: Python :: 2',
	      'Programming Language :: Python :: 2.7',
	      'Programming Language :: Python :: 3',
	      'Programming Language :: Python :: 3.4', 
	      ]
       )


