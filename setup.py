#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

#%% Load packages
from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from setuptools import Extension 
from setuptools.command.test import test as TestCommand

from codecs import open  # To use a consistent encoding
from os import path
import os,sys
import numpy as np
import platform

here = path.abspath(path.dirname(__file__))


#%% Hack to remove option for c++ code
try:
	# see http://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
	from setuptools.py31compat import get_path, get_config_vars

	(opt,) = get_config_vars('OPT')

	#print('OPT %s' % opt)

	if not opt is None:
		opt = " ".join( flag for flag in opt.split() if flag != '-Wstrict-prototypes' )
		os.environ['OPT'] = opt
except:
	import setuptools
	print('old version of setuptools: %s'  % setuptools.__version__ )
	pass
#print('OPT %s' % opt)


#%% Test suite

class OATest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        print('## oapackage test: load package' )
        import oapackage
        print('## oapackage test: oalib version %s' % oapackage.version() )
        
        oapackage.unittest(verbose=1)
        if 0:
            ii=0
            al=oapackage.exampleArray(ii, 0)
            Deff=al.Defficiency()
            print('## oapackage test: example array %d: Deff %.3f' % (ii, Deff))
        errno=0
        #errno = pytest.main(self.pytest_args)
        sys.exit(errno)


#%% Define sources of the package
oadev=0
srcs=[ 'arraytools.cpp', 'arrayproperties.cpp', 'pareto.cpp', 'nonroot.cpp','mathtools.cpp', 'oaoptions.cpp', 'tools.cpp',  'md5.cpp','strength.cpp']
srcs=srcs+[ 'Deff.cpp' ]
#srcs=srcs+[ 'lmc.h', 'Deff.h', 'mathtools.h', 'tools.h', 'arraytools.h' ]

srcs=srcs+[ 'lmc.cpp', 'extend.cpp']	# code used for extension
srcs=[ 'src/' + ff for ff in srcs]
if os.path.exists('dev/oadevelop.cpp') and 0:
  oadev=1
  print('Building development code')
  srcs=[ 'dev/oadevelop.cpp']+srcs

sources =   srcs + ['src/bitarray/bit_array.cpp']
swig_opts=[]

if oadev:
  #sources = ['oalib_wrap.cxx'] + srcs + ['bitarray/bit_array.cpp']

  sources = ['oalib.i'] + sources
  swig_opts+=['-modern', '-DSWIGCODE', '-DFULLPACKAGE', '-DOADEV', '-c++', '-w503,401,362' , '-Isrc/', '-Idev/'] # , '-o oalib_wrap_dev.cxx']
else:
  if 0:
    sources += ['oalib_wrap.cxx'] 
  else:
    sources = ['oalib.i'] + sorted(sources)
    swig_opts+=['-modern', '-DSWIGCODE', '-DFULLPACKAGE',  '-c++', '-w503,401,362,302,389,446,509,305' , '-Isrc/']

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
                           include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries, swig_opts=swig_opts
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
    

oalib_module.extra_compile_args = ['-DNOOMP', '-DSWIGCODE', '-DFULLPACKAGE'] # '-DHAVE_BOOST'

if platform.system()=='Windows':
	oalib_module.extra_compile_args.append('-DWIN32')
	
if os.name=='nt':
  oalib_module.extra_compile_args += [];  

    # for cygwin/mingw ?
  #oalib_module.extra_compile_args += ['-fpermissive', '-std=gnu++11' ];  
else:
  oalib_module.extra_compile_args += ['-O2', '-Wno-unknown-pragmas', '-Wno-sign-compare', '-Wno-return-type' , '-Wno-unused-variable','-Wno-unused-result','-fPIC'];


if platform.node()=='marmot' or  platform.node()=='goffer' or platform.node()=='pte':
  # openmp version of code
  oalib_module.extra_compile_args+=['-fopenmp', '-DDOOPENMP']
  oalib_module.extra_link_args+=['-fopenmp']
      
if oadev:
  oalib_module.extra_compile_args.append('-DOADEV')

print('find_packages: %s' % find_packages() )
#print('swig_opts: %s' % str(swig_opts) )

data_files=[]
#data_files+=[('', ['oalib/oahelper.py'])]
#data_files+=[('oalib', ['oalib/test.oa'])]
#data_files+=[ ('', ['scripts/example_python_testing.py'])]

scripts=['scripts/example_python_testing.py']

packages=['oapackage']

# fix from:
# http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module

#from distutils.command.build import build
from setuptools.command.install import install

#class CustomBuild(build):
#    def run(self):
#        self.run_command('build_ext')
#        build.run(self)


class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        install.run(self)
        #self.run_command('install')
        #self.do_egg_install()


#setup(
#    cmdclass={'build': CustomBuild, 'install': CustomInstall},
#    py_modules=['module1'],
#    ext_modules=[module1]
#)

# PyPi does not support markdown....
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md', 'rt').read()
    
setup (name = 'OApackage',
      #cmdclass = {'test': OATest },
      cmdclass = {'test': OATest, 'install': CustomInstall},
       version = '2.0.24',
       author      = "Pieter Eendebak",
       description = "Package to generate and analyse orthogonal arrays and optimal designs",
       long_description=long_description,
       author_email='pieter.eendebak@gmail.com',
	license="BSD",
       url='http://www.pietereendebak.nl/oapackage/index.html',
#       description = """Python interface to Orthogonal Array library""",
	keywords=[ "orthogonal arrays, design of experiments"],
       ext_modules = [oalib_module],
       py_modules = ['oalib'],	
      # packages=find_packages(exclude=['oahelper']),
       packages=packages,
       #packages=['oalib'],
       #packages=find_packages(exclude=['oahelper']),
        data_files = data_files,
        test_suite = "oapackage.unittest",
    scripts=scripts,
    tests_require=['numpy'],
       zip_safe=False,
	requires=['numpy', 'matplotlib'],
	classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Science/Research', 
	      'Programming Language :: Python :: 2',
	      'Programming Language :: Python :: 2.7',
	      'Programming Language :: Python :: 3',
	      'Programming Language :: Python :: 3.3',
	      'Programming Language :: Python :: 3.4', 
	      ]
       )


