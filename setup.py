#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

#%% Load packages
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from setuptools import Extension
from setuptools.command.test import test as TestCommand

from codecs import open  # To use a consistent encoding
from os import path
import os
import sys
import platform

try:
    import numpy as np
except ImportError:
    raise RuntimeError(
        "numpy cannot be imported. numpy must be installed "
        "prior to installing OApackage")

npinclude = np.get_include()

here = path.abspath(path.dirname(__file__))

#%%


def checkZlib(verbose=0):
    """ Check whether zlib is available 

    Code adapted from http://stackoverflow.com/questions/28843765/setup-py-check-if-non-python-library-dependency-exists    
    """
    ret_val = True
    try:
        import os
        import distutils.ccompiler
        import distutils.sysconfig
        import tempfile
        from distutils.errors import CompileError, LinkError

        # create test file...
        c_code = """
            #include <zlib.h>
            #include <stdio.h>
        
            int main(int argc, char* argv[])
            {
                printf("Hello zlib test...\\n");
        
                return 0;
            }
            """
        tmp_dir = tempfile.mkdtemp(prefix='helloztest_')
        bin_file_name = os.path.join(tmp_dir, 'helloz')
        file_name = bin_file_name + '.c'

        #print('tmp_dir: %s'  % tmp_dir)

        with open(file_name, 'w') as fp:
            fp.write(c_code)

        # compile and link code
        compiler = distutils.ccompiler.new_compiler()
        assert isinstance(compiler, distutils.ccompiler.CCompiler)
        distutils.sysconfig.customize_compiler(compiler)

        libraries = ['z']
        ret_val = True
        try:
            if verbose:
                print('checkZlib: compile and link')
            compiler.link_executable(
                compiler.compile([file_name]), bin_file_name, libraries=libraries, )
        except CompileError as e:
            if verbose:
                print('helloz compile error')
            ret_val = False
        except LinkError as e:
            if verbose:
                print('helloz link error')
            ret_val = False
        except Exception as e:
            if verbose:
                print('helloz  error')
                print(e)
            ret_val = False
    except Exception as e:
        ret_val = False

    return ret_val


#%%

import re


def get_version_info(verbose=0):
    """ Extract version information from source code """

    GIT_REVISION = None
    try:
        if os.path.exists('src/version.h'):
            with open('src/version.h') as f:
                ln = f.readline()
                # print(ln)
                m = re.search('.* "(.*)"', ln)
                FULLVERSION = (m.group(1))
        else:
            FULLVERSION = '0.0'
    except Exception as E:
        FULLVERSION = '0.0'
    if verbose:
        print('get_version_info: %s' % FULLVERSION)
    return FULLVERSION, GIT_REVISION

# print(get_version_info())


#%% Hack to remove option for c++ code
try:
        # see
        # http://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
    from setuptools.py31compat import get_path, get_config_vars

    (opt,) = get_config_vars('OPT')

    if not opt is None:
        opt = " ".join(flag for flag in opt.split()
                       if flag != '-Wstrict-prototypes')
        os.environ['OPT'] = opt
except:
    import setuptools
    print('old version of setuptools: %s' % setuptools.__version__)
    pass


#%% Test suite

class OATest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        # New setuptools don't need this anymore, thus the try block.
        try:
            self.test_args = []
            self.test_suite = True
        except AttributeError:
            pass

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        print('## oapackage test: load package')
        import oapackage
        print('## oapackage test: oalib version %s' % oapackage.version())
        print('## oapackage test: package compile options\n%s\n' % oapackage.oalib.compile_information())

        oapackage.unittest(verbose=1)
        if 0:
            ii = 0
            al = oapackage.exampleArray(ii, 0)
            Deff = al.Defficiency()
            print('## oapackage test: example array %d: Deff %.3f' % (ii, Deff))
        errno = 0
        #errno = pytest.main(self.pytest_args)
        sys.exit(errno)


#%% Define sources of the package
oadev = 0
srcs = ['arraytools.cpp', 'arrayproperties.cpp', 'pareto.cpp', 'nonroot.cpp',
        'mathtools.cpp', 'oaoptions.cpp', 'tools.cpp', 'md5.cpp', 'strength.cpp']
srcs = srcs + ['Deff.cpp', 'evenodd.cpp']
srcs = srcs + ['conference.cpp']
#srcs=srcs+[ 'lmc.h', 'Deff.h', 'mathtools.h', 'tools.h', 'arraytools.h' ]

srcs = srcs + ['lmc.cpp', 'extend.cpp']  # code used for extension
srcs = ['src/' + ff for ff in srcs]
if os.path.exists('dev/oadevelop.cpp') and 1:
    oadev = 1
    print('Building development code')
    srcs = ['dev/oadevelop.cpp'] + srcs

sources = srcs + ['src/bitarray/bit_array.cpp']
swig_opts = []

compile_options = []

if oadev:
    sources = ['oalib.i'] + sources
    swig_opts += ['-modern', '-c++', '-w503,401,362,509,389',
                  '-Isrc/', '-Idev/'] 
    compile_options += ['-DSWIGCODE', '-DFULLPACKAGE', '-DOADEV', '-Idev/']
    swig_opts += ['-DSWIGCODE', '-DFULLPACKAGE', '-DOADEV']
else:
    sources = ['oalib.i'] + sorted(sources)
    swig_opts += ['-modern', '-c++',
                  '-w503,401,362,302,389,446,509,305', '-Isrc/']
    compile_options += ['-DSWIGCODE', '-DFULLPACKAGE']
    swig_opts += ['-DSWIGCODE', '-DFULLPACKAGE']

if 1:
    swig_opts += ['-Isrc/nauty/']
    compile_options += ['-Isrc/nauty/']

    sources += ['src/graphtools.cpp']

    # nauty/gtools.c
    for f in 'nauty/nauty.c nauty/nautinv.c nauty/nautil.c nauty/naurng.c nauty/naugraph.c nauty/schreier.c nauty/naugroup.c'.split(' '):
        sources += ['src/' + f]

if platform.system() == 'Windows':
    compile_options += ['-DWIN32', '-D_WIN32']
    swig_opts += ['-DWIN32', '-D_WIN32']


if 'VSC_SCRATCH' in os.environ.keys():
    # we are running on the VSC cluster
    zlibdir = os.environ.get('EBROOTZLIB')  # SOFTROOTZLIB

    libraries = ['z']
    library_dirs = [zlibdir + '/lib']
    include_dirs = ['.', 'src', npinclude, zlibdir + '/include']
    oalib_module = Extension('_oalib',
                             sources=sources,
                             include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries, swig_opts=swig_opts
                             )
else:
    libraries = []
    library_dirs = []
    include_dirs = ['.', 'src', npinclude]

    oalib_module = Extension('_oalib', sources=sources,
                             include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries, swig_opts=swig_opts
                             )
    progs = ['oainfo', 'oasplit', 'oacat']
    progs = []
    pm = []
    for ii, p in enumerate(progs):
        prog_module = Extension(p, sources=sources + ['utils/%s.cpp' % p],
                                include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries)
        pm.append(prog_module)

compile_options += ['-DNOOMP']
swig_opts += ['-DNOOMP']

# ['-DNOOMP', '-DSWIGCODE', '-DFULLPACKAGE'] # '-DHAVE_BOOST'
oalib_module.extra_compile_args = compile_options

if checkZlib(verbose=0):
    if platform.system() == 'Windows':
        pass
    else:
        zlibflag = '-DUSEZLIB'
        oalib_module.extra_compile_args += [zlibflag]
        swig_opts += [zlibflag]
        oalib_module.extra_link_args += ['-lz']
else:
    zlibflag = '-DNOZLIB'
    oalib_module.extra_compile_args += [zlibflag]
    swig_opts += [zlibflag]

if os.name == 'nt':
    oalib_module.extra_compile_args += []

    # for cygwin/mingw ?
    #oalib_module.extra_compile_args += ['-fpermissive', '-std=gnu++11' ];
else:
    oalib_module.extra_compile_args += ['-O3', '-Wno-unknown-pragmas', '-Wno-sign-compare',
                                        '-Wno-return-type', '-Wno-unused-variable', '-Wno-unused-result', '-fPIC']
    oalib_module.extra_compile_args += ['-Wno-date-time',]
    #swig_opts += [ '-Wno-delete-non-virtual-dtor' ]

if platform.node() == 'marmot' or platform.node() == 'goffer' or platform.node() == 'pte':
    # openmp version of code
    oalib_module.extra_compile_args += ['-fopenmp', '-DDOOPENMP']
    oalib_module.extra_link_args += ['-fopenmp']

print('find_packages: %s' % find_packages())
#print('swig_opts: %s' % str(swig_opts) )

data_files = []
scripts = ['scripts/example_python_testing.py']
packages = ['oapackage']

# fix from:
# http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module

from distutils.command.build import build
from setuptools.command.install import install


# see: http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module
class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)
        # self.run_command('install')
        # self.do_egg_install()

class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        install.run(self)
        # self.run_command('install')
        #self.do_egg_install()

# PyPi does not support markdown....
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md', 'rt').read()

version = get_version_info()[0]
print('OApackage: version %s'  % version)

setup(name='OApackage',
      #cmdclass = {'test': OATest },
      cmdclass={'test': OATest, 'install': CustomInstall, 'build': CustomBuild},
      version=version,
      author="Pieter Eendebak",
      description="Package to generate and analyse orthogonal arrays and optimal designs",
      long_description=long_description,
      author_email='pieter.eendebak@gmail.com',
      license="BSD",
      url='http://www.pietereendebak.nl/oapackage/index.html',
      keywords=["orthogonal arrays, design of experiments"],
      ext_modules=[oalib_module],
      py_modules=['oalib'],
      # packages=find_packages(exclude=['oahelper']),
      packages=packages,
      data_files=data_files,
      test_suite="oapackage.unittest",
      scripts=scripts,
      # nose and coverage are only for tests, but we'd like to encourage
      # people to run tests!
      tests_require=['numpy', 'nose>=1.3', 'coverage>=4.0'],
      zip_safe=False,
      install_requires=['numpy>=1.10'],
      requires=['numpy', 'matplotlib'],
      classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5'
                   ]
      )
