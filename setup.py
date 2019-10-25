#!/usr/bin/env python

"""
setup.py file for OApackage
"""

# %% Load packages
from setuptools import setup, find_packages
from setuptools import Extension
from setuptools.command.test import test as TestCommand
from distutils.command.build import build as distutils_build
from setuptools.command.install import install as setuptools_install
import setuptools.command.build_ext

from os import path
import os
import sys
import logging
import platform
import subprocess
import re

try:
    import numpy as np
except ImportError:
    raise RuntimeError(
        "numpy cannot be imported. numpy must be installed "
        "prior to installing OApackage")

npinclude = np.get_include()
setup_directory = path.abspath(path.dirname(__file__))

is_python3 = sys.version_info >= (3, 4)

# %%


def checkZlib(verbose=0):
    """ Check whether zlib is available

    Code adapted from http://stackoverflow.com/questions/28843765/setup-py-check-if-non-python-library-dependency-exists
    """
    ret_val = True
    try:
        import distutils.ccompiler
        import distutils.sysconfig
        import tempfile
        from distutils.errors import CompileError, LinkError

        # create test file...
        c_code = """
            #include <zlib.h>
            #include <stdio.h>

            int main(int argc, char* argv[]) {
                printf("Hello zlib test...\\n");
                return 0;
            }
            """
        tmp_dir = tempfile.mkdtemp(prefix='helloztest_')
        bin_file_name = os.path.join(tmp_dir, 'helloz')
        file_name = bin_file_name + '.c'

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
                print('checkZlib: compile error in %s, zlib not available' % file_name)
            ret_val = False
        except LinkError as e:
            if verbose:
                print('checkZlib: link error in %s, zlib not available' % file_name)
            ret_val = False
        except Exception as e:
            if verbose:
                print('checkZlib: unknown error in %s, zlib not available' % file_name)
                logging.exception(e)
            ret_val = False
    except Exception as e:
        ret_val = False

    return ret_val


# %%

def get_version_info(verbose=0):
    """ Extract version information from source code """

    GIT_REVISION = None
    try:
        if os.path.exists('src/version.h'):
            with open('src/version.h') as f:
                ln = f.readline()
                m = re.search('.* "(.*)"', ln)
                FULLVERSION = (m.group(1))
        else:
            FULLVERSION = '0.0'
    except Exception as E:
        FULLVERSION = '0.0'
    if verbose:
        print('get_version_info: %s' % FULLVERSION)
    return FULLVERSION, GIT_REVISION


try:
    from distutils.version import LooseVersion
    from distutils.spawn import find_executable

    def get_swig_executable(swig_minimum_version='3.0', verbose=0):
        """ Get SWIG executable """
        # stolen from https://github.com/FEniCS/ffc/blob/master/setup.py

        swig_executable = None
        swig_version = None
        swig_valid = False
        for executable in ["swig", "swig3.0"]:
            swig_executable = find_executable(executable)
            if swig_executable is not None:
                # Check that SWIG version is ok
                output = subprocess.check_output([swig_executable, "-version"]).decode('utf-8')
                swig_version = re.findall(r"SWIG Version ([0-9.]+)", output)[0]
                if LooseVersion(swig_version) >= LooseVersion(swig_minimum_version):
                    swig_valid = True
                    break
        if verbose:
            print("Found SWIG: %s (version %s)" % (swig_executable, swig_version))
        return swig_executable, swig_version, swig_valid
    swig_executable, swig_version, swig_valid = get_swig_executable()
    print('swig_version %s, swig_executable %s' % (swig_version, swig_executable))
except BaseException:
    def get_swig_executable():
        return None, None, False
    swig_valid = False

# %% Trick to remove option for c++ code compilation
try:
    # see
    # http://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
    from setuptools.py31compat import get_config_vars

    (opt,) = get_config_vars('OPT')

    if not opt is None:
        opt = " ".join(flag for flag in opt.split()
                       if flag != '-Wstrict-prototypes')
        os.environ['OPT'] = opt
except BaseException:
    import setuptools
    print('old version of setuptools: %s' % setuptools.__version__)


# %% Test suite

class OATest(TestCommand):
    """ Run a limited set of tests for the package """

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
        print('## oapackage test: load package')
        # import here, cause outside the eggs aren't loaded
        import oapackage
        import oapackage.oahelper
        import oapackage.scanf
        print('## oapackage test: oalib version %s' % oapackage.version())
        print('## oapackage test: package compile options\n%s\n' % oapackage.oalib.compile_information())

        oapackage.oalib.test_array_manipulation(verbose=0)
        oapackage.oalib.test_conference_candidate_generators(verbose=0)
       
        errno = 0
        sys.exit(errno)


# %% Define sources of the package
oadev = 0
srcs = ['arraytools.cpp', 'arrayproperties.cpp', 'pareto.cpp', 'nonroot.cpp',
        'mathtools.cpp', 'oaoptions.cpp', 'tools.cpp', 'md5.cpp', 'strength.cpp', 'graphtools.cpp',
        'conference.cpp', 'unittests.cpp','Deff.cpp', 'evenodd.cpp']


srcs = srcs + ['lmc.cpp', 'extend.cpp']  # code used for extension
srcs = ['src/' + ff for ff in srcs]
if os.path.exists('dev/oadevelop.cpp'):
    oadev = 1
    print('Building development code')
    srcs = ['dev/oadevelop.cpp'] + srcs


sources = srcs + ['src/bitarray/bit_array.cpp']
oaheaders = [cppfile.replace('.cpp', '.h') for cppfile in srcs] + [os.path.join('src', 'version.h')]


for nauty_file in 'nauty.c nautinv.c nautil.c naurng.c naugraph.c schreier.c naugroup.c'.split(' '):
    sources += [os.path.join('src', 'nauty', nauty_file)]

nautyheaders = [
    os.path.join(
        'src',
        'nauty',
        headerfile) for headerfile in [
            'gtools.h',
            'naugroup.h',
            'nautinv.h',
            'naurng.h',
            'naugraph.h',
            'nausparse.h',
            'nautil.h',
            'nauty.h',
        'schreier.h']]
sources = sources

swig_opts = []
compile_options = []

sources = ['oalib.i'] + sorted(sources)
if oadev:
    swig_opts += ['-modern', '-c++', '-w503,401,362,509,389',
                  '-Isrc/', '-Idev/']
    compile_options += ['-DSWIGCODE', '-DFULLPACKAGE', '-DOADEV', '-Idev/']
    swig_opts += ['-DSWIGCODE', '-DFULLPACKAGE', '-DOADEV']
else:
    swig_opts += ['-modern', '-c++',
                  '-w503,401,362,302,389,446,509,305', '-Isrc/']
    compile_options += ['-DSWIGCODE', '-DFULLPACKAGE']
    swig_opts += ['-DSWIGCODE', '-DFULLPACKAGE']

# add nauty files
swig_opts += ['-Isrc/nauty/']
compile_options += ['-Isrc/nauty/']


if platform.system() == 'Windows':
    compile_options += ['-DWIN32', '-D_WIN32']
    swig_opts += ['-DWIN32', '-D_WIN32']

rtd = os.environ.get('READTHEDOCS', False)
print('Readthedocs environment: %s' % (rtd,))

if 'VSC_SCRATCH' in os.environ.keys():
    # we are running on the VSC cluster
    zlibdir = os.environ.get('EBROOTZLIB')  # SOFTROOTZLIB

    libraries = ['z']
    library_dirs = [zlibdir + '/lib']
    include_dirs = ['.', 'src', npinclude, zlibdir + '/include']
else:
    libraries = []
    library_dirs = []
    include_dirs = ['.', 'src', npinclude]

oalib_module = Extension('_oalib', sources=sources,
                         include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries, swig_opts=swig_opts
                         )

compile_options += ['-DNOOMP']
swig_opts += ['-DNOOMP']

oalib_module.extra_compile_args = compile_options

if checkZlib(verbose=1):
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
else:
    oalib_module.extra_compile_args += ['-O3', '-Wno-unknown-pragmas', '-Wno-sign-compare',
                                        '-Wno-return-type', '-Wno-unused-variable', '-Wno-unused-result', '-fPIC']
    oalib_module.extra_compile_args += ['-Wno-date-time', ]

if platform.node() == 'marmot' or platform.node() == 'goffer' or platform.node() == 'woelmuis':
    # openmp version of code
    oalib_module.extra_compile_args += ['-fopenmp', '-DDOOPENMP']
    oalib_module.extra_link_args += ['-fopenmp']

print('find_packages: %s' % find_packages())

data_files = []
scripts = ['misc/scripts/example_oapackage_python.py']
packages = find_packages()

# fix from:
# http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module

if rtd and 0:
    ext_modules = []  # do not build on RTD, this generates a time-out error
    swigcmd = '%s -python -modern -c++ -w503,401,362,302,389,446,509,305 -Isrc/ -DSWIGCODE -DFULLPACKAGE -Isrc/nauty/ -DWIN32 -D_WIN32 -DNOOMP -DNOZLIB -o oalib_wrap.cpp oalib.i' % swig_executable
    print('RTD: run swig command: %s' % (swigcmd,))
    output = subprocess.check_output(swigcmd.split(' '))
    print('swig output:')
    print(output)
else:
    if not swig_valid:
        raise Exception('could not find a recent version if SWIG')

    ext_modules = [oalib_module]

# see: http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module


class CustomBuild(distutils_build):

    def run(self):
        self.run_command('build_ext')
        distutils_build.run(self)


class CustomInstall(setuptools_install):

    def run(self):
        self.run_command('build_ext')
        setuptools_install.run(self)


class BuildExtSwig3(setuptools.command.build_ext.build_ext):
    def find_swig(self):
        swig_executable, _, _ = get_swig_executable()
        return swig_executable


def readme():
    with open('README.md') as f:
        return f.read()


long_description = readme()

version = get_version_info()[0]
print('OApackage: version %s' % version)

if is_python3:
    python27_requirements = []
else:
    python27_requirements = ['mock; python_version <"3.0"', 'backports.functools_lru_cache;python_version<"2.9"']

setup(name='OApackage',
      cmdclass={'test': OATest, 'install': CustomInstall, 'build': CustomBuild, 'build_ext': BuildExtSwig3},
      version=version,
      author="Pieter Eendebak",
      description="Package to generate and analyse orthogonal arrays, conference designs and optimal designs",
      long_description=long_description,
      long_description_content_type='text/markdown',
      author_email='pieter.eendebak@gmail.com',
      license="BSD",
      url='http://www.pietereendebak.nl/oapackage/index.html',
      keywords=["orthogonal arrays, design of experiments, conference designs, isomorphism testing"],
      ext_modules=ext_modules,
      py_modules=['oalib'],
      packages=packages,
      data_files=data_files,
      scripts=scripts,
      tests_require=['numpy', 'nose>=1.3', 'coverage>=4.0', 'mock', 'python-dateutil'] + python27_requirements,
      zip_safe=False,
      install_requires=['numpy>=1.13', 'python-dateutil'] + python27_requirements,
      extras_require={
          'GUI': ["qtpy", 'matplotlib'],
          'documentation': ['sphinx']
      },
      requires=['numpy', 'matplotlib'],
      classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   "Programming Language :: Python :: 3.8",
                   'License :: OSI Approved :: BSD License'
                   ]
      )
