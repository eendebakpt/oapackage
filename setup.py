#!/usr/bin/env python

"""
setup.py file for OApackage
"""

import logging
import os
import platform
import re
import subprocess
import sys

# from distutils.command.build import build as distutils_build
from os import path

import setuptools.command.build_ext
from setuptools import Extension, find_packages, setup

if sys.version_info.minor > 10 or sys.version_info.major > 3:
    from setuptools.command.build import build as setuptools_build
else:
    from distutils.command.build import build as setuptools_build
from setuptools.command.install import install as setuptools_install

try:
    import numpy as np
except ImportError:
    raise RuntimeError("numpy cannot be imported. numpy must be installed prior to installing OApackage")

npinclude = np.get_include()
setup_directory = path.abspath(path.dirname(__file__))

# %%


def checkZlib(verbose=0):
    """Check whether zlib is available

    Code adapted from http://stackoverflow.com/questions/28843765/setup-py-check-if-non-python-library-dependency-exists
    """
    ret_val = True
    return True
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
        tmp_dir = tempfile.mkdtemp(prefix="helloztest_")
        bin_file_name = os.path.join(tmp_dir, "helloz")
        file_name = bin_file_name + ".c"

        with open(file_name, "w") as fp:
            fp.write(c_code)

        # compile and link code
        compiler = distutils.ccompiler.new_compiler()
        assert isinstance(compiler, distutils.ccompiler.CCompiler)
        distutils.sysconfig.customize_compiler(compiler)

        libraries = ["z"]
        ret_val = True
        try:
            if verbose:
                print("checkZlib: compile and link")
            compiler.link_executable(
                compiler.compile([file_name]),
                bin_file_name,
                libraries=libraries,
            )
        except CompileError:
            if verbose:
                print("checkZlib: compile error in %s, zlib not available" % file_name)
            ret_val = False
        except LinkError:
            if verbose:
                print("checkZlib: link error in %s, zlib not available" % file_name)
            ret_val = False
        except Exception as e:
            if verbose:
                print("checkZlib: unknown error in %s, zlib not available" % file_name)
                logging.exception(e)
            ret_val = False
    except Exception:
        ret_val = False

    return ret_val


# %%


def get_version_info(verbose=0):
    """Extract version information from source code"""

    GIT_REVISION = None
    if os.path.exists("src/version.h"):
        with open("src/version.h") as f:
            ln = f.readline()
            m = re.search('.* "(.*)"', ln)
            FULLVERSION = m.group(1)
    else:
        FULLVERSION = "0.0"
    if verbose:
        print("get_version_info: %s" % FULLVERSION)
    return FULLVERSION, GIT_REVISION


try:
    from shutil import which as find_executable
    # from distutils.spawn import find_executable

    try:
        from packaging.version import Version
    except ImportError:
        from distutils.version import LooseVersion as Version

    def get_swig_executable(swig_minimum_version="4.0", verbose=0):
        """Get SWIG executable"""
        # stolen from https://github.com/FEniCS/ffc/blob/master/setup.py

        swig_executable = None
        swig_version = None
        swig_valid = False
        for executable in ["swig", "swig3.0"]:
            swig_executable = find_executable(executable)
            if swig_executable is not None:
                # Check that SWIG version is ok
                output = subprocess.check_output([swig_executable, "-version"]).decode("utf-8")
                swig_version = re.findall(r"SWIG Version ([0-9.]+)", output)[0]
                if Version(swig_version) >= Version(swig_minimum_version):
                    swig_valid = True
                    break
        if verbose:
            print(f"Found SWIG: {swig_executable} (version {swig_version})")
        return swig_executable, swig_version, swig_valid

    swig_executable, swig_version, swig_valid = get_swig_executable()
    print(f"swig_version {swig_version}, swig_executable {swig_executable}")
except BaseException:

    def get_swig_executable():
        return None, None, False

    swig_valid = False


# %% Define sources of the package
oadev = 0
srcs = [
    "arraytools.cpp",
    "arrayproperties.cpp",
    "pareto.cpp",
    "nonroot.cpp",
    "mathtools.cpp",
    "oaoptions.cpp",
    "tools.cpp",
    "md5.cpp",
    "strength.cpp",
    "graphtools.cpp",
    "conference.cpp",
    "unittests.cpp",
    "Deff.cpp",
    "evenodd.cpp",
    "lmc.cpp",
    "extend.cpp",
]

srcs = ["src/" + ff for ff in srcs]
if os.path.exists("dev/oadevelop.cpp"):
    oadev = 1
    print("Building development code")
    srcs = ["dev/oadevelop.cpp"] + srcs


sources = srcs + ["src/bitarray/bit_array.cpp"]
oaheaders = [cppfile.replace(".cpp", ".h") for cppfile in srcs] + [os.path.join("src", "version.h")]


for nauty_file in "nauty.c nautinv.c nautil.c naurng.c naugraph.c schreier.c naugroup.c".split(" "):
    sources += [os.path.join("src", "nauty", nauty_file)]

nautyheaders = [
    os.path.join("src", "nauty", headerfile)
    for headerfile in [
        "gtools.h",
        "naugroup.h",
        "nautinv.h",
        "naurng.h",
        "naugraph.h",
        "nausparse.h",
        "nautil.h",
        "nauty.h",
        "schreier.h",
    ]
]
sources = sources

swig_opts = []
compile_options = []

sources = ["oalib.i"] + sorted(sources)
if oadev:
    swig_opts += ["-c++", "-doxygen", "-w503,401,362,509,389", "-Isrc/", "-Idev/"]
    compile_options += ["-DSWIGCODE", "-DFULLPACKAGE", "-DOADEV", "-Idev/"]
    swig_opts += ["-DSWIGCODE", "-DFULLPACKAGE", "-DOADEV"]
else:
    swig_opts += ["-c++", "-doxygen", "-w503,401,362,302,389,446,509,305", "-Isrc/"]
    compile_options += ["-DSWIGCODE", "-DFULLPACKAGE"]
    swig_opts += ["-DSWIGCODE", "-DFULLPACKAGE"]

# add nauty files
swig_opts += ["-Isrc/nauty/"]
compile_options += ["-Isrc/nauty/"]

# specifically for clang
# compile_options += ["-std=c++11"]


if platform.system() == "Windows":
    compile_options += ["-DWIN32", "-D_WIN32"]
    swig_opts += ["-DWIN32", "-D_WIN32"]

rtd = os.environ.get("READTHEDOCS", False)
print(f"Readthedocs environment: {rtd}")

if "VSC_SCRATCH" in os.environ.keys():
    # we are running on the VSC cluster
    zlibdir = os.environ.get("EBROOTZLIB")  # SOFTROOTZLIB

    libraries = ["z"]
    library_dirs = [zlibdir + "/lib"]
    include_dirs = [".", "src", npinclude, zlibdir + "/include"]
else:
    libraries = []
    library_dirs = []
    include_dirs = [".", "src", npinclude]

oalib_module = Extension(
    "_oalib",
    sources=sources,
    include_dirs=include_dirs,
    library_dirs=library_dirs,
    libraries=libraries,
    swig_opts=swig_opts,
)

compile_options += ["-DNOOMP"]
swig_opts += ["-DNOOMP"]

oalib_module.extra_compile_args = compile_options

if checkZlib(verbose=1):
    if platform.system() == "Windows":
        pass
    else:
        zlibflag = "-DUSEZLIB"
        oalib_module.extra_compile_args += [zlibflag]
        swig_opts += [zlibflag]
        oalib_module.extra_link_args += ["-lz"]
else:
    zlibflag = "-DNOZLIB"
    oalib_module.extra_compile_args += [zlibflag]
    swig_opts += [zlibflag]

if os.name == "nt":
    oalib_module.extra_compile_args += []
else:
    oalib_module.extra_compile_args += [
        "-O3",
        "-Wno-unknown-pragmas",
        "-Wno-sign-compare",
        "-Wno-return-type",
        "-Wno-unused-variable",
        "-Wno-unused-result",
        "-fPIC",
    ]
    oalib_module.extra_compile_args += [
        "-Wno-date-time",
    ]

if platform.node() == "marmot" or platform.node() == "goffer" or platform.node() == "woelmuis":
    # openmp version of code
    oalib_module.extra_compile_args += ["-fopenmp", "-DDOOPENMP"]
    oalib_module.extra_link_args += ["-fopenmp"]

print("find_packages: %s" % find_packages())

data_files = []
scripts = ["misc/scripts/example_oapackage_python.py"]
packages = find_packages()

# fix from:
# http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module

if not swig_valid:
    raise Exception("could not find a recent version if SWIG")

ext_modules = [oalib_module]

# see: http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module


class CustomBuild(setuptools_build):
    def run(self):
        self.run_command("build_ext")
        setuptools_build.run(self)


class CustomInstall(setuptools_install):
    def run(self):
        self.run_command("build_ext")
        setuptools_install.run(self)


class BuildExtSwig3(setuptools.command.build_ext.build_ext):
    def find_swig(self):
        swig_executable, _, _ = get_swig_executable()
        return swig_executable


def readme():
    with open("README.md") as f:
        return f.read()


long_description = readme()

version = get_version_info()[0]
# print("OApackage: version %s" % version)


setup(
    name="OApackage",
    cmdclass={"install": CustomInstall, "build": CustomBuild, "build_ext": BuildExtSwig3},
    version=version,
    author="Pieter Eendebak",
    description="Package to generate and analyse orthogonal arrays, conference designs and optimal designs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author_email="pieter.eendebak@gmail.com",
    license="BSD",
    url="http://www.pietereendebak.nl/oapackage/index.html",
    keywords=["orthogonal arrays, design of experiments, conference designs, isomorphism testing"],
    ext_modules=ext_modules,
    py_modules=["oalib"],
    packages=packages,
    data_files=data_files,
    scripts=scripts,
    tests_require=[
        "numpy>=1.26",
        "nose",
        "coverage",
        "matplotlib",
        "mock",
        "python-dateutil",
        "types-python-dateutil",
    ],
    zip_safe=False,
    install_requires=["numpy>=1.24", "python-dateutil", "matplotlib"],
    extras_require={
        "doc": ["sphinx", "sphinxcontrib.bibtex", "sphinxcontrib.napoleon", "breathe"],
    },
    requires=["numpy", "matplotlib"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: BSD License",
    ],
)
