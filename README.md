Orthogonal Array Package
========================

The Orthogonal Array package contains functionality to generate and analyse orthogonal arrays, optimal designs and conference designs.
Features include generation of complete series of orthogonal arrays, 
reduction of arrays to normal form and calculation of properties such as the strength or D-efficiency of an array.
For more information about the package see the
documentation at [oapackage.readthedocs.io](https://oapackage.readthedocs.io/en/latest/). A large collection of results generated
with the package can be found at <http://pietereendebak.nl/oapackage/>.

Usage
-------

The package can be used from Python:
``` python
>>> import oapackage
>>> al=oapackage.exampleArray(0)
>>> al.showarraycompact()
00
00
01
01
10
10
11
11
>>> print('D-efficiency %f, rank %d' % (al.Defficiency(), al.rank()) )
D-efficiency 1.000000, rank 2
>>> print('Generalized wordlength pattern: %s' % str(al.GWLP()))
Generalized wordlength pattern: (1.0, 0.0, 0.0)
```

For more examples see the Jupyter notebooks in the
[docs/examples](docs/examples/). 

Acknowledgements
----------------

If you use this code or any of the results, please cite this program as follows:

* [OApackage: A Python package for generation and analysis of orthogonal arrays, optimal designs and conference designs](https://doi.org/10.21105/joss.01097), P.T. Eendebak, A.R. Vazquez, Journal of Open Source Software, 2019
* [Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays](https://doi.org/10.1002/jcd.20236), E.D. Schoen, P.T. Eendebak, M.V.M. Nguyen, Volume 18, Issue 2, pages 123-140, 2010
* [Two-Level Designs to Estimate All Main Effects and Two-Factor Interactions](https://doi.org/10.1080/00401706.2016.1142903), Pieter T. Eendebak, Eric D. Schoen, Technometrics Vol. 59 , Iss. 1, 2017
* [A classification criterion for definitive screening designs](https://projecteuclid.org/euclid.aos/1547197252), E.D. Schoen, P.T. Eendebak, P. Goos, Ann. Statist. 47, no. 2, 2019.

The code was written by:

* Pieter Eendebak <pieter.eendebak@gmail.com>
* Alan Vazquez-Alcocer
* Vincent Brouerius van Nidek

Ideas contributed by:

* Eric Schoen <eric.schoen@tno.nl>
* Alan Vazquez-Alcocer <alanrvazquez@gmail.com>

See the file LICENSE for copyright details.

Installation
------------

[![PyPI version](https://badge.fury.io/py/OApackage.svg)](https://badge.fury.io/py/OApackage)
[![Build status](https://ci.appveyor.com/api/projects/status/f6ia9br95soimf9u/branch/master?svg=true)](https://ci.appveyor.com/project/eendebakpt/oapackage-4lws8)
[![Build Status](https://travis-ci.org/eendebakpt/oapackage.svg?branch=master)](https://travis-ci.org/eendebakpt/oapackage)

The Python interface to the package is available on the [Python Package index](https://pypi.python.org/pypi/OApackage/).
Installation can be done using the following command:
``` console
$ pip install OApackage 
```
(or `pip install OApackage --user` if you do not have admin rights). To compile the package you need Python, Numpy and Swig 3.x.

The command line tools have been tested using Linux, Windows XP/Win7/Win10 and Raspberry Pi.
The program uses a `cmake` build system. From the command line type:
```
$ mkdir -p build; cd build
$ cmake ..
$ make
$ make install
````

Contributing, unit testing and support
--------------------------------------

See the file [CONTRIBUTING.md](https://github.com/eendebakpt/oapackage/blob/master/CONTRIBUTING.md) on GitHub.
