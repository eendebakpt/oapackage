Orthogonal Array Package
========================

(version 2.3.3)

The code allows to work with orthogonal arrays. Features include generation of complete series of orthogonal arrays, 
reduction of arrays to normal form and calculation of properties such as the strength or D-efficiency of an array.
For more information about the package see the page <http://pietereendebak.nl/oapackage/>.

Usage
-------

The package can be used from Python:
``` python
>>> import oapackage
>>> al=oapackage.exampleArray(0)
>>> al.showarray()
array: 
  0   0
  0   0
  0   1
  0   1
  1   0
  1   0
  1   1
  1   1
>>> print('D-efficiency %f, rank %d' % (al.Defficiency(), al.rank()) )
D-efficiency 1.000000, rank 2
>>> print('Generalized wordlength pattern: %s' % str(al.GWLP()))
Generalized wordlength pattern: (1.0, 0.0, 0.0)
```

For for examples see the Ipython notebooks in the
[examples](examples/).

Acknowledgements
------------

If you use this code or any of the results, please cite this program as follows:

Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays, P.T. Eendebak, E.D. Schoen, M.V.M. Nguyen, Volume 18, Issue 2, pages 123-140, 2010.

The code was written by:

* Pieter Eendebak <pieter.eendebak@gmail.com>
* Vincent Brouerius van Nidek
* Alan Vazquez-Alcocer

Ideas contributed by:

* Eric Schoen <eric.schoen@tno.nl>
* Alan Vazquez-Alcocer <alanrvazquez@gmail.com>

See the file LICENSE for copyright details.

Installation
------------

[![PyPI version](https://badge.fury.io/py/OApackage.svg)](https://badge.fury.io/py/OApackage)
[![Build status](https://ci.appveyor.com/api/projects/status/f6ia9br95soimf9u?svg=true)](https://ci.appveyor.com/project/eendebakpt/oapackage-4lws8)

The Python interface to the package is available on the [Python Package index](https://pypi.python.org/pypi/OApackage/).
Installation can be done using the following command::

> pip install OApackage --user

To compile the package you need Python, Numpy and Swig 3.x.

The binary tools have been tested using Linux, Windows XP/Win7 and Raspberry Pi.
The program uses a cmake build system. From the commandline type::

> mkdir -p build; cd build
> cmake ..
> make
> make install


Data format
-----------

Arrays are stored in plain files or binary files. For text files the first line contains the number of columns, the number of rows and the number of arrays (or -1 if the number of arrays is not specified). Then for each array a single line with the index of the array, followed by N lines containing the array. The binary format is suitable for storing a very high numbers of arrays and supports random access. Also see the file `FORMAT.txt`.

To examine the result files one can use the oacat tool from the package, or the standard UNIX utilities cat, less, head and tail.






