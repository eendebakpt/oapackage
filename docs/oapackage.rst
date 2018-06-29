============================
The Orthogonal Array package
============================

:Author: P.T. Eendebak [1]_
:Date:   July, 2015

.. role:: math(raw)
   :format: html latex
..

.. role:: raw-latex(raw)
   :format: latex
..

|image|

Introduction
============

Orthogonal arrays are an important tool in the design of
experiments :cite:`Hedayat1999`. The Orthogonal Array
package contains functionality the generate orthogonal arrays and to
analyse their properties. The algorithms and methods in this package
have been described in :cite:`Eendebak2009`. A large
collection of arrays can be found on the Orthogonal Array
packagewebsite :cite:`EendebakOA` or on the website of Neil
Sloane :cite:`Sloanewebsite`.


Example usage
-------------

The Orthogonal Array packagecan be used to work with arrays and
calculate statistical properties. For example to calculate the
:math:`D`-efficiency and rank of a design we can use:


.. code-block:: python3
   :caption: Calculate D-efficiency 

   >>> al=oapackage.exampleArray(0)
   >>> al.showarray()
   array: 0 0 0 0 0 1 0 1 1 0 1 0 1 1 1 1
   >>> print(’D-efficiency D-efficiency 1.000000, rank 2
   >>> print(’Generalized wordlength pattern: %s' % al.GWLP()
   Generalized wordlength pattern: (1.0, 0.0, 0.0)

Interfaces
----------

The Orthogonal Array package has several interfaces. First of all there
are command line tools for manipulating arrays and generating. All
functions of the package can be used from either C++ or Python. For a
restricted set of functionality also Matlab or R can be used.

Compilation and installation
----------------------------

The program has been tested using Linux and Windows (XP, Windows 7 and
Windows 8). The Python interface is available from the Python Package
Index https://pypi.python.org/pypi/OApackage/. The package can be
installed from the command line using pip:

> pip install –user OApackage

The R interface to the package is available from CRAN
http://cran.r-project.org/web/packages/oapackage/index.html.

The command line tools uses a cmake build system. From the command line
type:

> mkdir -p build; cd build > cmake .. > make > make install

This creates the command line utilities and a C++ library. To compile
the Python interface using Linux use

> python setup.py build > python setup.py install –user

The Python interface requires Numpy :cite:`NumPy`,
Matplotlib :cite:`Matplotlib` and Swig. The code has been
tested with Python 2.7, 3.5 and 3.6.

Using Windows start Cygwin or the Visual Studio command prompt. From the
package source directory run:

> python setup.py bdist\_wininst

This creates a binary installer package.

License
-------

The code is available under a BSD style license, see the file LICENSE
for details. If you use this code or any of the results, please cite
this program as follows:

*Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays*,
P.T. Eendebak, E.D. Schoen, M.V.M. Nguyen, Volume 18, Issue 2, pages
123-140, 2010.

Acknowledgements
----------------

The code and ideas for this package have been contributed by Eric
Schoen, Ruben Snepvangers, Vincent Brouerius van Nidek, Alan
Vazquez-Alcocer and Pieter Thijs Eendebak.

