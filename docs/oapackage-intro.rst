Introduction
============

Orthogonal arrays are an important tool in the design of
experiments :cite:`Hedayat1999`. The Orthogonal Array
package contains functionality to generate and analyse orthogonal arrays, optimal designs and conference designs. The algorithms and methods in this package
have been described in :cite:`Eendebak2009`. A large
collection of arrays can be found on the Orthogonal Array
packagewebsite :cite:`EendebakOA` or on the website of Neil
Sloane :cite:`Sloanewebsite`.


Example usage
-------------

The Orthogonal Array package can be used to work with arrays and
calculate statistical properties. For example to calculate the
:math:`D`-efficiency and rank of a design we can use:


.. code-block:: python
   :caption: Calculate D-efficiency 

   >>> import oapackage
   >>> al=oapackage.exampleArray(0)
   >>> al.showarray()
   array: 0 0 0 0 0 1 0 1 1 0 1 0 1 1 1 1
   >>> print('D-efficiency D-efficiency %f, rank %d’ % (al.Defficiency(), al.rank()) )
   D-efficiency 1.000000, rank 2
   >>> print('Generalized wordlength pattern: %s' % al.GWLP() )
   Generalized wordlength pattern: (1.0, 0.0, 0.0)

Interfaces
----------

The Orthogonal Array package has several interfaces. First of all there
are command line tools for manipulating arrays and generating designs. All
functions of the package can be used from either C++ or Python.

For the generation of D-optimal designs also Matlab or R can be used, see
the documentation `README.Matlab.md <https://github.com/eendebakpt/oapackage/README.Matlab.md>`_
and `README.R.md <https://github.com/eendebakpt/oapackage/README.R.md>`_.

License 
-------

The code is available under a BSD style license, see the file `LICENSE <https://github.com/eendebakpt/oapackage/blob/master/LICENSE>`_
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

Installation
------------

The program has been tested using Linux and Windows (XP, Windows 7 and
Windows 8, Windows 10). The Python interface is available from the `Python Package
Index <https://pypi.python.org/pypi/OApackage/>`_. The package can be
installed from the command line using pip:

.. code-block:: console

  $ pip install OApackage

The source code for the package is available on https://github.com/eendebakpt/oapackage.
The command line tools use a cmake build system. From the command line
type:

.. code-block:: console

  $ mkdir -p build
  $ cd build
  $ cmake .. 
  $ make
  $ make install

This creates the command line utilities and a C++ library.


To compile the Python interface use

.. code-block:: console

  $ python setup.py build 
  $ python setup.py install --user

The Python interface requires Numpy :cite:`NumPy`,
Matplotlib :cite:`Matplotlib` and Swig. The code has been
tested with Python 2.7, 3.5 and 3.6.

The R interface to the optimal design functionality of the package is available from CRAN
http://cran.r-project.org/web/packages/oapackage/index.html.
For the Matlab and Octave interface to the optimal design functionality see the 
file `README.Matlab.md <https://github.com/eendebakpt/oapackage/blob/master/README.Matlab.md>`_.



