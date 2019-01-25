Introduction
============

Orthogonal arrays, optimal designs and conference designs are important tools for the design of
experiments :cite:`Elster1995` :cite:`Hedayat1999` :cite:`Wu2009`. The Orthogonal Array
package contains functionality to generate and analyse these types of designs. The algorithms and methods in this package
are described in :cite:`Eendebak2009` and :cite:`EendebakSO`. A large
collection of results generated with the package can be found on the Orthogonal Array
package website :cite:`EendebakOApackageResults`. An alternative collection of orthogonal arrays can be found on the website of Neil
Sloane :cite:`Sloanewebsite`.


Example usage
-------------

The Orthogonal Array package can be used to generate and manipulate designs. Additionally, it can
calculate some of their statistical properties. The following example shows how to generate an orthogonal array with
8 runs and 2 factors, and calculate three relevant statistical properties:

.. admonition::  Calculate D-efficiency

  .. doctest:: 
  
   >>> import oapackage
   >>> array=oapackage.exampleArray(0) # define an orthogonal array 
   >>> array.showarray() 
   array:
     0   0
     0   0
     0   1
     0   1
     1   0
     1   0
     1   1
     1   1
   >>> D = array.Defficiency() # calculate the D-efficiency for estimating the [main-effects model?]
   >>> array_rank = array.rank() # calculate the rank of the design
   >>> print('D-efficiency %f, rank %d' % (D, array_rank) )
   D-efficiency 1.000000, rank 2
   >>> gwlp = array.GWLP() # calculate the generalized word length pattern
   >>> print('Generalized wordlength pattern: %s' % (gwlp,) )
   Generalized wordlength pattern: (1.0, 0.0, 0.0)

Interfaces
----------

The Orthogonal Array package has interfaces in C++ and Python. These package also includes command line 
tools for generating, manipulating and analyzing all the types of designs. In the documentation here you will find references to
both the Python and the C++ interface.

For the generation of optimal designs :cite:`EendebakSO`, the Orthogonal Array package has also a Matlab and R interfaces; see
the documentation `README.Matlab.md <https://github.com/eendebakpt/oapackage/README.Matlab.md>`_
and `README.R.md <https://github.com/eendebakpt/oapackage/README.R.md>`_.

License 
-------

The code is available under a BSD style license; see the file `LICENSE <https://github.com/eendebakpt/oapackage/blob/master/LICENSE>`_
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

The Python interface requires Numpy :cite:`NumPy2012`,
Matplotlib :cite:`Matplotlib` and Swig. The code has been
tested with Python 2.7, 3.5, 3.6 and 3.7.

The R interface to the optimal design functionality of the package is available from
`CRAN <http://cran.r-project.org/web/packages/oapackage/index.html>`_.
For the Matlab and Octave interface to the optimal design functionality see the 
file `README.Matlab.md <https://github.com/eendebakpt/oapackage/blob/master/README.Matlab.md>`_.



