Introduction
============

Orthogonal arrays, optimal designs and conference designs are important tools for the design of
experiments :cite:`Elster1995` :cite:`Hedayat1999` :cite:`Wu2009`. The Orthogonal Array
package contains functionality to generate and analyse these types of designs.
To generate the arrays and designs, the package uses the exhaustive enumeration algorithm of :cite:`Eendebak2009` and
the optimization algorithm of :cite:`EendebakSO`.
To analyze the arrays and designs, the package includes a wide variety of relevant statistical and combinatorial
criteria.
A large collection of orthogonal arrays, optimal designs and conference designs generated with the package are available in the Orthogonal Array package website :cite:`EendebakOApackageResults`.


Example usage
-------------

The Orthogonal Array package can be used to generate and manipulate arrays and designs. Additionally, it can
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
   >>> D = array.Defficiency() # calculate the D-efficiency for estimating the interaction effects model
   >>> array_rank = array.rank() # calculate the rank of the design
   >>> print('D-efficiency %f, rank %d' % (D, array_rank) )
   D-efficiency 1.000000, rank 2
   >>> gwlp = array.GWLP() # calculate the generalized word length pattern
   >>> print('Generalized wordlength pattern: %s' % (gwlp,) )
   Generalized wordlength pattern: (1.0, 0.0, 0.0)

The statistical properties of the arrays and designs are introduced in :ref:`Properties of designs`.

Interfaces
----------

The Orthogonal Array package has interfaces in C++ and Python for generating, manipulating and analyzing all the types of arrays and designs. In this documentation, you will find references to
both the Python and the C++ interface. The package also includes several command line tools.

For the generation of optimal designs :cite:`EendebakSO`, the Orthogonal Array package has also a Matlab interface; see
the documentation `README.Matlab.md <https://github.com/eendebakpt/oapackage/README.Matlab.md>`_.

License 
-------

The code is available under a BSD style license; see the file `LICENSE <https://github.com/eendebakpt/oapackage/blob/master/LICENSE>`_
for details. If you use this code or any of the results, please cite
this program as follows:

* `OApackage: A Python package for generation and analysis of orthogonal arrays, optimal designs and conference designs <https://doi.org/10.21105/joss.01097>`_, P.T. Eendebak, A.R. Vazquez, Journal of Open Source Software, 2019
* *Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays*, E.D. Schoen, P.T. Eendebak, M.V.M. Nguyen, Volume 18, Issue 2, pages 123-140, 2010.

Acknowledgements
----------------

The code and ideas for this package have been contributed by Eric
Schoen, Ruben Snepvangers, Vincent Brouerius van Nidek, Alan
Roberto Vazquez and Pieter Thijs Eendebak.

Installation
------------

The packge is continously tested on Linux and Windows. The Python interface is available from the `Python Package
Index <https://pypi.python.org/pypi/OApackage/>`_. The package can be
installed from the command line using pip:

.. code-block:: console

  $ pip install OApackage

The source code for the package is available on https://github.com/eendebakpt/oapackage.
The command line tools use a cmake build system. From the command line,
type the following:

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
tested with Python 3.6, 3.7 and 3.8.

The R interface to the optimal design functionality of the package is available from
`CRAN <http://cran.r-project.org/web/packages/oapackage/index.html>`_.
For the Matlab and Octave interface to the optimal design functionality see the 
file `README.Matlab.md <https://github.com/eendebakpt/oapackage/blob/master/README.Matlab.md>`_.

Related sites of orthogonal arrays
----------------------------------

There are several related sites available online which include collections
of orthogonal arrays. For instance, the website of Neil Sloane :cite:`Sloanewebsite`,
the website of Hongquan Xu :cite:`HongquanXuOnline`, the SAS website managed
by Warren Kuhfeld :cite:`WK19`, and the R package _DoE.base_ :cite:`DoEbase` include lists
and surveys of attractive orthogonal arrays gathered from different sources. 


