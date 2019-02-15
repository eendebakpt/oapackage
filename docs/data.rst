


Data representation
===================

All designs handled by the OApackage are integer valued. The designs (whether these are orthogonal arrays, optimal designs or conferences designs)
are stored in an :meth:`array_link` object.

An orthogonal array (OA) of strength :math:`{t}`, :math:`{N}` runs and
:math:`{n}` factors at :math:`{s}` levels is an :math:`{N}\times {n}`
array of symbols :math:`0,
\ldots,({s}-1)`, such that for every subset of :math:`{t}` columns,
every :math:`{t}`-tuple occurs equally
often :cite:`Rao1947`. The set of all strength-:math:`{t}` OAs with 
:math:`{N}` runs and :math:`{n}` factors at :math:`{s}` levels 
is denoted by :math:`{\operatorname{OA}({N}; {t}; {s}^{n})}`. 
The OAs are represented by arrays (data in memory is stored in column-major form).

A D-optimal design :cite:`Donev2007` (:math:`X`) is an :math:`{N}\times {n}` array 
that maximizes the :math:`{(\operatorname{det}({X}^{T}_{M}{X}^{\phantom{T}}_{M})^{1/p})/N}`,
for a given :math:`{N}\times {p}` model matrix :math:`{X}_{M}` (for details see :ref:`Model matrices`).
An orthogonal array is called D-optimal if it provides the largest determinant among all comparable orthogonal arrays.

For :math:`{N}` even, a conference design :math:`C` is 
an :math:`{N}\times {n}` array which satisfies :math:`{C}^{T}C = (n-1) I_{n}`,
with :math:`{C}_{ii} = 0` and :math:`{C}_{ij} \in \{-1,1\}`, for 
:math:`{i} \neq {j}` and :math:`{i}, {j} = 1, \ldots, n`.
See :cite:`Xiao2012`, :cite:`Wu2009`.


Data structures
---------------

The package contains several data structures. Here we describe the main
structures and their use.

  :meth:`~oalib.array_link`
    The structure containing an orthogonal array is called the
    :meth:`~oalib.array_link` structure. Lists of arrays are stored in the
    :meth:`~oalib.arraylist_t` object, which is implemented as a :code:`std::deque` container.

  :meth:`~oalib.arrayfile_t`
    This is an object that allows for reading and writing of arrays to
    disk.

  :meth:`~oalib.arraydata_t`
    The structure describing a certain class of orthogonal arrays or
    optimal designs.

  :meth:`~oalib.conference_t`
    The structure describing a certain class of conference designs.

  :meth:`~oalib.array_transformation_t`
    This describes a transformation of an orthogonal array, which 
    includes row and column permutations, as well as permutations 
    to the symbols in one or more columns.

  :meth:`~oalib.conference_transformation_t`
    This describes a transformation of conference design or double conference design, which includes 
    row and column permutations, as well as sign switches to the elements
    in one or more rows and columns.

Representing arrays
-------------------

The structure containing an orthogonal array is called the
:cpp:class:`array_link` structure. It consists of a specified number of rows and
columns, the data (integer valued) and an index. In the Python interface the :meth:`array_link` object can be indexed just as
normal arrays. 

It is also possible to convert to a Numpy array. The
:class:`~oalib.array_link` object implements the Python array interface, so most
operations from packages such as Numpy work on the :meth:`~oalib.array_link`
object.

.. testsetup::
   
   import sys
   # doctest: +NORMALIZE_WHITESPACE


.. admonition:: Array representation and indexing in Python

  .. doctest:: 
  
    >>> import oapackage; import numpy as np
    >>> al=oapackage.exampleArray(0)
    >>> al.showarray() # doctest: +NORMALIZE_WHITESPACE
    array:
      0   0
      0   0
      0   1
      0   1
      1   0
      1   0
      1   1
      1   1
      
    >>> al[2,1] 
    1
    >>> X=np.array(al)
    >>> X 
    array([[0, 0],
           [0, 0],
           [0, 1],
           [0, 1],
           [1, 0],
           [1, 0],
           [1, 1],
           [1, 1]], dtype=int16)

The C++ class is :cpp:class:`array_link`.
    
Classes of arrays
-----------------

The :cpp:class:`arraydata_t` object represents data about a class of orthogonal
arrays, e.g. the class :math:`{\operatorname{OA}(N; t; s^k)}`.
The :cpp:class:`conference_t` object represents data about a class of conference designs.


Array transformations
---------------------

Transformations of (orthogonal) arrays consist of row, column and 
level permutations. A transformation is represented by 
the :cpp:class:`array_transformation_t` object.

For a given transformation the column permutations are applied first,
then the level permutations and finally the row permutations. The level
and column permutations are not commutative.

The conference transformations also allow for row sign switches and are
described by the class :cpp:class:`conference_transformation_t`.

Reading and writing arrays
--------------------------

Reading and writing arrays to disk can be done with the :cpp:class`arrayfile_t`
class. 

.. admonition:: Write an array or a list of arrays to disk

  .. doctest:: 

   >>> import oapackage
   >>> list_of_arrays = [oapackage.exampleArray(24), oapackage.exampleArray(25)]
   >>> _ = oapackage.writearrayfile('test.oa', list_of_arrays)
   >>> oapackage.oainfo('test.oa')
   file test.oa: 64 rows, 16 columns, 2 arrays, mode text, nbits 0
   >>> al=oapackage.exampleArray()
   >>> af=oapackage.arrayfile_t('test.oa', al.n_rows, al.n_columns)
   >>> af.append_array(al)
   >>> print(af)
   file test.oa: 8 rows, 2 columns, 1 arrays, mode text, nbits 8
   >>> af.closefile()

The arrays can be written in text or binary format. For more details on
the file format, see Section :ref:`File formats`.

The Python interface is :meth:`oalib.arrayfile_t` and the C++ interface is  :cpp:class`arrayfile_t`.

File formats
------------

The Orthogonal Array package stores arrays in a custom file
format. There is a text format which is easily readable by humans and a
binary format which is faster to process and memory efficient.

Plain text array files
~~~~~~~~~~~~~~~~~~~~~~

Arrays are stored in plain text files with extension ``.oa``. The first line
contains the number of columns, the number of rows and the number of
arrays (or -1 if the number of arrays is not specified). Then, for each
array, a single line with the index of the array, followed by N lines
containing the array.

A typical example of a text file is the following:

.. code-block:: c

  5 8 1
  1
  0 0 0 0 0
  0 0 0 1 1
  0 1 1 0 0
  0 1 1 1 1
  1 0 1 0 1
  1 0 1 1 0 
  1 1 0 0 1 
  1 1 0 1 0
  -1

This file contains exactly 1 array with 8 rows and 5 columns.

Binary array files
~~~~~~~~~~~~~~~~~~

Every binary file starts with a header, which has the following format:

.. code-block:: c

  [INT32] 65 (magic identifier) 
  [INT32] b: Format: number of bits per number. Currently supported are 1 and 8
  [INT32] N: number of rows 
  [INT32] k: kumber of columns 
  [INT32] Number of arrays (can be -1 if unknown)
  [INT32] Binary format number: 1001: normal, 1002: binary diff, 1003: binary diff zero
  [INT32] Reserved integer
  [INT32] Reserved integer

The format of the remainder of the binary file depends on the binary format specified.
For the normal binary format the format is as follows. For each array (the
number is specified in the header):

.. code-block:: c

  [INT32] Index
  [Nxk elements] The elements contain b bits

If the number of bits per number is 1 (e.g. a 2-level array), then the
data is padded with zeros to a multiple of 64 bits. The data of the
array is stored in column-major order. The binary file format allows for
random access reading and writing. The ``binary diff`` and ``binary diff
zero`` formats are special formats.

A binary array file can be compressed using gzip. Most tools in the
Orthogonal Array package can read these compressed files transparently.
Writing to compressed array files is not supported at the moment.

Data files
~~~~~~~~~~

The analysis tool (``oaanalyse``) writes data to disk in binary format.
The format consists of a binary header:

::

  [FLOAT64] Magic number 30397995;
  [FLOAT64] Magic number 12224883;
  [FLOAT64] nc: Number of rows
  [FLOAT64] nr: Number of columns

After the header there follow ``nc*nr`` ``[FLOAT64]`` values.

MD5 sums
~~~~~~~~

To check data integrity on disk the packages includes functions to
generate MD5 sums of designs. 

.. admonition:: Calculate md5 sum of a design

  .. doctest:: 

     >>> import oapackage; al=oapackage.exampleArray(0)
     >>> al.md5()
     '6454c492239a8e01e3c01a864583abf2'

The C++ functions :cpp:func:`array_link::md5` and :cpp:func:`md5`.

Command line interface
----------------------

Several command line tools are included in the Orthogonal Array package. For each tool,
help can be obtained from the command line by using the switch ``-h``.
The tools include the following:

`oainfo`
    This program reads Orthogonal Array package data files and reports
    the contents of the files. For example:

    .. code-block:: console
    
        $ oainfo result-8.2-2-2-2.oa
        Orthogonal Array package 1.8.7
        oainfo: reading 1 file(s)
        file result-8.2-2-2.oa: 8 rows, 3 columns, 2 arrays, mode text, nbits 0
        $

`oacat`
    Shows the contents of a file with orthogonal arrays for a data file.

`oacheck`
    Checks or reduces an array to canonical form.

`oaextendsingle`
    Extends a set of arrays in LMC form with one or more columns.

`oacat`
    Shows the contents of an array file or data file.

    Usage: ``oacat [OPTIONS] [FILES]``

`oajoin`
    Reads one or more files from disk and join all the array files into a
    single list.

`oasplit`
    Takes a single array file as input and splits the arrays into a
    specified number of output files.

`oapareto`
    Calculates the set of Pareto optimal arrays in a file with arrays.

`oaanalyse`
    Calculates various statistical properties of arrays in a file. 
    The properties are described in section :ref:`Properties of designs`.


.. figure:: images/oaimage-18_2-3-3-3-3-3-n17.png
   :alt: alternate text
   :align: center

   Orthogonal array in :math:`\mathrm{OA}(18, 2 3^a, 2)`. 
