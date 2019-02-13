Normal form of arrays
=====================


The Orthogonal Array package contains functions to reduce
designs to canonical form with respect to some ordering. The
default ordering for arrays is the lexicographic ordering in
columns :cite:`Eendebak2009`. Alternative orderings include the
delete-one-factor projection ordering introduced
in :cite:`EendebakDOF` or the even-odd ordering.
For a given ordering of a set of arrays, the minimal element of all arrays in an
isomorphism class defines a unique representative of that isomorphism class. 

Specialized packages such as Nauty :cite:`nautyI`, :cite:`nautyII` can also reduce
arrays to their canonical form using state-of-the-art, graph-isomorphism methods.
However, these methods do not take into account the special structure of the arrays
and so, they cannot be tailored to create normal forms of a specific form.


                       
Reduction to LMC normal form
----------------------------

The OApackage implements theory and methods from the article `Complete enumeration of pure-level and mixed-level orthogonal arrays, Schoen et al. <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcd.20236>`_ to
reduce orthogonal arrays to their LMC normal form. The C++ function to perform
the reduction is :cpp:func:`reduceLMCform`:

.. doctest::

  >>> import oapackage
  >>> oapackage.set_srand(1)
  >>> array = oapackage.exampleArray(1, 0).selectFirstColumns(3)
  >>> array = array.randomperm()
  >>> print('input array:'); array.transposed().showarraycompact()
  input array:
  1100010111010001
  0101100110100011
  0000001111101101
  >>> reduced_array = oapackage.reduceLMCform(array)
  >>> print('reduced array:'); reduced_array.transposed().showarraycompact()
  reduced array:
  0000000011111111
  0000111100001111
  0001011101110001

It is also possible to check whether an array is in normal form
with the :cpp:func:`LMCcheck` method:

.. doctest::
   
    >>> import oapackage
    >>> array = oapackage.exampleArray(1)
    >>> lmc_type = oapackage.LMCcheck(array)
    >>> if lmc_type == oapackage.LMC_MORE:
    ...      print('array is in minimal form')
    ... elif lmc_type == oapackage.LMC_LESS:
    ...      print('array is not in minimal form')
    array is in minimal form

.. comment
    The :cpp:func:`LMCcheck` method can also check on other normal form such as the J4 ordering or J5 ordering. 
    For example: oaextend.setAlgorithm(oapackage.MODE_J4)
    
Reduction to delete-one-factor projection form
----------------------------------------------

The article `A canonical form for non-regular arrays based on generalized word length pattern values of delete-one-factor projections <https://econpapers.repec.org/paper/antwpaper/2014007.htm>`_
:cite:`EendebakDOF` describes a canonical form of an orthogonal array based on delete-one-factor projections. 
The C++ interface to delete-one-factor projection form is :cpp:func:`reduceDOPform`.
The reduction method works well for large arrays with a large variation in the projection values.
    

.. comment
    .. doxygenfunction:: reduceDOPform(const array_link&)


An example on how to use this reduction is shown in :ref:`Example code for delete-one-factor projections`, which can be found
in the example notebooks section.
    

Reduction using graph isomorphisms
----------------------------------

The function :py:meth:`~oalib.reduceOAnauty` reduces an orthogonal array to Nauty canonical form. To reduce general graphs to Nauty canonical form, the OApackage includes the function :py:meth:`~oalib.reduceGraphNauty`.


.. admonition:: Reduce a design to normal form using Nauty
 
  .. testsetup::
     
     import oapackage
     oapackage.set_srand(1)
     
  .. doctest::
    
    >>> oapackage.set_srand(1)
    >>> al = oapackage.exampleArray(0).randomperm()
    >>> al.showarray()
    array: 
      0   0
      0   1
      1   1
      0   1
      1   0
      0   0
      1   0
      1   1
    >>> transformation=oapackage.reduceOAnauty(al, 0)
    >>> transformation.show()
    array transformation: N 8
    column permutation: {0,1}
    level perms:
    {0,1}
    {0,1}
    row permutation: {0,5,1,3,4,6,2,7}
    >>> alr=transformation.apply(al)
    >>> alr.showarray()
    array: 
      0   0
      0   0
      0   1
      0   1
      1   0
      1   0
      1   1
      1   1

.. _LMC0:

Normal forms for conference designs
-----------------------------------

For conference designs a convenient normal form is the LMC0 ordering (sometimes also called L0 ordering) :cite:`Schoen2018dsd`.

.. admonition:: LMC0 ordering

 The LMC0 ordering for conference designs is defined in three steps:
  
  Definition LMC0 i: *Order of elements*
    The LMC0 order of the factor levels -1, 0 and +1 is 0 < +1 < -1.

  Definition LMC0 ii: *Order of columns*
    A column a is smaller than a column b according to LMC0 ordering, notated as a < b, if either
    of the following conditions hold:
    
       1. If we replace the values -1 by +1 in both columns, then the first element where the
          columns differ is smaller in a than in b according to Definition 1.
       2. The zeros in column a are in the same position as the zeros in column b, and the first element where the
          columns differ is smaller in a than in b according to Definition i.
          
  Definition LMC0 iii: *Order of designs*
    Conference design A is smaller than conference design B according to LMC0 ordering, notated
    as A < B, if the first column where the designs differ is smaller in A than in B.

The definition implies that the ordering of designs is column-by-column and that the position of zeros in the columns is dominant over the values +1, -1.
To check whether a design is in LMC0 form we can use :cpp:func:`LMC0check`.

.. admonition:: Conference design in normal form
 
  .. testsetup::
     
     import oapackage
     
  .. doctest::
    
    >>> array = oapackage.exampleArray(53,1)
    exampleArray 53: third array in C(12,4)
    >>> array.showarray()
    array:
      0   1   1   1
      1   0  -1   1
      1   1   0  -1
      1   1   1  -1
      1   1   1  -1
      1   1  -1   1
      1   1  -1   1
      1  -1   1   0
      1  -1   1   1
      1  -1   1   1
      1  -1  -1  -1
      1  -1  -1  -1
    >>> oapackage.LMC0check(array) == oapackage.LMC_LESS
    True

  