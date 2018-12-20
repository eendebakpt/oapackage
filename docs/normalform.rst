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
reduce orthogonal arrays to their LMC normal form. The C++ function to perform the reduction is:

.. doxygenfunction:: reduceLMCform(const array_link&)


Reduction to delete-one-factor projection form
----------------------------------------------

The article `A canonical form for non-regular arrays based on generalized word length pattern values of delete-one-factor projections <https://econpapers.repec.org/paper/antwpaper/2014007.htm>`_
:cite:`EendebakDOF` describes a canonical form of an orthogonal array based on delete-one-factor projections. 
The C++ interface to delete-one-factor projection form is:

.. doxygenfunction:: reduceDOPform(const array_link&)


An example on how to use this reduction is shown in :ref:`Example code for delete-one-factor projections`, which can be found
in the example notebooks section.
    

Reduction using graph isomorphisms
----------------------------------

The function :py:meth:`~oalib.reduceOAnauty` reduces an orthogonal array to Nauty canonical form. To reduce general graphs to Nauty canonical form, the OApackage includes the function :py:meth:`~oalib.reduceGraphNauty`.


.. admonition:: Reduce a design to normal form using Nauty
 
   .. code-block:: python
    
    >>> al = oapackage.exampleArray(0).randomperm()
    >>> al.showarray()
    array: 
      1   1
      0   0
      0   1
      1   0
      1   0
      0   1
      1   1
      0   0
    >>> transformation=oapackage.reduceOAnauty(al, 0)
    >>> transformation.show()
    array transformation: N 8
    column permutation: {0,1}
    level perms:
    {0,1}
    {0,1}
    row permutation: {1,7,2,5,3,4,0,6}
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

