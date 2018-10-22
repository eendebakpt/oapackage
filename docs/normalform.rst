Normal form of designs
======================


If we introduce an ordering a set of arrays, then for each
isomorphism class of arrays the minimal element defines a unique
canonical form. The Orthogonal Array package contains functions to reduce
designs to canonical form with respect to some ordering. The
default ordering for arrays is the lexicographic ordering in
columns :cite:`Eendebak2009`. An alternative ordering is the
delete-one-factor projection ordering as described
in :cite:`EendebakDOF` or the even-odd ordering.

Another approach to generation of canonical forms for designs is to use
graph-isomorphism packages such as
Nauty :cite:`nautyI`, :cite:`nautyII`. An advantage of this approach is 
that the graph isomorphism packages contain state-of-the-art methods for reduction to canonical form.
A disadvantage is that these methods are unaware of the special structure of designs and cannot be tailored
to create normal forms of a specific form.

                       
Reduction to LMC normal form
----------------------------

Based on theory from `Complete enumeration of pure-level and mixed-level orthogonal arrays, Schoen et al. <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcd.20236>`_ we can reduce
orthogonal arrays to LMC normal form. The C++ function to perform reduction is

.. doxygenfunction:: reduceLMCform(const array_link&)


Reduction to delete-one-factor projection form
----------------------------------------------

The canonical form is described in `A canonical form for non-regular arrays based on generalized word length pattern values of delete-one-factor projections <https://econpapers.repec.org/paper/antwpaper/2014007.htm>`_
:cite:`EendebakDOF`.

An example with the methods is :ref:`Example code for delete-one-factor projections` which can be found
in the example notebooks section.


The C++ interface to delete-one-factor projection form is:

.. doxygenfunction:: reduceDOPform(const array_link&)

.. comment
    .. code-block:: c++
       :caption: C++ interface to delete-one-factor projection form
    
        /// reduce an array to canonical form using delete-1-factor ordering
        array_link reduceDOPform(const array_link &array);
    

Reduction using graph isomorphisms
----------------------------------

To reduce a general graph to Nauty canonical form one can use :py:meth:`~oalib.reduceGraphNauty`. For orthogonal arrays we can
encode the array structure as a graph. The reduction can then be done
with :py:meth:`~oalib.reduceOAnauty`.


.. code-block:: python
   :caption: Reduce a design to normal form using Nauty
   
   >>> al = oapackage.exampleArray(0).randomperm()
   >>> al.showarray()
   array: 0 1 1 1 1 0 0 0 0 0 1 1 1 0 0 1
   >>> t=oapackage.reduceOAnauty(al, verbose=0)
   >>> t.show()
   array transformation: N 8 column permutation: 0,1 level perms: 0,1 0,1 row permutation: 3,4,0,7,2,6,1,5
   >>> alr=t.apply(al)
   >>> alr.showarray()
   array: 0 0 0 0 0 1 0 1 1 0 1 0 1 1 1 1

