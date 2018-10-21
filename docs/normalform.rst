Normal form of arrays
======================


The Orthogonal Array package provides functions to reduce
the arrays to their canonical form, with respect to some ordering. The
default ordering is the lexicographically minimum in columns (LMC) ordering 
:cite:`Eendebak2009`. Alternative orderings include the
delete-one-factor projection ordering introduced
in :cite:`EendebakDOF` and the even-odd ordering.

Specialized packages such as Nauty :cite:`nautyI`, :cite:`nautyII` can also reduce  arrays to their canonical form using state-of-the-art, graph-isomophism methods. However, these methods do not take into account the special structure of the arrays and so, they  cannot be tailored to create normal forms of a specific form.

                       
Reduction to LMC normal form
----------------------------

The OApackage implements theory and methods from the article `Complete enumeration of pure‐level and mixed‐level orthogonal arrays, Schoen et al., <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcd.20236>`_ to reduce
orthogonal arrays to their LMC normal form. The C++ function to perform the reduction is:

.. doxygenfunction:: reduceLMCform(const array_link&)

.. comment
    .. admonition:: C++ block
    
        .. doxygenfunction:: reduceLMCform(const array_link&)

    .. sidebar:: Sidebar Title
        :subtitle: Optional Sidebar Subtitle
    
        Subsequent indented lines comprise
        the body of the sidebar, and are
        interpreted as body elements.


    
    .. topic:: C++ code
    
        .. doxygenfunction:: reduceLMCform(const array_link&)
    
    .. code-block:: c++
       :caption: Reduction to normal form
    
    
        /// Reduce an array to canonical form using LMC ordering.
        array_link reduceLMCform(const array_link &al);

Reduction to delete-one-factor projection form
----------------------------------------------

The article `A canonical form for non-regular arrays based on generalized word length pattern values of delete-one-factor projections, <https://econpapers.repec.org/paper/antwpaper/2014007.htm>`_
:cite:`EendebakDOF` describes a canonical form of an orthogonal array based on delete-one-factor projections. The C++ function function to perform the reduction in terms of this form is:

.. doxygenfunction:: reduceDOPform(const array_link&)

.. comment
    .. code-block:: c++
       :caption: C++ interface to delete-one-factor projection form
    
        /// reduce an array to canonical form using delete-1-factor ordering
        array_link reduceDOPform(const array_link &al);
    
An example on how to use this reduction is shown in :ref:`Example code for delete-one-factor projections`, which can be found
in the example notebooks section.

Reduction using graph isomorphisms
----------------------------------

The function :py:meth:`~oalib.reduceOAnauty` reduces an orthogonal array to Nauty canonical form. To reduce general graphs to Nauty canonical form, the OApackage includes the function :py:meth:`~oalib.reduceGraphNauty`.

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

