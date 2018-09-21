Generation of designs
=====================

The package can be used to generate several classes of designs. Generated
designs are available on the website http://www.pietereendebak.nl/oapackage/index.html

The main functions:

.. autosummary::

        oapackage.Doptim.Doptimize
        oapackage.oalib.generateConferenceExtensions
        
        
Generation of orthogonal arrays
-------------------------------

A list of arrays in LMC form can be extended to a list of arrays in LMC
form with one additional column. Details for the algorithm are described
in :cite:`Eendebak2009`.

The main functions for array extension are the following:

.. doxygenfunction:: extend_arraylist(const arraylist_t&, arraydata_t&, OAextend const&)
.. doxygenfunction:: extend_arraylist(const arraylist_t&, const arraydata_t&)

                     
.. comment
    .. code-block:: c++
       :caption: C++ interface
       
        /// extend a list of arrays
        arraylist_t & extend_arraylist(arraylist_t & alist, arraydata_t &fullad, 
                    OAextend const &oaextendoptions);

Here :meth:`~oapackage.oalib.arraydata_t` is the structure describing the type of arrays and
:meth:`~oapackage.oalib.OAextend` contains various options for the algorithm.

A typical session could be:

.. code-block:: python
   :caption: Extend an array
   
   >>> N=8; ncols=3;
   >>> arrayclass=oapackage.arraydata_t(2, N, 2, ncols)
   >>> al=arrayclass.create_root() 
   >>> al.showarray()
   array: 0 0 0 0 0 1 0 1 1 0 1 0 1 1 1 1
   >>> 
   >>> alist=oapackage.extend_array(al, arrayclass)
   >>> for al in alist:
   ... al.showarray()
   array: 0 0 0 0 0 0 0 1 1 0 1 1 1 0 1 1 0 1 1 1 0 1 1 0
   array: 0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1

Even-odd
--------

The even-odd arrays are a special class of orthognal arrays with at least one of the odd :math:`J`-characteristics unequal to zero.
More information on this class of designs will be written later.

Conference designs
------------------

An :math:`n\times k` conference design is an :math:`N\times k` matrix
with entries 0, -1, +1 such that i) in each column the symbol 0 occurs
exactly one time and ii) all columns are orthogonal to each other. A
more detailed description is given
in :cite:`wiki:ConferenceMatrix`.

.. code-block:: python
 :caption: Generate conference designs with 8 rows
                    
 >>> import oapackage
 >>> ctype=oapackage.conference_t(N=8, k=8)
 >>> al = ctype.create_root_three()
 >>> al.showarray() array: 0 1 1 1 0 -1 1 1 0 1 1 1 1 1 -1 1 -1 1 1 -1 1 1 -1 -1
 >>> l4=oapackage.extend_conference ([al], ctype, verbose=0)
 >>> l5=oapackage.extend_conference ( l4, ctype,verbose=0) 
 >>> l6=oapackage.extend_conference ( l5, ctype, verbose=0)
 >>>
 >>> print('number of non-isomorphic conference designs: number of conference designs: %d'  % len(l6) )
 non-isomorphic conference designs: 11


The full interface for conference designs is available
in the :ref:`Interface for conference designs`.

.. comment 
    .. doxygenfile:: conference.h

Calculation of D-optimal designs
--------------------------------

D-optimal designs can be calculated with the function :py:meth:`oapackage.Doptim.Doptimize`.
This function uses a coordinate exchange algorithm to generate designs
with good properties for the :math:`D`-efficiency.

An example script with Python to generate optimal designs with 40 runs
and 7 factors is shown below.

.. code-block:: python
 :caption: Doptimize
 
 >>> N=40; s=2; k=7;
 >>> arrayclass=oapackage.arraydata_t(s, N, 0, k) 
 >>> print('We generate optimal designs with: %s' % arrayclass)
 We generate optimal designs with: arrayclass: N 40, k 7, strength 0, s 2,2,2,2,2,2,2, order 0.
 >>> alpha=[1,2,0] 
 >>> method=oapackage.DOPTIM_UPDATE 
 >>> scores, dds, designs, ngenerated = oapackage.Doptimize(arrayclass, nrestarts=40, optimfunc=alpha, selectpareto=True)
 Doptim: optimization class 40.2-2-2-2-2-2-2
 Doptimize: iteration 0/40
 Doptimize: iteration 39/40 Doptim: done (8 arrays, 0.6 [s]) 
 >>> print('Generated %d designs, the best D-efficiency is %.4f’ % (len(designs), dds[:,0].max() ))
 Generated 8 designs, the best D-efficiency is 0.9098

The parameters of the function are documented in the code.

To calculate properties of designs we can use the following functions.
For :math:`D`-efficiencies we can use

.. doxygenfunction:: array_link::Defficiencies(int)
    :no-link:
    :outline:

.. comment
    .. code-block:: c++
    
        std::vector<double> array_link::Defficiencies ( int verbose ) const;

to calculate the :math:`D`-, :math:`D_s`- and :math:`D_1`-efficiency.
For details see :cite:`EendebakSO`.

The projective estimation capacity (PEC) sequence
from :cite:`loeppky2004ranking` can be calculated with:

.. doxygenfunction:: PECsequence(const array_link&, int)
    :no-link:
    :outline:
.. doxygenfunction:: array_link::PECsequence()
    :no-link:
    :outline:

.. figure:: images/motivating-40-d-2-2-2-2-2-2-2-scatterplot-ndata2.png

   Scatterplot for the :math:`D`-efficiency and :math:`D_s`-efficiency
   for generated designs in :math:`{\operatorname{OA}(40; 2; 2^7)}`. The
   Pareto optimal designs are colored, while the non-Pareto optimal
   designs are grey. For reference the strength-3 orthogonal array with
   highest D-efficiency is also included in the plot.
