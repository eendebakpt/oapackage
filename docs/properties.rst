Properties of designs
=====================

Statistical properties of an array
----------------------------------

Most properties of an array can be calculated using the \|array\_link\|
object. The interface is listed below.

::

    struct array_link
    {
        ... 
        
    public:
        // statistical properties of the array

        /// calculate D-efficiency
        double Defficiency() const;
        /// calculate main effect robustness (or Ds optimality)
        double DsEfficiency(int verbose=0) const;
        /// calculate A-efficiency
        double Aefficiency() const;
        /// calculate E-efficiency
        double Eefficiency() const;
        /// calculate rank of array
        int rank() const;
        /// calculate generalized wordlength pattern
        std::vector<double> GWLP() const;
        /// return true if the array is a foldover array
        bool foldover() const;
        /// calculate centered L2 discrepancy
        double CL2discrepancy() const;
        /// Calculate the projective estimation capacity sequence
        std::vector<double> PECsequence() const;
    }

The :math:`D`-efficiency, :math:`A`-efficiency and :math:`E`-efficiency
are calculated by calculating the SVD of the second order interaction
matrix. The efficiencies can then be calculated using the eigenvalues of
the SVD. For the definition of the :math:`D`-, :math:`A`- and
:math:`E`-efficiency see Definition :ref:`DAE`. For the
rank of a matrix the LU decomposition of the matrix is calculated using
the Eigen package :cite:`eigenweb`.

.. topic:: D-efficiency and average VIF

   Let :math:`X` be an :math:`N\times k` :math:`2`-factor
   array with second order model :math:`{F(X)}`. Then we define the :math:`{{\color{darkblue}D}}`-efficiency and the average variance inflation factor as

   .. math::
       :name: DAE
    
       {{\color{darkblue}D}(X)} = \left( \det {F(X)}^T {F(X)}\right)^{1/m} / N , 
       \label{formula:Defficiency} \\
       {\mathrm{VIF}(X)} = N \operatorname{tr}\left( \frac{1}{ {F(X)}^T {F(X)}} \right) /m . \label{formula:VIF}
       
   The matrix :math:`{F(X)}^T {F(X)}` is called the information matrix. Let :math:`\lambda_1, \ldots, \lambda_m` be the eigenvalues of the information matrix. Then the :math:`{{\color{darkblue}E}}`-efficiency of a matrix is [definition:Eefficiency] defined as

   .. math::
       :name: Eefficiency
       {{\color{darkblue}E}(X)} = \min_j \lambda_j .
       \label{formula:E-efficiency}

Note that in terms of the eigenvalues we have
:math:`{{\color{darkblue}D}(X)} = (\prod_j \lambda_j)^{1/m} / N` and
:math:`{\mathrm{VIF}(X)} = N (\sum_j \lambda_j^{-1})/m`.

The :math:`D_s`-effiency is the main effect robustness, see the appendix
in :cite:`Schoen2010` for more details.


GWLP and J-characteristics
--------------------------

From an \|array\_link\| object we can calculate the generalized
worldlength patterns :cite`Xu2001`, :math:`F`-values and
:math:`J`-characteristics.

.. code-block:: python
 :caption: Calculate GWLP and :math:`F`-values 
   
 >>> al=oapackage.exampleArray(1)
 >>> al.showarray() array: 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 0 1 0 1 0 0 1 1 0 0 0 1 1 1 1 0 1 1 1 1 1 0 0 1 1 1 0 1 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0 0 1 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0
 >>> g=al.GWLP() >>> print(’GWLP: GWLP: (1.0, 0.0, 0.0, 1.0, 1.0, 0.0)
 >>> print('F3-value: %s' % str(al.Fvalues(3)))
 F3-value: (4, 6)
 >>> print('F4-value: %s' % str(al.Fvalues(4)))
 F4-value: (1, 4)
 >>> print('J3-characteristics: %s’ % str(al.Jcharacteristics(3)))
 J3-characteristics: (8, 8, 0, 0, 0, 8, 0, 8, 0, 0)



MD5 sums
--------

To check data structures on disk the packages includes functions to
generate MD5 sums of designs. 

.. code-block:: python
 :caption: Calculate md5 sum of a design

 >>> import oapackage; al=oapackage.exampleArray(0)
 >>> al.md5()
 '6454c492239a8e01e3c01a864583abf2'

The C++ functions are:

.. code-block:: c

    /// calculate md5 sum of a data block in memory
    std::string md5(void *data, int numbytes);
    /// calculate md5 sum of a file on disk
    std::string md5(const std::string filename);

