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

   Let :math:`{{\color{darkblue}X}}` be an :math:`N\times k` :math:`2`-factor
   array with second order model :math:`{F({{\color{darkblue}X}})}`. Then we define the :math:`{{\color{darkblue}D}}`-efficiency and the average variance inflation factor as

   .. math::
       :name: DAE
    
       {{\color{darkblue}D}({{\color{darkblue}X}})} = \left( \det {F({{\color{darkblue}X}})}^T {F({{\color{darkblue}X}})}\right)^{1/m} / N , 
       \label{formula:Defficiency} \\
       {\mathrm{VIF}({{\color{darkblue}X}})} = N \operatorname{tr}\left( \frac{1}{ {F({{\color{darkblue}X}})}^T {F({{\color{darkblue}X}})}} \right) /m . \label{formula:VIF}
       
   The matrix :math:`{F({{\color{darkblue}X}})}^T {F({{\color{darkblue}X}})}` is called the information matrix. Let :math:`\lambda_1, \ldots, \lambda_m` be the eigenvalues of the information matrix. Then the :math:`{{\color{darkblue}E}}`-efficiency of a matrix is [definition:Eefficiency] defined as

   .. math::
       :name: Eefficiency
       {{\color{darkblue}E}({{\color{darkblue}X}})} = \min_j \lambda_j .
       \label{formula:E-efficiency}

Note that in terms of the eigenvalues we have
:math:`{{\color{darkblue}D}(X)} = (\prod_j \lambda_j)^{1/m} / N` and
:math:`{\mathrm{VIF}(X)} = N (\sum_j \lambda_j^{-1})/m`.

The :math:`D_s`-effiency is the main effect robustness, see the appendix
in :cite:`Schoen2010` for more details.
