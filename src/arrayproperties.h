/*! \file arrayproperties.h
 \brief Contains functions to calculate properties of arrays

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#ifndef ARRAYPROPERTIES_H
#define ARRAYPROPERTIES_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "arraytools.h"
#include "oaoptions.h"
#include "tools.h"

#define stringify(name) #name

/// Calculate D-efficiency and VIF-efficiency and E-efficiency values using SVD
void DAEefficiencyWithSVD (const Eigen::MatrixXd &x, double &Deff, double &vif, double &Eeff, int &rank, int verbose);

/** Calculate the rank of the second order interaction matrix of an orthogonal array
 *
 * The model is the intercept, main effects and interaction effects
 * The rank, D-efficiency, VIF-efficiency and E-efficiency are appended to the second argument
 *
 * The vector ret is filled with the rank, Defficiency, VIF efficiency and Eefficiency
 */
int array_rank_D_B (const array_link &al, std::vector< double > *ret = 0, int verbose = 0);

/// Calculate D-efficiency for a 2-level array using symmetric eigenvalue decomposition
double Defficiency (const array_link &al, int verbose = 0);

std::vector< double > Defficiencies (const array_link &al, const arraydata_t &arrayclass, int verbose = 0,
                                     int addDs0 = 0);

/// Calculate VIF-efficiency of matrix
double VIFefficiency (const array_link &al, int verbose = 0);

/// Calculate A-efficiency of matrix
double Aefficiency (const array_link &al, int verbose = 0);

/// Calculate E-efficiency of matrix (1 over the VIF-efficiency)
double Eefficiency (const array_link &al, int verbose = 0);

/// calculate various A-efficiencies
std::vector< double > Aefficiencies (const array_link &al, int verbose = 0);

#ifdef FULLPACKAGE
/// Return the D-efficiencies for the projection designs
std::vector< double > projDeff (const array_link &al, int kp, int verbose);

/// Return the projection estimation capacity sequence of a design
std::vector< double > PECsequence (const array_link &al, int verbose = 0);
#endif

/// Return the distance distribution of a design
std::vector< double > distance_distribution (const array_link &al);

/// Calculate J-characteristics of matrix (the values are signed)
std::vector< int > Jcharacteristics (const array_link &al, int jj = 4, int verbose = 0);

/** @brief Calculate GWLP (generalized wordlength pattern)
 *
 *  The method used for calculation is from Xu and Wu (2001), "Generalized minimum aberration for asymmetrical
 * fractional factorial desings"
 *  The non-symmetric arrays see "Algorithmic Construction of Efficient Fractional Factorial Designs With Large Run
 * Sizes", Xu
 *
 * \param al Array to calculate the GWLP value for
 * \param verbose Verbosity level
 * \param truncate If True then round values near zero to solve double precision errors
 */
std::vector< double > GWLP (const array_link &al, int verbose = 0, int truncate = 1);

/** @brief Calculate GWLP (generalized wordlength pattern)
*
*  The method used for calculation is from Xu and Wu (2001), "Generalized minimum aberration for asymmetrical
* fractional factorial desings"
*  The non-symmetric arrays see "Algorithmic Construction of Efficient Fractional Factorial Designs With Large Run
* Sizes", Xu
*
* \param al Array to calculate the GWLP value for
* \param verbose Verbosity level
* \param truncate If True then round values near zero to solve double precision errors
*
*/
std::vector< double > GWLPmixed (const array_link &al, int verbose = 0, int truncate = 1);

// SWIG has some issues with typedefs, so we use a define
#define GWLPvalue mvalue_t< double >

typedef mvalue_t< double > DOFvalue;

/// calculate delete-one-factor GWLP (generalized wordlength pattern) projections
std::vector< GWLPvalue > projectionGWLPs (const array_link &al);

std::vector< GWLPvalue > sortGWLP (std::vector< GWLPvalue >);

/// calculate delete-one-factor GWLP (generalized wordlength pattern) projection values
std::vector< double > projectionGWLPvalues (const array_link &al);

/** calculate centered L2-discrepancy of a design
 *
 * The method is from "A connection between uniformity and aberration in regular fractions of two-level factorials",
 * Fang and Mukerjee, 2000
 */
double CL2discrepancy (const array_link &al);

/** Calculate second order interaction model for 2-level array
*
* \param array Array to calculate second order interaction model from
* \returns Array interaction effects
*/
array_link array2secondorder (const array_link &al);

/** calculate second order interaction model for 2-level array
 *
 * \param array Array to calculate second order interaction model from
 * \returns Array with intercept, main effects and interaction effects
 */
array_link array2xf (const array_link &al);

/** calculate second order interaction model for 2-level array
*
* \param array Array to calculate second order interaction model from
* \returns Array with intercept, main effects and interaction effects
*/
Eigen::MatrixXd array2xfeigen (const array_link &al);

/// return rank of an array based on FullPivHouseholderQR
int arrayrankFullPivQR (const array_link &al, double threshold = -1);

/// return rank of an array based on ColPivHouseholderQR
int arrayrankColPivQR (const array_link &al, double threshold = -1);

/// return rank of an array based on arrayrankFullPivLU
int arrayrankFullPivLU (const array_link &al, double threshold = -1);

/// return rank of an array based on Eigen::JacobiSVD
int arrayrankSVD (const array_link &al, double threshold = -1);

/// calculate the rank of an array
int arrayrank (const array_link &al);

/// Return rank of an array. Information about the different methods for rank calculation is printed to stdout
int arrayrankInfo (const Eigen::MatrixXd &, int verbose = 1);

/// Return rank of an array. Information about the different methods for rank calculation is printed to stdout
int arrayrankInfo (const array_link &al, int verbose = 1);

/// convert array_link to Eigen matrix
Eigen::MatrixXd arraylink2eigen (const array_link &al);

/** Structure to efficiently calculate the rank of the second order interaction matrix of many arrays sharing a common
 *subarray
 *
 * The input arrays are assumed to be of the form A_i = [A_0 X_i]
 *
 *
 **/
class rankStructure {
      public:
        typedef Eigen::FullPivHouseholderQR< Eigen::MatrixXd > EigenDecomp;

      public:
        array_link alsub;
        int r;
        /// verbosity level
        int verbose;
        /// number of columns of subarray in cache
        int ks;
        /// number of columns to subtract from array when updating cache
        int nsub;
        /// used for debugging
        int id;

      private:
        /// decomposition of subarray
        EigenDecomp decomp;
        Eigen::MatrixXd Qi;

        /// internal structure
        long ncalc, nupdate;

      public:
        /// constructor
        rankStructure (const array_link &al, int nsub = 3, int verbose = 0) {
                this->verbose = verbose;
                ks = 0;
                this->nsub = nsub;
                this->id = -1;
                ncalc = 0;
                nupdate = 0;
                updateStructure (al);
        }
        /// constructor
        rankStructure (int nsub = 3, int id = -1) {
                this->verbose = 0;
                ks = 0;
                this->nsub = nsub;
                this->id = id;
                ncalc = 0;
                nupdate = 0;

                array_link al = exampleArray (1);
                updateStructure (al);
        }

        void info () const {
                printf ("	rankStructure: submatrix %dx%d, rank %d, rank of xf %d\n", alsub.n_rows,
                        alsub.n_columns, this->alsub.rank (), (int)decomp.rank ());
        }

        /// update the structure cache with a new array
        void updateStructure (const array_link &al) {
                this->alsub = al;
                this->ks = al.n_columns;
                Eigen::MatrixXd A = array2xf (al).getEigenMatrix ();
                decomp.compute (A);

                this->Qi = decomp.matrixQ ().inverse ();

                nupdate++;

                if (this->verbose >= 1 && nupdate % 30 == 0) {
                        printfd ("updateStructure: ncalc %d, nupdate %d\n", ncalc, nupdate);
                }
        }

        /// helper function
        EigenDecomp::PermutationType matrixP () const { return decomp.colsPermutation (); }

        /// calculate the rank of an array directly, uses special threshold
        int rankdirect (const Eigen::MatrixXd &A) const {
                EigenDecomp decomp (A);
                decomp.setThreshold (1e-12);

                if (0) {
                        printf ("rankdirect: threshold: %e, rank %d\n", decomp.threshold (), (int)decomp.rank ());
                }
                int rank = decomp.rank ();
                return rank;
        }

        /// calculate the rank of the second order interaction matrix of an array directly
        int rankxfdirect (const array_link &al) const {
                Eigen::MatrixXd mymatrix = arraylink2eigen (array2xf (al)); // array2xf;
                return rankdirect (mymatrix);
        }

        /// calculate the rank of the second order interaction matrix of an array using the cache system
        int rankxf (const array_link &al);
};

/// Return the condition number of a matrix
double conditionNumber (const array_link &M);

#ifdef FULLPACKAGE

#include "pareto.h"

enum paretomethod_t { PARETOFUNCTION_DEFAULT, PARETOFUNCTION_J5 };

/** Calculate the Pareto optimal arrays from a list of array files

    Pareto optimality is calculated according to (rank; A3,A4; F4)
*/
void calculateParetoEvenOdd (const std::vector< std::string > infiles, const char *outfile, int verbose = 1,
                             arrayfilemode_t afmode = ABINARY, int nrows = -1, int ncols = -1,
                             paretomethod_t paretomethod = PARETOFUNCTION_DEFAULT);

// Calculate the Pareto optimal desings from a list of arrays (rank; A3,A4; F4)
Pareto< mvalue_t< long >, long > parsePareto (const arraylist_t &arraylist, int verbose,
                                              paretomethod_t paretomethod = PARETOFUNCTION_DEFAULT);

/** calculate A3 and A4 value for array
 *
 * \param al Array for which to calculate A3 and A4
 * \returns Object with A3 and A4
 */
mvalue_t< long > A3A4 (const array_link &al);

/// calculate F4 value for 2-level array
inline mvalue_t< long > F4 (const array_link &al, int verbose = 1) {
        jstruct_t js (al, 4);
        std::vector< int > FF = js.calculateF ();
        if (verbose >= 3) {
                myprintf ("  parseArrayPareto: F (high to low): ");
                display_vector (FF);
                myprintf ("\n");
        }

        mvalue_t< long > v (FF, mvalue_t< long >::LOW);
        return v;
}

template < class IndexType >
/** Add array to list of Pareto optimal arrays
 *
 * The values to be optimized are:
 *
 * 1) Rank (higher is better)
 * 2) A3, A4 (lower is better)
 * 3) F4 (lower is better, sum of elements is constant)
 *
 * Valid for 2-level arrays of strength at least 3
 *
 * */
inline typename Pareto< mvalue_t< long >, IndexType >::pValue calculateArrayParetoRankFA (const array_link &al,
                                                                                          int verbose) {
        int N = al.n_rows;
        int model_rank = arrayrankFullPivLU (array2secondorder (al), 1e-12) + 1 +
                         al.n_columns; // valid for 2-level arrays of strength at least 3
        mvalue_t< long > a3a4_values = A3A4 (al);
        mvalue_t< long > f4 = F4 (al);

        // add the 3 values to the combined value
        typename Pareto< mvalue_t< long >, IndexType >::pValue p;
        p.push_back (model_rank);
        p.push_back (a3a4_values);
        p.push_back (f4);

        if (verbose >= 2) {
                if (verbose >= 3) {
                        std::vector< double > gwlp = al.GWLP ();
                        myprintf ("parseArrayPareto: A4 (scaled) %ld, %f\n", a3a4_values.v[1], gwlp[4]);
                }

                myprintf ("  parseArrayPareto: rank %d, verbose %d\n", al.rank (), verbose);
        }

        return p;
}

/// add Jmax criterium to Pareto set
template < class IndexType >
void addJmax (const array_link &al, typename Pareto< mvalue_t< long >, IndexType >::pValue &p, int verbose = 1) {
        std::vector< int > j5 = al.Jcharacteristics (5);
        int j5max = vectormax (j5, 0);

        int v1 = (j5max == al.n_rows);
        int v2 = 1 - v1;

        if (verbose >= 3) {
                printf ("calculateArrayParetoJ5: j5max %d: %d %d\n", j5max, v1, v2);
        }

        p.push_back (v1);
        p.push_back (v2);
}

template < class IndexType >
typename Pareto< mvalue_t< long >, IndexType >::pValue calculateArrayParetoJ5 (const array_link &al, int verbose) {
        typename Pareto< mvalue_t< long >, IndexType >::pValue p =
            calculateArrayParetoRankFA< IndexType > (al, verbose);
        addJmax< IndexType > (al, p, verbose);

        return p;
}

template < class IndexType >
/** Add array to list of Pareto optimal arrays
 *
 * The values to be optimized are:
 *
 * 1) Rank (higher is better)
 * 2) A3, A4 (lower is better)
 * 3) F4 (?? is better, sum of elements is constant)
 *
 * */
inline void parseArrayPareto (const array_link &al, IndexType i, Pareto< mvalue_t< long >, IndexType > &pset,
                              int verbose) {
        typename Pareto< mvalue_t< long >, IndexType >::pValue p =
            calculateArrayParetoRankFA< IndexType > (al, verbose);

        // add the new tuple to the Pareto set
        pset.addvalue (p, i);
}

#endif // FULLPACKAGE

/// convert C value to D-efficiency value
inline double Cvalue2Dvalue (double C, int ka) {
        double ma = 1 + ka + ka * (ka - 1) / 2;
        double A = pow (C, 1. / ma);

        return A;
}

/// convert D-efficiency value to C value
inline double Dvalue2Cvalue (double A, int ka) {
        int ma = 1 + ka + ka * (ka - 1) / 2;
        double C = pow (A, ma);

        return C;
}

/// Return index of an array
inline int get_oaindex (const array_t *s, const colindex_t strength, const colindex_t N) {
        int oaindex = N;
        for (colindex_t z = 0; z < strength; z++) {
                oaindex /= s[z];
        }

        return oaindex;
}

#endif // ARRAYPROPERTIES_H
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;
