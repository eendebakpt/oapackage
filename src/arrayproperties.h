/*! \file arrayproperties.h
 \brief Contains functions to calculate properties of arrays

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "arraytools.h"
#include "oaoptions.h"
#include "tools.h"

/// Calculate D-efficiency and VIF-efficiency and E-efficiency values using SVD
void DAEefficiencyWithSVD (const Eigen::MatrixXd &secondorder_interaction_matrix, double &Deff, double &vif, double &Eeff, int &rank, int verbose);

/** Calculate the rank of the second order interaction matrix of an orthogonal array
 *
 * The model is the intercept, main effects and interaction effects
 * The rank, D-efficiency, VIF-efficiency and E-efficiency are appended to the second argument
 *
 * The return vector is filled with the rank, Defficiency, VIF efficiency and Eefficiency
 */
int array2rank_Deff_Beff (const array_link &al, std::vector< double > *ret = 0, int verbose = 0);

/// Calculate D-efficiency for a 2-level array using symmetric eigenvalue decomposition
double Defficiency (const array_link &orthogonal_array, int verbose = 0);

/** Calculate efficiencies for an array
 *
 * \param array Array to use in calculation
 * \param arrayclass Specification of the array class
 * \param verbose Verbosity level
 * \param addDs0 If True, then add the Ds0-efficiency to the output
 * \returns Vector with the calculate D-efficiency, the main effect robustness (or Ds-optimality) and D1-efficiency for an orthogonal array
 */
std::vector< double > Defficiencies (const array_link &array, const arraydata_t &arrayclass, int verbose = 0,
                                     int addDs0 = 0);

/// Calculate VIF-efficiency of matrix
double VIFefficiency (const array_link &orthogonal_array, int verbose = 0);

/// Calculate A-efficiency of matrix
double Aefficiency (const array_link &orthogonal_array, int verbose = 0);

/// Calculate E-efficiency of matrix (1 over the VIF-efficiency)
double Eefficiency (const array_link &orthogonal_array, int verbose = 0);

/// calculate various A-efficiencies
std::vector< double > Aefficiencies (const array_link &orthogonal_array, int verbose = 0);

/** Calculate D-efficiencies for all projection designs
*
* \param array Design to calculate D-efficiencies for
* \param number_of_factors Number of factors into which to project
* \param verbose Verbosity level
* \returns Vector with calculated D-efficiencies
*/
std::vector< double > projDeff (const array_link &array, int number_of_factors, int verbose = 0);

/** Calculate the projection estimation capacity sequence for a design
*
* \param array Input array
* \param verbose Verbosity level
* \returns Vector with the caculated PEC sequence
*
* The PECk of a design is the fraction of estimable second-order models in k factors.
* The vector (PEC1, PEC2, ..., ) is called the projection estimation capacity sequence.
* See "Ranking Non-regular Designs", J.L. Loeppky, 2004.
*
*/
std::vector< double > PECsequence (const array_link &array, int verbose = 0);

/**Calculate the projection information capacity sequence for a design.
*
* \param array Input array
* \param verbose Verbosity level
* \returns Vector with the caculated PIC sequence
* 
* The PICk of a design is the average D-efficiency of estimable second-order models in k factors. The vector
* (PIC1, PIC2, ..., ) is called the PIC sequence.
*
*/
std::vector< double > PICsequence(const array_link &array, int verbose = 0);

/** Return the distance distribution of a design
 *
 * The distance distribution is described in "Generalized minimum aberration for asymmetrical fractional factorial designs", Wu and Xu, 2001
 */
std::vector< double > distance_distribution (const array_link &array);

/** Calculate Jk-characteristics of a matrix
 *
 * The calcualted Jk-values are signed.
 * 
 * \param array Array to calculate Jk-characteristics for
 * \param number_of_columns Number of columns
 * \param verbose Verbosity level
 * \returns Vector with calculated Jk-characteristics
 */
std::vector< int > Jcharacteristics (const array_link &array, int number_of_columns = 4, int verbose = 0);

/** @brief Calculate GWLP (generalized wordlength pattern)
 *
 * The method used for calculation is from Xu and Wu (2001), "Generalized minimum aberration for asymmetrical
 * fractional factorial desings". For non-symmetric arrays see "Algorithmic Construction of Efficient Fractional Factorial Designs With Large Run
 * Sizes", Xu, Technometrics, 2009.
 *
 * \param array Array to calculate the GWLP value for
 * \param verbose Verbosity level
 * \param truncate If True then round values near zero to solve double precision errors
 * \returns Vector with calculated generalized wordlength pattern
 * 
 * A more detailed description of the generalized wordlength pattern can also be found in the documentation at https://oapackage.readthedocs.io/.
 */
std::vector< double > GWLP (const array_link &array, int verbose = 0, int truncate = 1);

/** @brief Calculate GWLP (generalized wordlength pattern) for mixed-level arrays
*
*  The method used for calculation is from "Algorithmic Construction of Efficient Fractional Factorial Designs With Large Run
* Sizes", Xu, Technometrics, 2009.
*
* \param array Array to calculate the GWLP value for
* \param verbose Verbosity level
* \param truncate If True then round values near zero to solve double precision errors
* \returns Vector with calculated generalized wordlength pattern
*
*/
std::vector< double > GWLPmixed (const array_link &array, int verbose = 0, int truncate = 1);

// SWIG has some issues with typedefs, so we use a define
#define GWLPvalue mvalue_t< double >

/// delete-one-factor projection value
typedef mvalue_t< double > DOFvalue;

/// calculate delete-one-factor GWLP (generalized wordlength pattern) projections
std::vector< GWLPvalue > projectionGWLPs (const array_link &al);

/// sort a list of GWLP values and return the sorted list
std::vector< GWLPvalue > sortGWLP (std::vector< GWLPvalue >);


/** Calculate centered L2-discrepancy of a design
 *
 * The method is from "A connection between uniformity and aberration in regular fractions of two-level factorials",
 * Fang and Mukerjee, 2000
 */
double CL2discrepancy (const array_link &array);

/** Calculate second order interaction model for 2-level array
*
* \param array Array to calculate second order interaction model from
* \returns Array interaction effects
*/
array_link array2secondorder (const array_link &array);

/** calculate second order interaction model for 2-level array
 *
 * \param array Array to calculate second order interaction model from
 * \returns Array with intercept, main effects and interaction effects
 */
array_link array2xf (const array_link &array);

enum model_matrix_t {
	/// only the intercept
	MODEL_CONSTANT,
	/// intercept and main effects
	MODEL_MAIN,
	/// intercept, main effects and second order interactions
	MODEL_INTERACTION,
	/// intercept, main effects and second order effects(interactions and quadratic effects)
	MODEL_SECONDORDER,
	/// invalid model
	MODEL_INVALID
};

/** Calculate model matrix for a conference design
 *
 * \param conference_design Conference design
 * \param mode Can be 'm' for main effects, 'i' for interaction effects or 'q' for quadratic effects
 * \param verbose Verbosity level
 * \returns Calculated model matrix
 */
array_link conference_design2modelmatrix(const array_link & conference_design, const char*mode, int verbose= 0);

/** Convert orthogonal array or conference design to model matrix
 *
 * The model matrix consists of the intercept, main effects and (optionally) the interaction effects and quadratic effects.
 * The order in the interaction effects is (c1, c2)=(0,0), (1,0), (2,0), (2,1), ... with c2<c1 for columns c1, c2.
 * The size of the model matrix calculated by this function is given by @ref array2modelmatrix_sizes.
 *
 * \param array Orthogonal array or conference design
 * \param mode Type of model matrix to calculate. Can be 'm' for main effects, 'i' for interaction effects or 'q' for quadratic effects
 * \param verbose Verbosity level
 * \returns Calculated model matrix
 * 
 * For conference designs the method @ref conference_design2modelmatrix is used. For orthogonal array the calculated is performed with @ref array2eigenModelMatrixMixed.
 * 
 */
Eigen::MatrixXd array2modelmatrix(const array_link &array, const char *mode, int verbose = 0);


/** Return the sizes of the model matrices calculated
 * 
 * \param array Orthogonal array or conference designs
 * \returns List with the sizes of the model matrix for: only intercept; intercept, main; intercept, main, and iteraction terms, intercept, main and full second order
 */
std::vector<int> array2modelmatrix_sizes(const array_link &array);

/** calculate second order interaction model for 2-level array
*
* \param array Array to calculate second order interaction model from
* \returns Array with intercept, main effects and interaction effects
*/
Eigen::MatrixXd array2xfeigen (const array_link &array);

/// return rank of an array based on Eigen::FullPivHouseholderQR
int arrayrankFullPivQR (const array_link &al, double threshold = -1);

/// return rank of an array based on Eigen::ColPivHouseholderQR
int arrayrankColPivQR (const array_link &al, double threshold = -1);

/// return rank of an array based on Eigen::FullPivLU
int arrayrankFullPivLU (const array_link &al, double threshold = -1);

/// return rank of an array based on Eigen::JacobiSVD
int arrayrankSVD (const array_link &al, double threshold = -1);

/// calculate the rank of an array
int arrayrank (const array_link &array);

/// Return rank of an array. Information about the different methods for rank calculation is printed to stdout
int arrayrankInfo (const Eigen::MatrixXd &, int verbose = 1);

/// Return rank of an array. Information about the different methods for rank calculation is printed to stdout
int arrayrankInfo (const array_link &array, int verbose = 1);

/// convert array_link to Eigen matrix
Eigen::MatrixXd arraylink2eigen (const array_link &array);

/** Structure to efficiently calculate the rank of the second order interaction matrix of many arrays
 *
 * The efficiency is obtained if the arrays share a common subarray. The theory is described in "Efficient rank calculation for matrices with a common submatrix", Eendebak, 2016
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

		/// print information about the rank structure
		void info() const;

        /// update the structure cache with a new array
		void updateStructure(const array_link &al);

        /// calculate the rank of an array directly, uses special threshold
		int rankdirect(const Eigen::MatrixXd &array) const;

        /// calculate the rank of the second order interaction matrix of an array directly
		int rankxfdirect(const array_link &array) const;

        /// calculate the rank of the second order interaction matrix of an array using the cache system
        int rankxf (const array_link &array);
};

/// Return the condition number of a matrix
double conditionNumber (const array_link &matrix);

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
/** Calculate properties of an array and create a Pareto element
 *
 * The values calculated are:
 *
 * 1) Rank (higher is better)
 * 2) A3, A4 (lower is better)
 * 3) F4 (lower is better, sum of elements is constant)
 *
 * Valid for 2-level arrays of strength at least 3
 *
 * */
typename Pareto< mvalue_t< long >, IndexType >::pValue calculateArrayParetoRankFA (const array_link &array,
                                                                                          int verbose) {
        int N = array.n_rows;
        int model_rank = arrayrankFullPivLU (array2secondorder (array), 1e-12) + 1 +
                         array.n_columns; // valid for 2-level arrays of strength at least 3
        mvalue_t< long > a3a4_values = A3A4 (array);
        mvalue_t< long > f4 = F4 (array);

        // add the 3 values to the combined value
        typename Pareto< mvalue_t< long >, IndexType >::pValue p;
        p.push_back (model_rank);
        p.push_back (a3a4_values);
        p.push_back (f4);

        if (verbose >= 2) {
                if (verbose >= 3) {
                        std::vector< double > gwlp = array.GWLP ();
                        myprintf ("parseArrayPareto: A4 (scaled) %ld, %f\n", a3a4_values.raw_values()[1], gwlp[4]);
                }

                myprintf ("  parseArrayPareto: rank %d, verbose %d\n", array.rank (), verbose);
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
/// Calculate Pareto element with J5 criterium
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
 * 3) F4 (lower is better, sum of elements is constant)
 *
 * */
inline void parseArrayPareto (const array_link &array, IndexType i, Pareto< mvalue_t< long >, IndexType > &pset,
                              int verbose) {
        typename Pareto< mvalue_t< long >, IndexType >::pValue p =
            calculateArrayParetoRankFA< IndexType > (array, verbose);

        // add the new tuple to the Pareto set
        pset.addvalue (p, i);
}

#endif // FULLPACKAGE

/// convert C value to D-efficiency value
inline double Cvalue2Dvalue (double Cvalue, int number_of_columns) {
        double ma = 1 + number_of_columns + number_of_columns * (number_of_columns - 1) / 2;
        double Defficiency = pow (Cvalue, 1. / ma);

        return Defficiency;
}

/// convert D-efficiency value to C value
inline double Dvalue2Cvalue (double Defficiency, int number_of_columns) {
        int ma = 1 + number_of_columns + number_of_columns * (number_of_columns - 1) / 2;
        double Cvalue = pow (Defficiency, ma);

        return Cvalue;
}



