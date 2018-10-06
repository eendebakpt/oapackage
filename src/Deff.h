/*! \file Deff.h
 *  \brief Contains functions to generate optimal designs
 *
 * For more information see "Two-Level Designs to Estimate All Main Effects and Two-Factor Interactions", P.T. Eendebak
 * and E.D. Schoen, 2017
 */

#ifndef DEFF_H
#define DEFF_H

#include "arrayproperties.h"
#include "arraytools.h"

/** calculate score from from set of efficiencies
 *
 * The score is the weighted sum of the efficiencies.
 *
 * \param efficiencies Vector with calculated efficiencies
 * \param alpha Weights for the efficiencies
 */
double scoreD (const std::vector< double > efficiencies, const std::vector< double > alpha);

/// different algorithms for the optimization routines
enum { DOPTIM_UPDATE, DOPTIM_SWAP, DOPTIM_FLIP, DOPTIM_AUTOMATIC, DOPTIM_NONE };

/** Optimize a design according to the optimization function specified.
 *
 * Arguments:
 * \param A0			Array to be optimized
 * \param arrayclass	Structure describing the design class
 * \param alpha			3x1 array with optimization parameters
 * \param verbose		Verbosity level
 * \param optimmethod	Optimization method to use
 * \param niter			Number of iterations
 * \param nabort		Number of iterations after which to abort when no improvements are found
 * \returns		Optimized designs
 */
array_link optimDeff (const array_link &A0, const arraydata_t &arrayclass, std::vector< double > alpha,
                      int verbose = 1, int optimmethod = DOPTIM_AUTOMATIC, int niter = 100000, int nabort = 0);

/** Structure containing results of the Doptimize function
 */
struct DoptimReturn {
        /// calculated efficiencies for the designs
        std::vector< std::vector< double > > dds;
        /// designs generated
        arraylist_t designs;
        /// number of restarts performed
        int nrestarts;
        int nimproved;
};

/** Generates optimal designs for the specified class of designs
 *
 * \param arrayclass Class of designs to optimize
 * \param nrestarts Number of restarts to perform
 * \param alpha Optimization parameters
 * \param verbose Verbosity level
 * \param method Method for optimization algorithm
 * \param niter Maximum number of iterations for each restart
 * \param maxtime Maximum calculation time. If this time is exceeded, the function is aborted
 * \param nabort Maximum number of iterations when no improvement is found
 *
 */
DoptimReturn Doptimize (const arraydata_t &arrayclass, int nrestarts, const std::vector< double > alpha, int verbose,
                        int method = DOPTIM_AUTOMATIC, int niter = 300000, double maxtime = 100000, int nabort = 5000);

/// function to generate optimal designs with mixed optimization approach
DoptimReturn DoptimizeMixed (const arraylist_t &sols, const arraydata_t &arrayclass, const std::vector< double > alpha,
                             int verbose = 1, int nabort = -1);

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
