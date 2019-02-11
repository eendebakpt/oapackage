/*! \file Deff.h
 *  \brief Contains functions to generate optimal designs
 *
 * For more information see "Two-Level Designs to Estimate All Main Effects and Two-Factor Interactions", P.T. Eendebak
 * and E.D. Schoen, 2017
 */

#pragma once
#include <vector>

#include "arrayproperties.h"
#include "arraytools.h"

/** Calculate score from a set of efficiencies
 *
 * The score is the weighted sum of the efficiencies.
 *
 * \param efficiencies Vector with calculated efficiencies
 * \param weights Weights for the efficiencies
 * \returns Weighted sum of the efficiencies
 */
double scoreD (const std::vector< double > efficiencies, const std::vector< double > weights);

/// Different methods for the optimization. The default method DOPTIM_SWAP is a coordinate-exchange algorithms
enum coordinate_exchange_method_t {
	/// replace a random element with a random value
	DOPTIM_UPDATE,
	/// swap two elements at random
	DOPTIM_SWAP,
	/// randomly flip an element between 0 and 1
	DOPTIM_FLIP,
	/// automatically select one of the methods
	DOPTIM_AUTOMATIC, 
	/// perform no optimization
	DOPTIM_NONE } ;

/** Structure containing results of the Doptimize function
 */
struct DoptimReturn {
        /// calculated efficiencies for the generated designs
        std::vector< std::vector< double > > dds;
        /// designs generated
        arraylist_t designs;
        /// number of restarts performed
        int nrestarts;
        int _nimproved;
};

/** Generates optimal designs for the specified class of designs
 *
 * The method uses a coordinate-exchange algorithm to optimze a target function defined by the optimziation paramaters.
 * The optimization is performed multiple times to prevent finding a design in a local minmum of the target function.
 *
 * The method is described in more detail in "Two-Level Designs to Estimate All Main Effects and Two-Factor Interactions",
 * Eendebak et al., 2015, Technometrics, https://doi.org/10.1080/00401706.2016.1142903.
 *
 * \param arrayclass Class of designs to optimize
 * \param nrestarts Number of restarts to perform
 * \param alpha Optimization parameters. The target function is alpha_1 D + alpha_2 D_s + alpha D_1
 * \param verbose Verbosity level
 * \param method Method for optimization algorithm
 * \param niter Maximum number of iterations for each restart
 * \param maxtime Maximum calculation time. If this time is exceeded, the function is aborted
 * \param nabort Maximum number of iterations when no improvement is found
 * \returns A structure with the generated optimal designs
 */
DoptimReturn Doptimize (const arraydata_t &arrayclass, int nrestarts, const std::vector< double > alpha, int verbose,
                        coordinate_exchange_method_t method = DOPTIM_AUTOMATIC, int niter = 300000, double maxtime = 100000, int nabort = 5000);

/** Function to generate optimal designs with mixed optimization approach
 *
 * This function is beta code. See Doptimize for detauls of the  parameters.
 *
 */
DoptimReturn DoptimizeMixed (const arraylist_t &sols, const arraydata_t &arrayclass, const std::vector< double > alpha,
                             int verbose = 1, int nabort = -1);

/** Optimize a design according to the optimization function specified.
 *
 * Arguments:
 * \param array			Array to be optimized
 * \param arrayclass	Structure describing the design class
 * \param alpha			3x1 array with optimization parameters
 * \param verbose		Verbosity level
 * \param optimmethod	Optimization method to use
 * \param niter			Number of iterations
 * \param nabort		Number of iterations after which to abort when no improvements are found
 * \returns		Optimized designs
 */
array_link optimDeff (const array_link &array, const arraydata_t &arrayclass, std::vector< double > alpha,
                      int verbose = 1, coordinate_exchange_method_t optimmethod = DOPTIM_AUTOMATIC, int niter = 100000, int nabort = 0);
