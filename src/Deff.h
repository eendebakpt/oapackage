/*! \file Deff.h
 *  \brief Contains functions to generate optimal designs
 *
 */

#ifndef DEFF_H
#define DEFF_H

#include "arraytools.h"
#include "arrayproperties.h"

/// calculate score from from set of efficiencies
double scoreD ( const std::vector<double> dd, const std::vector<double> alpha );

/// different algorithms for the optimization routines
enum {DOPTIM_UPDATE, DOPTIM_SWAP, DOPTIM_FLIP, DOPTIM_AUTOMATIC, DOPTIM_NONE};

/** Optimize a design according to the optimization function specified.
 *
 * Arguments:
 * 	arrayclass: structure describing the design class
 * 	alpha: (3x1 array)
 * 	verbose: output level
 */
array_link  optimDeff ( const array_link &A0,  const arraydata_t &arrayclass, std::vector<double> alpha, int verbose=1, int optimmethod = DOPTIM_AUTOMATIC, int niter=100000, int nabort=0 );

/// debugging function
array_link  optimDeff2level ( const array_link &A0,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose=1, int optimmethod= DOPTIM_AUTOMATIC, int niter=100000, int nabort = 0 );


/** @brief Structure containing results of the Doptimize function
 *
 */
struct DoptimReturn {
    std::vector<std::vector<double> > dds;	/// scores generated
    arraylist_t designs;	/// designs generated
    int nrestarts;	/// final number of restarts performed
    int nimproved;
};


/// function to generate optimal designs
DoptimReturn Doptimize ( const arraydata_t &arrayclass, int nrestarts, const std::vector<double> alpha, int verbose, int method = DOPTIM_AUTOMATIC, int niter = 300000, double maxtime = 100000, int nabort=5000 );

/// helper function
DoptimReturn DoptimizeDebug ( const arraydata_t &arrayclass, int nrestarts, const std::vector<double> alpha, int verbose, int method = DOPTIM_AUTOMATIC, int niter = 300000, double maxtime = 100000, int nabort=5000 );

DoptimReturn DoptimizeMixed(const arraylist_t &sols, const arraydata_t &arrayclass, const std::vector<double> alpha, int verbose=1, int nabort=-1);



#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
