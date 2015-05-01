/*! \file Deff.h
 *  \brief Contains functions to generate optimal designs
 * 
 */

#ifndef DEFF_H
#define DEFF_H

#include "arraytools.h"
#include "arrayproperties.h"

/// calculate score from from set of efficiencies
double scoreD(const std::vector<double> dd, const std::vector<double> alpha);

/// different algorithms for the optimization routines
enum {DOPTIM_SWAP, DOPTIM_UPDATE, DOPTIM_AUTOMATIC, DOPTIM_FLIP, DOPTIM_NONE};

/** Optimize a design according to the optimization function specified.
 * 
 * Arguments:
 * 	arrayclass: structure describing the design class
 * 	alpha: (3x1 array)
 * 	verbose: output level
 */
array_link  optimDeff(const array_link &A0,  const arraydata_t &arrayclass, std::vector<double> alpha, int verbose=1, int optimmethod = DOPTIM_AUTOMATIC, int niter=10000, int nabort=2500);


//typedef std::pair< std::vector<std::vector<double> >, arraylist_t > DoptimReturn;

/** @brief Structure containing results of the Doptimize function
 * 
 */
struct DoptimReturn {
  std::vector<std::vector<double> > dds;	/// scores generated
  arraylist_t designs;	/// designs generated
  int nrestarts;	/// final number of restarts performed
};


/// function to generate optimal designs 
DoptimReturn Doptimize(const arraydata_t &arrayclass, int nrestarts, int niter, std::vector<double> alpha, int verbose, int method, double maxtime = 100000, int nabort=5000 );


#endif 