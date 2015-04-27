/*! \file Deff.h
 *  \brief Contains functions to optimize designs
 * 
 */

#ifndef DEFF_H
#define DEFF_H


#include "arraytools.h"
#include "arrayproperties.h"

//double tmp[] = { 1,1,0};
//std::vector<int> v( tmp, tmp+3 );
//std::vector<double> alpha0(tmp, tmp+3);

/// calculate score from from set of efficiencies
double scoreD(const std::vector<double> dd, const std::vector<double> alpha);

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

struct DoptimReturn {
  std::vector<std::vector<double> > dds;	/// scores generated
  arraylist_t designs;	/// designs generated
  int nrestarts;	/// final number of restarts performed
};


/// optimize designs 
DoptimReturn Doptimize(const arraydata_t &arrayclass, int nrestarts, int niter, std::vector<double> alpha, int verbose, int method, double maxtime = 100000, int nabort=5000 );

extern "C" {

// wrapper function for R
//double DoptimizeR(int *N, int *k, int *nrestarts, int *niter, double *alpha1, double *alpha2, double *alpha3, int *verbose, int *method, double *maxtime , int *nabort );

}
#endif 