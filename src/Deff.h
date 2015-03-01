/*! \file Deff.h
 *  \brief Contains functions to optimize designs
 * 
 */

#pragma once

#include "arraytools.h"
#include "arrayproperties.h"

//double tmp[] = { 1,1,0};
//std::vector<int> v( tmp, tmp+3 );
//std::vector<double> alpha0(tmp, tmp+3);

double scoreD(const std::vector<double> dd, const std::vector<double> alpha);

array_link  optimDeff(const array_link &A0,  arraydata_t &arrayclass, std::vector<double> alpha, int verbose=1, int method=-1, int niter=10000, int nabort=2500);

