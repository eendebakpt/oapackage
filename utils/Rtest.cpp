/** \file Rtest.cpp

 C++ program: Rtest

 Dummy file for testing R code

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "arraytools.h"
#include "arrayproperties.h"
#include "Deff.h"
#include "tools.h"

#include "R_ext/Print.h"


int main ( int argc, char* argv[] )
{

	int rr=40;
	int cc=7;
	int verbose=1;

	arraydata_t arrayclass ( 2, rr, 0, cc );
	int niter=10000;
	std::vector<double> alpha ( 3 );
	alpha[0]=1;
	alpha[1]=1;
	double t0=get_time_ms();
	array_link al = arrayclass.randomarray ( 0 );
	std::vector<double> dd = Defficiencies ( al, arrayclass, 0, 0 );
	int method=1;

	Rprintf("Rprintf test...\n"); 
	
	DoptimReturn drr = Doptimize(arrayclass, 20, alpha, verbose, method);
	return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
