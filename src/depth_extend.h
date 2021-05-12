/** \file depth_extend.h

 This file contains development code.

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#ifndef OADEVELOP_H
#define OADEVELOP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <deque>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <numeric>

#include "lmc.h"
#include "arrayproperties.h"
#include "extend.h"
#include "nonroot.h"

#include "pareto.h"
#include "evenodd.h"

#ifndef NOOMP
// set this to disable openmp
#ifdef OADEBUG
#define NOOMP
#endif
#endif

#ifdef _OPENMP
#else
#ifndef NOOMP
#define NOOMP
#endif
#endif

#ifdef NOOMP
#else
#ifndef DOOPENMP
#define DOOPENMP
#endif
#endif


/// extend array using dynamic pruning based on D-efficiency
int extend_array_dynamic ( const arraylist_t &alist, arraydata_t &adx, OAextend &oaextend, arraylist_t &earrays, std::vector<std::vector<double> > &edata, int kfinal, double Afinal, int directcheck, int verbose=1 );
/// extend array using dynamic pruning based on D-efficiency
int extend_array_dynamic ( const array_link &al, arraydata_t &adx, OAextend &oaextend, arraylist_t &earrays, std::vector<std::vector<double> > &edata, dextend_t &dextend, int kfinal, double Afinal, int rs, int verbose=1 );


#endif
