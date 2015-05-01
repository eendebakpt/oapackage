/** \file nonroot.h

\brief This file contains definitions are functions related to the non-root stage of the algorithm

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See COPYING file that comes with this distribution
*/

#pragma once

#include <list>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <set>

#include "mathtools.h"
#include "tools.h"
#include "arrayproperties.h"

#include "lmc.h"


/// default reduction function for non-root stage
lmc_t LMCreduce_non_root ( const array_t * original, const arraydata_t* ad, dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, const LMC_static_struct_t &tmpStatic ) ;
lmc_t LMCreduce_non_root_j4 ( const array_t * original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic );

lmc_t LMCreduce_non_root_2level ( const array_t * original, const arraydata_t* ad, dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, const LMC_static_struct_t &tmpStatic ) ;

lmc_t LMC_check_col_rowsymm ( const array_t *arraycol, const arraydata_t *ad, const symmdata &sd, int col, int dverbose=0 );
