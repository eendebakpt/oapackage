/** \file nonroot.h

\brief This file contains definitions are functions related to the non-root stage of the LMC algorithm

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See COPYING file that comes with this distribution
*/

#pragma once

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>

#include "arrayproperties.h"
#include "mathtools.h"
#include "tools.h"

#include "lmc.h"

/// default reduction function for non-root stage
lmc_t LMCreduce_non_root (const array_t *original, const arraydata_t *ad, dyndata_t *dyndata,
                          LMCreduction_t *reduction, const OAextend &oaextend, const LMC_static_struct_t &tmpStatic);

/// default reduction function for non-root stage with J4 algorithm
lmc_t LMCreduce_non_root_j4 (const array_t *original, const arraydata_t *ad, const dyndata_t *dyndata,
                             LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic);

/// specialized reduction function for 2-level arrays
lmc_t LMCreduce_non_root_2level (const array_t *original, const arraydata_t *ad, dyndata_t *dyndata,
                                 LMCreduction_t *reduction, const OAextend &oaextend,
                                 const LMC_static_struct_t &tmpStatic);

/// check column for row symmetry exchanges
lmc_t LMC_check_col_rowsymm (const array_t *arraycol, const arraydata_t *ad, const symmdata &sd, int col,
                             int dverbose = 0);
