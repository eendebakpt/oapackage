/*! \file evenodd.cpp
 *  \brief Contains functions to calculate even-odd designs
 *
 */

#include <printfheader.h>

#include "arraytools.h"
#include "arrayproperties.h"

#ifdef DOOPENMP
#include "omp.h"
#endif


#include "evenodd.h"

#ifndef myprintf
#define myprintf printf
#endif

#ifdef MAINMEX
#define MATLABUPDATE
#else
#ifdef MATLAB_MEX
#define MATLABUPDATE   mexEvalString ("drawnow" );
#else
#define MATLABUPDATE
#endif
#endif



// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
