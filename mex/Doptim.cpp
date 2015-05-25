/*! D-optimization
 *
 * \author Pieter Eendebak <pieter.eendebak@gmail.com>
 */


#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <sys/timeb.h>


#include "Deff.h"
#include "arraytools.h"
#include "arrayproperties.h"

#include "mex.h"

#ifdef DOOPENMP
#include "omp.h"
#endif


/// Get optional int scalar from Matlab parameters
inline int mxGetOptionalInt ( int n_in, const mxArray **prhs, int pos, int defaultvalue )
{
	int value = defaultvalue;
	if ( n_in > pos )
		value = static_cast < int > ( mxGetScalar ( prhs[pos] ) );
	return value;
}
/// Get optional double scalar from Matlab parameters
inline double mxGetOptionalDouble ( int n_in, const mxArray **prhs, int pos, double defaultvalue )
{
	double value = defaultvalue;
	if ( n_in > pos )
		value = mxGetScalar ( prhs[pos] );
	return value;
}


void mexFunction ( int n_out, mxArray *plhs[], int n_in, const mxArray *prhs[] )
{

#ifdef DOOPENMP
#include "omp.h"
#endif

	/* Parse input */
	if ( n_in < 3 ) {
		mexPrintf ( "Usage: Doptim(N, levels, [nrestarts], [alpha], [verbose], [maxtime])\n" );
		mexPrintf ( "   calculates D-optimal designs, " );
		mexPrintf ( "written by Pieter Eendebak <pieter.eendebak@gmail.com>\n\n" );
		mexPrintf ( "The function optimizes designs according to the optimization function:\n" );
		mexPrintf ( "\n	F= alpha(1)*D + alpha(2) * Ds + alpha(3) * D1\n\n" );
		mexPrintf ( "Input arguments:\n" );
		mexPrintf ( "   N (integer): number of runs\n" );
		mexPrintf ( "   s (array): levels of the factors\n" );
		mexPrintf ( "   nrestarts (integer, default is 40): number of restarts\n" );
		mexPrintf ( "   alpha (3x1 array, default is [1,1,0]): weights for the various efficiency measures\n" );
		mexPrintf ( "   verbose (integer, default is 1): Determines the output level. Higher is more output\n" );
		mexPrintf ( "   maxtime (double, default: 150): maximum running time in seconds\n" );
		mexPrintf ( "Output arguments:\n" );
		mexPrintf ( "   A (Nxk or Nxkxm array): the generated design(s). m is equal to the number of restarts specified.\n" );
		mexPrintf ( "   d (3x1 or 3xm array): computed values of D-efficiency, Ds and D1\n" );
		mexPrintf ( "\n" );
		mexPrintf ( "For more details see the webpage http://www.pietereendebak.nl/oapackage.\n" );
		mexPrintf ( "\n" );

		mexPrintf ( "Example: [A, d] = Doptim(40, [2,2,2,2,2,2], 10, [1,2,0])\n\n" );
		//mexPrintf ( "Example: [A, d] = Doptim(18, [3,3,2,2], 10, [1,2,0])\n\n" );

		mexPrintf ( "   built: %s %s\n", __DATE__, __TIME__ );
		mexPrintf ( "error: Doptim: requires at least 3 input arguments!\n" );
		return;
	}


	int N = mxGetOptionalInt ( n_in, prhs, 0, 40 );

	int nrestarts = mxGetOptionalInt ( n_in, prhs, 2, 40 );
	int verbose = mxGetOptionalInt ( n_in, prhs, 4, 1 );
	double maxtime = mxGetOptionalDouble ( n_in, prhs, 5, 150 );
	int dverbose = mxGetOptionalInt ( n_in, prhs, 6, 1 );
	int nout=1;

	// error checking
	int k = mxGetM ( prhs[1] ) *mxGetN ( prhs[1] );

	if ( N<1 || N>2000 )
		mexErrMsgTxt ( "Doptim: parameter N should be between 1 and 2000!" );

	if ( nrestarts<1 ||nrestarts>1000000 )
		mexErrMsgTxt ( "Doptim: parameter nrestarts should be between 1 and 200.000!" );

	if ( k<1 || k>32 )
		mexErrMsgTxt ( "Doptim: number of factors (size of levels array) should be between 1 and 32!" );

	std::vector<int> s ( k );
	double *x = mxGetPr ( prhs[1] );
	for ( int i=0; i<k; i++ ) {
		s[i] = x[i];
		if ( s[i]<2 || s[i]>20 )
			mexErrMsgTxt ( "Doptim: factor levels should be between 2 and 20!" );
	}

	std::vector<double> alpha ( 3 );
	if ( n_in<4 ) {
		alpha[0]=1;
		alpha[1]=1;
		alpha[2]=0;
	} else {
		int na = mxGetM ( prhs[3] ) *mxGetN ( prhs[3] );
		if ( na!=3 ) {
			mexErrMsgTxt ( "Doptim: parameter alpha should be 3x1 double array!" );

		}
		x = mxGetPr ( prhs[3] );
		for ( int i=0; i<3; i++ )
			alpha[i] = x[i];
	}

	if ( verbose ) {
		mexPrintf ( "Doptim: N %d, k %d, nrestarts %d\n", N, k, nrestarts );
	}

	//const mwSize *dims = mxGetDimensions(prhs[0]);
	arraydata_t arrayclass ( s, N, 0, k );

	int niter=12000;
	int nabort=3000;
	int method=1; 

	double t0 = get_time_ms();
	bool abort = false;


	DoptimReturn rr = Doptimize(arrayclass, nrestarts, niter, alpha, verbose, method, maxtime, nabort );
	//std::vector<std::vector<double> > dds = rr.first; 	arraylist_t AA = rr.second;
	std::vector<std::vector<double> > dds = rr.dds; 	arraylist_t AA = rr.designs;
	
	
	if ( verbose>=2 ) {
		mexPrintf ( "generated %d designs, selecting best %d design(s)\n", AA.size(), nout );
	}

	nout = std::min ( nout, ( int ) AA.size() );
	// select best arrays


		if ( verbose>=1 && nout>0) {
		mexPrintf ( "Doptim: generated design with D-efficiency %f\n",  AA[0].Defficiency());
	}


	// sort according to values
	std::vector<double> sval ( AA.size() );
	for ( size_t i=0; i<AA.size(); i++ ) {
		sval[i]=-scoreD ( dds[i], alpha );
	}

	indexsort sorter ( sval );

	AA=sorter.sorted ( AA );
	dds=sorter.sorted ( dds );

	// reduce to requested output size
	AA.resize ( nout );
	dds.resize ( nout );

	int nn = AA.size();
	plhs[1] = mxCreateNumericMatrix ( 3, nn, mxDOUBLE_CLASS, mxREAL );
	mwSize dims[3];
	dims[0]=N;
	dims[1]=k;
	dims[2]=nn;
	plhs[0] = mxCreateNumericArray ( 3, dims, mxDOUBLE_CLASS, mxREAL );
	double *out = mxGetPr ( plhs[1] );

	for ( int i=0; i<nn; i++ ) {
		for ( int j=0; j<3; j++ )
			out[3*i+j]=dds[i][j];
	}

	out = mxGetPr ( plhs[0] );

	for ( int i=0; i<nn; i++ ) {
		for ( int j=0; j<N*k; j++ )
			out[N*k*i+j]=AA[i].array[j];
	}



//	mexPrintf("go!\n");
	return;
}

#include "mexstandalone.h"
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
