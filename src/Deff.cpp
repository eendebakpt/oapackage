/*! \file Deff.cpp
 *  \brief Contains functions to optimize designs
 *
 */

#include <stdio.h>

#include "arraytools.h"
#include "arrayproperties.h"

#ifdef DOOPENMP
#include "omp.h" 
#endif


#include "Deff.h"


#ifndef myprintf
#define myprintf printf
#endif


DoptimReturn Doptimize(const arraydata_t &arrayclass, int nrestartsmax, int niter, std::vector<double> alpha, int verbose, int method, double maxtime , int nabort )
{
		double t0 = get_time_ms();
	std::vector<std::vector<double> > dds;
	arraylist_t AA;


	
	bool abort=false;
	
	int nrestarts=0;
	
	#ifdef DOOPENMP
	#pragma omp parallel for num_threads(4) schedule(dynamic,1)
#endif
	for ( int i=0; i<nrestartsmax; i++ ) {
		if ( abort )
			continue;

#ifdef DOOPENMP
		#pragma omp critical
#endif
		{
			if ( verbose && i%25==0 ) {
				myprintf ( "Doptim: iteration %d/%d\n", i, nrestartsmax );
#ifdef MAINMEX
#else
#ifdef MATLAB_MEX
				mexEvalString ( "drawnow" );
#endif
#endif
			}
		}

		array_link al = arrayclass.randomarray ( 1 );


		array_link  A = optimDeff ( al,  arrayclass, alpha, verbose>=2, method, niter,  nabort );
		std::vector<double> dd = A.Defficiencies();
		if ( verbose>=2 ) {
			myprintf ( "Doptim: iteration %d/%d: %f %f %f\n", i, nrestartsmax, dd[0], dd[1], dd[2] );
		}

#ifdef DOOPENMP
		#pragma omp critical
#endif
		{
			AA.push_back ( A );
			dds.push_back ( dd );
			nrestarts++;
		}

		if ( ( get_time_ms()-t0 ) > maxtime ) {
			if ( verbose )
				myprintf ( "max running time exceeded, aborting\n" );
			#pragma omp critical
			{
				abort=true;
			}
			//break;
		}
	}

	// loop is complete

	DoptimReturn a = {dds, AA, nrestarts};
	//a.first;
	
	return a; // DoptimReturn( dds, AA);

}

double scoreD ( const std::vector<double> dd, const std::vector<double> alpha )
{
	double v=0;
	for ( size_t i=0; i<dd.size();  i++ )
		v+= dd[i]*alpha[i];
	return v;

}

array_link  optimDeff ( const array_link &A0,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niter, int nabort )
{

	int nx=0;
	int N = arrayclass.N;
	int k = arrayclass.ncols;
	//#get factor levels
	std::vector<int> s = arrayclass.getS();

	if ( optimmethod==DOPTIM_UPDATE ) {
		if ( arrayclass.is2level() )
			optimmethod=DOPTIM_FLIP;
	}
	array_link A = A0;
	symmetry_group sg=symmetry_group ( s );

	std::vector<int> gidx = sg.gidx;

	std::vector<double> dd0 = A0.Defficiencies();

	if ( verbose ) {
		printf ( "optimDeff: initial D-efficiency %.4f\n",  dd0[0] );
	}

	// initialize score
	double d = scoreD ( dd0, alpha );
//   gsize=tuple(sg.gsize)
//   gstart=tuple(sg.gstart)

	int lc=0;	// index of last change to array


//#pragma omp for
	for ( int ii=0; ii<niter; ii++ ) {

		// select random row and column
		int r = fastrandK ( N );		
		int c = fastrandK ( k );
		int r2 = fastrandK ( N );
		
		r=rand()%N; c=rand()%k; r2=rand()%N;
		
		//int c2 = fastrandK(k);
		// make sure column is proper column group
		int c2 = sg.gstart[sg.gidx[c]] + fastrandK ( sg.gsize[gidx[c]] );

		// get values
		array_t o = A._at ( r,c );
		array_t o2 = A._at ( r2,c2 ); // no extra error checking

		// update
		
					if (optimmethod==DOPTIM_SWAP && o==o2)
				continue;

		switch ( optimmethod ) {
		case DOPTIM_SWAP: // swap
			
			
			A._setvalue ( r,c,o2 );
			A._setvalue ( r2,c2,o );
			break;
		case DOPTIM_UPDATE: { // random update
			int val = fastrandK ( s[c] );
			A._setvalue ( r,c,val );
			break;
		}
		case DOPTIM_FLIP: // flip
			A._setvalue ( r,c,1-o );
			break;
			case DOPTIM_NONE:
			break;
		}
		// evaluate
		std::vector<double> dd = A.Defficiencies();
		nx++;
		double dn = scoreD ( dd, alpha );


		// switch back if necessary

		if ( dn>=d ) {
			if ( dn>d )  {
				lc=ii;
				if ( verbose>=3 )
					printf ( "optimDeff: ii %d: %.6f -> %.6f\n", ii, d, dn );
				d=dn;
			}
		} else {
			// restore to original
			switch ( optimmethod ) {
			case DOPTIM_SWAP:
				A._setvalue ( r,c,o );
				A._setvalue ( r2,c2,o2 );
				break;
			case DOPTIM_FLIP:
			case DOPTIM_UPDATE:
				A._setvalue ( r,c,o );
				break;
			case DOPTIM_NONE:
			break;
			}
		}

		// check progress
		if ( ( ii-lc ) >nabort ) {
			if ( verbose>=2 )
				printf ( "optimDeff: early abort ii %d, lc %d\n", ii, lc );
			//abort=true;
			break;
		}

	}

	std::vector<double> dd = A.Defficiencies();
	double dn = scoreD ( dd, alpha );

		if ( verbose ) {
		printf ( "optimDeff: final score %.4f, final D-efficiency %.4f\n",  dn, dd0[0] );
	}

	//printf("nx %d\n", nx);
	return A;
//      return std::pair<array_link, std::vector<double> >(A, dd);
}

#include <algorithm> 

extern "C" {
	
double DoptimizeR(int *_N, int *_k, int *nrestarts, int *_niter, double *alpha1, double *alpha2, double *alpha3, int *_verbose, int *_method, double *maxtime , int *nabort, double *output )
{

int niter=*_niter;
int method=*_method;
int N = *_N;
int k = *_k;
int verbose = *_verbose;

output[0]=1;
output[1]=2;
output[2]=3;
output[3]=4;

if (verbose>=2)
	printf("DoptimizeR: N %d, k %d, nrestarts %d, niter %d, alpha1 %f\n", N, k, *nrestarts, niter, *alpha1);

 arraydata_t arrayclass(2, N, 0, k);
std::vector<double> alpha(3);
alpha[0]=std::max(*alpha1,0.);
alpha[1]=std::max(*alpha2, 0.);
alpha[2]=std::max(*alpha3,0.);

DoptimReturn rr = Doptimize(arrayclass, *nrestarts, niter, alpha,  verbose,  method, *maxtime,  *nabort);

	std::vector<std::vector<double> > dds = rr.dds; 	arraylist_t AA = rr.designs;

	// sort according to values
	std::vector<double> sval ( AA.size() );
	for ( size_t i=0; i<AA.size(); i++ ) {
		sval[i]=-scoreD ( dds[i], alpha );
	}

	indexsort sorter ( sval );

	AA=sorter.sorted ( AA );
	dds=sorter.sorted ( dds );

array_link best = AA[0];

std::copy(best.array, best.array+N*k, output);

if (verbose>=2) {
	printf("DoptimizeR: done\n");
}
//best.showarray();

return best.Defficiency();

}

} // extern "C"

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
