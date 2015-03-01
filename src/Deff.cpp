/*! \file Deff.cpp
 *  \brief Contains functions to optimize designs
 *
 */

#include "arraytools.h"
#include "arrayproperties.h"


#include "Deff.h"



double scoreD ( const std::vector<double> dd, const std::vector<double> alpha )
{
	double v=0;
	for ( size_t i=0; i<dd.size();  i++ )
		v+= dd[i]*alpha[i];
	return v;

}

array_link  optimDeff ( const array_link &A0,  arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niter, int nabort )
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

	//printf("nx %d\n", nx);
	return A;
//      return std::pair<array_link, std::vector<double> >(A, dd);
}


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
