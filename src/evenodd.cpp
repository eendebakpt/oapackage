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

/// initialize the new list of extension columns
arraylist_t depth_extend_sub_t::initialize ( const arraylist_t& alist, const arraydata_t &adf, const OAextend &oaextend )
{
	myassert ( alist.size() == lmctype.size() , "depth_extend_t: update" );

	int ncolsx=0;

	valididx.clear();

	if ( alist.size() >0 ) {
		ncolsx=alist[0].n_columns;
	} else {
		arraylist_t v;
		return v;
	}
	arraydata_t ad ( &adf, ncolsx );

	if ( verbose>=2 )
		printfd ( "initialize: %d \n", alist.size() );

	#pragma omp parallel for schedule(dynamic,1)	// FIXME: implement this
	for ( int k=0; k<(int)alist.size(); ++k ) {
		LMCreduction_t reduction ( &ad );
		// needed to make the code thread-safe
		reduction.initStatic();

		reduction.reset();
		reduction.setArray ( alist[k] );
		lmc_t rx = LMC_EQUAL;
		reduction.updateSDpointer ( alist[k] );
		/*
		if ( 0 ) { // FIXME: enable this?
			reduction.updateSDpointer ( alist[k] );
			int col=alist[k].n_columns-1;
			rx = LMC_check_col_rowsymm ( alist[k].array+alist[k].n_rows* ( alist[k].n_columns-1 ), &ad,  *reduction.sd.get(),col, 0 );
			printf ( "\n%s: rx %d\n", __FUNCTION__, rx );

			printf ( "col: " );
			print_perm ( alist[k].array+alist[k].n_rows* ( alist[k].n_columns-1 ), ad.N );
		} */

		lmc_t lmc =  LMCcheck ( alist[k].array, ad, oaextend, reduction );
		int lc=reduction.lastcol;
//#pragma omp critical
		{
			this->lmctype[k]=lmc;
			this->lastcol[k]=lc;
		}
		if ( verbose>=2 ) {
			printf ( "   depth_extend_sub_t.initialize: initialize: array %d lmc %d\n", k, ( int ) lmc );
			//  printf("   -> lastcol %d md5: %s\n", this->lastcol[k], alist[k].md5().c_str() );
			//oaextend.info(); reduction->show();

		}

		if ( k>25 && 0 ) {
			printf ( "file %s: line %d: exit\n", __FILE__, __LINE__ );
			exit ( 0 );
		}
		reduction.releaseStatic();
	}


	for ( size_t k=0; k<alist.size(); ++k ) {
		if ( verbose>=1 ) {
			printf ( "  depth_extend_sub_t.initialize: array %ld: lmc %d, lastcol %d, ncolumns %d\n", k, ( int ) this->lmctype[k], this->lastcol[k], ncolsx );
		}
		bool b1= ( this->lastcol[k]>=ncolsx || this->lastcol[k]==-1 );
		bool b2 = lmctype[k]>=LMC_EQUAL;
		if ( b1!=b2 ) {
			printf ( "oadevelop:initialize: huh? b1 %d b2 %d, lmctype[k] %d, lastcol[k] %d, col %d\n", b1, b2, lmctype[k], lastcol[k], ncolsx );
			oaextend.info();
			//writearrayfile("/home/eendebakpt/tmp/tmp.oa", alist[k]);
			//exit(0);

		}

		if ( this->lastcol[k]>=ncolsx-1 || this->lastcol[k]==-1 ) {
			valididx.push_back ( k );
		} else {
			if ( lmctype[k]>=LMC_EQUAL ) {
				//printf("oadevelop.h: initialize: huh?\n");
				//exit(0);
			}
		}
	}

	if ( verbose )
		printf ( "depth_extend_sub_t: initialize: selected %ld valid extension indices of %ld arrays\n", valididx.size(), alist.size() );
	arraylist_t v = ::selectArrays ( alist, valididx );

	// reset valididx
	for ( size_t i=0; i<valididx.size(); i++ )
		valididx[i]=i;

 
	if ( verbose>=3 )
		printf ( "   v %ld\n", v.size() );
	return v;
}


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
