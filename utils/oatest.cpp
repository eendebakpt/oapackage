/** \file oatest.cpp

 C++ program: oatest

 oatest: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "oadevelop.h"
#include "lmc.h"


//void mydebug(array_link al, arraydata_t &adata, OAextend &oaextend, int);
//void mydebug2(array_link al, arraydata_t &adata, OAextend &oaextend);

void mydebug10 ( array_link al, arraydata_t &adata, OAextend &oaextend, int dverbose=0 )
{
	double t0,t,dt;
	printf ( "mydebug2 (J5)! al %d %d\n", al.n_rows, al.n_columns );
	oaextend.info();

	al.showarray();
	al.showproperties();
	al.selectFirstColumns ( 6 ).showproperties();

	lmc_t r;
	LMCreduction_t reduction ( &adata );
	t0 =get_time_ms();
	LMCreduction_t reductionx=reduction;
	r = LMCcheck ( al, adata, oaextend, reductionx ) ;
	dt =get_time_ms()-t0;
	reductionx.symms.makeColpermsUnique();
	reductionx.symms.showColperms();
	reductionx.symms.showSymmetries();
	printf ( "original lmc_t: %d, time: %.3f [ms]\n", ( int ) r, 1e3*dt );


	printf ( "--------------------------------------------\n" );

	// pre-compute
	LMCreduction_t reductionsub = calculateSymmetryGroups ( al.deleteColumn ( -1 ), adata,  oaextend );
	reductionsub.symms.colperms[3].clear();
	reductionsub.symms.colperms[4].clear();
	reductionsub.symms.showColperms ( 1 );
	reductionsub.symms.showColcombs ( 1 );

	printf ( "--------------------------------------------\n" );
	if ( 0 ) {
		arraydata_t adfix = adata;
		// split column groups
		std::vector<int> splits; // = symmetrygroup2splits(sg, ad.ncols, verbose);
		// for(int k=0; k<=i+1; k++) splits.push_back(k);
		splits.push_back ( 0 ); // splits.push_back(5);
		adfix.set_colgroups ( splits );

		std::vector<int> w;
		w.push_back ( 0 );
		w.push_back ( 1 );
		w.push_back ( 2 );
		w.push_back ( 4 );
		w.push_back ( 11 );
		std::vector<colindex_t> ww = comb2perm<colindex_t> ( w, adata.ncols );
		print_perm ( ww );
		array_link alx = al.selectColumns ( ww );
		reduction.setArray ( al );

		// FIXME: allow non-properinitialization of array
		//reduction.mode=LMC_REDUCE_INIT;

		r = LMCcheckj5 ( alx, adfix, reduction, oaextend, 1 );
		printf ( "xxx r %d\n", r );
		return;
	}

	lmc_t rn;
	t0 =get_time_ms();
	int niter=1;
	for ( int ix=0; ix<niter; ix++ ) {
		copy_array ( al.array, reduction.array, adata.N, adata.ncols ); // hack?
		rn = LMCcheckSymmetryMethod ( al, adata, oaextend, reduction, reductionsub, dverbose ) ;
	}
	dt =get_time_ms()-t0;
	printf ( "new lmt_t %d: time: %.3f [ms]\n", rn, 1e3*dt );

	if ( r!=rn ) {
		printf ( "ERROR: orig lmc_t %d, new lmc_t %d\n", r, rn );
	}

}

void mydebug5 ( array_link al, arraydata_t &adata, OAextend &oaextend, int dverbose=1 )
{
	double t0,dt;
	oaextend.setAlgorithm ( MODE_J5ORDERXFAST, &adata );

	lmc_t r=LMC_MORE;
	printf ( "mydebug5! al %d %d\n", al.n_rows, al.n_columns );
	LMCreduction_t reduction ( &adata );

	t0=get_time_ms();

	// pre-compute
	//oaextend.setAlgorithm(MODE_LMC_SYMMETRY);
	LMCreduction_t reductionsub = calculateSymmetryGroups ( al.deleteColumn ( -1 ), adata,  oaextend );

	reductionsub.symms.showSymmetries();
	reductionsub.symms.showColcombs();

	dt =get_time_ms()-t0;
	printf ( "### pre-compute: time: %.3f [ms]\n", 1e3*dt );

	int niter=250;
	lmc_t rn = r;
	t0 =get_time_ms();
	for ( int ix=0; ix<niter; ix++ ) {
		copy_array ( al.array, reduction.array, adata.N, adata.ncols ); // hack?
		reduction.reset();
		reduction.updateSDpointer ( al );
		rn = LMCcheckSymmetryMethod ( al, adata, oaextend, reduction, reductionsub, dverbose ) ;
	}
	dt =get_time_ms()-t0;
	printf ( "### new lmt_t %d: time: %.3f [ms]\n", rn, 1e3*dt/niter );


	t0 =get_time_ms();
	LMCreduction_t reductionx=reduction;
	niter=4;
	for ( int ix=0; ix<niter; ix++ ) {
		reductionx.init_state=COPY;
		r = LMCcheck ( al, adata, oaextend, reductionx ) ;
	}
	dt =get_time_ms()-t0;
	if ( niter>0 ) {
		printf ( "### original lmc_t: %d, time: %.3f [ms]\n\n", ( int ) r, 1e3*dt/niter );
	}
}


void mydebug ( array_link al, arraydata_t &adata, OAextend &oaextend, int dverbose=1 )
{
	double t0, dt, t;
	oaextend.setAlgorithm ( MODE_J5ORDERXFAST, &adata );

	lmc_t r=LMC_MORE;
	int niter=1;
	if ( dverbose<0 ) {
		niter=10;
		dverbose=0 ;
	}

	printf ( "mydebug! al %d %d\n", al.n_rows, al.n_columns );
	LMCreduction_t reduction ( &adata );
	LMCreduction_t reductionx = reduction;

	if ( 0 ) {
		t0 =get_time_ms();
		reductionx=reduction;
		//oaextend.setAlgorithm(MODE_J4);
		for ( int ix=0; ix<niter; ix++ ) {
			reductionx.init_state=COPY;
			r = LMCcheck ( al, adata, oaextend, reductionx ) ;
		}
		dt =get_time_ms()-t0;
		printf ( "### original lmc_t: %d, time: %.3f [ms]\n\n", ( int ) r, 1e3*dt );
	}

	if ( 0 ) {
		t0 =get_time_ms();
		reductionx=reduction;
		oaextend.setAlgorithm ( MODE_LMC_SYMMETRY );
		for ( int ix=0; ix<niter; ix++ ) {
			reductionx.init_state=COPY;
			r = LMCcheck ( al, adata, oaextend, reductionx ) ;
		}
		dt =get_time_ms()-t0;
		printf ( "### symmetry lmc_t: %d, time: %.3f [ms]\n\n", ( int ) r, 1e3*dt );

	}

	// FIXME: create safe combination class (using std::vector)

	t0=get_time_ms();

	// pre-compute
	//oaextend.setAlgorithm(MODE_LMC_SYMMETRY);

	LMCreduction_t reductionsub = calculateSymmetryGroups ( al.deleteColumn ( -1 ), adata,  oaextend );

	reductionsub.symms.showSymmetries();

	dt =get_time_ms()-t0;
	printf ( "### pre-compute: time: %.3f [ms]\n", 1e3*dt );

	if ( 0 ) {
		int nc=6;
		printf ( "## symmetrices with %d cols:\n", nc );
		symmetryset xx= reductionsub.symms.symmetries[nc];
		for ( symmetryset::const_iterator it = xx.begin(); it != xx.end(); it++ ) {
			it->show();
		}
	}


	//reductionsub.showColperms(1);
	if ( dverbose ) {
		reductionsub.symms.showColcombs ( 1 );
		reductionsub.symms.showSymmetries ( 1 );
	}


	printf ( "running LMCcheckSymmetryMethod: " );
	if ( dverbose ) {
		adata.show();
		reduction.symms.show();
		reductionsub.symms.show();
	}

	lmc_t rn = r;
	t0 =get_time_ms();
	for ( int ix=0; ix<niter; ix++ ) {
		copy_array ( al.array, reduction.array, adata.N, adata.ncols ); // hack?
		reduction.updateSDpointer ( al );

		rn = LMCcheckSymmetryMethod ( al, adata, oaextend, reduction, reductionsub, dverbose ) ;
	}
	dt =get_time_ms()-t0;
	printf ( "### new lmt_t %d: time: %.3f [ms]\n", rn, 1e3*dt );


	if ( r!=rn ) {
		printf ( "###\nERROR: orig lmc_t %d, new lmc_t %d\n", r, rn );
		reductionx.mode=OA_REDUCE;
		r = LMCcheck ( al, adata, oaextend, reductionx ) ;

		reductionx.transformation->show();

	} else {
		printf ( "check: good!\n\n" );

	}
}

void mydebug1 ( array_link al, arraydata_t &adata, OAextend &oaextend, int dverbose=1 )
{
	oaextend.setAlgorithm ( MODE_J5ORDERX, &adata );

	lmc_t r=LMC_MORE;
	int niter=1;

	if ( dverbose<0 ) {
		niter=5;
		dverbose=0 ;
	}

	printf ( "mydebug1! al %d %d\n", al.n_rows, al.n_columns );
	LMCreduction_t reduction ( &adata );
	reduction.init_state=COPY;

	double t0 =get_time_ms();
	LMCreduction_t reductionx=reduction;
	//oaextend.setAlgorithm(MODE_J4);
	for ( int ix=0; ix<niter; ix++ ) {
		r = LMCcheck ( al, adata, oaextend, reductionx ) ;
	}
	double dt =get_time_ms()-t0;
	printf ( "original lmc_t: %d, time: %.3f [ms]\n", ( int ) r, 1e3*dt );

}

void mydebug2 ( array_link al, arraydata_t &adata, OAextend &oaextend, int dverbose=1 )
{
	lmc_t r=LMC_MORE;
	int niter=1;
	if ( dverbose<0 ) {
		niter=10;
		dverbose=0 ;
	}
	printf ( "mydebug2! al %d %d\n", al.n_rows, al.n_columns );
	oaextend.setAlgorithm ( MODE_J5ORDERX, &adata );

	double t0=get_time_ms(), dt;

	LMCreduction_t reduction ( &adata );
	//  oaextend.setAlgorithm(MODE_J5ORDERX, &adata);

	// pre-compute
	LMCreduction_t reductionsub = calculateSymmetryGroups ( al.deleteColumn ( -1 ), adata,  oaextend );

	dt =get_time_ms()-t0;
	printf ( "pre-compute: time: %.3f [ms]\n", 1e3*dt );

	//reductionsub.showColperms(1);
	if ( dverbose ) {
		reductionsub.symms.showColcombs ( 1 );
		reductionsub.symms.showSymmetries ( 1 );
	}

	printf ( "running LMCcheckXX: " );
	adata.show();
	lmc_t rn = r;
	t0 =get_time_ms();
	for ( int ix=0; ix<niter; ix++ ) {
		copy_array ( al.array, reduction.array, adata.N, adata.ncols ); // hack?
		rn = LMCcheckSymmetryMethod ( al, adata, oaextend, reduction, reductionsub, dverbose ) ;
	}
	dt =get_time_ms()-t0;
	printf ( "new lmt_t %d: time: %.3f [ms]\n", rn, 1e3*dt );

}

void mydebug2d ( array_link al, arraydata_t &adata, OAextend &oaextend, int dverbose=1 )
{
	lmc_t r=LMC_MORE;
	int niter=1;
	if ( dverbose<0 ) {
		niter=10;
		dverbose=0 ;
	}
	printf ( "mydebug2! al %d %d\n", al.n_rows, al.n_columns );
	oaextend.setAlgorithm ( MODE_LMC_SYMMETRY, &adata );

	double t0=get_time_ms(), dt;

	LMCreduction_t reduction ( &adata );

	adata.show();
	lmc_t rn = r;
	t0 =get_time_ms();
	for ( int ix=0; ix<niter; ix++ ) {
		copy_array ( al.array, reduction.array, adata.N, adata.ncols ); // hack?
		reduction.init_state=COPY;
		rn = LMCcheck ( al, adata, oaextend, reduction ) ;
	}
	dt =get_time_ms()-t0;
	printf ( "lmc_symmetry direct lmt_t %d: time: %.3f [ms]\n", rn, 1e3*dt );

}

#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>

#include "Deff.h"
#include "arraytools.h"
#include "strength.h"

using namespace Eigen;


//std::vector<double> Defficiencies (const array_link &al, const arraydata_t & arrayclass, int verbose ) ;

array_link  finalcheck ( const array_link &al,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niterx, int nabortx )
{
	int nabort=al.n_rows*al.n_columns+2;
	int niter= nabort+5;
	array_link  al2 = optimDeff2level ( al, arrayclass, alpha, verbose>=3,  DOPTIM_UPDATE,  niter, nabort );


	std::vector<double> dd0 = al.Defficiencies();
	double d0 = scoreD ( dd0, alpha );

	std::vector<double> dd = al2.Defficiencies();
	double d = scoreD ( dd, alpha );

	if ( d>d0 ) {
		printf ( "finalcheck: %f -> %f\n", d0, d );
	}

	return al2;
}


//TODO: use optimDeff2level by default
//TODO: make enum for selection strategy?
//TODO: IMPLEMENT flip method as well
//TODO: Run valgrind...

int main ( int argc, char* argv[] )
{
	AnyOption opt;
	/* parse command line options */
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setOption ( "output", 'o' );
	opt.setOption ( "rand", 'r' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "ii", 'i' );
	opt.setOption ( "dverbose", 'd' );
	opt.setOption ( "rows" );
	opt.setOption ( "cols" );
	opt.setOption ( "nrestarts" );
	opt.setOption ( "niter" );
	opt.setOption ( "mdebug", 'm' );
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.addUsage ( "Orthonal Array: oatest: testing platform" );
	opt.addUsage ( "Usage: oatest [OPTIONS] [FILE]" );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.processCommandArgs ( argc, argv );


			int randvalseed = opt.getIntValue ( 'r', 1 );

				srand ( randvalseed );

				
	print_copyright();
	//cout << system_uname();
	setloglevel ( NORMAL );

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	setloglevel ( SYSTEM );

		if (1) {
		int ii = opt.getIntValue ( 'i', 1 );
		array_link al = exampleArray(ii, 1);
		
		al.reduceDOP();
		//al=al.randomperm();
		al=al.randomcolperm();
		
		arraydata_t arrayclass = arraylink2arraydata(al);
		array_transformation_t trans ( &arrayclass );
		trans.randomizecolperm();
		trans.randomizerowperm();
		trans.randomize();
		al=trans.apply(al);

	
		
		array_link alr = al.reduceDOP();
		
		
		array_transformation_t at = reductionDOP(al, 1);
		array_link alr2 = at.apply(al);
		
//alr2 = at.inverse().apply(al);

		al.showarraycompact();
		printf("----\n");
		alr.showarraycompact();
		printf("----\n");
		alr2.showarraycompact();
		
		if (alr==alr2)
			printf("arrays equal!\n");
		else
			printf("arrays unequal!\n");
	return 0;	
	}
	
	if (1) {
		int ii = opt.getIntValue ( 'i', 1 );
		array_link al = exampleArray(ii, 1);
		
		double A= al.Aefficiency();
		std::vector<double> aa = Aefficiencies(al, 1);
		
		printf(" A : %f -> %f\n" , A, aa[0] );
		printf(" A A1 A2: %f, %f, %f\n" , aa[0], aa[1], aa[2] );
		
		al.showarray();
		
	return 0;	
	}
	if (0)
	{
		int r = opt.getIntValue ( 'r', 8 );
		int niter = opt.getIntValue ( 'i', 100 );

		int race=0;
		
		#pragma omp parallel for // num_threads(8) // schedule(dynamic,1)
		for ( int i=0; i<r; i++ ) {
			array_link al = exampleArray ( 12 );
			//al = exampleArray ( 1 );
			al = exampleArray ( 2 );
			//al = exampleArray ( 12 );
			//al = al.selectFirstColumns(9);
//		int strength  = 3;
			
			int r=race;
			race++;
			if( r+1!=race )  {
				printf("race condition\n");
			}
			
			int N = al.n_rows;
			//printf ( "iter loop %d: strength %d\n", i, ad.strength ) ;
			printf ( "iter loop %d:\n", i ) ;
			for ( int ij=0; ij<1; ij++ ) {
				for ( int j=0; j<niter; j++ ) {
#ifdef _OPENMP
					printf(" loop %d: tid %d\n", j, omp_get_thread_num() );
#endif
					OAextend oaextend;
			arraydata_t ad = arraylink2arraydata ( al, 0, al.strength() );

					arraydata_t adlocal ( &ad, al.n_columns );
					LMCreduction_t reduction ( &adlocal ); // TODO: place outside loop (only if needed)
					// make sure code is thread safe
					reduction.initStatic();

					{
						OAextend oaextendx=oaextend;
						//oaextendx.setAlgorithm ( MODE_J5ORDERXFAST, &adlocal );
						oaextendx.setAlgorithm ( MODE_ORIGINAL, &adlocal );
						oaextendx.extendarraymode=OAextend::NONE;
						reduction.init_state=COPY;

						int extensioncol = al.n_columns-1;
						arraylist_t extensions0;

//printf("here\n");
//	        setloglevel(NORMAL);
						extend_array ( al.array,  &adlocal, extensioncol, extensions0, oaextendx );
//printf("here 2\n");
					}

					//lmc_t lmc =  LMCcheck ( al, ( adlocal ), oaextend, reduction ); // omp good
				}



			}

		}
		exit ( 0 );
	}

	if (0)
	{
		int r = opt.getIntValue ( 'r', 3 );

		array_link al = exampleArray ( 12 );
		int strength  = 3;
		arraydata_t ad = arraylink2arraydata ( al, 0, strength );
		int N = al.n_rows;
		for ( int ij=0; ij<2; ij++ ) {
			for ( strength=1; strength<5; strength++ ) {
//int s = strength_check ( al,  strength );
				ad = arraylink2arraydata ( al, 0, strength );
				int s = strength_check ( ad, al, 1 );
				printf ( "iter loop %d: strength %d: %d\n", ij, strength, s );
			}
			strength=r;
			ad = arraylink2arraydata ( al, 0, strength );
			for ( int j=0; j<12000; j++ ) {
				//strength_check ( al,  strength );
				strength_check ( ad, al );
			}
		}
		exit ( 0 );

	}
	{

		array_link al = exampleArray ( 12 );
		int N = al.n_rows;
		colindex_t firstcolcomb[5];
		init_perm ( firstcolcomb, 5 );
		firstcolcomb[4]=8;
		//retx= jj45split ( al.array, N, jj, firstcolcomb, ad, oaextend, tmpStatic, reduction );
		for ( int ij=0; ij<20; ij++ ) {
			printf ( "iter loop %d\n", ij );
			for ( int j=0; j<200000; j++ ) {
				//jj45_t xx = jj45val ( al.array, N, 5, firstcolcomb, -1, 1 );
				fastj(al.array, N, 5, firstcolcomb );
			}
		}
		exit ( 0 );
	}

	if ( 0 ) {
		const int maxval=5;
		int nelements = 8;
		array_t array[8];
		array[0]=maxval-1;
		array[1]=2;
		array[3]=2;
		array[4]=0;
		array[5]=4;
		array[6]=1;

		int elements[maxval];
		for ( int i=0; i<10000; i++ ) {
			std::fill ( elements, elements+maxval, 0 );
			for ( int j=0; j<30000; j++ ) {
				countelements ( array, nelements, maxval, elements );
			}
		}
		exit ( 0 );
	}

	int r = opt.getIntValue ( 'r', 0 );
	int dverbose = opt.getIntValue ( 'd', 1 );
	int md = opt.getIntValue ( 'm', 0 );
	int aidx = opt.getIntValue ( 'i', 0 );
	int rr = opt.getIntValue ( "rows", 40 );
	int cc = opt.getIntValue ( "cols", 7 );

	srand ( r );
	int verbose = opt.getIntValue ( 'v', 1 );
	setloglevel ( verbose );


	if ( 1 ) {
		array_link A ( 0,8,0 );
		printf ( "should return an error\n  " );
		A.Defficiencies();

		A = array_link ( 1,8,0 );
		printf ( "should return an error\n  " );
		A.at ( 0,0 ) =-2;
		A.Defficiencies();
		return 0;
	}


	if ( 0 ) {
		double dt, mean;
		arraydata_t arrayclass ( 2, rr, 0, cc );
		int nrestarts=opt.getIntValue ( "nrestarts", 10 );
		int niter=opt.getIntValue ( "niter", 160000 );
		std::vector<double> alpha ( 3 );
		alpha[0]=1;
		alpha[1]=0.5;
		alpha[2]=0;
		array_link al = arrayclass.randomarray();
		int ro = aidx;

		double t0=get_time_ms();
		for ( int j=0; j<nrestarts; j++ ) {
			seedfastrand ( ro+100*j );
			srand ( ro+100*j );
			array_link  al2 = optimDeff ( al, arrayclass, alpha, 0,  DOPTIM_UPDATE,  niter, 0 );
		}
		printf ( "dt %.1f [ms]\n", 1e3*get_time_ms ( t0 ) );
		return 0;
	}

	if ( 1 ) {
		//arraylist_t sols = readarrayfile("/home/eendebakpt/misc/oa/oacode/testdata/design-D-nearzero.oa");
		arraylist_t sols = readarrayfile ( "/home/eendebakpt/misc/oa/oacode/testdata/design-D-nearzero2.oa" );
		array_link al = sols[0];
		arraydata_t arrayclass = arraylink2arraydata ( al );
		arrayclass=arraydata_t ( 2, 70, 0, 7 );

		if ( aidx>10 )
			al=arrayclass.randomarray ( 1 );
		double D = Defficiency ( al,2 );

		printf ( "---------\n" );
		std::vector<double> dd = Defficiencies ( al,  arrayclass, 2, 0 );


		double t0=get_time_ms();
		for ( int i=0; i<1000; i++ ) {
			double D = Defficiency ( al,0 );
		}
		printf ( "time %.1f [ms]\n", 1e3*get_time_ms ( t0 ) );
		t0=get_time_ms();
		for ( int i=0; i<1000; i++ ) {
			Defficiencies ( al,arrayclass,0,0 );
		}
		printf ( "time %.1f [ms]\n", 1e3*get_time_ms ( t0 ) );

		//throw; printf("throw done!\n");
		return 0;
	}


	if ( 0 ) {
		int N = 10;
		int k =3;
		const int nn = N*k;
		std::vector<int> updatepos = permutation<int> ( nn );
		if ( r==0 ) {
			std::random_shuffle ( updatepos.begin(), updatepos.end() );
		}
		int updateidx = 0;

		int niter=nn;

//#pragma omp for
		for ( int ii=0; ii<niter; ii++ ) {
			// select random row and column
			int r = updatepos[updateidx] % N;
			int c = updatepos[updateidx] / N;
			printf ( "ii %d: pos %3d: %3d %3d\n", ii, updatepos[ii], r, c );

			updateidx++;
		}

		return 0;
	}

	if ( 0 ) {
		double dt, mean;
		arraydata_t arrayclass ( 2, rr, 0, cc );
		int nrestarts=opt.getIntValue ( "nrestarts", 10 );
		int niter=opt.getIntValue ( "niter", 60000 );
		std::vector<double> alpha ( 3 );
		alpha[0]=1;
		alpha[1]=0;
		alpha[2]=0;
		int nabort=0;
		array_link al = arrayclass.randomarray();
		array_link al3;
		al = arrayclass.randomarray();
		int ro = aidx;

		int ni=10;
		std::vector<double> mm ( ni );
		std::vector<double> m ( nrestarts );

//					seedfastrand(ro); srand(ro); nabort=md; array_link  al2 = optimDeff ( al, arrayclass, alpha, verbose,  DOPTIM_UPDATE,  niter, nabort ); return 0;


		double t0=get_time_ms();
		for ( int j=0; j<ni; j++ ) {
			seedfastrand ( ro+100*j );
			srand ( ro+100*j );
			al = arrayclass.randomarray();
			for ( int i=0; i<nrestarts; i++ ) {
				seedfastrand ( i+ro +100*j );
				srand ( i+ro+100*j );
				array_link  al2 = optimDeff ( al, arrayclass, alpha, dverbose,  DOPTIM_UPDATE,  niter, nabort );
				//DoptimReturn rr=Doptimize(arrayclass, nrestarts, niter, alpha, 1, 0, 1000, nabort);
				//array_link  alx = finalcheck(al2, arrayclass, alpha, verbose,  DOPTIM_UPDATE,  niter, nabort);

				std::vector<double> dd = al2.Defficiencies();
				double d = scoreD ( dd, alpha );
				m[i]=d;

				if ( 1 ) {
					srand ( i+ro+200*j );
					//al2 = optimDeff ( al2, arrayclass, alpha, dverbose,  DOPTIM_UPDATE,  niter, 3000 );
					al3 = optimDeff ( al2, arrayclass, alpha, dverbose,  DOPTIM_SWAP,  niter, 6300 );
					dd = al3.Defficiencies();
					d = scoreD ( dd, alpha );
					if ( d> m[i]+1e-8 ) {
						printf ( "check! %.8f (%.6f %.6f) -> %.8f (%.6f %.6f)\n", m[i], al2.Defficiency(), al2.DsEfficiency() , d, al3.Defficiency(), al3.DsEfficiency() );
					}
					//m[i]=d;

				}
			}
			mm[j]=*std::max_element ( m.begin(), m.end() ) ;
			printf ( "  mm[%d] %.4f\n", j, mm[j] );

		}
		dt=get_time_ms()-t0;
		mean = accumulate ( mm.begin(), mm.end(), 0.0 ) / mm.size();
		printf ( "time %.1f [ms]: mean %.8f\n", dt*1e3, mean );



		t0=get_time_ms();
		for ( int j=0; j<ni; j++ ) {
			seedfastrand ( ro+ 100*j );
			srand ( ro+ 100*j );
			al = arrayclass.randomarray();

			for ( int i=0; i<nrestarts; i++ ) {
// al = arrayclass.randomarray();
				seedfastrand ( i+ro +100*j );
				srand ( i+ro+100*j );

				nabort=al.n_columns*al.n_rows+1;
				array_link  al2 = optimDeff2level ( al, arrayclass, alpha, dverbose,  DOPTIM_UPDATE,  niter, nabort );
				//DoptimReturn rr=Doptimize(arrayclass, nrestarts, niter, alpha, 1, 0, 1000, nabort);
				std::vector<double> dd = al2.Defficiencies();
				double d = scoreD ( dd, alpha );

				m[i]=d;
			}
			mm[j]=*std::max_element ( m.begin(), m.end() ) ;
			printf ( "  mm[%d] %.4f\n", j, mm[j] );

		}
		dt=get_time_ms()-t0;
		mean = accumulate ( mm.begin(), mm.end(), 0.0 ) / mm.size();
		printf ( "time %.1f [ms]: mean %.8f\n", dt*1e3, mean );
		return 0;

	}










	// ################################


	if ( 1 ) {
		arraydata_t adata ( 2, rr, 0, cc );
		int niter=10000;
		std::vector<double> alpha ( 3 );
		alpha[0]=1;
		alpha[1]=1;
		double t0=get_time_ms();
		arraylist_t sols;
		for ( int i=0; i<20; i++ ) {
			array_link al = adata.randomarray ( 0 );
			sols.push_back ( al );
		}
		DoptimReturn a = DoptimizeMixed ( sols, adata, alpha, verbose );
		DoptimReturn a2 = DoptimizeMixed ( a.designs, adata, alpha, verbose );
		DoptimReturn a3 = DoptimizeMixed ( a2.designs, adata, alpha, verbose );

		// TODO: incorporate into python code, C++ solution for sorting designs
		// balance iteration with niter
		return 0;
	}

	if ( 0 ) {
		arraydata_t adata ( 2, 64, 0, 7 );
		int niter=10000;
		std::vector<double> alpha ( 3 );
		alpha[0]=1;
		alpha[1]=1;
		double t0=get_time_ms();
		for ( int i=0; i<10000+10000*aidx; i++ ) {
			array_link al = adata.randomarray ( 0 );
			//std::vector<double> dd = al.Defficiencies();

			arraydata_t arrayclass = arraylink2arraydata ( al );
			std::vector<double> dd = Defficiencies ( al, arrayclass, 0, 0 );

			//int dmethod = DOPTIM_SWAP;
			//dmethod = DOPTIM_NONE;
			//array_link  alx = optimDeff ( al, adata, alpha, 2, dmethod, niter, 3000 );
		}
		printf ( "dt %.3f [s]\n", get_time_ms()-t0 );

		return 0;
	}

	if ( 0 ) {

		arraydata_t adata ( 2, 80, 0, 7 );
		int niter=10000;
		std::vector<double> alpha ( 3 );
		alpha[0]=1;
		alpha[1]=1;
		double t0=get_time_ms();
		for ( int jjj=0; jjj<10; jjj++ ) {
			for ( int i=0; i<4800; i++ ) {
				array_link al = adata.randomarray ( 0 );
				std::vector<double> dd = al.Defficiencies();
			}
		}
		printf ( "dt %.3f [s]\n", get_time_ms()-t0 );
	}
	if ( 1 ) {
		arraydata_t adata ( 2, 80, 0, 7 );
		int niter=10000;
		std::vector<double> alpha ( 3 );
		alpha[0]=1;
		alpha[1]=1;
		double t0=get_time_ms();
		for ( int i=0; i<10; i++ ) {
			array_link al = adata.randomarray ( 0 );
			std::vector<double> dd = al.Defficiencies();
			int dmethod = DOPTIM_SWAP;
			//dmethod = DOPTIM_NONE;
			array_link  alx = optimDeff ( al, adata, alpha, 2, dmethod, niter, 3000 );
		}
		printf ( "dt %.3f [s]\n", get_time_ms()-t0 );

		return 0;
	}

	if ( 0 ) {
// test PEC sequence
		srand ( get_time_ms() );
		array_link al ( rr,cc,-1 );
		arraydata_t arrayclass ( 2, rr, 1, cc );
		arrayclass.show();
		for ( int i=0; i<al.n_columns*al.n_rows; i++ )
			al.array[i]=rand() %2;

		double t0= ( get_time_ms() );

		std::vector<double> pec = PECsequence ( al );
		printf ( "PEC: " );
		display_vector ( pec );
		printf ( " \n" );
		printf ( "dt: %.1f [s]\n" , get_time_ms()-t0 );
		return 0;
	}


	if ( 0 ) {
		array_link al = exampleArray ( aidx );
		al.showarray();
		return 0;
	}
	{
		for ( int k=0; k<aidx; k++ )
			fastrand();

		const int N = 9;
		arraydata_t arrayclass ( 3, N, 2, 2 );
		arrayclass.show();
		array_link al0=arrayclass.randomarray ( 0 )	;
		//al0 = array_link(3,2,0);al0.at(0,0) =0; al0.at(1,0)=0; al0.at(2,0)=1;	al0.at(0,1) =0; al0.at(1,1)=1; al0.at(2,1)=1;
		al0.showarray();


		printf ( "D %f\n", al0.Defficiency() );
		printf ( "-----\n" );
		MatrixFloat mm= al0.getModelMatrix ( 2 );
		std::cout << mm << std::endl;

		printf ( "-----\n" );
		mm = array2eigenModelMatrix ( al0 );
		std::cout << mm << std::endl;
		exit ( 0 );
	}


	{
		array_link al = exampleArray ( aidx );
		al.show();
		std::cout << al.getModelMatrix ( 2 ) << std::endl;
		printf ( "------\n" );
		std::cout << array2eigenModelMatrix ( al ) << std::endl;
		exit ( 0 );
		al.showarray();
		al.show();

		MatrixFloat ww;
		ww = array2eigenModelMatrix ( al );
		for ( int ii=0; ii<10; ii++ ) {
			al.Defficiencies();
			//detXtX(ww);
		}
		exit ( 0 );

		std::vector<double> dd;
		for ( int ii=0; ii<10000; ii++ ) {
			al = exampleArray ( aidx );
			dd=al.Defficiencies();
		}
		printf ( "dd " );
		display_vector ( dd );
		std::cout << std::endl;
		exit ( 0 );


	}
	{
		array_link al = exampleArray ( aidx );
		arraydata_t arrayclass=arraylink2arraydata ( al );

		arrayclass.randomarray ( 1 ).showarray();
		return 0;

		al = exampleArray ( aidx );
		al.showarray();

		MatrixFloat ME = array2eigenME ( al );
		std::cout << "---- ME ----\n";
		std::cout << ME << std::endl;


		//std::pair<Eigen::MatrixXd,Eigen::MatrixXd> mm3 = array2eigenModelMatrix2 ( al, 2 );
		//std::cout <<"gr\n" << mm3.first << std::endl;

//		std::pair<Eigen::MatrixXd,Eigen::MatrixXd> mmx = array2eigenModelMatrix2 (al );

		std::vector<double> dd = al.Defficiencies ( 2 );
		printf ( "Defficiencies: %f %f %f\n", dd[0], dd[1], dd[2] );
		printf ( "Defficiency: %f\n", al.Defficiency() );
		//std::cout << "mmx\n"<< mmx.first << std::endl ;

		if ( verbose>=2 ) {

			MatrixFloat mm = al.getModelMatrix ( 2, 1 );
			std::cout << "model matrix\n" << mm << std::endl;
		}

		return 0 ;

	}

	if ( 0 ) {
		const char *fname = "/home/eendebakpt/tmp/test.oa";
		arraylist_t ll = readarrayfile ( fname );
		printf ( "read %ld arrays\n", ll.size() );
		array_link al= ll[0];

		al.show();
		al.showarray();
		printf ( "go!\n" );
		array_link r = reduceDOPform ( al,1 );
		printf ( "done!\n" );
//	sleep(1);
		std::cout.flush();
		exit ( 0 );
	}


	if ( 0 ) {
		const char *fname = "/home/eendebakpt/tmp/x.oa";
		const char *fm = "r+b";
		FILE *nfid = fopen ( fname, fm );
		printf ( "## opened file with mode %s", fm );
		printf ( "  " );
		printfd ( "file is at position %d\n", ftell ( nfid ) );

		char buf[4]= {1,2,3,4};
		int p = fseek ( nfid, 0, SEEK_END );
		printf ( "  p %d, nfid %ld\n", p, ( long ) nfid );
		printf ( " fseek errno: %s\n", strerror ( errno ) );
		int rr=fread ( buf, 1, 4, nfid );
		printf ( "rr %d, something went wrong with fread()? %s\n", rr, strerror ( errno ) );
		int pp = ftell ( nfid );
		fseek ( nfid, 0, SEEK_CUR );
		rr=fwrite ( buf, 1, 4, nfid );
		printf ( "rr %d, something went wrong with fwrite()! %s\n", rr, strerror ( errno ) );
		printf ( "  written rr %d\n", rr );
//afile->seek(0);
		exit ( 0 );
	}


	// load array

	const char *afile="/home/eendebakpt/misc/oa/oacode/testdata/test-547-21.oa";
	//const char *afile="/home/eendebakpt/misc/oa/oacode/build/xx32-11.oa";
	// afile="/home/eendebakpt/misc/oa/oacode/build/xx32-6.oa";
	//afile="/home/eendebakpt/misc/oa/oacode/build/xx12.oa";
	// afile="/home/eendebakpt/misc/oa/oacode/build/xl.oa";
	afile="/home/eendebakpt/misc/oa/oacode/build/x3.oa";
//     afile="/home/eendebakpt/misc/oa/oacode/build/x3r.oa";

	if ( opt.getArgc() >0 ) {
		afile = opt.getArgv ( 0 );
	}
	printf ( "using array file %s\n", afile );

	arraylist_t lst= readarrayfile ( afile );

	if ( 0 ) {
		int nn=lst.size();
		for ( int ii=0; ii<nn; ii++ ) {
			lst.push_back ( lst[ii] );
		}
		nn=lst.size();

//    int ii=420; lst.erase(lst.begin(), lst.begin()+ii); lst.erase(lst.begin()+1, lst.end());

//#define DOOPENMP 1

		std::vector<int> rr ( lst.size() );
#ifdef _OPENMP
		printf ( "openmp: num arrays %ld, omp_get_num_threads() %d\n", lst.size(), omp_get_num_threads() );
#endif
		int lmccount=0;

#ifdef _OPENMP
		omp_set_num_threads ( 4 );
#endif

		#pragma omp parallel for
		for ( int j=0; j< ( int ) lst.size(); j++ ) {
			#pragma omp critical
			{
#ifdef _OPENMP
				printf ( "oatest: j %d: thread %d\n", j, omp_get_thread_num() );
#endif
			}
			array_link ee = lst[j];
			int extcol=ee.n_columns;
			arraydata_t adlocal=arraylink2arraydata ( ee,0,3 );
			LMCreduction_t reduction ( &adlocal ); // TODO: place outside loop (only if needed)
			LMCreduction_t tmp = reduction;

			OAextend oaextend;
			oaextend.setAlgorithm ( MODE_J5ORDERXFAST, &adlocal );

			//getGlobalStatic;
			#pragma omp critical
			{
#ifdef DOOPENMP
				int jg = j % omp_get_num_threads();
				jg=j % 8;

				printf ( "  %s: j %d, reduction.staticdata %ld\n", __FUNCTION__, ( int ) j, ( long ) reduction.staticdata );

				// #pragma omp critical
				{
					//           reduction.initStatic();
					reduction.staticdata=new LMC_static_struct_t();
					reduction.staticdata->update ( &adlocal );
				}

				printf ( "openmp: col %d: j %d/%ld, jg %d\n", extcol, j, ( long ) nn, jg );


#else
				// still make sure code is thread safe
				reduction.initStatic();

#endif

			}

			double t0, dt;
			lmc_t lmc;

			reduction.init_state=COPY;

			t0=get_time_ms();
			printfd ( "at critical point: LMCcheck... j %d\n", ( int ) j );
			//#pragma omp critical
			lmc =  LMCcheck ( ee, ( adlocal ), oaextend, reduction );
			dt=get_time_ms()-t0;

			#pragma omp critical
			{
				if ( lmc!=LMC_LESS )
					lmccount++;
				rr[j]=lmc;
				reduction.releaseStatic();
			}


		}


		printf ( "lmccount %d\n", lmccount );


		exit ( 0 );

	}


	if ( 0 ) {
		printf ( "FIXME: caching of J45 values?\n" );
		printf ( "split acts on root part??\n" );
		printf ( "fixme: should be perm instead of comb? ! ?\n" );
		printf ( "we should be able to skip column selection in root stage\n" );

		printf ( "FIXME: calculation of colperms group -> should be at most a factor 2 extra (since array is 1 column less!\n" );
		printf ( "FIXME: calculation of colperms group -> eliminate overhead in sort and unique\n" );
		printf ( "FIXME: LMCreduce -> for non root initialized arrays there are some funny comparisons!!!\n" );
		printf ( " 	i) reduction.array is not initialized at all\n" );
		printf ( " 	ii) check_root_update does not use dyndata\n" );
		printf ( " 	iii) create unit tests\n" );
		printf ( "FIXME: calculation of colperms group -> share between multiple arrays (LATER)\n" );
	}

	int strength = lst[0].strength();
	//strength=1;
	arraydata_t adata = arraylink2arraydata ( lst[0], 0, strength );
	LMCreduction_t reduction ( &adata );


	if ( 0 ) {
		OAextend oaextend;
		oaextend.setAlgorithm ( MODE_J5ORDERX, &adata );
		oaextend.j5structure=J5_45;
		double t0 =get_time_ms();
		for ( size_t i=0; i<lst.size(); i++ ) {
			printf ( "array %ld\n", i );
			array_link al = lst[i];
			lmc_t r = LMCcheck ( al, adata, oaextend, reduction ) ;
			//lmc_t r = LMCcheckj5(al, adata, reduction, oaextend) ;

			printf ( "lmc_t %d\n" , ( int ) r );
			reduction.symms.showColperms();
		}
		double dt =get_time_ms()-t0;
		printf ( "time: %.3f [ms]\n", 1e3*dt );
	}

	array_link al=lst[0];
//    mydebug(al, adata, oaextend);
	OAextend oaextend;
//    oaextend.setAlgorithm(MODE_J5ORDERX, &adata);
//   oaextend.j5structure=J5_45;


	switch ( md ) {
	case 0:
		mydebug ( al, adata, oaextend, dverbose );
		break;
	case 1:
		mydebug1 ( al, adata, oaextend, dverbose );
		break;
	case 2:
		mydebug2 ( al, adata, oaextend, dverbose );
		break;
	case 3:
		mydebug2d ( al, adata, oaextend, dverbose );
		break;
	case 5:
		mydebug5 ( al, adata, oaextend, dverbose );
		break;
	case 10:
		mydebug10 ( al, adata, oaextend, dverbose );
		break;
	}


	arraysymmetry::rowpermpool.reset();

	return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
