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
	 if (niter>0){
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

	if (0) {
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

	if (0) {
int nc=6;
printf("## symmetrices with %d cols:\n", nc);	
	symmetryset xx= reductionsub.symms.symmetries[nc];
	                    for( symmetryset::const_iterator it = xx.begin(); it != xx.end(); it++) {
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

using namespace Eigen;



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
	opt.setOption ( "mdebug", 'm' );
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.addUsage ( "Orthonal Array: oatest: testing platform" );
	opt.addUsage ( "Usage: oatest [OPTIONS] [FILE]" );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.processCommandArgs ( argc, argv );


	print_copyright();
	//cout << system_uname();
	setloglevel ( NORMAL );

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	int r = opt.getIntValue ( 'r', 0 );
	int dverbose = opt.getIntValue ( 'd', 1 );
	int md = opt.getIntValue ( 'm', 0 );
	int aidx = opt.getIntValue ( 'i', 0 );
    int rr = opt.getIntValue("rows", 40);
    int cc = opt.getIntValue("cols", 7);

	srand ( r );
	int verbose = opt.getIntValue ( 'v', NORMAL );
	setloglevel ( verbose );

	
	if(1)
	    {
// test PEC sequence
      srand(get_time_ms() );
    array_link al(rr,cc,-1);
    arraydata_t arrayclass(2, rr, 1, cc);
    arrayclass.show();
for(int i=0; i<al.n_columns*al.n_rows; i++) al.array[i]=rand()%2;

      double t0=(get_time_ms() );

    std::vector<double> pec = PECsequence(al);
    printf("PEC: " ); display_vector(pec);  printf(" \n" );
printf("dt: %.1f [s]\n" , get_time_ms()-t0);
return 0;
    }

    
	if(0)
	{
			array_link al = exampleArray(aidx);
		al.showarray();
		return 0;
	}
	{
		for(int k=0; k<aidx; k++) fastrand();
		
		const int N = 9;
		arraydata_t arrayclass(3, N, 2, 2);
		arrayclass.show();
	array_link al0=arrayclass.randomarray(0)	;
	 //al0 = array_link(3,2,0);al0.at(0,0) =0; al0.at(1,0)=0; al0.at(2,0)=1;	al0.at(0,1) =0; al0.at(1,1)=1; al0.at(2,1)=1;
	al0.showarray();
	
	
	printf("D %f\n", al0.Defficiency() );
	printf("-----\n");
	Eigen::MatrixXd mm= al0.getModelMatrix(2);
	std::cout << mm << std::endl;
	
	printf("-----\n");
	mm = array2eigenModelMatrix(al0);
	std::cout << mm << std::endl;
	exit(0);
	}
	
	
	{
		array_link al = exampleArray(aidx);
		al.show();
		std::cout << al.getModelMatrix(2) << std::endl;
printf("------\n");
		std::cout << array2eigenModelMatrix(al) << std::endl;
			exit(0);
		al.showarray();
		al.show();
		
		Eigen::MatrixXd ww;
		  ww = array2eigenModelMatrix(al);
		for(int ii=0; ii<10; ii++) {
			al.Defficiencies();
		//detXtX(ww);
		}
		exit(0);
		
		std::vector<double> dd;
		for(int ii=0; ii<10000; ii++) {
			al = exampleArray(aidx);
			 dd=al.Defficiencies();
		}
		printf("dd "); display_vector(dd); std::cout << std::endl;
		exit(0);
		
		
	}
	{
		array_link al = exampleArray(aidx);
		arraydata_t arrayclass=arraylink2arraydata(al);
		
		arrayclass.randomarray(1).showarray();
		return 0;
		
		 al = exampleArray(aidx);
	al.showarray();
	
	Eigen::MatrixXd ME = array2eigenME(al);
	std::cout << "---- ME ----\n";
	std::cout << ME << std::endl;
	

	//std::pair<Eigen::MatrixXd,Eigen::MatrixXd> mm3 = array2eigenModelMatrix2 ( al, 2 );
	//std::cout <<"gr\n" << mm3.first << std::endl;

//		std::pair<Eigen::MatrixXd,Eigen::MatrixXd> mmx = array2eigenModelMatrix2 (al );

	std::vector<double> dd = al.Defficiencies(2);
	printf("Defficiencies: %f %f %f\n", dd[0], dd[1], dd[2]);
	printf("Defficiency: %f\n", al.Defficiency() );
	//std::cout << "mmx\n"<< mmx.first << std::endl ;

	if (verbose>=2) {
	
Eigen::MatrixXd mm = al.getModelMatrix(2, 1);
	std::cout << "model matrix\n" << mm << std::endl;
	}
	
	return 0 ;
	
}
	
	if (0) {
		const char *fname = "/home/eendebakpt/tmp/test.oa";
	arraylist_t ll = readarrayfile(fname);
	printf("read %zu arrays\n", ll.size());
	array_link al= ll[0];
	
	al.show();
	al.showarray();
	printf("go!\n"); 
	array_link r = reduceDOPform(al,1);
	printf("done!\n");
//	sleep(1);
	std::cout.flush();
	exit(0);
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
		printf ( "openmp: num arrays %zu, omp_get_num_threads() %d\n", lst.size(), omp_get_num_threads() );
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
			printf ( "array %zu\n", i );
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
