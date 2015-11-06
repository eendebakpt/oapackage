/** \file oa_depth_extend.cpp

 C++ program: oa_depth_extend

 oa_depth_extend: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

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
#include "strength.h"
#include "extend.h"

#include "oadevelop.h"
//#include "lmc.h"

#ifdef OADEBUG
#else
#define DOOPENMP
#endif

#ifdef DOOPENMP
#include "omp.h"
#endif


/**
 * @brief Read in files with arrays and join them into a single file
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] )
{
	AnyOption opt;

	/* parse command line options */
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setOption ( "output", 'o' );
//	opt.setFlag ( "sort", 's' );
	opt.setFlag ( "split", 's' );
	opt.setFlag ( "extend", 'e' );
	opt.setFlag ( "prune", 'p' );
	opt.setOption ( "j5structure", 'j' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "kfinal", 'k' );
	opt.setOption ( "hack", 'H' );
	opt.setOption ( "method", 'm' );
	opt.setOption ( "writearrays", 'w' );
	opt.setOption ( "maxarrays" );
	opt.setOption ( "discardJ5" );
	opt.setOption ( "maxk", 'K' );
	opt.setOption ( "Avalue", 'A' );
	opt.setFlag ( "coptions" );	/* print compile time options */
	opt.setOption ( "logtime", 'Q' );
	opt.setOption ( "format", 'f' );
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.addUsage ( "Orthonal Array: oa_depth_extend: testing platform" );
	opt.addUsage ( "Usage: oa_depth_extend [OPTIONS] [FILES]" );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.addUsage ( " -f [FORMAT]					Output format (default: TEXT, or BINARY) " );
	opt.addUsage ( " -o [FILE] --output [FILE]	Output prefix (default: test-depth) " );
	opt.addUsage ( " -A [Afinal] " );
	opt.addUsage ( " -v [verbose level] " );
	opt.addUsage ( " --split 				Split calculation " );
	opt.addUsage ( " -p prune, -e extend " );
	opt.addUsage ( " -K --maxK [n] " );

	opt.addUsage ( " -Q [time] 			Logging time (in seconds)" );
	opt.addUsage ( " --j5structure [VALUE] 			..." );
	std::string ss = printfstring ( " -m [MODE]			Algorithm (" ) + algorithm_t_list() + ")\n" ;
	opt.addUsage ( ss.c_str() );
	opt.processCommandArgs ( argc, argv );

	print_copyright();
	setloglevel ( NORMAL );

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	if ( opt.getValue ( "coptions" ) != 0 ) {
		opt.printUsage();
		print_options ( cout );
		exit ( 0 );
	}

	logstream ( SYSTEM ) << "#time start: "<< currenttime() << std::endl;
	double time0=get_time_ms();

	const char *oaconfigfile = "oaconfig.txt";

	if ( opt.getValue ( "oaconfig" ) !=NULL )
		oaconfigfile = opt.getValue ( 'c' );

	arrayfile::arrayfilemode_t mode = arrayfile::ABINARY_DIFF;

#ifdef _OPENMP
	//omp_set_num_threads(4);
	printf ( "openmp: num threads %d, max num threads %d\n", omp_get_num_threads(), omp_get_max_threads() );
	//omp_set_dynamic(1);

	if ( 0 ) {
		//printf ( "openmp: omp_get_dynamic() %d, omp_get_nested() %d, omp_get_max_active_levels() %d\n", omp_get_dynamic(), omp_get_nested(), omp_get_max_active_levels() );
		omp_set_nested ( 1 );
		//omp_set_dynamic(4);
		printf ( "openmp: num threads %d\n", omp_get_max_threads() );
//		printf ( "openmp: omp_get_dynamic() %d, omp_get_nested() %d, omp_get_max_active_levels() %d\n", omp_get_dynamic(), omp_get_nested(), omp_get_max_active_levels() );
	}
	if ( omp_get_nested() !=0 )
		printf ( "note: omp_get_nested()=%d, make sure to set OMP_THREAD_LIMIT!\n", omp_get_nested() );
#endif

	if ( opt.getValue ( "format" ) !=0 ) {
		std::string format = opt.getValue ( 'f' );
		mode = arrayfile_t::parseModeString ( format );
	}

	std::string outputprefix = "test-depth";
	if ( opt.getValue ( "output" ) !=0 )
		outputprefix = opt.getValue ( 'o' );

	int verbose = opt.getIntValue ( 'v', 2 );
	double Afinal = opt.getDoubleValue ( 'A', 0 );
	int kfinal = opt.getIntValue ( 'k', 7 );
	int splitcalc = opt.getFlag ( "split" );
	int discardJ5 = opt.getIntValue ( "discardJ5", -1 ); // discard arrays with J5==N with size larger then specfied number of columns

	int maxk = opt.getIntValue ( "maxk", 1000 );
	int writedeptharrays = opt.getIntValue ( 'w', 1 );
	long maxarrays = opt.getIntValue ( "maxarrays", 2047483000 );

	int hack = opt.getIntValue ( 'H', 0 );

#ifdef OADEBUG
	if ( hack ) {
		printfd ( "setting hack to on\n" );
		globalHackOption ( 0, 1 );
	} else {
		globalHackOption ( 0, 0 );
	}
	printfd ( "###### hack: %d\n", globalHackOption ( 0 ) );
#endif

	int doextend = opt.getFlag ( 'e' );
	int dopruning = opt.getFlag ( 'p' );
	j5structure_t j5structure = ( j5structure_t ) opt.getIntValue ( 'j', J5_45 );
	int method = opt.getIntValue ( 'm', MODE_J5ORDERX );

	arraydata_t *adfull = readConfigFile ( oaconfigfile );
	if ( adfull==0 ) {
		fprintf ( stderr, "oa_depth_extend: could not read config file" );
		exit ( 1 );
	}
	if ( maxk<adfull->ncols ) {
		if ( maxk<=5 ) {
			printf ( "oa_depth_extend: maxk should be >= 5\n" );
			return 1;
		}
		arraydata_t *adfullr = new arraydata_t ( adfull, maxk );
		delete adfull;
		adfull = adfullr;
	}
	arraylist_t *arraylist = new arraylist_t;
	if ( opt.getArgc() ==0 ) {
		cout << "  adding root to arraylist\n";
		array_link al ( adfull->N, adfull->strength, -1 );
		al.create_root ( *adfull );
		arraylist->push_back ( al );
	} else {
		/* read in the arrays */
		cout << "Reading " << opt.getArgc() << " file(s)" << endl;
		for ( int i = 0 ; i < opt.getArgc() ; i++ ) {
			int n = readarrayfile ( opt.getArgv ( i ), arraylist );
			cout << "file " <<  opt.getArgv ( i ) << ":   read " << n << " array(s)" << endl;
		}
	}

	arraylist_t earrays;


	OAextend oaextendx;
	oaextendx.setAlgorithm ( ( algorithm_t ) method, adfull );
	oaextendx.j5structure=j5structure;

	oaextendx.checkarrays=0;
	oaextendx.use_row_symmetry=0;
	oaextendx.extendarraymode=1;
	oaextendx.init_column_previous=INITCOLUMN_J5;
	//oaextendx.init_column_previous=INITCOLUMN_ZERO;
	oaextendx.nLMC=500000;
	oaextendx.info();

	if ( verbose )
		printf ( "oa_depth_extend: initialization phase (splitcalc %d)\n", splitcalc );
	depth_extend_t dextend ( adfull, 10, discardJ5 );
	dextend.oaextend = oaextendx;
	dextend.logtime= opt.getDoubleValue ( "logtime", 240 );
	dextend.setNarraysMax ( maxarrays );

	int initcolumn=adfull->strength+1;
	if ( arraylist->size() >0 )
		initcolumn=arraylist->at ( 0 ).n_columns;

	dextend.arraywriter = new arraywriter_t();
	dextend.arraywriter->initArrayFiles ( *dextend.ad, initcolumn+1, outputprefix, mode );
	dextend.arraywriter->writearrays=writedeptharrays;
	dextend.counter = new counter_t ( adfull->N );
	dextend.loglevelcol=7;

	setloglevel ( SYSTEM );
	//setloglevel (QUIET );

	double t0=get_time_ms();

	setdtsymm ( 0 );
	setdtlmc ( 0 );

	depth_extensions_storage_t *ds = 0;

	if ( splitcalc ) {
		printf ( "splitting calculation\n" );
		ds = new depth_extensions_storage_t( );
		ds->resize ( arraylist->size() );
	} else {
		//exit(0);
	}
	// loop over all arrays

	//printf("loop start\n"); double t0x=get_time_ms(); // HACK
	//dextend.oaextend.extendarraymode=OAextend::NONE; // HACK

	//FIXME2: enable this parallel loop?, num_threads(4)
	#pragma omp parallel for  schedule(dynamic,1)
	for ( int ai=0; ai< ( int ) arraylist->size(); ai++ ) {
		const array_link &al = arraylist->at ( ai );

		if ( al.n_columns>=maxk )
			continue;
//printfd("dextend.setposition ( al.n_columns=%d, ai=%d, arraylist->size(), 0, 0 )\n", al.n_columns, ai );

		//printfd ( "  openmp: num threads %d, omp_get_dynamic() %d, omp_get_nested() %d, omp_get_max_active_levels() %d\n", omp_get_max_threads(), omp_get_dynamic(), omp_get_nested(), omp_get_max_active_levels() );

		myassert ( adfull->N==al.n_rows, "oa_depth_extend: nrows array, nrows config\n" );

		if ( verbose>=3 || ( verbose>=2 && ai%40==0 ) )
			printf ( "oa_depth_extend: array %d/%ld (%d %d): time %.1f [s]\n", ai, arraylist->size(), al.n_rows, al.n_columns, get_time_ms()-t0 );
		ff();


		depth_extend_t dextendloop ( dextend );
		dextendloop.setposition ( al.n_columns, ai, arraylist->size(), 0, 0 );

		if ( 0 ) {
			double s=0;
			//for(int ix=0; ix<2000; ix++) s+=al.strength(); printf("s %d\n", s);
			for ( int jx=0; jx<1000; jx++ ) {
				for ( int ix=0; ix<20000; ix++ )
					s+=sqrt ( ix );
			}
			printf ( "s %f\n", s );
			continue; // HACK
		}
		if ( 0 ) {
			// HACK
			int extensioncol = al.n_columns;
			arraylist_t extensions0;
			OAextend oaextendx = dextend.oaextend;
			for ( int ni=0; ni<3000; ni++ ) {

				extend_data_t *es = new extend_data_t ( adfull, 10 );
				delete es;
//	extend_array ( al.array,  adfull, extensioncol, extensions0, oaextendx );
			}
			continue;
		}
		depth_extend_array ( al, dextendloop, *adfull, verbose, ds, ai );

		//ds->columnextensionsList[i];
		// 		printfd("#### omp_get_thread_num() %d: num extensions %d\n", omp_get_thread_num(), ds->columnextensionsList[ai].size() );

		//printfd ( "array %zu: thread %d/%d\n", ai , omp_get_thread_num(), omp_get_num_threads() );


		if ( ai%800==0 ) {
			printf ( "array %d: ", ai );
			dextend.counter->showcountscompact();
		}
	}

	//printf("time loop %.3f [s]\n", get_time_ms()-t0x);	exit(0); // HACK

	if ( ds!=0 ) {
		if ( verbose ) {
			printf ( "oa_depth_extend: calling processDepth\n" );
			std::cout << "#time step: " << printfstring ( "%.1f", get_time_ms()-time0 ) << " [s]" << std::endl;
		}

		for ( size_t ai=0; ai<arraylist->size(); ai++ ) {
			printf ( "ai %ld: %ld arrays for extension\n", ai, ds->goodarrayslist[ai].size() );
		}
		dextend.counter->showcounts ( "after init", adfull->strength, adfull->ncols );

		//printf("## adfull: %ld\n",  (long) adfull );
//adfull->show();

		//#pragma omp parallel for schedule(dynamic,1)
		for ( size_t ai=0; ai<arraylist->size(); ai++ ) {
			const array_link &al = arraylist->at ( ai );
			if ( verbose>=3 || ( verbose>=2 && ai%40==0 ) )
				printf ( "oa_depth_extend process: array %ld/%ld (%d %d): time %.1f [s]\n", ai, arraylist->size(), al.n_rows, al.n_columns, get_time_ms()-t0 );
			ff();
			int extensioncol = al.n_columns+1;

			depth_extend_t dextendloop ( dextend );
			dextendloop.extension_column_list=ds->columnextensionsList[ai];
			dextendloop.setposition ( al.n_columns, ai, arraylist->size(), 0, 0 );
			dextendloop.loglevelcol=al.n_columns;

//printf("ai %zu: dextendloop.extension_column_list.size(), =dextendsub.valididx.size(): %zu %zu\n", ai, dextendloop.extension_column_list.size(), ds->dextendsubList[ai].valididx.size() );
//printf("ai %zu: ->ad: %ld\n", ai , (long) dextendloop.ad );
//printf("ai %zu: ->ad.s: %ld\n", ai, (long) dextendloop.ad->s );
//adfull->show();
//dextendloop.ad->show();

			/*
			for (int kk=0; kk<ds->goodarrayslist[ai].size(); kk++) {
								int ps = ds->goodarrayslist[ai][kk].row_symmetry_group().permsize();
			printfd("## array %zu: group size %d\n", kk, ps);
			}
			*/

//printf("## oa_depth_extend: ds->goodarrayslist[ai][0].n_columns %d, extensioncol %d\n", ds->goodarrayslist[ai][0].n_columns, extensioncol );

//ds->depthalglist[ai]=DEPTH_DIRECT; // HACK
			processDepth ( ds->goodarrayslist[ai], ds->depthalglist[ai], dextendloop, ds->dextendsubList[ai], extensioncol, verbose );

//		dextendloop.showprogress(1, extensioncol, 1);

		}
	}

	dextend.counter->showcounts ( *adfull );

	//dextend.arraywriter->closeafiles();
	delete dextend.arraywriter;
	delete dextend.counter;

	if ( ds!=0 ) {
		delete ds;
	}
	delete adfull;
	delete arraylist;

	cleanGlobalStatic();
	arraysymmetry::rowpermpool.reset();


	logstream ( SYSTEM ) << "#time end: "<< currenttime() << std::endl;
	if ( verbose )
		std::cout << "#time total: " << printfstring ( "%.1f", get_time_ms()-time0 ) << " [s]" << std::endl;

	return 0;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
