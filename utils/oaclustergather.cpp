/** \file oaclustergather.cpp

 C++ program: oaclustergather

 Dedicated tool to gather data on cluster

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2013

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <numeric>

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"
#include "arrayproperties.h"

#include "oadevelop.h"
#include "pareto.h"

std::string splitDir ( std::vector<int> ii )
{
	std::string s ;
	for ( size_t i=0; i<ii.size(); i++ ) {
		s+= printfstring ( "sp%d-split-%d/", i, ii[i] );
	}
	return s;
}
std::string splitFile ( std::vector<int> ii )
{
	std::string s ;
	for ( size_t i=0; i<ii.size(); i++ ) {
		s+= printfstring ( "sp%d-split-%d", i, ii[i] );
		if ( i<ii.size()-1 )
			s+="-";
	}
	return s;
}
std::string splitTag ( std::vector<int> ii )
{
	std::string s ;
	if (ii.size()==0){
			s+= "base";
			return s;
	}
	for ( size_t i=0; i<ii.size(); i++ ) {
		s+= printfstring ( "%d", ii[i] );
		if ( i<ii.size()-1 )
			s+=".";
	}
	return s;
}

inline std::vector<int> tovec ( std::vector<int> ii, int i )
{
	ii.push_back ( i );
	return ii;
}

inline std::vector<int> tovec ( int i )
{
	std::vector<int> ii;
	ii.push_back ( i );
	return ii;
}
inline std::vector<int> tovec ( int i, int j )
{
	std::vector<int> ii;
	ii.push_back ( i );
	ii.push_back ( j );
	return ii;
}
inline std::vector<int> tovec ( int i, int j, int k )
{
	std::vector<int> ii;
	ii.push_back ( i );
	ii.push_back ( j );
	ii.push_back ( k );
	return ii;
}

bool readNumbersFile ( const char *numbersfile, std::vector<long> &na, std::vector<long> &npareto, int kmin=-1, int kmax= -1 )
{
	FILE *fid = fopen ( numbersfile, "rt" );
	//printf("read numbers file %s\n", numbersfile);
	if ( fid==0 ) {
		//printf ( "error: read numbers file %s: error\n" );
		return 0;
	}
	long np, n, k;
	//int k;

	while ( 1 ) {
//		for ( int k=kmin; k<=kmax; k++ ) {
		int r = fscanf ( fid, "k %ld: %ld %ld\n", &k, &n, &np );
		//printf("read numbers file: %s: r %d: %ld %ld %ld\n", numbersfile, r, k, n, np);
		if ( r<3 )
			break;
		if ( k>= ( int ) na.size() )
			na.resize ( k+1 );
		if ( k>= ( int ) npareto.size() )
			npareto.resize ( k+1 );
		na[k]=n;
		npareto[k]=np;
	}

	fclose ( fid );

	return 1;
}

void writeNumbersFile ( const char *numbersfile, std::vector<long> na, std::vector<long> npareto, int kmin=-1, int kmax= -1 )
{

	if ( kmin<0 )
		kmin=0;
	if ( kmax<0 )
		kmax=na.size() +1;

	FILE *fid = fopen ( numbersfile, "wt" );
	for ( int k=kmin; k<=kmax; k++ ) {
		fprintf ( fid, "k %d: %ld %ld\n", k, na[k], npareto[k] );
	}

	fclose ( fid );
}

const std::string filesep = "/";

/**
 * @brief Analyse arrays using Pareto optimality
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] )
{

	/* parse command line options */
	AnyOption opt;
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setFlag ( "coptions" );	/* print compile time options */
	opt.setOption ( "output", 'o' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "basedir", 'b' );
	opt.setOption ( "split0" );
	opt.setOption ( "split1" );
	opt.setOption ( "split2" );
	opt.setOption ( "numbersfile" );
	opt.setOption ( "config", 'c' );
	opt.setOption ( "nsplit0" );
	opt.setOption ( "nsplit1" );
	opt.setOption ( "nsplit2" );
	opt.setOption ( "nparetodiff" );	// allow difference in pareto files for .oa files and numbers file
	opt.setOption ( "kmin" );
	opt.setOption ( "kmax" );

	opt.addUsage ( "Orthonal Array Cluster Gather: special tool" );
	opt.addUsage ( "Usage: oaclustergather [OPTIONS] " );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.addUsage ( " -v --verbose  			Verbose level (default: 1) " );
	opt.addUsage ( " -b --basedir  			Base calculation dir " );
	opt.addUsage ( " --numbersfile [FILENAME]  Output name of number of arrays " );
	opt.addUsage ( " --nsplit0 [NUMBER]		Number of split files at level 0" );
	opt.addUsage ( " --nsplit1 [NUMBER]		Number of split files at level 1" );
	opt.addUsage ( " --nsplit2 [NUMBER]		Number of split files at level 2" );
	opt.addUsage ( " --kmin [NUMBER] --kmax [NUMBER]	Min and max number of columns (inclusive) to process" );
	opt.addUsage ( " -o [FILE] --output [FILE]	Output prefix for filtered arrays (default: no output) " );
	opt.processCommandArgs ( argc, argv );

	print_copyright();

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 || ( opt.getValue ( "coptions" ) != 0 ) ) {
		opt.printUsage();
		if ( opt.getValue ( "coptions" ) != 0 ) {
			print_options ( cout );
		}
		exit ( 0 );
	}

	const char *outputprefix = 0;
	if ( opt.getValue ( "output" ) !=0 )
		outputprefix = opt.getValue ( 'o' );

	const std::string basedir = opt.getStringValue ( "basedir", "" );
	const char *configfile = opt.getStringValue ( "config", "oaconfig.txt" );
	int verbose = opt.getIntValue ( 'v', 1 );
	int split0 = opt.getIntValue ( "split0", -1 );
	int split1 = opt.getIntValue ( "split1", -1 );
	int nsplit0 = opt.getIntValue ( "nsplit0", -1 );
	int nsplit1 = opt.getIntValue ( "nsplit1", -1 );
	int allowparetodiff = opt.getIntValue ( "nparetodiff", 0 );
	int nsplit2 = opt.getIntValue ( "nsplit2", -1 );
	int kmin = opt.getIntValue ( "kmin", 9 );
	int kmax = opt.getIntValue ( "kmax", 24 );
	const char *numbersfile = opt.getStringValue ( "numbersfile", 0 );

	arraydata_t *adata = readConfigFile ( configfile );
	assert ( adata!=0 );

	if ( verbose ) {
		printf ( "oaclustergather: basedir %s, split0 %d, nsplit1 %d, kmin %d, kmax %d, verbose %d\n", basedir.c_str(), split0, nsplit1, kmin, kmax, verbose );
	}
	if ( verbose ) {
		std::cout << "#time start: "<< currenttime() << std::endl;
	}
	double time0=get_time_ms();

	/* Loop over all subdirectories for all number of columns */

	std::vector<long> na ( kmax+1 );	  /// total number of arrays
	std::vector<long> npareto ( kmax+1 ); /// total number of pareto arrays

	int cleanrun=1; /// indicates whether all necessary files have been found


	std::vector<int> nsplit;
	nsplit.push_back ( nsplit0 );
	nsplit.push_back ( nsplit1 );
	nsplit.push_back ( nsplit2 );

	std::vector<int> lvls;
	if ( split0>=0 ) {
		lvls.push_back ( split0 );
		if ( split1>=0 ) {
			lvls.push_back ( split1 );
		}
	}
	if ( verbose ) {
		std::cout << "oaclustergather: levels "<< splitTag ( lvls ) << std::endl;
	}

	//printf("xxx nsplit0 %d\n", nsplit0 );

	// loop over all columns
	for ( int k=kmin; k<=kmax; k++ ) {

		if ( verbose>=2 )
			printf ( "## oaclustergather: %d columns (time %.1f [s])\n", k, get_time_ms()-time0 );
		Pareto<mvalue_t<long>,array_link> pset;

		arraydata_t adata0 ( adata, k );

		if ( nsplit0==-1 || 1 ) {
			int level = lvls.size();

			//printf("------------------- nsplit[level] %d\n", nsplit[level]);
			std::string splittag = splitTag ( lvls );

			if ( verbose )
				printf ( "\n## oaclustergather: %d columns, gathering results for stage %d: split %s (time %.1f [s])\n", k, level, splittag.c_str(), get_time_ms()-time0 );
			for ( int jj=0; jj<nsplit[level]; jj++ ) {
				std::string subdir = splitDir ( tovec ( lvls, jj ) );

				std::string  nfilesub0= "numbers-" + splitTag ( tovec ( lvls, jj ) ) + ".txt";
				std::string  nfilesub= basedir + filesep + subdir + filesep + nfilesub0;
				std::vector<long> nasub ( kmax+1 );
				std::vector<long> nparetosub ( kmax+1 );
				bool b = readNumbersFile ( nfilesub.c_str(), nasub, nparetosub, kmin, kmax );

				if (verbose>=2)
				{
					printf("   --> read numbers file %s\n", nfilesub.c_str() );
				}
				if ( jj==0 ) {
					//printf("numbers: nasub[%d]=%d\n", k, nasub[k]);
				}

				std::string subfile0 = splitFile ( tovec ( lvls, jj ) ) + printfstring ( "-extend-%s.oa", adata0.idstr().c_str() );
				std::string subfilepareto0 = splitFile ( tovec ( lvls, jj ) ) + printfstring ( "-pareto-%s.oa", adata0.idstr().c_str() );

				std::string afile = basedir + filesep + subdir + filesep + subfile0;

				bool paretofile = 0;
				if (verbose>=2)
					printf ( "  check %s\n", subfilepareto0.c_str() );
				if ( oa_file_exists ( basedir + filesep + subdir + subfilepareto0 ) ) {
					if (verbose>=2) {
					printf ( "  switching to Pareto file %s->%s!\n", subfile0.c_str(), subfilepareto0.c_str() );
					}
					afile =  basedir + filesep + subdir + subfilepareto0;
					paretofile = 1;
				}
				int nn = nArrays ( afile.c_str() );
				// printf ( "file %s: %d arrays\n", afile.c_str(), nn );

				if ( verbose>=2 )
					printf ( "### file %s: %d arrays\n", base_name ( afile ).c_str(), nn );

				if ( nn>=0 )
					if ( paretofile ) {
						if ( verbose>=2 ) {
							printf ( "  k %d: adding %ld arrays\n", k, nasub[k] );
						}
						na[k]+=nasub[k];
						// check
						if ( nparetosub[k]!=nn )  {
							printfd ( "\n### warning: column %d: nparetosub[%d] %d, nn %d (number of arrays in .oa file)\n\n", jj, k, nparetosub[k], nn );
							if (! allowparetodiff ) {
							cleanrun=0;
							}
						}
					} else {
						na[k] += nn;
					}
				else {
					if ( cleanrun || verbose>=2 ) {	// only report error if we are in a clean run
						fprintf ( stderr, "file %s: error (nn %d)\n", subfile0.c_str(), nn );
						fprintf ( stderr, "   error: file %s\n", afile.c_str() );
					}
					cleanrun=0;
					exit ( 0 );
				}

				const arraylist_t arraylist = readarrayfile ( afile.c_str(), 0 );

				#pragma omp parallel for
				for ( int i=0; i< ( int ) arraylist.size(); i++ ) {
					if ( verbose>=3 || ( ( i%5000==0 ) && verbose>=2 ) ) {
						printf ( "oaclustergather: file %d, array %d/%ld\n", jj, i, arraylist.size() );
						printf ( "  " );
						pset.show ( 1 );
//#ifdef OPENMP
						//printf(" openmp: %d\n",  omp_get_num_threads() );
//#endif

					}

					const array_link &al = arraylist.at ( i );

					//parseArrayPareto ( al, al, pset, verbose );
					Pareto<mvalue_t<long>, array_link >::pValue p = calculateArrayPareto<array_link> ( al, verbose>=3 );
					//typename

					#pragma omp critical
					{
						// add the new tuple to the Pareto set
						pset.addvalue ( p, al );
					}
				}

				if ( verbose>=2 || ( ( jj%20==0 || ( jj==nsplit1-1 ) ) && verbose>=1 ) ) {
					printf ( "oaclustergather: file %d/%d, %ld arrays: %d Pareto values, %d Pareto elements\n", jj, nsplit1, arraylist.size(), pset.number(), pset.numberindices() );
					//printf ( "  " ); pset.show ( 1 );
				}
			}
		} else {
			printf("ERROR: code should not be used any more...\n");
			exit(1);			
		}

		if ( verbose ) {
			printf ( "final pareto set %d cols: ", k );
			pset.show ( 2 );
		}

		// write pareto set to disk

		arraylist_t pp = pset.allindicesdeque();
		npareto[k] =  pset.numberindices();

		if ( cleanrun ) {
			std::string cdir =  splitDir ( lvls ); //printfstring ( "sp0-split-%d/", split0 );
			std::string xxx =  splitFile ( lvls );
			if ( lvls.size() >0 ) {
			} else {
				xxx =  "results-j5evenodd";
			}
			std::string pfile = printfstring ( "%s-pareto-%s.oa", xxx.c_str(), adata0.idstr().c_str() );
			std::string afile = basedir + "/" + cdir + pfile;

			if ( verbose )
				printf ( "  writing pareto file %s\n", pfile.c_str() );
			writearrayfile ( afile.c_str(), &pp, ABINARY, adata->N, k );
		}
	}

	if ( verbose )
		printf ( "totals (cleanrun %d):\n", cleanrun );
	for ( int k=kmin; k<=kmax; k++ ) {
		if ( verbose ) {
			printf ( "%d columns: %ld arrays, %ld pareto\n", k, na[k], npareto[k] );
		}
	}


	if ( numbersfile && cleanrun ) {
		if ( verbose )
			printf ( "writing numbers file %s\n", numbersfile );
		writeNumbersFile ( numbersfile, na, npareto, kmin, kmax );
	}

	/* free allocated structures */
	delete adata;

	if ( verbose ) {
		std::cout << "#time end: "<< currenttime() << std::endl;
		std::cout << "#time total: " << printfstring ( "%.1f", get_time_ms()-time0 ) << " [s]" << std::endl;
	}
	return 0;
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
