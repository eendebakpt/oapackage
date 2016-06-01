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
#include <map>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "evenodd.h"
#include "oadevelop.h"
#include "lmc.h"

#include "conference.h"


#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>

#include "Deff.h"
#include "arraytools.h"
#include "strength.h"

using namespace Eigen;


//std::vector<double> Defficiencies (const array_link &al, const arraydata_t & arrayclass, int verbose ) ;

array_link finalcheck ( const array_link &al,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niterx, int nabortx )
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



#include "graphtools.h"


template < class Type>
/// return the condition number of a matrix
double conditionNumber ( const array_link M )
{
	MatrixXd A = arraylink2eigen ( M );
	JacobiSVD<Matrix<Type,-1,-1> > svd ( A );
	double cond = svd.singularValues() ( 0 ) / svd.singularValues() ( svd.singularValues().size()-1 );
	return cond;
}






/// convert 2-level design to second order interaction matrix
inline void array2eigenxf ( const array_link &al, Eigen::MatrixXd &mymatrix )
{
	int k = al.n_columns;
	int n = al.n_rows;
	int m = 1 + k + k* ( k-1 ) /2;

	mymatrix = Eigen::MatrixXd::Zero ( n,m );

	// init first column
	int ww=0;
	for ( int r=0; r<n; ++r ) {
		mymatrix ( r, ww ) = 1;
	}

	// init array
	ww=1;
	for ( int c=0; c<k; ++c ) {
		int ci = c*n;
		for ( int r=0; r<n; ++r ) {
			mymatrix ( r, ww+c ) = al.array[r+ci];
		}
	}

	// init interactions
	ww=k+1;
	for ( int c=0; c<k; ++c ) {
		int ci = c*n;
		for ( int c2=0; c2<c; ++c2 ) {
			int ci2 = c2*n;

			const array_t * p1 = al.array+ci;
			const array_t * p2 = al.array+ci2;
			for ( int r=0; r<n; ++r ) {
				mymatrix ( r, ww ) = ( *p1+*p2 ) %2;
				p1++;
				p2++;
			}
			ww++;
		}
	}

	mymatrix.array() *= 2;
	mymatrix.array() -= 1;
}



int main ( int argc, char* argv[] )
{
	AnyOption opt;
	/* parse command line options */
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setOption ( "output", 'o' );
	opt.setOption ( "input", 'I' );
	opt.setOption ( "rand", 'r' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "ii", 'i' );
	opt.setOption ( "jj" );
	opt.setOption ( "xx", 'x' );
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


	double t0=get_time_ms(), dt=0;
	int randvalseed = opt.getIntValue ( 'r', 1 );
	int ix = opt.getIntValue ( 'i', 1 );
	int r = opt.getIntValue ( 'r', 20 );
	int jj = opt.getIntValue ( "jj", 5 );

	int xx = opt.getIntValue ( 'x', 0 );
	int niter  = opt.getIntValue ( "niter", 10 );
	int verbose  = opt.getIntValue ( "verbose", 1 );

	char *input = opt.getValue ( 'I' );
	if ( input==0 )
		input="test.oa";

	srand ( randvalseed );
	if ( randvalseed==-1 ) {
		randvalseed=time ( NULL );
		printf ( "random seed %d\n", randvalseed );
		srand ( randvalseed );
	}


	print_copyright();
	//cout << system_uname();
	setloglevel ( NORMAL );

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	setloglevel ( SYSTEM );


	if ( 1 ) {

// FIXME: make unittest of this one
		//arraylist_t ll = readarrayfile ( "/home/eendebakpt/oatmp/conf/dconferencej1j3-36-8.oa" );
		//arraylist_t ll = readarrayfile ( "/home/eendebakpt/oatmp/conf/dconferencej1j3-28-4.oa" );
		arraylist_t ll = readarrayfile ( "/home/eendebakpt/oatmp/conf/dconferencej1j3-32-8.oa" );
		
		//ll[52]=ll[51]; 
		//ll[51]=ll[52]; 
		
		array_link alx = ll[0];
		alx.showproperties();
		int N = alx.n_rows;

		int fi=51;
		printf("first:\n"); ll[fi].showarray(); printf("next :\n"); ll[fi+1].showarray(); printf("  first diff index: %d\n", ll[fi].firstColumnDifference(ll[fi+1]) );
		

		conference_t ct(N, 2*N);
		CandidateGenerator cgenerator(array_link() , ct);
		cgenerator.verbose=1;
		cgenerator.show();

double t0;
if (0) {
t0=get_time_ms();		

				for(size_t i=0; i<ll.size(); i++) {
			const array_link &al = ll[i];
			std::vector<cperm> cc1 = generateDoubleConferenceExtensionsInflate (al, ct, 1, 1, 1);		
		}
				printf ( "   dtt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );
}

		t0=get_time_ms();
		std::vector<cperm> cc1;
		for(size_t i=0; i<ll.size(); i++) {
			printf("--- i %d -------\n", (int) i);
			const array_link &al = ll[i];
			int nc1=-1;
			if (1) {
				printfd("   _________________________\n");
			 cc1 = generateDoubleConferenceExtensionsInflate (al, ct, 1, 1, 1);
			 nc1=cc1.size();
				printfd("   _________________________\n");
			}
			//cperm tmp= cc1[0]; printf("size cc1: %d\n", (int)tmp.size() );
			
			cgenerator.verbose=2;
			std::vector<cperm> cc2 = cgenerator.generateCandidates(al);
			//printf("size cc2: %d\n", cc2[0].size());
			
			
			
			int nc2=cc2.size();
			
		printf("%d: number of candidates: %d/%d\n", (int)i, nc1, nc2);	
		cgenerator.show();
		
		if(i==52 && 0) {
			printf("cc1: \n");
			showCandidates(cc1);
			printf("cc2: \n");
			showCandidates(cc2);
			break;
		}
		}
				printf ( "   dtt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );

		exit ( 0 );
	}

	{
		int filterip=1;
		int filtersymm=1;
		int filterj3 = 1;
		array_link al = exampleArray ( r, 1 );
		al.show();

		al=al.selectFirstColumns ( xx );
		int kstart=xx-1;
		//al=al.selectFirstColumns(4);		int kstart=3;
		array_link als = al.selectFirstColumns ( kstart );

		//als.show();

		printf ( "find extensions of:\n" );
		al.showarray();
		al.row_symmetry_group().show();

		int N = al.n_rows;

		conference_t				ct = conference_t ( N, N );
		ct.j1zero=1;
		ct.j3zero=1;
		ct.itype=CONFERENCE_RESTRICTED_ISOMORPHISM;
		ct.ctype=conference_t::DCONFERENCE;


		//for(int z=0; z<10; z++)
		std::vector<cperm> cci = generateDoubleConferenceExtensionsInflate ( al, ct, verbose, 1, 1 );


		//printf ( "no symm:\n" );
		//cc = generateDoubleConferenceExtensions ( als, ct, verbose, 0, filterip, filterj3, 0 );

		if ( ix )
			exit ( 0 );
		printf ( "## full array (with symm):\n" );

		t0=get_time_ms();
		std::vector<cperm> cc3 = generateDoubleConferenceExtensions ( al, ct, verbose, 0, filterip, filterj3, 1 );
		for ( size_t i=0; i<cc3.size(); i++ ) {
			printf ( "  %d: ", ( int ) i );
			print_cperm ( cc3[i] );
			printf ( "\n" );
		}
		printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );
		exit ( 0 );

		printf ( "## full array (no symm):\n" );
		t0=get_time_ms();
		cc3 = generateDoubleConferenceExtensions ( al, ct, verbose, 0, filterip, filterj3, 0 );
		printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );


		//	printf ( "extend_conference: extended array %d/%d to %d arrays\n", ( int ) i, ( int ) lst.size(), nn );
		exit ( 0 );
	}



	{


		//		long imax = std::numeric_limits<long>::max(); printf("max for long %ld\n" , imax);
		//		 imax = std::numeric_limits<int>::max();
		//	printf("max for int %ld\n" , imax);
		//	exit(0);

		arraylist_t lst= readarrayfile ( input );
		printf ( "read %d arrays\n", ( int ) lst.size() );
		std::sort ( lst.begin(), lst.end() );



		Jcounter jcounter = calculateJstatistics ( input, jj, 1 );
		//Jcounter jcounter2 = calculateJstatistics ( input, jj, 1 ); jcounter += jcounter2;

		printf ( "--- results ---\n" );
		jcounter.show();
		jcounter.showPerformance();

		writeStatisticsFile ( "numbers-J.txt", jcounter, 1 );

		const char *numbersfile = "numbers-J.txt";
		Jcounter jc = readStatisticsFile ( numbersfile, verbose );
		jc.show();
		exit ( 0 );


		int jj=0;
		if ( xx ) {
			Pareto<mvalue_t<long>,array_link> pset;
			for ( int i=0; i<niter; i++ )
				addArraysToPareto ( pset, calculateArrayParetoJ5Cache<array_link>, lst, jj, verbose );
		}

		exit ( 0 );
	}




	/*

	if ( 1 ) {

		arraydata_t ad;

		conference_t ctype ( 8, 3 );
		ctype.itype=CONFERENCE_RESTRICTED_ISOMORPHISM;
		ctype.ctype=conference_t::DCONFERENCE;
		arraylist_t lst= readarrayfile ( "test.oa" );
		int verbose=1;

		for ( int i=0; i< ( int ) lst.size() ; i++ ) {
			array_transformation_t t = reduceOAnauty ( lst[i]+1, 2 );

			array_link A = t.apply ( lst[i]+1 ) + ( -1 );
			printf ( "array %d\n", i );
			lst[i].showarray();
			printf ( "array %d reduced\n", i );
			A.showarray();
		}

		arraylist_t lst2 = addConstant ( lst, 0 );

		arraylist_t outlist = selectConferenceIsomorpismClasses ( lst2, verbose, ctype.itype );
		outlist = addConstant ( outlist, 0 );
		writearrayfile ( "test2.oa", outlist );
		exit ( 0 );
	}
	*/

	if ( 1 ) {

		array_link al = exampleArray ( 5 );
		array_link alx = al;
		alx.randomperm();

		array_transformation_t t1 = reduceOAnauty ( al, 1 );
		//t1.show();	return 0;

		array_link alr1 = t1.apply ( al );


		array_transformation_t t2 = reduceOAnauty ( alx, 1 );
		array_link alr2 = t2.apply ( alx );


		printf ( "reduced:\n" );
		alr1.showarray();
		printf ( "random reduced:\n" );
		alr2.showarray();

		if ( alr1 != alr2 )
			printf ( "error: reductions unequal!\n" );

		return 0;
	}
	if ( 0 ) {

		arraylist_t ll = readarrayfile ( "/home/eendebakpt/tmp/sp0-split-10/sp0-split-10-pareto-64.2-2-2-2-2-2-2-2-2-2.oa" );
		array_link al=ll[0];

		int r0= ( al ).rank();
		int r=array2xf ( al ).rank();
		printf ( "rank: %d %d\n",  r0,r );

		arraydata_t arrayclass=arraylink2arraydata ( al, 1 );

		OAextend oaextend=OAextend();
		oaextend.checkarrays=0;
		oaextend.setAlgorithm ( MODE_J5ORDERXFAST, &arrayclass );
		setloglevel ( NORMAL );

		arrayclass.show();
		oaextend.info();

		printf ( "extend!\n" );
		al.show();

		int current_col=al.n_columns;
		arraylist_t extensions;
		int nr_extensions = extend_array ( al.array, &arrayclass, current_col,extensions, oaextend );

		//arraylist_t ww=extend_array(al, arrayclass, oaextend);

		return 0;
	}


	{
		arraylist_t lst = readarrayfile ( input );
		arraylist_t lstgood  = selectConferenceIsomorpismClasses ( lst, verbose );

		return 0;
	}


	return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
