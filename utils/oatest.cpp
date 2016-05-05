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


	int randvalseed = opt.getIntValue ( 'r', 1 );
	int ix = opt.getIntValue ( 'i', 11 );
	int r = opt.getIntValue ( 'r', 0 );
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

	{

		const char *numbersfile = "numbers-J.txt";
		Jcounter jc = readStatisticsFile ( numbersfile, verbose );
		jc.show();
		exit ( 0 );

//		long imax = std::numeric_limits<long>::max(); printf("max for long %ld\n" , imax);
//		 imax = std::numeric_limits<int>::max();
//	printf("max for int %ld\n" , imax);
//	exit(0);

		arraylist_t lst= readarrayfile ( input );
		printf ( "read %d arrays\n", ( int ) lst.size() );
		std::sort ( lst.begin(), lst.end() );



		Jcounter jcounter = calculateJstatistics ( input, jj, 1 );
		Jcounter jcounter2 = calculateJstatistics ( input, jj, 1 );
		jcounter += jcounter2;

		printf ( "--- results ---\n" );
		jcounter.show();
		jcounter.showPerformance();

		writeStatisticsFile ( "numbers-J.txt", jcounter, 1 );


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
