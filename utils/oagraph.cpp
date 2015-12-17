/** \file oagraph.cpp

 C++ program: oagraph

 oagraph: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "graphtools.h"






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
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.addUsage ( "Orthogonal Array: oagraph: testing platform" );
	opt.addUsage ( "Usage: oagraph [OPTIONS] " );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.processCommandArgs ( argc, argv );


	int randvalseed = opt.getIntValue ( 'r', 1 );

	srand ( randvalseed );

	int verbose = opt.getIntValue ( 'v', 1 );
	int ix = opt.getIntValue ( 'i', 2 );

	print_copyright();
	setloglevel ( NORMAL );

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}


	setloglevel ( SYSTEM );

	// select a design
	array_link al = exampleArray ( ix, 2 );
	arraydata_t arrayclass = arraylink2arraydata ( al );

	array_transformation_t trans ( arrayclass );
	trans.reset();

	trans.randomize();
	array_link alr = trans.apply ( al );
	alr=al;

	if (0) {
	std::pair<array_link, std::vector<int> > Gc = array2graph ( alr,  verbose );
	std::vector<int> tr = nauty::reduceNauty ( Gc.first, Gc.second );
	printf ( "canon: " ); display_vector ( tr ); printf ( "\n" );
	array_transformation_t ttm = oagraph2transformation ( tr, arrayclass, verbose );
}
	
	array_transformation_t ttm = reduceOAnauty ( alr, 1 );
	ttm.show();

	array_link alm = ttm.apply ( alr );

	al.showarraycompact();
	printf ( "---\n" );
	alr.showarraycompact();
	printf ( "---\n" );
	alm.showarraycompact();



//return 0;

	array_link alx = al.randomperm();
	alx=alr;
	array_transformation_t ttx = reduceOAnauty ( alx, 0 );
	array_link alxm = ttx.apply ( alx );

	printf ( "--- alxm \n" );
	alxm.showarraycompact();

	if ( alxm!=alm ) {
		printf ( "error with nauty reduction...\n" );
	}
	// TODO: extract arraytransformation_t from canon
	// TODO: print minimal arrays and make unit test using 2 random transformations
	// TODO: clean up interfaces



	/*
	// convert design (with isomorphism type) to a colored graph
	nauty::NyGraph G = array2graph(al, 1);

	if (verbose)
	printf("creatd graph with %d vertices and ? edges\n", G.num_nodes );

	*/
	if ( verbose )
		printf ( "done\n" );


	return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
