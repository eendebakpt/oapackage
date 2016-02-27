/** \file oaconference.cpp

 C++ program: oaconference

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
#include "extend.h"

#include "oadevelop.h"
#include "lmc.h"

#include "conference.h"


int main ( int argc, char* argv[] )
{
	AnyOption opt;
	/* parse command line options */
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setOption ( "output", 'o' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "ii", 'i' );
	opt.setOption ( "rows", 'N' );
	opt.setOption ( "cols" );
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.setOption ( "input", 'i' );
	opt.setOption ( "output", 'o' );

	opt.addUsage ( "Orthonal Array: oaconference: testing platform" );
	opt.addUsage ( "Usage: oaconference [OPTIONS] [FILE]" );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.addUsage ( " -N --rows  			Number of rows " );
	opt.addUsage ( " -v --verbose [INT] 			Verbosity level " );
	opt.addUsage ( " -i --input [FILENAME]		Input file to use" );
	opt.processCommandArgs ( argc, argv );


	print_copyright();
	//cout << system_uname();
	setloglevel ( NORMAL );

	int verbose = opt.getIntValue ( 'v', 1 );
	int N = opt.getIntValue ( 'N', 10 );

	const std::string output = opt.getStringValue ( 'o', "" );
	const std::string input = opt.getStringValue ( 'i', "" );


	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	setloglevel ( SYSTEM );

	int k=3;
	conference_t ctype ( N, k );

	if ( verbose )
		printf ( "create conference matrices of size %dx%d\n", ctype.N, ctype.ncols );

	arraylist_t kk;
	if ( input.length() >1 ) {
		kk = readarrayfile ( input.c_str() )	 ;
		//al=kk[0];
	}
else {
	array_link al = ctype.create_root();
kk.push_back(al);	
}
	//printf ( "initial array:\n" );
	//al.showarray();

	arraylist_t outlist = extend_conference ( kk, ctype,  verbose );

	if (0) {
	array_link al=kk[0];
	int extcol=al.n_columns;
	conference_extend_t ce = extend_conference_matrix ( al, ctype, extcol, verbose );
	outlist = ce.getarrays ( al ) ;
	}
	
	if ( output.length() >1 ) {
		printf ( "oaconference: write %d arrays to file %s...\n", ( int ) outlist.size(), output.c_str() );
		writearrayfile ( output.c_str(),outlist );
	}
	printf ( "done...\n" );

	return 0;

}


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
