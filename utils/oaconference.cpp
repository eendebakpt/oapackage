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
	opt.setOption ( "select", 's' );
	opt.setOption ( "cols" );
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.setOption ( "input", 'i' );
	opt.setOption ( "output", 'o' );
	opt.setOption ( "ctype" );

	opt.addUsage ( "Orthonal Array: oaconference: testing platform" );
	opt.addUsage ( "Usage: oaconference [OPTIONS] [FILE]" );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.addUsage ( " -N --rows  			Number of rows " );
	opt.addUsage ( " -v --verbose [INT] 			Verbosity level " );
	opt.addUsage ( " -i --input [FILENAME]		Input file to use" );
	opt.addUsage ( " -o --output [FILEBASE]		Output file to use" );
	opt.addUsage ( " --ctype [TYPE]				Zero for normal type" );
	opt.processCommandArgs ( argc, argv );


	print_copyright();
	//cout << system_uname();
	setloglevel ( NORMAL );

	int verbose = opt.getIntValue ( 'v', 1 );
	int N = opt.getIntValue ( 'N', 10 );
	int select = opt.getIntValue ( 's', 1 );
	const conference_t::conference_type ctx = (conference_t::conference_type)opt.getIntValue ( "ctype", 0 );

	const std::string output = opt.getStringValue ( 'o', "" );
	const std::string input = opt.getStringValue ( 'i', "" );


	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	setloglevel ( SYSTEM );

	int kstart=-1;
	int kmax=N;
	conference_t ctype ( N, N );
	ctype.ctype=ctx;


	arraylist_t kk;
	if ( input.length() >1 ) {
		kk = readarrayfile ( input.c_str() )	 ;
		//al=kk[0];

		if ( kk.size() >0 ) {
			ctype = 	conference_t ( N, kstart );
		}
		kstart=kk[0].n_columns;
		kmax=kk[0].n_columns+1;
	} else {
		array_link al = ctype.create_root();
		kk.push_back ( al );
		kstart=2;
	}

	if ( verbose )
		printf ( "oaconference: extend %d conference matrices of size %dx%d\n", ( int ) kk.size(), ctype.N, ctype.ncols );

	for ( int extcol=kstart; extcol<kmax; extcol++ ) {
		printf ( "oaconference: extend column %d (max number of columns %d)\n", extcol, kmax);
		arraylist_t outlist = extend_conference ( kk, ctype,  verbose );

		if ( select )
			outlist = selectConferenceIsomorpismClasses ( outlist, verbose );

		sort ( outlist.begin(), outlist.end(), compareLMC0 );

		if ( output.length() >1 ) {
			std::string outfile = output + printfstring ( "-%d-%d", ctype.N, extcol+1 )  + ".oa";
			printf ( "oaconference: write %d arrays to file %s...\n", ( int ) outlist.size(), outfile.c_str() );
			writearrayfile ( outfile.c_str(),outlist );
		}

		//printf("generated columns: %d\n", outlist[0].n_columns);
				printf ( "oaconference: extend column %d: generated %d arrays\n", extcol, (int)outlist.size() );

		// loop
		kk=outlist;
	}
	printf ( "done...\n" );

	return 0;

}


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
