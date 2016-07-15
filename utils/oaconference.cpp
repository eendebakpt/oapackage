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


#include "conference.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <string>
using namespace std;

void print_all_permutations ( const string& s )
{
	string s1 = s;
	sort ( s1.begin(), s1.end() );
	int n=0;
	do {
		cout << s1 << endl;
		n++;
	} while ( next_permutation ( s1.begin(), s1.end() ) );
	printf ( "%d/%ld perms (len %ld)\n", n, factorial<long> ( s1.size() ), ( long ) s1.size() );
}

void testx()
{
	print_all_permutations ( "001111112222222" );
	exit ( 0 );
}

int main ( int argc, char* argv[] )
{
	AnyOption opt;
	/* parse command line options */
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setOption ( "output", 'o' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "ii", 'i' );
	opt.setOption ( "kmax", 'k' );
	opt.setOption ( "itype" );
	opt.setOption ( "j1zero" );
	opt.setOption ( "j3zero" );
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
	opt.addUsage ( " -k --kmax  			Max number of columns " );
	opt.addUsage ( " -v --verbose [INT] 			Verbosity level " );
	opt.addUsage ( " -i --input [FILENAME]		Input file to use" );
	opt.addUsage ( " -o --output [FILEBASE]		Output file to use" );
	opt.addUsage ( printfstring ( " --ctype [TYPE]				Zero for normal type, %d for diagonal, %d for double conference type ",conference_t::CONFERENCE_DIAGONAL, conference_t::DCONFERENCE ).c_str() );
	opt.addUsage ( printfstring ( " --itype [TYPE]				Matrix isomorphism type (CONFERENCE_ISOMORPHISM %d)", CONFERENCE_ISOMORPHISM ).c_str() );
	opt.addUsage ( " --j3zero [INT]		Restrict designs to J3=0" );
	opt.processCommandArgs ( argc, argv );


	//testx();

	print_copyright();
	//cout << system_uname();
	setloglevel ( NORMAL );

	int verbose = opt.getIntValue ( 'v', 1 );
	int N = opt.getIntValue ( 'N', 10 );
	int kmax = opt.getIntValue ( 'k', N );
	int select = opt.getIntValue ( 's', 1 );
	const conference_t::conference_type ctx = ( conference_t::conference_type ) opt.getIntValue ( "ctype", 0 );
	const matrix_isomorphism_t itype = ( matrix_isomorphism_t ) opt.getIntValue ( "itype", CONFERENCE_ISOMORPHISM );
	const int j1zero = opt.getIntValue ( "j1zero", 0 );
	const int j3zero = opt.getIntValue ( "j3zero", 0 );

	const std::string output = opt.getStringValue ( 'o', "" );
	const std::string input = opt.getStringValue ( 'i', "" );


	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}

	setloglevel ( SYSTEM );

	int kstart=-1;
	conference_t ctype ( N, N );

	arraylist_t inputarrays;
	if ( input.length() >1 ) {
		inputarrays = readarrayfile ( input.c_str() )	 ;
		//al=kk[0];

		if ( inputarrays.size() >0 ) {
			ctype = conference_t ( N, kstart );
		}
		ctype.ctype=ctx;
		ctype.itype=itype;
		ctype.j3zero=j3zero;
		ctype.j1zero=j1zero;
		kstart=inputarrays[0].n_columns;
		kmax=inputarrays[0].n_columns+1;
	} else {
		ctype.ctype=ctx;
		ctype.itype=itype;
		ctype.j3zero=j3zero;
		ctype.j1zero=j1zero;

		ctype.addRootArrays ( inputarrays );
		kstart=inputarrays[0].n_columns;
	}

	if(ctype.ctype==conference_t::DCONFERENCE) {
		kmax=std::min( (int)(ceil(N/2)+2), kmax);
	}
	
	if ( verbose ) {
		printf ( "oaconference: extend %d conference matrices to size %dx%d (itype %d (CONFERENCE_ISOMORPHISM %d, CONFERENCE_RESTRICTED_ISOMORPHISM %d) )\n", ( int ) inputarrays.size(), ctype.N, ctype.ncols, itype, CONFERENCE_ISOMORPHISM, CONFERENCE_RESTRICTED_ISOMORPHISM );
		printf ( "oaconference: ctype %d (DCONFERENCE %d, CONFERENCE_NORMAL %d) )\n", ctype.ctype, conference_t::DCONFERENCE, conference_t::CONFERENCE_NORMAL );
	}
	if ( verbose>=3 ) {
		printf ( "--- initial set of arrays ---\n" );
		showArrayList ( inputarrays );
	}

	double t0 = get_time_ms();
	
	for ( int extcol=kstart; extcol<kmax; extcol++ ) {
		printf ( "oaconference: extend column %d (max number of columns %d, time %.1f [s])\n", extcol, kmax, get_time_ms(t0) );

		arraylist_t outlist;

		switch ( ctype.ctype ) {
		case conference_t::DCONFERENCE: {
			outlist = extend_double_conference ( inputarrays, ctype,  verbose );
			sort ( outlist.begin(), outlist.end(), compareLMC0 );
			break;
		}
		case conference_t::CONFERENCE_NORMAL:
		case conference_t::CONFERENCE_DIAGONAL:
			switch ( ctype.itype ) {
			case CONFERENCE_RESTRICTED_ISOMORPHISM:
				outlist = extend_conference_restricted ( inputarrays, ctype,  verbose );
				break;

			case CONFERENCE_ISOMORPHISM:
				outlist = extend_conference ( inputarrays, ctype,  verbose, select );
				break;
			default
					:
				printfd ( "isomorphism type not implemented" );
				break;
			}
			break;
		default
				:
			printfd ( "not implemented: itype %d\n", ctype.itype );
			break;
		}

		if ( select ) {
			outlist = selectConferenceIsomorpismClasses ( outlist, verbose, ctype.itype );
		}
		sort ( outlist.begin(), outlist.end(), compareLMC0 );

		if ( ctype.ctype==conference_t::DCONFERENCE ) {
			for ( size_t i=0; i<outlist.size(); i++ ) {
				if ( ! isConferenceFoldover ( outlist[i] ) )
					printfd ( "#### found an even-odd conference matrix!!!\n" );
			}
		}

		if ( output.length() >=1 ) {
			std::string outfile = output + printfstring ( "-%d-%d", ctype.N, extcol+1 )  + ".oa";
			printf ( "oaconference: write %d arrays to file %s...\n", ( int ) outlist.size(), outfile.c_str() );

			if ( outlist.size() < 20000 )
				writearrayfile ( outfile.c_str(),outlist );
			else {
				//writearrayfile ( outfile.c_str(),outlist, ABINARY);
				writearrayfile ( outfile.c_str(),outlist, ABINARY_DIFFZERO );
			}
		}

		printf ( "oaconference: extend column %d: generated %d non-isomorphic arrays (%.1f [s])\n", extcol, ( int ) outlist.size(),  get_time_ms(t0) );

		// loop
		inputarrays=outlist;
	}
	printf ( "done...\n" );

	return 0;

}


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
