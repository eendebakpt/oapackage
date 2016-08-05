/** \file oajoin.cpp

 C++ program: oajoin

 oajoin can merge several files with arrays into a single file

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"

/**
 * @brief Read in files with arrays and join them into a single file
 * @param argc
 * @param argv[]
 * @return
 */
int main(int argc, char* argv[])
{
    AnyOption opt;

    /* parse command line options */
    opt.setFlag(  "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption("output", 'o');
    opt.setFlag("sort", 's');
    opt.setFlag("latex", 'l');
    opt.setOption("format", 'f');
    opt.setOption("verbose", 'v');

//    	printf("Orthogonal Arrays %s\n", version().c_str());

    opt.addUsage( "Orthonal Array Join: join several array files into a single file" );
    opt.addUsage( "Usage: oajoin [OPTIONS] [FILES]" );
    opt.addUsage( "" );
    opt.addUsage( " -h --help  			Prints this help " );
    opt.addUsage( " -s --sort			Sort the arrays " );
    opt.addUsage( " -l --latex			Output with LaTeX format " );
    opt.addUsage( " -f [FORMAT]					Output format (TEXT, BINARY (default), D (binary difference) ) " );
    opt.addUsage( " -o [FILE] --output [FILE]	Output prefix (default: standard output) " );
    opt.processCommandArgs(argc, argv);
    
    print_copyright();

    /* parse options */
    if( opt.getFlag( "help" ) || opt.getFlag( 'h' ) || opt.getArgc()==0 ) {
        opt.printUsage();
        exit(0);
    }

    const std::string outputprefix = opt.getStringValue('o', "");
    bool sortarrays = opt.getFlag('s');
    int verbose =  opt.getIntValue("verbose", 2);
    std::string format = opt.getStringValue('f', "BINARY");
    arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString(format);

    /* read in the arrays */
    arraylist_t *arraylist = new arraylist_t;
	if (verbose)
	  std::cout << "oajoin: reading " << opt.getArgc() << " file(s)" << endl;

	rowindex_t nrows; int ncols, nbits;
	
    for( int i = 0 ; i < opt.getArgc() ; i++ ) {
        if (verbose)
		  cout << "file " <<  opt.getArgv( i ) << endl ;
        int n = readarrayfile( opt.getArgv( i ), arraylist, 0, &ncols, &nrows, &nbits);
        if (verbose)
		  cout << "   read " << n << " array(s)" << endl;
    }

    /* perform operations on arrays */
    if (sortarrays) {
        if (verbose)
			std::cout << "Sorting arrays" << endl;
        sort(arraylist->begin(), arraylist->end());
    }

    /* write arrays to file */
    string outfile = "";
    if (outputprefix!="") {
        outfile += outputprefix;
        outfile += ".oa";
    }
    else {
        outfile += "<standard output>";
    }

    if (verbose)
	  std::cout << "oajoin: writing " << arraylist->size() << " arrays to file " << outfile << std::endl;

    if(opt.getFlag('l')!=0) {
        writearrayfile(outfile.c_str(), arraylist, arrayfile::ALATEX);
    }
    else {
      if (arraylist->size()==0) {
	if (verbose>=2)
	  printf("mode %d, nr %d, nc %d\n", mode, nrows, ncols);
      writearrayfile(outfile.c_str(), arraylist, mode, nrows, ncols);
      }
	else {
      writearrayfile(outfile.c_str(), arraylist, mode);
	}
    }

    /* free allocated structures */
    //free_sols(*arraylist);
    delete arraylist;

    return 0;
}
