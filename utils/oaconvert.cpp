/** \file oaconvert.cpp

 C++ program: oaconvert

 oaconvert can convert oa files

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2012

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
    opt.setOption("format", 'f');
    opt.setOption("verbose", 'v');

    opt.addUsage( "Orthonal Array Convert: convert array file into a different format" );
    opt.addUsage( "Usage: oaconvert [OPTIONS] [INPUTFILE] [OUTPUTFILE]" );
    opt.addUsage( "" );
    opt.addUsage( " -h --help  			Prints this help " );
    //opt.addUsage( " -s --sort			Sort the arrays " );
    opt.addUsage( " -f [FORMAT]					Output format (default: TEXT, or BINARY) " );
    opt.processCommandArgs(argc, argv);

    int verbose =  opt.getIntValue("verbose", 2);

    if (verbose>0)
        print_copyright();

    /* parse options */
    if( opt.getFlag( "help" ) || opt.getFlag( 'h' ) || opt.getArgc()==0 ) {
        opt.printUsage();
        exit(0);
    }

    const std::string outputprefix = opt.getStringValue('o', "");
    std::string format = opt.getStringValue('f', "BINARY");

    arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString(format);

    if (opt.getArgc() != 2) {
        opt.printUsage();
        exit(0);

    }

    std::string infile = opt.getArgv( 0 );
    std::string outfile = opt.getArgv( 1 );

    arrayfile_t af(infile.c_str());
    int nc= af.ncols;
    int nr=af.nrows;
	array_link al(nr, nc, -1);

	if (!af.isopen() ) {
	  fprintf(stderr, "oaconvert: could not open input file, aborting...\n");
	  return 1;
	}
	if (verbose)
	  printf("oaconvert: output mode %d, nr %d nc %d\n", mode, nr, nc);
	
    if (0) {
	bool sortarrays = opt.getFlag('s');
        /* read in the arrays */
        arraylist_t *arraylist = new arraylist_t;
        if (verbose)
            std::cout << "oaconvert: reading " << opt.getArgc() << " file(s)" << endl;
        int n = readarrayfile( infile.c_str(), arraylist);


        /* perform operations on arrays */
        if (sortarrays) {
            if (verbose)
                std::cout << "Sorting arrays" << endl;
            sort(arraylist->begin(), arraylist->end());
        }

        /* write arrays to file */

        if (verbose)
            std::cout << "oaconvert: writing " << arraylist->size() << " arrays to file " << outfile << std::endl;


        writearrayfile(outfile.c_str(), arraylist, mode, nr, nc);

        /* free allocated structures */
        delete arraylist;
    } else {
		// streaming mode
		
        int narrays=af.narrays;
        int nb=8; // TODO: we have not read any arrays so far, so nb is hard to predict
        if(mode==ABINARY_DIFFZERO)
		  nb=1;

        arrayfile_t afout(outfile.c_str(), nr, nc, narrays, mode, nb);

        if (af.narrays<0) {
            narrays = arrayfile_t::NARRAYS_MAX;
            printf("warning: untested code! (number of arrays undefined)\n");
        }

        long i;
        int index;
        for (i = 0; i < narrays; i++)
        {
            if ( (i%10000==0 && verbose) || (verbose>=3) ) {
                log_print(QUIET, "oaconvert: loading arrays: %d/%d\n", i, af.narrays);
            }

            //al.show();
            int g = af.read_array(al);
            if (g<0) {
			  if (verbose)
                printf("   oaconvert: read_array returned index %d, end of file\n", g);
                break;
            }

            afout.append_array(al);
        }

        afout.finisharrayfile();
    }

    return 0;
}
