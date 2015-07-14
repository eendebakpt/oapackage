/** \file oainfo.cpp

 C++ program: oainfo

 Print information about OA data file

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "arraytools.h"
#include "anyoption.h"

/**
 * @brief Read in files with arrays and join them into a single file
 * @param argc
 * @param argv[]
 * @return
 */
int main(int argc, char* argv[])
{

    /* parse command line options */
    AnyOption opt;
    opt.setFlag(  "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption("format", 'f');
    opt.setOption("verbose", 'v');

    opt.addUsage( "Orthonal Array Info: print information about files" );
    opt.addUsage( "Usage: oainfo [OPTIONS] [FILES]" );
    opt.addUsage( "" );
    opt.addUsage( " -h --help  			Prints this help " );
    opt.addUsage( " -v --verbose  			Verbose level " );
    opt.processCommandArgs(argc, argv);


    int verbose =  opt.getIntValue("verbose", 2);
    if (verbose>=2) {
        print_copyright_light();
    }

    /* parse options */
    if( opt.getFlag( "help" ) || opt.getFlag( 'h' ) || opt.getArgc()==0 ) {
        opt.printUsage();
        exit(0);
    }

    /* read in the arrays */
    if (verbose>=2)
        cout << "oainfo: reading " << opt.getArgc() << " file(s)" << endl;

    for( int i = 0 ; i < opt.getArgc() ; i++ ) {
        if (verbose>=3)
            std::cout << "file " <<  opt.getArgv( i ) << endl ;
        const char *fname = opt.getArgv( i );
	std::string fnames(fname);
        arrayfile_t afile(fnames, 3*(verbose>=3) );

        if (afile.isopen() ) {
            std::cout << afile.showstr() << std::endl;
        } else {
            // try to read as binary file
            if (verbose>=3) {
                printf("oainfo: opening as data file...\n");
            }
            int nr;
            int nc;
            bool valid=false;
            FILE *fid = fopen(fname, "rb");
            if (fid>0) {
                //printf("fid %d\n", fid);
                if (verbose>=3) {
                    printf("oainfo: fid %ld\n", (long)fid);
                }

                valid= readbinheader(fid, nr, nc);
                if (verbose>=4)
                    printf("oainfo: readbinheader valid: %d\n", valid);
                fclose(fid);
            }
            if (valid==true) {
                std::cout << printfstring("file %s: binary data file with %d rows, %d columns\n", fname, nr, nc);
            }
            else {
                std::cout << printfstring("file %s: invalid file\n", fname);

            }
        }
    }

    return 0;
}
