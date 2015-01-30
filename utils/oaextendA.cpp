/** \file oaextendA.cpp

 C++ program: oaextendA

 oaextendA: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "oadevelop.h"



/**
 * @brief Testing code
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
    opt.setFlag("directcheck", 'd');
    opt.setOption("verbose", 'v');
    opt.setOption("kfinal", 'k');
    opt.setOption("Avalue", 'A');
    opt.setOption("format", 'f');
    opt.setOption("oaconfig", 'c'); /* file that specifies the design */

    opt.addUsage( "Orthonal Array: oaextendA: testing platform" );
    opt.addUsage( "Usage: oaextendA [OPTIONS] [FILES]" );
    opt.addUsage( "" );
    opt.addUsage( " -h --help  			Prints this help " );
    opt.addUsage( " -v --verbose [INTEGER]  		Output level " );
    opt.addUsage( " -c --oaconfig [CONFIGFILE]  	Config file " );
    opt.addUsage( " -f [FORMAT]			Output format (default: TEXT, or BINARY) " );
    opt.addUsage( " -o [FILE] --output [FILE]		Output prefix (default: standard output) " );
    opt.addUsage( " -A [Dfinal]				D-efficiency threshold " );
    opt.addUsage( " -k [kfinal] 			Number of columns for threshold" );
    opt.processCommandArgs(argc, argv);

    print_copyright();
    //cout << system_uname();

    /* parse options */
    if ( opt.getFlag( "help" ) || opt.getFlag( 'h' ) || opt.getArgc()==0 ) {
        opt.printUsage();
        exit(0);
    }

    int verbose = opt.getIntValue('v', 1);

    logstream(QUIET) << "#time start: "<< currenttime() << std::endl;

    if (verbose<=1)
            setloglevel(QUIET);

        if (verbose>=3)
            setloglevel(DEBUG);

    const char *oaconfigfile = opt.getStringValue('c', "oaconfig.txt");
    const char *resultprefix;

    arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString(opt.getStringValue('f', "B"));

    const char *outputprefix = opt.getStringValue("output", "test-extendA.oa");
    double Afinal = opt.getDoubleValue('A', 0);
    int kfinal = opt.getIntValue('k', 7);

    int directcheck = opt.getFlag( "directcheck" );

    arraydata_t *ad;
    ad = readConfigFile(oaconfigfile);

    /* read in the arrays */
    arraylist_t arraylist;
    if (verbose)
      cout << "oaextendA: reading " << opt.getArgc() << " file(s)" << endl;

    for ( int i = 0 ; i < opt.getArgc() ; i++ ) {
        int n = readarrayfile( opt.getArgv( i ), &arraylist);
        cout << "file " <<  opt.getArgv( i ) << ":   read " << n << " array(s)" << endl;
    }

    OAextend oaextend;
    oaextend.singleExtendTime=60;	// report progress every minute
    
    arraylist_t earrays;
    if (verbose)
      printf("oaextendA: Afinal %.4f, kfinal %d\n", Afinal, kfinal);
  
    std::vector<std::vector<double> > edata; /// structure that contains the D-efficiencies and loss factors
    extend_array_dynamic(arraylist, *ad, oaextend, earrays, edata, kfinal,  Afinal, directcheck, verbose);

    int nc;
    if (arraylist.size()>0)
        nc = arraylist.at(0).n_columns+1;
    else
        nc=0;

    // write arrays to disk
    std::string outfile = outputprefix;
    printf("writing array file %s: %d arrays with %d columns\n", outfile.c_str(), (int)earrays.size(), nc);
    writearrayfile(outfile.c_str(), &earrays, arrayfile::ABINARY, ad->N, nc);

    // write loss factors to disk
    std::string datafile = replaceString(outfile, ".oa", "-lossfactor.bin");
    if(verbose>=2)
      printf("oaextendA: writing %d loss factors to %s\n", (int)edata.size(), datafile.c_str() );
    vectorvector2binfile(datafile, edata, 1, 2);
    delete ad;


    logstream(QUIET) << "#time end: "<< currenttime() << std::endl;

    return 0;
}
