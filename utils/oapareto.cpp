/** \file oapareto.cpp

 C++ program: oapareto

 oapareto: Calculate a set of Pareto optimal arrays from a list of arrays

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <numeric>

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"
#include "arrayproperties.h"

#include "pareto.h"


/**
 * @brief Analyse arrays
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] ) {
    AnyOption opt;

    /* parse command line options */
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "output", 'o' );
    opt.setOption ( "verbose", 'v' );
    opt.setOption ( "format", 'f' );
    opt.setOption ( "paretomethod" );
    opt.setFlag ( "test", 't' );

    opt.addUsage ( "Orthonal Array Analyse: convert set of arrays to set of Pareto optimal designs" );
    opt.addUsage ( "Usage: oapareto [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.addUsage ( " -v --verbose  			Verbose level (default: 1) " );
    opt.addUsage ( " -o [FILE] --output [FILE]	Output file for filtered arrays (default: no output) " );
    opt.addUsage ( " -f [FORMAT]				Output format (TEXT, BINARY (default), D (binary difference) ) " );
    opt.addUsage ( printfstring ( " --paretomethod [method]	Method to be used (PARETOFUNCTION_DEFAULT %d, PARETOFUNCTION_J5 %d ", PARETOFUNCTION_DEFAULT, PARETOFUNCTION_J5 ).c_str() );

    opt.processCommandArgs ( argc, argv );

    print_copyright();

    /* parse options */
    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() ==0 ) {
        opt.printUsage();
        exit ( 0 );
    }

    const char *outputprefix = 0;
    if ( opt.getValue ( "output" ) !=0 )
        outputprefix = opt.getValue ( 'o' );

    int verbose = opt.getIntValue ( 'v', 1 );
    paretomethod_t paretoj5 = ( paretomethod_t ) opt.getIntValue ( "paretomethod", PARETOFUNCTION_DEFAULT );
    std::string format = opt.getStringValue ( 'f', "BINARY" );
    arrayfile::arrayfilemode_t afmode = arrayfile_t::parseModeString ( format );


    /* read in the arrays */
    if ( verbose )
        std::cout << "oapareto: reading " << opt.getArgc() << " file(s)" << std::endl;

    int nr=0;
    int nc=0;
    int nn;

    initncombscache ( 20 );

    std::vector<std::string> infiles;
    for ( int i = 0 ; i < opt.getArgc() ; i++ ) {
        if ( verbose>=3 )
            std::cout << "file " <<  opt.getArgv ( i ) << std::endl ;
        infiles.push_back ( opt.getArgv ( i ) );
        arrayfileinfo ( opt.getArgv ( i ), nn, nr, nc );
        if ( verbose>=3 )
            std::cout << "   has " << nn << printfstring ( " array(s) (size %d %d)", nr, nc ) << std::endl;

    }
    calculateParetoEvenOdd ( infiles, outputprefix, verbose, afmode, nr, nc, paretoj5 );

    return 0;
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
