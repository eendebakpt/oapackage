/** \file oaconvert.cpp

 C++ program: oaconvert

 oaconvert can convert oa files

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2012

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arraytools.h"
#include "tools.h"

/**
 * @brief Read in files with arrays and join them into a single file
 * @param argc
 * @param argv[]
 * @return
 */
int main (int argc, char *argv[]) {
        AnyOption opt;

        /* parse command line options */
        opt.setFlag ("help", 'h'); /* a flag (takes no argument), supporting long and short form */
        opt.setOption ("output", 'o');
        opt.setOption ("format", 'f');
        opt.setOption ("verbose", 'v');

        opt.addUsage ("Orthonal Array Convert: convert array file into a different format");
        opt.addUsage ("Usage: oaconvert [OPTIONS] [INPUTFILE] [OUTPUTFILE]");
        opt.addUsage ("");
        opt.addUsage (" -h --help  			Prints this help ");
        opt.addUsage (" -f [FORMAT]					Output format (default: TEXT, or BINARY) ");
        opt.processCommandArgs (argc, argv);

        int verbose = opt.getIntValue ("verbose", 2);

        if (verbose > 0)
                print_copyright ();

        /* parse options */
        if (opt.getFlag ("help") || opt.getFlag ('h') || opt.getArgc () == 0) {
                opt.printUsage ();
                exit (0);
        }

        const std::string outputprefix = opt.getStringValue ('o', "");
        std::string format = opt.getStringValue ('f', "BINARY");

        arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString (format);

        if (opt.getArgc () != 2) {
                opt.printUsage ();
                exit (0);
        }

        std::string infile = opt.getArgv (0);
        std::string outfile = opt.getArgv (1);

		convert_array_file(infile, outfile, mode, verbose);

        return 0;
}
