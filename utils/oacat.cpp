/** \file oacat.cpp

 C++ program: oacat

 oacat: print to contents of an OA file or data file to stdout

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arraytools.h"
#include "tools.h"

/**
 * @brief Read in files with arrays and print to stdout
 * @param argc
 * @param argv[]
 * @return
 */
int main (int argc, char *argv[]) {
        AnyOption opt;

        /* parse command line options */
        opt.setFlag ("help", 'h'); 
        opt.setOption ("verbose", 'v');
        opt.setOption ("md5", 'm');

        opt.addUsage ("Orthonal Array cat: read array or data file and print to stdout");
        opt.addUsage ("Usage: oacat [OPTIONS] [FILES]");
        opt.addUsage ("");
        opt.addUsage (" -h --help  			Prints this help ");
        opt.addUsage (" -m --md5 [INT] 			If 1, then display md5 sum for each array. Default is 0.");
        opt.processCommandArgs (argc, argv);

        /* parse options */
        if (opt.getFlag ("help") || opt.getFlag ('h') || opt.getArgc () == 0) {
                print_copyright ();
                opt.printUsage ();
                exit (0);
        }

        int verbose = opt.getIntValue ("verbose", 2);
        int domd5 = opt.getIntValue ("md5", 0);

        if (verbose >= 2)
                print_copyright ();

        /* read in the files */
        if (verbose)
                std::cout << "oacat: reading " << opt.getArgc () << " file(s)" << endl;

        for (int i = 0; i < opt.getArgc (); i++) {
                const char *fname = opt.getArgv (i);
                arrayfile_t afile (opt.getArgv (i), 0);

                if (afile.isopen ()) {
                        if (verbose)
                                cout << "file " << opt.getArgv (i) << endl;
                        arraylist_t arraylist = readarrayfile (opt.getArgv (i), 0);
                        if (verbose >= 2) {
                                cout << "   read " << arraylist.size () << " array(s)" << endl;
                        }
                        for (size_t j = 0; j < arraylist.size (); j++) {
                                if (domd5) {
                                        std::string m = arraylist[j].md5 ();
                                        printf ("array %d: %s\n", (int)j, m.c_str ());
                                } else {
                                        arraylist[j].showarray ();
                                }
                        }
                } else {
                        // try to read as binary data file
                        int number_rows;
                        int number_columns;
                        bool valid_binary_data = false;
                        FILE *fid = fopen (fname, "rb");
                        if (fid != 0) {
                                valid_binary_data = readbinheader (fid, number_rows, number_columns);
                                if (valid_binary_data) {
                                        if (verbose) {
                                                printf ("data file %s: %d %d\n", fname, number_rows, number_columns);
                                        }
                                        double *dd = new double[number_columns];
                                        for (size_t r = 0; r < (size_t)number_rows; r++) {
                                                fread (dd, sizeof (double), number_columns, fid);
                                                for (size_t column = 0; column < (size_t)number_columns; column++) {
                                                        printf ("%f ", dd[column]);
                                                }
                                                printf ("\n");
                                        }
                                        delete[] dd;
                                } else {
                                        if (verbose) {
                                                printf ("file %s: no array file or data file\n", fname);
                                        }
                                }
                                fclose (fid);
                        }
                }
        }

        return 0;
}
