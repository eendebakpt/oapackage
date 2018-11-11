/** \file oa_depth_extend.cpp

 C++ program: oa_depth_extend

 oa_depth_extend: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arrayproperties.h"
#include "arraytools.h"
#include "extend.h"
#include "strength.h"
#include "tools.h"

#include "oadevelop.h"

#ifdef DOOPENMP
#include "omp.h"
#endif

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
        opt.setFlag ("split", 's');
        opt.setFlag ("extend", 'e');
        opt.setFlag ("prune", 'p');
        opt.setOption ("j5structure", 'j');
        opt.setOption ("verbose", 'v');
        opt.setOption ("kfinal", 'k');
        opt.setOption ("hack", 'H');
        opt.setOption ("method", 'm');
        opt.setOption ("writearrays", 'w');
        opt.setOption ("maxarrays");
        opt.setOption ("discardJ5");
        opt.setOption ("maxk", 'K');
        opt.setOption ("Avalue", 'A');
        opt.setFlag ("coptions"); /* print compile time options */
        opt.setOption ("logtime", 'Q');
        opt.setOption ("format", 'f');
        opt.setOption ("oaconfig", 'c'); /* file that specifies the design */

        opt.addUsage ("Orthonal Array: oa_depth_extend: testing platform");
        opt.addUsage ("Usage: oa_depth_extend [OPTIONS] [FILES]");
        opt.addUsage ("");
        opt.addUsage (" -h --help  			Prints this help ");
        opt.addUsage (" -f [FORMAT]					Output format (default: TEXT, or BINARY) ");
        opt.addUsage (" -o [FILE] --output [FILE]	Output prefix (default: test-depth) ");
        opt.addUsage (" -A [Afinal] ");
        opt.addUsage (" -v [verbose level] ");
        opt.addUsage (" --split 				Split calculation ");
        opt.addUsage (" -p prune, -e extend ");
        opt.addUsage (" -K --maxK [n] ");

        opt.addUsage (" -Q [time] 			Logging time (in seconds)");
        opt.addUsage (
            printfstring (
                " --j5structure [VALUE] 			Integer, can be J5_ORIGINAL (%d) or J5_45 (%d)",
                J5_ORIGINAL, J5_45)
                .c_str ());
        std::string ss = printfstring (" -m [MODE]			Algorithm (") + algorithm_t_list () + ")\n";
        opt.addUsage (ss.c_str ());
        opt.processCommandArgs (argc, argv);

        print_copyright ();
        setloglevel (NORMAL);

        /* parse options */
        if (opt.getFlag ("help") || opt.getFlag ('h') || opt.getArgc () < 0) {
                opt.printUsage ();
                exit (0);
        }

        if (opt.getValue ("coptions") != 0) {
                opt.printUsage ();
                print_options (cout);
                exit (0);
        }

        logstream (SYSTEM) << "#time start: " << currenttime () << std::endl;
        double time0 = get_time_ms ();

        const char *oaconfigfile = "oaconfig.txt";

        if (opt.getValue ("oaconfig") != NULL)
                oaconfigfile = opt.getValue ('c');

        arrayfile::arrayfilemode_t mode = arrayfile::ABINARY_DIFF;

        if (opt.getValue ("format") != 0) {
                std::string format = opt.getValue ('f');
                mode = arrayfile_t::parseModeString (format);
        }

        std::string outputprefix = "test-depth";
        if (opt.getValue ("output") != 0)
                outputprefix = opt.getValue ('o');

        int verbose = opt.getIntValue ('v', 2);
        double Afinal = opt.getDoubleValue ('A', 0);
        int kfinal = opt.getIntValue ('k', 7);
        int splitcalc = opt.getFlag ("split");
        int discardJ5 = opt.getIntValue (
            "discardJ5", -1); // discard arrays with J5==N with size larger then specfied number of columns

        int maxk = opt.getIntValue ("maxk", 1000);
        int writedeptharrays = opt.getIntValue ('w', 1);
        long maxarrays = opt.getIntValue ("maxarrays", 2047483000);

        int hack = opt.getIntValue ('H', 0);

        if (discardJ5 > 0) {
                printf ("  discardJ5: %d\n", discardJ5);
        }
#ifdef _OPENMP
        printf ("  openmp: num threads %d, max num threads %d\n", omp_get_num_threads (), omp_get_max_threads ());

        if (omp_get_nested () != 0)
                printf ("note: omp_get_nested()=%d, make sure to set OMP_THREAD_LIMIT!\n", omp_get_nested ());
#endif

        int doextend = opt.getFlag ('e');
        int dopruning = opt.getFlag ('p');
        j5structure_t j5structure = (j5structure_t)opt.getIntValue ('j', J5_45);
        int method = opt.getIntValue ('m', MODE_J5ORDERX);

        arraydata_t *adfull = readConfigFile (oaconfigfile);
        if (adfull == 0) {
                fprintf (stderr, "oa_depth_extend: could not read config file");
                exit (1);
        }
        if (maxk < adfull->ncols) {
                if (maxk <= 5) {
                        printf ("oa_depth_extend: maxk should be >= 5\n");
                        return 1;
                }
                arraydata_t *adfullr = new arraydata_t (adfull, maxk);
                delete adfull;
                adfull = adfullr;
        }
        arraylist_t *arraylist = new arraylist_t;
        if (opt.getArgc () == 0) {
                cout << "  adding root to arraylist\n";
                array_link al (adfull->N, adfull->strength, -1);
                al.create_root (*adfull);
                arraylist->push_back (al);
        } else {
                /* read in the arrays */
                cout << "Reading " << opt.getArgc () << " file(s)" << endl;
                for (int i = 0; i < opt.getArgc (); i++) {
                        int n = readarrayfile (opt.getArgv (i), arraylist);
                        cout << "file " << opt.getArgv (i) << ":   read " << n << " array(s)" << endl;
                }
        }

        arraylist_t earrays;

        OAextend oaextendx;
        oaextendx.setAlgorithm ((algorithm_t)method, adfull);
        oaextendx.j5structure = j5structure;

        oaextendx.checkarrays = 0;
        oaextendx.use_row_symmetry = 0;
        oaextendx.extendarraymode = 1;
        oaextendx.init_column_previous = INITCOLUMN_J5;
        // oaextendx.init_column_previous=INITCOLUMN_ZERO;
        oaextendx.nLMC = 500000;
        oaextendx.info ();

        if (verbose)
                printf ("oa_depth_extend: initialization phase (splitcalc %d)\n", splitcalc);
        depth_extend_t dextend (adfull, 10, discardJ5);
        dextend.oaextend = oaextendx;
        dextend.logtime = opt.getDoubleValue ("logtime", 240);
        dextend.setNarraysMax (maxarrays);

        int initcolumn = adfull->strength + 1;
        if (arraylist->size () > 0)
                initcolumn = arraylist->at (0).n_columns;

        dextend.arraywriter = new arraywriter_t ();
        dextend.arraywriter->initArrayFiles (*dextend.ad, initcolumn + 1, outputprefix, mode);
        dextend.arraywriter->writearrays = writedeptharrays;
        dextend.counter = new counter_t (adfull->N);
        dextend.loglevelcol = 7;

        setloglevel (SYSTEM);

        double t0 = get_time_ms ();

        depth_extensions_storage_t *ds = 0;

        if (splitcalc) {
                printf ("splitting calculation\n");
                ds = new depth_extensions_storage_t ();
                ds->resize (arraylist->size ());
        } else {
        }

// loop over all arrays
#pragma omp parallel for schedule(dynamic, 1)
        for (int ai = 0; ai < (int)arraylist->size (); ai++) {
                const array_link &al = arraylist->at (ai);

                if (al.n_columns >= maxk)
                        continue;

                myassert (adfull->N == al.n_rows, "oa_depth_extend: nrows array, nrows config\n");

                if (verbose >= 3 || (verbose >= 2 && ai % 40 == 0))
                        printf ("oa_depth_extend: array %d/%ld (%d %d): time %.1f [s]\n", ai, arraylist->size (),
                                al.n_rows, al.n_columns, get_time_ms () - t0);
                flush_stdout ();

                depth_extend_t dextendloop (dextend);
                dextendloop.setposition (al.n_columns, ai, arraylist->size (), 0, 0);

                depth_extend_array (al, dextendloop, *adfull, verbose, ds, ai);

                if (ai % 800 == 0) {
                        printf ("array %d: ", ai);
                        dextend.counter->showcountscompact ();
                }
        }

        if (ds != 0) {
                if (verbose) {
                        printf ("oa_depth_extend: calling processDepth\n");
                        std::cout << "#time step: " << printfstring ("%.1f", get_time_ms () - time0) << " [s]"
                                  << std::endl;
                }

                for (size_t ai = 0; ai < arraylist->size (); ai++) {
                        printf ("ai %ld: %ld arrays for extension\n", ai, ds->goodarrayslist[ai].size ());
                }
                dextend.counter->showcounts ("after init", adfull->strength, adfull->ncols);

                for (size_t ai = 0; ai < arraylist->size (); ai++) {
                        const array_link &al = arraylist->at (ai);
                        if (verbose >= 3 || (verbose >= 2 && ai % 40 == 0))
                                printf ("oa_depth_extend process: array %ld/%ld (%d %d): time %.1f [s]\n", ai,
                                        arraylist->size (), al.n_rows, al.n_columns, get_time_ms () - t0);
                        flush_stdout ();
                        int extensioncol = al.n_columns + 1;

                        depth_extend_t dextendloop (dextend);
                        dextendloop.extension_column_list = ds->columnextensionsList[ai];
                        dextendloop.setposition (al.n_columns, ai, arraylist->size (), 0, 0);
                        dextendloop.loglevelcol = al.n_columns;

                        processDepth (ds->goodarrayslist[ai], ds->depthalglist[ai], dextendloop,
                                      ds->dextendsubList[ai], extensioncol, verbose);
                }
        }

        if (dextend.discardJ5 >= 0) {
                printf ("dextend.discardJ5: %d, number %ld\n", dextend.discardJ5, dextend.discardJ5number);
        }
        dextend.counter->showcounts (*adfull);

        delete dextend.arraywriter;
        delete dextend.counter;

        if (ds != 0) {
                delete ds;
        }
        delete adfull;
        delete arraylist;

        cleanGlobalStatic ();
        arraysymmetry::rowpermpool.reset ();

        logstream (SYSTEM) << "#time end: " << currenttime () << std::endl;
        if (verbose)
                std::cout << "#time total: " << printfstring ("%.1f", get_time_ms () - time0) << " [s]" << std::endl;

        return 0;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
