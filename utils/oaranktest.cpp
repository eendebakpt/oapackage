/** \file oaranktest.cpp

C++ program: oaranktest

oaranktest: tool for testing speed of rank calculations

Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2016

Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <algorithm>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arrayproperties.h"
#include "arraytools.h"
#include "tools.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

#include "Deff.h"
#include "arraytools.h"
#include "strength.h"

using namespace Eigen;

int main (int argc, char *argv[]) {
        AnyOption opt;
        /* parse command line options */
        opt.setFlag ("help", 'h'); 
        opt.setOption ("output", 'o');
        opt.setOption ("input", 'I');
        opt.setOption ("rand", 'r');
        opt.setOption ("niter", 'n');
        opt.setOption ("verbose", 'v');
        opt.setOption ("ii", 'i');

        opt.addUsage ("Orthonal Array: oaranktest: testing platform");
        opt.addUsage ("Usage: oatest [OPTIONS] [FILE]");
        opt.addUsage ("");
        opt.addUsage (" -h --help  			Prints this help ");
        opt.processCommandArgs (argc, argv);

        double t0 = get_time_ms (), dt = 0;
        int randvalseed = opt.getIntValue ('r', 1);
        int xx = opt.getIntValue ('x', 0);
        int niter = opt.getIntValue ("niter", 2);
        int verbose = opt.getIntValue ("verbose", 1);

        const char *input = opt.getValue ('I');
        if (input == 0)
                input = "test.oa";

        srand (randvalseed);
        if (randvalseed == -1) {
                randvalseed = time (NULL);
                printf ("random seed %d\n", randvalseed);
                srand (randvalseed);
        }

        print_copyright ();

        /* parse options */
        if (opt.getFlag ("help") || opt.getFlag ('h') || opt.getArgc () < 0) {
                opt.printUsage ();
                exit (0);
        }

        /* read data from file */

        double t00 = get_time_ms ();
        t0 = get_time_ms ();
        arraylist_t ll = readarrayfile (input);
        std::sort (ll.begin (), ll.end ());
        dt = get_time_ms () - t0;

        int s = ll.size ();
        for (int j = 0; j < niter; j++) {
                for (int i = 0; i < s; i++) {
                        ll.push_back (ll[i]);
                }
        }

        if (verbose)
                printf ("oaranktest: %ld arrays (reading %.3f [s])\n", ll.size (), dt);

        const long nn = ll.size ();
        std::vector< int > rr (nn);
        t0 = get_time_ms ();
        for (size_t i = 0; i < ll.size (); i++) {
                array_link A = ll[i];
                array_link B = array2xf (A);
                int r1 = arrayrankColPivQR (B);
                rr[i] = r1;

                if (verbose >= 3) {
                        B.show ();
                        int r3 = arrayrankFullPivLU (B);
                        int r2 = arrayrankSVD (B);
                }
        }
        const array_link al0 = ll[0];
        printf ("warm-up complete\n");

        t0 = get_time_ms ();
        for (size_t i = 0; i < ll.size (); i++) {
                array_link A = ll[i];
                array_link B = array2secondorder (A);
        }
        dt = get_time_ms () - t0;
        printf ("oaranktest: conversion to second order (%.3f [s], %.3f Marrays/s)\n", dt, nn / dt);

        if (1) {
                t0 = get_time_ms ();
                for (size_t i = 0; i < ll.size (); i++) {
                        array_link A = ll[i];
                        array_link B = array2xf (A);
                        int r = arrayrankSVD (B);
                        if (r != rr[i]) {
                                printfd ("error: i %d, r %d rr[i] %d\n", i, r, rr[i]);
                        }
                        myassert (r == rr[i], "arrayrankSVD");
                }
                dt = get_time_ms () - t0;
                printf ("oaranktest: rank SVD (%.3f [s], %.3f Marrays/s)\n", dt, nn / dt);
        }

        if (1) {

                t0 = get_time_ms ();
                for (size_t i = 0; i < ll.size (); i++) {
                        array_link A = ll[i];
                        array_link B = array2xf (A);
                        int r = arrayrankFullPivLU (B);

                        if (r != rr[i]) {
                                printfd ("error: i %d, r %d rr[i] %d\n", i, r, rr[i]);
                        }
                        myassert (r == rr[i], "FullPivLU");
                }
                dt = get_time_ms () - t0;
                printf ("oaranktest: rank FullPivLU (%.3f [s], %.3f Marrays/s)\n", dt, nn / dt);
        }

        if (verbose >= 2) {
                array2xf (al0).showarray ();
        }
        if (1) {
                t0 = get_time_ms ();
                for (size_t i = 0; i < ll.size (); i++) {
                        array_link A = ll[i];
                        array_link B = array2xf (A);
                        int r = arrayrankColPivQR (B);
                        if (r != rr[i]) {
                                printfd ("error: i %d, r %d rr[i] %d\n", i, r, rr[i]);
                        }
                        myassert (r == rr[i], "arrayrankColPivQR");
                }
                dt = get_time_ms () - t0;
                printf ("oaranktest: rank ColPiv (%.3f [s], %.3f Marrays/s)\n", dt, nn / dt);
        }

        for (int nsub = 2; nsub < 5; nsub++) {
                rankStructure rank_calculator (al0.selectFirstColumns (al0.n_columns - nsub), nsub, 0);
                if (verbose >= 2) {
                        rank_calculator.alsub.show ();
                        printf ("subrank: %d\n", arrayrankFullPivLU (array2xf (rank_calculator.alsub)));
                        printf ("---\n\n");
                }
                t0 = get_time_ms ();
                for (size_t i = 0; i < ll.size (); i++) {
                        array_link A = ll[i];
                        int r = rank_calculator.rankxf (A);
                        if (r != rr[i]) {
                                printfd ("  i %d, r %d rr[i] %d\n", i, r, rr[i]);
                        }
                        assert (r == rr[i]);
                }
                dt = get_time_ms () - t0;
                printf ("oaranktest: rank ColPivHouseholderQR-cache (nsub %d, %.3f [s], %.3f Marrays/s)\n", nsub, dt,
                        nn / dt);
        }

        return 0;
}

