
/** \file oatest.cpp

C++ program: oatest

oatest: tool for testing new algorithms

Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

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
#include "extend.h"
#include "graphtools.h"
#include "graphtools.h"
#include "tools.h"

#include "evenodd.h"
#include "lmc.h"

#include "conference.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

#include "Deff.h"
#include "arraytools.h"
#include "strength.h"

#ifdef HAVE_BOOST
#include <boost/filesystem.hpp>
#include <string>
#endif

using namespace Eigen;


#include "graphtools.h"

array_link array2xf2 (const array_link &al) {
        const int k = al.n_columns;
        const int n = al.n_rows;
        const int m = 1 + k + k * (k - 1) / 2;
        array_link out (n, m, array_link::INDEX_DEFAULT);

        // init first column
        int ww = 0;
        for (int r = 0; r < n; ++r) {
                out.array[r] = 1;
        }

        // init array
        ww = 1;
        for (int c = 0; c < k; ++c) {
                int ci = c * n;
                array_t *pout = out.array + (ww + c) * out.n_rows;
                for (int r = 0; r < n; ++r) {
                        pout[r] = 2 * al.array[r + ci] - 1;
                }
        }

        // init interactions
        ww = k + 1;
        for (int c = 0; c < k; ++c) {
                int ci = c * n + n;
                for (int c2 = 0; c2 < c; ++c2) {
                        int ci2 = c2 * n + n;

                        const array_t *p1 = out.array + ci;
                        const array_t *p2 = out.array + ci2;
                        array_t *pout = out.array + ww * out.n_rows;

                        for (int r = 0; r < n; ++r) {
                                pout[r] = -(p1[r] * p2[r]);
                        }
                        ww++;
                }
        }
        return out;
}

/// show information about Pareto criteria for conference matrix
void paretoInfo (const array_link &alx) {
        std::vector< int > j5 = alx.Jcharacteristics (5);
        int j5max = vectormax (j5, 0);

        int v1 = (j5max == alx.n_rows);
        int v2 = 1 - v1;

        int N = alx.n_rows;
        int rank = array2xf (alx).rank ();
        std::vector< int > F4 = alx.Fvalues (4);
        std::vector< double > gwlp = alx.GWLP ();
        printf ("pareto data: %d ; ", rank);
        printf (" %d ", (int)(N * N * gwlp[4]));
        printf (" ; ");
        display_vector (F4);
        printf (" ; ");
        printf (" %d ; %d", v1, v2);
        printf ("\n");
}

/// check whether an array is contained in a Pareto set
int arrayInPareto (const Pareto< mvalue_t< long >, array_link > &pset, const array_link &al, int verbose = 1) {

        std::vector< array_link > llx = pset.allindices ();
        arraylist_t ll (llx.begin (), llx.end ());
        int jx = arrayInList (al, ll, 0);
        if (verbose)
                myprintf ("arrayInPareto: index in pareto list %d\n", jx);
        return jx;
}

/// check composition operator. returns 0 if test id good
int checkConferenceComposition (const array_link &al, int verbose = 0) {
        conference_transformation_t T1 (al);
        // T1.randomizecolperm();
        T1.randomizecolflips ();
        // T1.randomizerowperm();
        T1.randomizerowflips ();
        T1.randomize ();

        conference_transformation_t T2 (al);
        // T2.randomize();
        T2.randomizerowperm ();
        // T2.randomizerowflips();
        // T2.randomizecolperm();
        // T2.randomizecolflips();

        conference_transformation_t T3 = T2 * T1;

        array_link al1 = T1.apply (al);
        array_link al1t2 = T2.apply (al1);
        array_link al3 = T3.apply (al);

        if (verbose) {
                printfd ("checkTransformationComposition: transforms\n");
                printf ("-- T1 \n");
                T1.show ();
                printf ("-- T2 \n");
                T2.show ();
                printf ("-- T3 \n");
                T3.show ();
                printfd ("checkTransformationComposition: arrays\n");
                al.showarray ();
                al1.showarray ();
                al1t2.showarray ();
                al3.showarray ();
        }

        myassert (al3 == al1t2, "unittest error: composition of conference transformations\n");

        return 0;
}

int main (int argc, char *argv[]) {
        AnyOption opt;
        /* parse command line options */
        opt.setFlag ("help", 'h'); /* a flag (takes no argument), supporting long and short form */
        opt.setOption ("output", 'o');
        opt.setOption ("input", 'I');
        opt.setOption ("rand", 'r');
        opt.setOption ("verbose", 'v');
        opt.setOption ("ii", 'i');
        opt.setOption ("jj");
        opt.setOption ("xx", 'x');
        opt.setOption ("dverbose", 'd');
        opt.setOption ("rows");
        opt.setOption ("cols");
        opt.setOption ("nrestarts");
        opt.setOption ("niter");
        opt.setOption ("mdebug", 'm');
        opt.setOption ("oaconfig", 'c'); /* file that specifies the design */

        opt.addUsage ("Orthonal Array: oatest: testing platform");
        opt.addUsage ("Usage: oatest [OPTIONS] [FILE]");
        opt.addUsage ("");
        opt.addUsage (" -h --help  			Prints this help ");
        opt.processCommandArgs (argc, argv);

        double t0 = get_time_ms (), dt = 0;
        int randvalseed = opt.getIntValue ('r', 1);
        int ix = opt.getIntValue ('i', 1);
        int r = opt.getIntValue ('r', 1);
        int jj = opt.getIntValue ("jj", 5);

        int xx = opt.getIntValue ('x', 0);
        int niter = opt.getIntValue ("niter", 10);
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


		try {
			array_link al = exampleArray(r);
			al.show();
			al.showarray();

			std::vector<int> sizes = array2modelmatrix_sizes(al);
			display_vector(sizes); myprintf("\n");
			MatrixFloat modelmatrix = array2modelmatrix(al, "i", 1);
			array_link modelmatrixx = modelmatrix;
			modelmatrixx.show();
			modelmatrixx.showarray();

			modelmatrix = array2modelmatrix(al, "main", 1);
			modelmatrixx = modelmatrix;
			modelmatrixx.show();
			modelmatrixx.showarray();

			exit(0);

		}
		catch (const std::exception &e) {
			std::cerr << e.what() << std::endl;
			throw;
		}

        {
                array_link al = exampleArray (r);

                array_transformation_t tt = reduceOAnauty (al);
                array_link alx = tt.apply (al);
                exit (0);
        }



        return 0;
}

// kate: indent-mode cstyle; indent-width 5; replace-tabs on;
