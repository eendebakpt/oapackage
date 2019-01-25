
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
#include "unittests.h"
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



		array_link array = exampleArray(0);
		lmc_t lmc_type = LMCcheck(array);


		array = array.randomperm();
		array.showarray();
		array_link reduced_array = reduceLMCform(array);
		reduced_array.showarray();
		exit(0);

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
