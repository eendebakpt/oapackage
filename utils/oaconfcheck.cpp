/** \file oatest.cpp

C++ program: oatest

oatest: tool for testing new algorithms

Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <map>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "evenodd.h"
<<<<<<< HEAD
//#include "oadevelop.h"
=======
>>>>>>> upstream/master
#include "lmc.h"

#include "conference.h"

int LMC0check_firstrow ( const array_link &al, std::vector<int> &colperm,  int ci ) {
    printf ( "not implemented...\n" );
    return LMC_LESS;

}

lmc_t compare_columns ( const array_link &al, std::vector<int> &colperm, int column ) {
    const int nrows=al.n_rows;

    //printf("compare_columns: source %d, target %d\n", column, colperm[column]);
    for ( int r=0; r<nrows; r++ ) {
        if ( al.atfast ( r, column ) < al.atfast ( r, colperm[column] ) )
            return LMC_MORE;
        if ( al.atfast ( r, column ) > al.atfast ( r, colperm[column] ) )
            return LMC_LESS;
    }
    return LMC_EQUAL;
}

lmc_t LMC0_columns ( const array_link &al, std::vector<int> &colperm, int column, int verbose=0 ) {
    const int ncols = al.n_columns;
    lmc_t r = LMC_NONSENSE;

    if ( verbose ) {
        printf ( "LMC0_columns: column %d, colperm ", column );
        display_vector ( colperm );
        printf ( "\n" );
    }
    for ( int c=column; c<ncols ; c++ ) {
        // select first column
        int col = colperm[c];
        colperm[c]=colperm[column];
        colperm[column]=col;

        // compare the current pair of columns
        r = compare_columns ( al, colperm, column );


        //printf("LMC0_columns: column comparison at %d: %d\n", column, r);
        if ( r==LMC_EQUAL ) {
            // columns are equal, so go one column deeper
            if ( verbose>=2 )
                printf ( "LMC0_columns: go to col %d\n", column+1 );
            r = LMC0_columns ( al, colperm, column+1 );
        }
        if ( r==LMC_LESS ) {
            // we already know the array is not in minimal form
            break;
        }

        // result is LMC_MORE, continue with the calculation

        // swap back column
        colperm[c]=col;
        colperm[column]=colperm[col];
    }
    return r;
}

/** Algorithm: check an array X for LMC0 form
 *
 * The tag [B] indicates branching
 *
 *
 *
 * 1. Select the first column [B]
 * 2. Apply row level permutations such that the first column only contains ones (and a single zero)
 * 3. Sort the rows of the design
 * 4. Select one of two possible sign permutations for the first row [B]
 * 5. Select the next column k [B]
 * 6. With the selected column:
 *    i. Apply the correct column level permutation to make the element X(1, k) equal to 1
 *    ii. Sort the rows using the ordering 0, 1, -1
 *    iii. Compare the column the the original column k of the design. If less, then abort with LMC_LESS, if more, continue to 5, if equal then go to 5 with k+1
 * 7. Return LMC_MORE
 *
 *
 * When selecting columns we use a std::vector to keep track of selected indices
 * When sorting the rows to do not sort the rows, but use a rowsort_t structure to keep track of the order of rows
 **/

/**
 * Dummy algorithm: only allow column permutations
 *
 *
 */
lmc_t LMC0check ( const array_link &al ) {

    const int ncols = al.n_columns;

    // initialize the column permutation
    std::vector<int> colperm ( ncols ); // create pointer for the column permutations
    init_perm ( colperm );

    // call the recursive part starting at column 0
    lmc_t r = LMC0_columns ( al, colperm, 0 );

    return r;
}



int main ( int argc, char* argv[] ) {
    AnyOption opt;
    /* parse command line options */
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "output", 'o' );
    opt.setOption ( "input", 'I' );
    opt.setOption ( "rand", 'r' );
    opt.setOption ( "verbose", 'v' );
    opt.setOption ( "ii", 'i' );
    opt.setOption ( "jj" );
    opt.setOption ( "rows" );
    opt.setOption ( "cols" );
    opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

    opt.addUsage ( "Orthonal Array: oaconfcheck: testing platform" );
    opt.addUsage ( "Usage: oaconfcheck [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.processCommandArgs ( argc, argv );


    double t0=get_time_ms(), dt=0;
    int randvalseed = opt.getIntValue ( 'r', 1 );
    int ix = opt.getIntValue ( 'i', 1 );
    int jj = opt.getIntValue ( "jj", 5 );

    int xx = opt.getIntValue ( 'x', 0 );
    int niter  = opt.getIntValue ( "niter", 10 );
    int verbose  = opt.getIntValue ( "verbose", 1 );

    char *input = opt.getValue ( 'I' );
    if ( input==0 )
        input="test.oa";

    srand ( randvalseed );
    if ( randvalseed==-1 ) {
        randvalseed=time ( NULL );
        printf ( "random seed %d\n", randvalseed );
        srand ( randvalseed );
    }


    print_copyright();

    /* parse options */
    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
        opt.printUsage();
        exit ( 0 );
    }

    setloglevel ( SYSTEM );

    /* test on known set of arrays */
    printf("### test code on known array\n");
    array_link al = exampleArray ( 2 );
    lmc_t r=LMC0check ( al );
    printf ( "minimal form check result: %d (should be %d)\n", r, LMC_MORE );

    for ( int i=0; i<3; i++ ) {
        printf ( "round %d: apply random transformation to the array\n", i );
        array_link altmp =  al.randomcolperm();
        lmc_t r=LMC0check ( altmp );

        if ( altmp==al ) {
            printf ( "  minimal form check: %d (should be %d)\n", r, LMC_MORE );
        } else        printf ( "  minimal form check: %d (should be %d)\n", r, LMC_LESS );

    }

    printf("### read from file\n");

    arraylist_t ll= readarrayfile ( input );

    for ( size_t i=0; i<ll.size(); i++ ) {
        array_link al = ll[i];

        //al = al.randomrowperm();
        //al = al.randomcolperm();

        lmc_t r = LMC0check ( al );
        printf ( "array %d: result %d\n", (int) i, ( int ) r );
    }


    return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
