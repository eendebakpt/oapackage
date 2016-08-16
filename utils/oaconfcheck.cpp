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
#include "oadevelop.h"
#include "lmc.h"

#include "conference.h"

int LMC0check_firstrow ( const array_link &al, std::vector<int> &colperm,  int ci ) {
    printf ( "not implemented...\n" );
    return LMC_LESS;

}

int LMC0check ( const array_link &al ) {
    /** Aproach: check an array X
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

    const int ncols = al.n_columns;

    std::vector<int> colperm ( ncols ); // create pointer for the column permutations
    //init_perm(colperm);
    
    int ci=0;
    for ( int c=ci; c<ncols ; c++ ) {
        // select first column
        int col = colperm[c];
        colperm[c]=colperm[ci];
        colperm[ci]=col;

        // go to part 2.
        int r = LMC0check_firstrow ( al, colperm, ci );

        if ( r==LMC_LESS )
            break;

        // swap back column
        colperm[c]=col;
        colperm[ci]=colperm[col];
    }

    return LMC_MORE;
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
    opt.setOption ( "xx", 'x' );
    opt.setOption ( "dverbose", 'd' );
    opt.setOption ( "rows" );
    opt.setOption ( "cols" );
    opt.setOption ( "nrestarts" );
    opt.setOption ( "niter" );
    opt.setOption ( "mdebug", 'm' );
    opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

    opt.addUsage ( "Orthonal Array: oaconfcheck: testing platform" );
    opt.addUsage ( "Usage: oaconfcheck [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.processCommandArgs ( argc, argv );


    double t0=get_time_ms(), dt=0;
    int randvalseed = opt.getIntValue ( 'r', 1 );
    int ix = opt.getIntValue ( 'i', 1 );
    int r = opt.getIntValue ( 'r', 20 );
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

    arraylist_t ll= readarrayfile ( input );

    for ( size_t i=0; i<ll.size(); i++ ) {
        printf ( "reducing array %d\n", ( int ) i );
        array_link al = ll[i];
        
        al = al.randomrowperm();
        al = al.randomcolperm();
        
        int results = LMC0check ( al );

        
    }


    return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
