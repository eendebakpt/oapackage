/** \file oatest.cpp

C++ program: oatest

oatest: tool for testing new algorithms

Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

Algorithms: Alan Vazquez <alanrvazquez@gmail.com>

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
#include "mathtools.h"

#include "evenodd.h"
#include "lmc.h"

#include "conference.h"


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
        input="pexample_two.oa";

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

    printf("### read from file\n");

    arraylist_t ll= readarrayfile ( input );

    ll=arraylist_t();
    ll.push_back(exampleArray(30,1));

    for ( size_t i=0; i<ll.size(); i++ ) {
        array_link al = ll[i];

        //al = al.randomrowperm();
        //al = al.randomcolperm();
        al.showarray();
        lmc_t r = LMC0check ( al, 2);
        printf ( "array %d: result %d\n (should be %d)\n", (int) i, r, LMC_MORE );

        return 0;
        /* Apply random transformation */
        if (0) {
        conference_transformation_t T1(al);
        T1.randomize();
        T1.show();
        array_link al1 = T1.apply ( al );
        al1.showarray();
        lmc_t a = LMC0check ( al1 );
        printf ( "array %d: result %d\n (should be %d, or something else)\n", (int) i, a, LMC_LESS );
        }
    }

    return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
