
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

#include "graphtools.h"
#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "evenodd.h"
#include "oadevelop.h"
#include "lmc.h"

#include "conference.h"


#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>

#include "Deff.h"
#include "arraytools.h"
#include "strength.h"

#ifdef HAVE_BOOST
#include <string>
#include <boost/filesystem.hpp>
#endif

using namespace Eigen;

void speedcheck_conf ( const char *input, int verbose ) {
    double t0;
    arraylist_t ll= readarrayfile ( input );
    printf ( "### speedcheck: read from file (%d arrays)\n", ( int ) ll.size() );

    t0=get_time_ms();
    for ( size_t i=0; i<ll.size(); i++ ) {
        array_link al = ll[i];

        //al = al.randomrowperm();
        //al = al.randomcolperm();
        if ( verbose>=2 )
            al.showarray(); // print array
        lmc_t r = LMC0check ( al );
        if ( verbose>=2 )
            printf ( "array %d: result %d\n (should be %d)\n", ( int ) i, r, ( int ) LMC_MORE );
        if ( 0 ) {
            /* Apply random transformation */
            conference_transformation_t T1 ( al );
            T1.randomize();
            T1.show();
            array_link al1 = T1.apply ( al );
            if ( verbose>=2 ) {
                al1.showarray(); // Show transformed array
            }
            if ( 0 ) {
                lmc_t a = LMC0check ( al1 );
                if ( verbose )
                    printf ( "array %d: result %d\n (should be %d most of the time)\n", ( int ) i, a, ( int ) LMC_LESS );
            }
        }
    }
    printf ( "dt lmc0  %.3f \n", get_time_ms()-t0 );

    {
        double t0=get_time_ms();
        arraylist_t  out= selectConferenceIsomorpismClasses ( ll, 0,CONFERENCE_ISOMORPHISM );


        for ( size_t i=0; i<ll.size(); i++ ) {
            array_link al = ll[i];

        }
        printf ( "dt nauty %.3f \n", get_time_ms()-t0 );
    }


    return ;

}


//std::vector<double> Defficiencies (const array_link &al, const arraydata_t & arrayclass, int verbose ) ;

array_link finalcheck ( const array_link &al,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niterx, int nabortx ) {
    int nabort=al.n_rows*al.n_columns+2;
    int niter= nabort+5;
    array_link  al2 = optimDeff2level ( al, arrayclass, alpha, verbose>=3,  DOPTIM_UPDATE,  niter, nabort );


    std::vector<double> dd0 = al.Defficiencies();
    double d0 = scoreD ( dd0, alpha );

    std::vector<double> dd = al2.Defficiencies();
    double d = scoreD ( dd, alpha );

    if ( d>d0 ) {
        printf ( "finalcheck: %f -> %f\n", d0, d );
    }

    return al2;
}



#include "graphtools.h"






array_link array2xf2 ( const array_link &al ) {
    const int k = al.n_columns;
    const int n = al.n_rows;
    const int m = 1 + k + k* ( k-1 ) /2;
    array_link out ( n, m, array_link::INDEX_DEFAULT );

    // init first column
    int ww=0;
    for ( int r=0; r<n; ++r ) {
        out.array[r]=1;
    }

    // init array
    ww=1;
    for ( int c=0; c<k; ++c ) {
        int ci = c*n;
        array_t *pout = out.array+ ( ww+c ) *out.n_rows;
        for ( int r=0; r<n; ++r ) {
            pout[r] = 2*al.array[r+ci]-1;
        }
    }

    // init interactions
    ww=k+1;
    for ( int c=0; c<k; ++c ) {
        int ci = c*n+n;
        for ( int c2=0; c2<c; ++c2 ) {
            int ci2 = c2*n+n;

            const array_t * p1 = out.array+ci;
            const array_t * p2 = out.array+ci2;
            array_t *pout = out.array+ww*out.n_rows;

            for ( int r=0; r<n; ++r ) {
                pout[r] = - ( p1[r]*p2[r] )  ;
            }
            ww++;
        }
    }
    return out;
}

array_link sortrows ( const array_link &al ) {
    const int nc = al.n_columns;

    std::vector<mvalue_t<int> > rr ( al.n_rows );
    for ( int i=0; i<al.n_rows; i++ ) {
        mvalue_t<int> &m = rr[i];
        m.v.resize ( nc );

        for ( int k=0; k<nc; k++ ) {
            m.v[k]= al.atfast ( i, k );
        }
    }
    indexsort sorter ( rr );
    sorter.show();

    array_link out ( al.n_rows, al.n_columns, 0 );
    for ( int r=0; r<al.n_rows; r++ ) {
        for ( int i=0; i<al.n_columns; i++ ) {
            int newrow= sorter.indices[r];
            out._setvalue ( r,i,  al.atfast ( newrow,i ) );
        }
    }
    return out;
}

void paretoInfo ( const array_link & alx ) {
    std::vector < int >j5 = alx.Jcharacteristics ( 5 );
    int j5max = vectormax ( j5, 0 );

    int v1 = ( j5max == alx.n_rows );
    int v2 = 1 - v1;

    int N = alx.n_rows;
    int rank = array2xf ( alx ).rank();
    std::vector<int> F4 = alx.Fvalues ( 4 );
    std::vector<double> gwlp = alx.GWLP();
    printf ( "pareto data: %d ; ", rank );
    printf ( " %d ", ( int ) ( N*N*gwlp[4] ) );
    printf ( " ; " );
    display_vector ( F4 );
    printf ( " ; " );
    printf ( " %d ; %d", v1, v2 );
    printf ( "\n" );

}

/// check whether an array is contained in a Pareto set
int arrayInPareto ( const  Pareto < mvalue_t < long >, array_link >  &pset, const array_link &al, int verbose=1 ) {

    std::vector<array_link> llx = pset.allindices();
    arraylist_t ll ( llx.begin(), llx.end() );
    int jx=arrayInList ( al, ll, 0 );
    if ( verbose )
        myprintf ( "arrayInPareto: index in pareto list %d\n", jx );
    return jx;
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

    opt.addUsage ( "Orthonal Array: oatest: testing platform" );
    opt.addUsage ( "Usage: oatest [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.processCommandArgs ( argc, argv );


    double t0=get_time_ms(), dt=0;
    int randvalseed = opt.getIntValue ( 'r', 1 );
    int ix = opt.getIntValue ( 'i', 1 );
    int r = opt.getIntValue ( 'r', 0 );
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

    {
        speedcheck_conf ( input, verbose );
        exit ( 0 );
    }
    if ( 0 ) {
        int ei=26;

        array_link al = exampleArray ( ei,1 );
        arrayrankInfo ( array2xf ( al ) );
        exit ( 0 );

        rankStructure rs;
        rs.verbose=r;
        int r = array2xf ( al ).rank();
        int rc = rs.rankxf ( al );
        printf ( "rank of array %d: %d %d\n", ei, r, rc );
        exit ( 0 );
    }

    if ( 0 ) {
        array_link al= exampleArray ( 22,1 );
        al.showarray();
        symmdata sd ( al );
        sd.show();
        exit ( 0 );
    }

    {
        /// test performance if different rank algorithms

        arraylist_t lst = readarrayfile ( input );
        rankStructure rs;
        rs.verbose = r;
        int r, rc;
        printf ( "test singular values\n" );
        for ( int i=0; i< ( int ) lst.size(); i++ ) {
            array_link al = lst[i];
            if ( verbose>=2 )
                printf ( "-\n" );
            // arrayrankInfo(array2xf(al));
            // arrayrankInfo(array2secondorder(al));



            r = arrayrank ( array2xf ( al ) );
            //rc = arrayrank( array2secondorder( al ) ) + 1 + al.n_columns;
            //rc = array2secondorder ( al ).rank() + 1 + al.n_columns;
            rc = rs.rankxf ( al );
            //printf ( "r %d, rc %d\n", r, rc );
            // myassert ( r==rc, "rank calculations" );
        }
        printfd ( "done\n" );

        exit ( 0 );
    }
    {
        /// test performance if different rank algorithms

        arraylist_t lst = readarrayfile ( input );
        rankStructure rs;
        rs.verbose = r;
        int r, rc;
        printf ( "test performance\n" );
        for ( int i=0; i< ( int ) lst.size(); i++ ) {
            array_link al = lst[i];
            if ( verbose>=2 )
                printf ( "-\n" );

            switch ( jj ) {
            case 0:
                r=arrayrankColPivQR ( array2xf ( al ) );
                break;
            case 1:
                r=arrayrankFullPivQR ( array2xf ( al ) );
                break;
            case 2:
                r=arrayrankSVD ( array2xf ( al ) );
                break;
            case 3:
                r=arrayrankFullPivLU ( array2xf ( al ) );
                break;

            }

        }
        printfd ( "done\n" );

        exit ( 0 );
    }
    {
        for ( int i=0; i<27; i++ ) {
            array_link al = exampleArray ( i,0 );
            if ( al.n_columns<5 )
                continue;
            al = exampleArray ( i,1 );

            rankStructure rs;
            rs.verbose=r;
            int r = array2xf ( al ).rank();
            int rc = rs.rankxf ( al );
            if ( verbose>=2 ) {
                printf ( "rank of example array %d: %d %d\n", i, r, rc );
                if ( verbose>=3 ) {
                    al.showproperties();
                }
            }
            myassert ( r==rc, "rank calculations" );
        }
        exit ( 0 );
    }


    {

    }
    {
        pareto_cb_cache paretofunction =
            calculateArrayParetoJ5Cache < array_link >;

        Pareto < mvalue_t < long >, array_link > pset;
        std::vector<std::string> alist;
        alist.push_back ( "test.oa" );
        alist.push_back ( "p567.oa" );
        alist.push_back ( "p568.oa" );

        for ( int i=0; i< ( int ) alist.size(); i++ ) {
            std::string psourcefile = alist[i];


            //arrayfile_t afile ( psourcefile.c_str(), 0 );
            printf ( "### source file %s\n", psourcefile.c_str() );

            arraylist_t arraylist = readarrayfile ( psourcefile.c_str() );
            int apos = arrayInList ( exampleArray ( 24 ), arraylist );
            printf ( "  exampleArray 24: %d\n", apos );
            apos = arrayInList ( exampleArray ( 26 ), arraylist );
            printf ( "  exampleArray 26: %d\n", apos );

            for ( size_t ij=0; ij<arraylist.size(); ij++ ) {
                array_link alx = arraylist[ij];
                printf ( "array %d: ", ( int ) ij );
                paretoInfo ( alx );
            }
            fflush ( stdout );

            //gfx::timsort(arraylist.begin(), arraylist.end());	// sorting the arrays makes the rank calculations with subrank re-use more efficient
            if ( verbose >= 2 )
                printf ( "oatest: read arrays in file %s: %d arrays \n",psourcefile.c_str(), ( int ) arraylist.size() );
            addArraysToPareto ( pset, paretofunction, arraylist, jj, verbose );
            printf ( "  example 24: " );
            int t = arrayInPareto ( pset, exampleArray ( 24 ), 1 );
            printf ( "  example 26: " );
            t = arrayInPareto ( pset, exampleArray ( 26 ), 1 );
            pset.show ( 2 );
            exit ( 0 );
        }

        printf ( "### final result\n" );
        printf ( "example 24: " );
        int t = arrayInPareto ( pset, exampleArray ( 24 ), 1 );
        printf ( "example 26: " );
        t = arrayInPareto ( pset, exampleArray ( 26 ), 1 );

        exit ( 0 );
    }
    printf ( "----\n" );

    {

        arraylist_t ll = readarrayfile ( "x.oa" );
        array_link al = ll[0];

        if ( 1 ) {
            printf ( "test conf foldover\n" );
            for ( int i=0; i< ( int ) ll.size(); i++ ) {
                isConferenceFoldover ( ll[i] );

            }
            exit ( 0 );
        }
        //int ctype=2;
        conference_t ctype ( al.n_rows, al.n_rows );
        ctype.itype=CONFERENCE_RESTRICTED_ISOMORPHISM;
        ctype.ctype=conference_t::DCONFERENCE;

        CandidateGeneratorDouble cgenerator ( array_link() , ctype );
        cgenerator.verbose=verbose;

        for ( int i=0; i< ( int ) ll.size(); i++ ) {
            std::vector<cperm> cc = cgenerator.generateDoubleConfCandidates ( ll[i] );
            printfd ( "generated %d\n", cc.size() );
        }
        exit ( 0 );
    }


    return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
