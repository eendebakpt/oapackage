/** \file oaconference.cpp

 C++ program: oaconference

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"


#include "conference.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <string>
using namespace std;

void print_all_permutations ( const string& s ) {
    string s1 = s;
    sort ( s1.begin(), s1.end() );
    int n=0;
    do {
        cout << s1 << endl;
        n++;
    } while ( next_permutation ( s1.begin(), s1.end() ) );
    printf ( "%d/%ld perms (len %ld)\n", n, factorial<long> ( s1.size() ), ( long ) s1.size() );
}

void testx() {
    print_all_permutations ( "001111112222222" );
    exit ( 0 );
}

template<class fwditer>
fwditer random_unique ( fwditer begin, fwditer end, size_t num_random ) {
    size_t left = std::distance ( begin, end );
    while ( num_random-- ) {
        fwditer r = begin;
        std::advance ( r, rand() %left );
        std::swap ( *begin, *r );
        ++begin;
        --left;
    }
    return begin;
}

/// generate a range of numbers
struct rangegenerator {
    rangegenerator ( int init ) : start ( init ) { }

    int operator() () {
        return start++;
    }

    int start;
};

#include <algorithm>

template <typename T, typename T2>
/// return subvector of a vector specified by indices
T2 subvector ( const T2& full, const T& ind ) {
    int num_indices = ind.size();
    T2 target ( num_indices );
    for ( int i = 0; i < num_indices; i++ ) {
        target[i] = full[ind[i]];
    }
    return target;
}

enum reduction_method {NONE, NAUTY, LMC0};

enum max_selection_method {SELECT_RANDOM, SELECT_FIRST};

/** select a subset of arrays
 * 
 * 
 * \param lst Input list of arrays
 * \param nmax The maximum size of the subset
 * \param method Method for selecting the subset
 */
arraylist_t selectArraysMax ( const arraylist_t &lst, int nmax, max_selection_method method = SELECT_RANDOM, int verbose=0 ) {
    nmax = std::min ( nmax, ( int ) lst.size() );

    if ( nmax<0 )
        return lst;

    if ( method==SELECT_FIRST ) {
        arraylist_t out ( lst.begin(), lst.begin() +nmax );
        return out;

    }
    std::vector<size_t> indices ( lst.size() );
    generate ( indices.begin(), indices.end(), rangegenerator ( 0 ) );
    std::random_shuffle ( indices.begin(), indices.end() ); // TODO: more efficient with Fischer-Yates shuffle
    indices.resize ( nmax );
    std::sort ( indices.begin(), indices.end() );

    if ( verbose ) {
        printfd ( "indices: generated " );
        display_vector ( indices );
        printf ( "\n" );
    }
    arraylist_t out = subvector ( lst, indices );

    return out;
}

int main ( int argc, char* argv[] ) {
    AnyOption opt;
    /* parse command line options */
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "output", 'o' );
    opt.setOption ( "verbose", 'v' );
    opt.setOption ( "randseed", 'r' );
    opt.setOption ( "kmax", 'k' );
    opt.setOption ( "itype" );
    opt.setOption ( "j1zero" );
    opt.setOption ( "j3zero" );
    opt.setOption ( "rows", 'N' );
    opt.setOption ( "format", 'f' );
    opt.setOption ( "nmax" );
    opt.setOption ( "nmaxmethod" );
    opt.setOption ( "select", 's' );
    opt.setOption ( "cols" );
    opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

    opt.setOption ( "input", 'i' );
    opt.setOption ( "output", 'o' );
    opt.setOption ( "ctype" );

    opt.setOption ( "reduction" );
    opt.setOption ( "debug" );

    opt.addUsage ( "Orthonal Array: oaconference: testing platform" );
    opt.addUsage ( "Usage: oaconference [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.addUsage ( " -N --rows  			Number of rows " );
    opt.addUsage ( " -k --kmax  			Max number of columns " );
    opt.addUsage ( " --nmax  			Max number of designs to sample for the next stage " );
    opt.addUsage ( printfstring ( " --nmaxmethod [INT]         Method (FIRST %d, RANDOM %d) ", SELECT_FIRST, SELECT_RANDOM ).c_str() );
    opt.addUsage ( " -v --verbose [INT] 			Verbosity level " );
    opt.addUsage ( " -i --input [FILENAME]		Input file to use" );
    opt.addUsage ( " -o --output [FILEBASE]		Output file to use" );
    opt.addUsage ( " -f [FORMAT]                                 Output format (TEXT, BINARY (default), D (binary difference) ) " );
    opt.addUsage ( printfstring ( " --ctype [TYPE]				Zero for normal type, %d for diagonal, %d for double conference type ",conference_t::CONFERENCE_DIAGONAL, conference_t::DCONFERENCE ).c_str() );
    opt.addUsage ( printfstring ( " --itype [TYPE]				Matrix isomorphism type (CONFERENCE_ISOMORPHISM %d, CONFERENCE_RESTRICTED_ISOMORPHISM %d)", CONFERENCE_ISOMORPHISM, CONFERENCE_RESTRICTED_ISOMORPHISM ).c_str() );
    opt.addUsage ( " --j1zero [INT]             Restrict designs to J1=0" );
    opt.addUsage ( " --j3zero [INT]		Restrict designs to J3=0" );
    opt.addUsage ( printfstring ( " --select [INT]             Method to reduce generated designs (None %d, NAUTY %d, LMC0 %d)", NONE, NAUTY, LMC0 ).c_str() );
    opt.processCommandArgs ( argc, argv );


    //testx();

    print_copyright();
    //cout << system_uname();
    setloglevel ( NORMAL );

    std::string format = opt.getStringValue ( 'f', "AUTO" );
    arrayfile::arrayfilemode_t file_mode_base = arrayfile_t::parseModeString ( format );
    
    int debug = opt.getIntValue ( "debug", 0 );
    int r = opt.getIntValue ( 'r', 1 );
    int verbose = opt.getIntValue ( 'v', 1 );
    int N = opt.getIntValue ( 'N', 10 );
    int kmax = opt.getIntValue ( 'k', N );
    int nmax = opt.getIntValue ( "nmax", -1 );
    max_selection_method nmaxmethod = ( max_selection_method ) opt.getIntValue ( "nmaxmethod", SELECT_FIRST );
    reduction_method select = ( reduction_method ) opt.getIntValue ( 's', 1 );
    const conference_t::conference_type ctx = ( conference_t::conference_type ) opt.getIntValue ( "ctype", 0 );
    const matrix_isomorphism_t itype = ( matrix_isomorphism_t ) opt.getIntValue ( "itype", CONFERENCE_ISOMORPHISM );
    const int j1zero = opt.getIntValue ( "j1zero", 0 );
    const int j3zero = opt.getIntValue ( "j3zero", 0 );
    const reduction_method reduction = ( reduction_method ) opt.getIntValue ( "reduction", NAUTY );

    const std::string output = opt.getStringValue ( 'o', "" );
    const std::string input = opt.getStringValue ( 'i', "" );

    if ( r==-1 ) {
        int randvalseed=time ( NULL );
        printf ( "random seed %d\n", randvalseed );
        srand ( randvalseed );
    } else
        srand ( r );

    /* parse options */
    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
        opt.printUsage();
        exit ( 0 );
    }

    setloglevel ( SYSTEM );

    int kstart=-1;
    conference_t ctype ( N, N, 0 );

    arraylist_t inputarrays;
    if ( input.length() >1 ) {
        inputarrays = readarrayfile ( input.c_str() )	 ;
        //al=kk[0];
        kstart=inputarrays[0].n_columns;

        if ( inputarrays.size() >0 ) {
            ctype = conference_t ( N, kstart, j1zero );
            assert ( inputarrays[0].n_rows==N );
        }
        ctype.ctype=ctx;
        ctype.itype=itype;
        ctype.j3zero=j3zero;
        ctype.j1zero=j1zero;
        //kmax=inputarrays[0].n_columns+1;

    } else {
        ctype.ctype=ctx;
        ctype.itype=itype;
        ctype.j3zero=j3zero;
        ctype.j1zero=j1zero;

        ctype.addRootArrays ( inputarrays );
        kstart=inputarrays[0].n_columns;
    }

    if ( ctype.ctype==conference_t::DCONFERENCE ) {
        if ( ctype.j1zero==1 )
            kmax=std::min ( ( int ) ( ceil ( N/2 ) +2 ), kmax );
    }

    if ( verbose ) {
        printf ( "oaconference: extend %d conference matrices to size %dx%d (itype %d (CONFERENCE_ISOMORPHISM %d, CONFERENCE_RESTRICTED_ISOMORPHISM %d) )\n", ( int ) inputarrays.size(), ctype.N, ctype.ncols, itype, CONFERENCE_ISOMORPHISM, CONFERENCE_RESTRICTED_ISOMORPHISM );
        printf ( "oaconference: ctype %d (DCONFERENCE %d, CONFERENCE_NORMAL %d) )\n", ctype.ctype, conference_t::DCONFERENCE, conference_t::CONFERENCE_NORMAL );
    }
    if ( verbose>=3 ) {
        printf ( "--- initial set of arrays ---\n" );
        showArrayList ( inputarrays );
    }

    double t0 = get_time_ms();

    for ( int extcol=kstart; extcol<kmax; extcol++ ) {
        printf ( "\n### oaconference: extend column %d (max number of columns %d, time %.1f [s])\n", extcol, kmax, get_time_ms ( t0 ) );

        arraylist_t outlist;

        switch ( ctype.ctype ) {
        case conference_t::DCONFERENCE: {
            outlist = extend_double_conference ( inputarrays, ctype,  verbose );
            sort ( outlist.begin(), outlist.end(), compareLMC0 );
            break;
        }
        case conference_t::CONFERENCE_NORMAL:
        case conference_t::CONFERENCE_DIAGONAL:
            switch ( ctype.itype ) {
            case CONFERENCE_RESTRICTED_ISOMORPHISM:
                outlist = extend_conference_restricted ( inputarrays, ctype,  verbose );
                break;

            case CONFERENCE_ISOMORPHISM:
                // TODO: make a version with symmetry inflation
                if ( debug ) {
                    printfd ( "debug option set: using extend_conference_plain\n" );
                    outlist = extend_conference_plain ( inputarrays, ctype,  verbose, select==NAUTY );
                } else
                    outlist = extend_conference ( inputarrays, ctype,  verbose, select==NAUTY );
                break;
            default
                    :
                printfd ( "isomorphism type not implemented" );
                break;
            }
            break;
        default
                :
            printfd ( "not implemented: itype %d\n", ctype.itype );
            break;
        }

        if ( select != NONE ) {
            switch ( select ) {
            case NAUTY:
                outlist = selectConferenceIsomorpismClasses ( outlist, verbose, ctype.itype );
                break;
            case LMC0:
                outlist = selectLMC0 ( outlist, verbose>=1, ctype );
                break;
            default:
                printfd ( "error: selection method %d not implemented\n", ( int ) select );
                break;
            }
        }
        sort ( outlist.begin(), outlist.end(), compareLMC0 );

        if ( ctype.ctype==conference_t::DCONFERENCE ) {
            long nevenodd=0;
            for ( size_t i=0; i<outlist.size(); i++ ) {
                if ( ! isConferenceFoldover ( outlist[i] ) ) {
                    if ( verbose>=2 )
                        printfd ( "#### found an even-odd conference matrix!!!\n" );
                    nevenodd++;
                }
            }
            printfd ( "## found %ld/%ld designs to be even-odd\n", nevenodd, outlist.size() );
        }

        if ( output.length() >=1 ) {
            arrayfilemode_t file_mode = file_mode_base;
            if ( file_mode==arrayfile::A_AUTOMATIC ) {
                if ( outlist.size() < 1000 ) {
                    // small files in text format for easy reading
                    file_mode= arrayfile::ATEXT;
                } else {
                    // larger files in binary format
                    if ( outlist.size() < 5000 ) {
                        file_mode= arrayfile::ABINARY;
                    } else {
                        file_mode= arrayfile::ABINARY_DIFF;
                    }
                }
            }
            //printfd("file_mode: %d, size %d\n", file_mode, outlist.size() );

            std::string outfile = output + printfstring ( "-%d-%d", ctype.N, extcol+1 )  + ".oa";
            if ( verbose )
                printf ( "oaconference: write %d arrays to file %s...\n", ( int ) outlist.size(), outfile.c_str() );

            writearrayfile ( outfile.c_str(),outlist, file_mode, N, extcol+1 );

        }

        printf ( "oaconference: extend column %d: generated %d non-isomorphic arrays (%.1f [s])\n", extcol, ( int ) outlist.size(),  get_time_ms ( t0 ) );

        if ( nmax>=0 ) {
            // select arrays
            outlist = selectArraysMax ( outlist, nmax, nmaxmethod );
        } else {
        }

        // loop
        inputarrays=outlist;
    }
    printf ( "done...\n" );

    return 0;

}


// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
