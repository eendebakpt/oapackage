/*! \file Deff.cpp
 *  \brief Contains functions to optimize designs
 *
 */

#include <printfheader.h>

#include "arraytools.h"
#include "arrayproperties.h"

#ifdef DOOPENMP
#include "omp.h"
#endif


#include "Deff.h"

#ifndef myprintf
#define myprintf printf
#endif

#ifdef MAINMEX
#define MATLABUPDATE
#else
#ifdef MATLAB_MEX
#define MATLABUPDATE   mexEvalString ("drawnow" );
#else
#define MATLABUPDATE
#endif
#endif


DoptimReturn Doptimize ( const arraydata_t &arrayclass, int nrestartsmax, std::vector<double> alpha, int verbose, int method, int niter, double maxtime , int nabort )
{
    if ( method==DOPTIM_AUTOMATIC )
        method = DOPTIM_UPDATE;

    double t0 = get_time_ms();
    std::vector<std::vector<double> > dds;
    arraylist_t AA;

    bool abort=false;
    int nrestarts=0;

#ifdef DOOPENMP
    #pragma omp parallel for num_threads(4) schedule(dynamic,1)
#endif
    for ( int i=0; i<nrestartsmax; i++ ) {
        if ( abort )
            continue;

#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            if ( verbose && ( i%500==0 || i==nrestartsmax-1 ) ) {
                myprintf ( "Doptimize: iteration %d/%d\n", i, nrestartsmax );
                MATLABUPDATE	// helper macro from matlab interface
            }
        }

        array_link al = arrayclass.randomarray ( 1 );


        array_link  A = optimDeff ( al,  arrayclass, alpha, verbose>=2, method, niter,  nabort );
        std::vector<double> dd = A.Defficiencies();
        if ( verbose>=2 ) {
#ifdef DOOPENMP
            #pragma omp critical
#endif
            {
                double score = scoreD(dd, alpha);
                myprintf ( "Doptimize: iteration %d/%d: %f %f %f: score %.3f\n", i, nrestartsmax, dd[0], dd[1], dd[2], score );
            }

        }

#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            AA.push_back ( A );
            dds.push_back ( dd );
            nrestarts++;
        }

        if ( ( get_time_ms()-t0 ) > maxtime ) {
            if ( verbose )
                myprintf ( "max running time exceeded, aborting\n" );
            #pragma omp critical
            {
                abort=true;
            }
        }
    }

    // loop is complete
    DoptimReturn a = {dds, AA, nrestarts, nrestarts};
    return a;
}


DoptimReturn DoptimizeDebug ( const arraydata_t &arrayclass, int nrestartsmax, std::vector<double> alpha, int verbose, int method, int niter, double maxtime , int nabort )
{
    if ( method==DOPTIM_AUTOMATIC )
        method = DOPTIM_UPDATE;

    double t0 = get_time_ms();
    std::vector<std::vector<double> > dds;
    arraylist_t AA;

    bool abort=false;
    int nrestarts=0;

#ifdef DOOPENMP
    #pragma omp parallel for num_threads(4) schedule(dynamic,1)
#endif
    for ( int i=0; i<nrestartsmax; i++ ) {
        if ( abort )
            continue;

#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            if ( verbose && ( i%500==0 || i==nrestartsmax-1 ) ) {
                myprintf ( "Doptimize: iteration %d/%d\n", i, nrestartsmax );
                MATLABUPDATE	// helper macro from matlab interface
            }
        }

        array_link al = arrayclass.randomarray ( 1 );


        array_link  A = optimDeff2level ( al,  arrayclass, alpha, verbose>=2, method, niter,  nabort );
        std::vector<double> dd = A.Defficiencies();
        if ( verbose>=2 ) {
#ifdef DOOPENMP
            #pragma omp critical
#endif
            {
                double score = scoreD(dd, alpha);
                myprintf ( "Doptimize: iteration %d/%d: %f %f %f: score %.3f\n", i, nrestartsmax, dd[0], dd[1], dd[2], score );
            }

        }

#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            AA.push_back ( A );
            dds.push_back ( dd );
            nrestarts++;
        }

        if ( ( get_time_ms()-t0 ) > maxtime ) {
            if ( verbose )
                myprintf ( "max running time exceeded, aborting\n" );
            #pragma omp critical
            {
                abort=true;
            }
        }
    }

    // loop is complete
    DoptimReturn a = {dds, AA, nrestarts, nrestarts};
    return a;
}

double scoreD ( const std::vector<double> dd, const std::vector<double> alpha )
{
    double v=0;
    for ( size_t i=0; i<dd.size();  i++ )
        v+= dd[i]*alpha[i];
    return v;

}

DoptimReturn DoptimizeMixed ( const arraylist_t &sols, const arraydata_t &arrayclass, const std::vector<double> alpha, int verbose, int nabort )
{
    const size_t nn = sols.size();
    double t0 = get_time_ms();
    std::vector<std::vector<double> > dds ( nn );
    arraylist_t AA ( nn );

    bool abort=false;
    int nimproved=0;

    int method1 = DOPTIM_SWAP;
    int method2 = DOPTIM_UPDATE;
    const int niter=600000;
    if ( nabort<0 )
        nabort = arrayclass.N*arrayclass.ncols*20+500;

    if ( verbose>=3 )
        myprintf ( "DoptimizeMixed: nabort %d\n", nabort );
#ifdef DOOPENMP
    #pragma omp parallel for num_threads(4) schedule(dynamic,1)
#endif
    for ( size_t i=0; i<nn; i++ ) {
        if ( abort )
            continue;

#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            if ( verbose && ( i%100==0 || i==nn-1 ) ) {
                myprintf ( "DoptimizeMixed: iteration %d/%d\n", ( int ) i, ( int ) nn );
                MATLABUPDATE	// helper macro from matlab interface
            }
        }

        const array_link &al = sols[i];
        double score0 = scoreD ( al.Defficiencies(), alpha );

        array_link  alu = optimDeff ( al,  arrayclass, alpha, verbose>=3, method1 , niter,  nabort );
        double score1 = scoreD ( alu.Defficiencies(), alpha );

        array_link  alu2 = optimDeff ( alu,  arrayclass, alpha, verbose>=3, method2 , niter,  0 );
        double score2 = scoreD ( alu2.Defficiencies(), alpha );


        //#pragma omp critical
        dds[i]=alu2.Defficiencies();
        AA[i]=alu2;

        if ( score2>score0 ) {
            if ( verbose>=2 )
                myprintf ( "DoptimizeMixed: array %d/%d: improve %.6f -> %.6f -> %.6f\n", ( int ) i, ( int ) nn, score0, score1, score2 );

            nimproved=nimproved+1;
        }

        //if ( verbose>=2 ) {
        //	myprintf ( "DoptimizeMixed: iteration %d/%d: %f %f %f\n", i, nn, dd[0], dd[1], dd[2] );
        //}
    }

    double dt = get_time_ms()-t0;
    if ( verbose ) {
        myprintf ( "DoptimizeMixed: improved %d/%ld arrays, %.2f [s]\n", nimproved, (long)nn, dt );
    }

    // loop is complete
    DoptimReturn a = {dds, AA, ( int ) nn, ( int ) nimproved};
    return a;
}


array_link  optimDeff ( const array_link &A0,  const arraydata_t &arrayclass, const std::vector<double> alpha, int verbose, int optimmethod, int niter, int nabort )
{
    const int N = arrayclass.N;
    const int k = arrayclass.ncols;


    if ( nabort<=0 ) {
        if ( arrayclass.is2level() ) {
            nabort= ( N*k ) *2.5+1;	// factor 2 to compensate for equal error switches
        } else {
            nabort= ( N*k ) * ( arrayclass.s[0] ) +2;
        }
    }

    int nx=0;
    //#get factor levels
    std::vector<int> s = arrayclass.getS();

    if ( optimmethod==DOPTIM_UPDATE ) {
        if ( arrayclass.is2level() )
            optimmethod=DOPTIM_FLIP;
    }
    array_link A = A0;
    symmetry_group sg=symmetry_group ( s );

    std::vector<int> gidx = sg.gidx;

    std::vector<double> dd0 = A0.Defficiencies();

    if ( verbose ) {
        myprintf ( "optimDeff: initial D-efficiency %.4f\n",  dd0[0] );
    }

    // initialize score
    double d = scoreD ( dd0, alpha );

    int lc=0;	// index of last change to array

    // initialize arary with random permutation
    const int nn = N*k;
    std::vector<int> updatepos = permutation<int> ( nn );
    my_random_shuffle ( updatepos.begin(), updatepos.end() );
    int updateidx = 0;


//#pragma omp for
    for ( int ii=0; ii<niter; ii++ ) {
        // select random row and column
        int r = updatepos[updateidx] % N;
        int c = updatepos[updateidx] / N;
        updateidx= ( updateidx+1 ) % ( nn );

        int r2 = fastrandK ( N );

        //r=rand()%N; c=rand()%k; r2=rand()%N;

        // make sure column is proper column group
        int c2 = sg.gstart[sg.gidx[c]] + fastrandK ( sg.gsize[gidx[c]] );

        // get values
        array_t o = A._at ( r,c );
        array_t o2 = A._at ( r2,c2 ); // no extra error checking

        // update

        if ( optimmethod==DOPTIM_SWAP && o==o2 )
            continue;

        switch ( optimmethod ) {
        case DOPTIM_SWAP: // swap
            A._setvalue ( r,c,o2 );
            A._setvalue ( r2,c2,o );
            break;
        case DOPTIM_UPDATE: { // random update
            int val = fastrandK ( s[c] );
            A._setvalue ( r,c,val );
            break;
        }
        case DOPTIM_FLIP: // flip
            A._setvalue ( r,c,1-o );
            break;
        case DOPTIM_NONE:
            break;
        default
                :
            myprintf ( "optimDeff: no such optimization method\n" );
            break;
        }
        // evaluate
        std::vector<double> dd = A.Defficiencies();
        nx++;
        double dn = scoreD ( dd, alpha );


        if ( verbose>=4 ) {
            myprintf ( "  optimDeff: switch %d, %d: %d with %d, %d: %d: d %.4f -> dn %.4f \n", r, c, r2, c2,o, o2, d, dn );
        }

        // switch back if necessary

        if ( dn>=d ) {
            if ( dn>d )  {
                lc=ii;
                //std::random_shuffle ( updatepos.begin(), updatepos.end() );	// update the sequence of positions to try

                if ( verbose>=3 )
                    myprintf ( "optimDeff: ii %d: %.6f -> %.6f\n", ii, d, dn );
                d=dn;
            } else {
                if ( verbose>=2 ) {
                    myprintf ( "optimDeff: equal ii %d: %.6f -> %.6f\n", ii, d, dn );
                }
            }
        } else {

            if ( ( dn+1e-10 ) >=d ) {
                if ( verbose>=2 ) {
                    myprintf ( "element within 1e-12! d %f, dn %f\n", d, dn );
                }
            }
            // restore to original
            switch ( optimmethod ) {
            case DOPTIM_SWAP:
                A._setvalue ( r,c,o );
                A._setvalue ( r2,c2,o2 );
                break;
            case DOPTIM_FLIP:
            case DOPTIM_UPDATE:
                A._setvalue ( r,c,o );
                break;
            case DOPTIM_NONE:
                break;
            }
        }

        // check progress
        if ( ( ii-lc ) >nabort ) {
            if ( verbose>=2 )
                myprintf ( "optimDeff: early abort ii %d, lc %d: %d\n", ii, lc, ( ii-lc ) );
            //abort=true;
            break;
        }

    }

    std::vector<double> dd = A.Defficiencies();
    double dn = scoreD ( dd, alpha );

    if ( verbose ) {
        myprintf ( "optimDeff: final score %.4f, final D-efficiency %.4f\n",  dn, dd[0] );
    }

    return A;
}


array_link  optimDeff2level ( const array_link &A0,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niter, int nabort )
{
    int nx=0;
    int N = arrayclass.N;
    int k = arrayclass.ncols;
    std::vector<int> s = arrayclass.getS();

    if ( ! arrayclass.is2level() ) {
        myprintf ( "optimDeff2level: error arrayclass is not 2-level\n" );
        return A0;
    }

    if ( nabort<=0 ) {
        nabort=N*k+2;
    }

    if ( optimmethod==DOPTIM_UPDATE ) {
        if ( arrayclass.is2level() )
            optimmethod=DOPTIM_FLIP;
    }
    array_link A = A0;

    std::vector<double> dd0 = A0.Defficiencies();

    if ( verbose ) {
        myprintf ( "optimDeff: initial D-efficiency %.4f\n",  dd0[0] );
    }

    // initialize arary with random permutation
    std::vector<int> updatepos = permutation<int> ( N*k );
    my_random_shuffle ( updatepos.begin(), updatepos.end() );
    int updateidx = 0;

    // initialize score
    double d = scoreD ( dd0, alpha );

    int lc=0;	// index of last change to array

    int nn = updatepos.size();

//#pragma omp for
    for ( int ii=0; ii<niter; ii++ ) {

        // select random row and column
        //int r = updatepos[updateidx] % N; int c = updatepos[updateidx] / N;
        int r = fastrandK ( N );
        int c = fastrandK ( k );
        updateidx= ( updateidx+1 ) % ( nn );

        int r2 = fastrandK ( N );
        int c2 = fastrandK ( k );

        // get values
        array_t o = A._at ( r,c );
        array_t o2 = A._at ( r2,c2 ); // no extra error checking

        // update

        if ( optimmethod==DOPTIM_SWAP && o==o2 )
            continue;

        switch ( optimmethod ) {
        case DOPTIM_SWAP: // swap
            A._setvalue ( r,c,o2 );
            A._setvalue ( r2,c2,o );
            break;
        case DOPTIM_UPDATE: { // random update
            int val = fastrandK ( s[c] );
            A._setvalue ( r,c,val );
            break;
        }
        case DOPTIM_FLIP: // flip
            A._setvalue ( r,c,1-o );
            break;
        case DOPTIM_NONE:
            break;
        default
                :
            myprintf ( "optimDeff: no such optimization method\n" );
            break;
        }
        // evaluate
        std::vector<double> dd = A.Defficiencies();
        nx++;
        double dn = scoreD ( dd, alpha );


        if ( verbose>=4 ) {
            myprintf ( "  optimDeff: switch %d, %d: %d with %d, %d: %d: d %.4f -> dn %.4f \n", r, c, r2, c2,o, o2, d, dn );
        }

        // switch back if necessary

        if ( dn>=d ) {
            if ( dn>d )  {
                lc=ii;
                if ( verbose>=3 )
                    myprintf ( "optimDeff: ii %d: %.6f -> %.6f\n", ii, d, dn );
                d=dn;
            } else {
                myprintf ( "optimDeff: equal ii %d: %.6f -> %.6f\n", ii, d, dn );
            }
        } else {
            // restore to original
            switch ( optimmethod ) {
            case DOPTIM_SWAP:
                A._setvalue ( r,c,o );
                A._setvalue ( r2,c2,o2 );
                break;
            case DOPTIM_FLIP:
            case DOPTIM_UPDATE:
                A._setvalue ( r,c,o );
                break;
            case DOPTIM_NONE:
                break;
            }
        }

        // check progress
        if ( ( ii-lc ) >nabort ) {
            if ( verbose>=2 )
                myprintf ( "optimDeff: early abort ii %d, lc %d: %d\n", ii, lc, ( ii-lc ) );
            //abort=true;
            break;
        }

    }

    std::vector<double> dd = A.Defficiencies();
    double dn = scoreD ( dd, alpha );

    if ( verbose ) {
        myprintf ( "optimDeff: final score %.4f, final D-efficiency %.4f\n",  dn, dd[0] );
    }

    return A;
}

#include <algorithm>
#include <cmath>

const double NaN = std::numeric_limits<double>::quiet_NaN();

extern "C" {

    void DefficienciesR(int *N, int *k, double *input,  double *D, double *Ds, double *D1 )
    {
        array_link al(*N, *k, array_link::INDEX_DEFAULT);
        std::copy ( input, input+(*N)*(*k), al.array );

        if (al.min()<0) {
            *D = NaN;
            *Ds = NaN;
            *D1=NaN;
            myprintf("DefficienciesR: error: design should have input elements from 0 to s-1\n");
        }

        std::vector<double> dd = al.Defficiencies();

        *D = dd[0];
        *Ds=dd[1];
        *D1=dd[2];
        return;
    }

    double DoptimizeR ( int *pN, int *pk, int *nrestarts, double *alpha1, double *alpha2, double *alpha3, int *_verbose, int *pointer_method, int *_niter, double *maxtime , int *nabort, double *output )
    {
        int niter=*_niter;
        int method=*pointer_method;
        int N = *pN;
        int k = *pk;
        int verbose = *_verbose;

        output[0]=1;
        output[1]=2;
        output[2]=3;
        output[3]=4;

        if ( verbose>=2 )
            myprintf ( "DoptimizeR: N %d, k %d, nrestarts %d, niter %d, alpha1 %f\n", N, k, *nrestarts, niter, *alpha1 );

        arraydata_t arrayclass ( 2, N, 0, k );
        std::vector<double> alpha ( 3 );
        alpha[0]=std::max ( *alpha1,0. );
        alpha[1]=std::max ( *alpha2, 0. );
        alpha[2]=std::max ( *alpha3,0. );

        DoptimReturn rr = Doptimize ( arrayclass, *nrestarts, alpha,  verbose,  method, niter, *maxtime,  *nabort );

        std::vector<std::vector<double> > dds = rr.dds;
        arraylist_t AA = rr.designs;

        // sort according to values
        std::vector<double> sval ( AA.size() );
        for ( size_t i=0; i<AA.size(); i++ ) {
            sval[i]=-scoreD ( dds[i], alpha );
        }

        indexsort sorter ( sval );

        AA=sorter.sorted ( AA );
        dds=sorter.sorted ( dds );

        array_link best = AA[0];
        std::vector<double> dd = best.Defficiencies();

        std::copy ( best.array, best.array+N*k, output );

        if (verbose>=1) {
            myprintf("Doptimize: generated design with D = %.6f, Ds = %.6f\n",dd[0], dd[1] );
        }
        if ( verbose>=2 ) {
            myprintf ( "DoptimizeR: done\n" );
        }

        return best.Defficiency();
    }

} // extern "C"

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
