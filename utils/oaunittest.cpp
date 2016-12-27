
/** \file oaunittest.cpp

C++ program: oaunittest

oaunittest: run some tests on the code

Author: Pieter Eendebak <pieter.eendebak@gmail.com>
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
#include "Deff.h"
#include "conference.h"

#include "graphtools.h"

#include "Eigen/Dense"

#ifdef HAVE_BOOST
#include <string>
#include <boost/filesystem.hpp>
#endif

template <class Type>
void permute ( Type *source, Type *target, const std::vector<int> p ) {
    for ( size_t i=0; i<p.size(); i++ )
        target[i]=source[p[i]];
}

/// check transformation inverse. return 0 if test is good
int checkConferenceInverse ( const array_link &al ) {
    conference_transformation_t T1 ( al );
    T1.randomize();

    conference_transformation_t T1i = T1.inverse();
    conference_transformation_t II = T1i * T1;

    myassert ( II.isIdentity(), "unittest error: inverse of conference matrix transformation\n" );
    myassert ( al == II.apply ( al ), "unittest error: inverse of conference matrix transformation\n" );

    return 0;
}

/// check composition operator. returns 0 if test id good
int checkConferenceComposition ( const array_link &al, int verbose=0 ) {
    conference_transformation_t T1 ( al );
    T1.randomize();

    conference_transformation_t T2 ( al );
    T2.randomize();

    conference_transformation_t T3 = T2 * T1;

    array_link al1 = T1.apply ( al );
    array_link al1t2 = T2.apply ( al1 );
    array_link al3 = T3.apply ( al );

    if ( verbose ) {
        printfd ( "checkTransformationComposition: transforms\n" );
        T1.show();
        T2.show();
        T3.show();
        printfd ( "checkTransformationComposition: arrays\n" );
        al.showarray();
        al1.showarray();
        al1t2.showarray();
        al3.showarray();
    }

    myassert ( al3==al1t2, "unittest error: composition of conference transformations\n" );

    return 0;
}

void checkGenerators(int verbose = 1)
{
    int N = 20;
    conference_t ct ( N, 4, 1 );
          ct.ctype=conference_t::DCONFERENCE;
          ct.j3zero=1;

          array_link al = exampleArray(35, 1);
          //al.showarray();

          CandidateGeneratorDouble cgenerator ( array_link(), ct );
          cgenerator.verbose=0;
          for ( int i=0; i<2; i++ ) {
        {
                const std::vector<cperm> &cl = cgenerator.generateCandidates ( al );
                    myassert ( cl.size()==3, "unittest error: inverse of array transformation\n" );

                    if (verbose>=2)
                    {
                        printfd ( "generated %d\n", cl.size() );
                    cgenerator.showCandidates ( 2 );
        }
               }
          }
          
          
          // generator for conference matrices
          {
              const int example_idx=39;
                     array_link al = exampleArray ( example_idx, 1 );
          conference_t ct ( al.n_rows, al.n_columns+4, 0 );

          if ( verbose>=2 ) {
              myprintf("checking generator on array:\n");
               al.showarray();
          }
          myassert ( al.is_conference() );
          myassert ( al.min() ==-1 );

          int filterj2=1;
          int filtersymminline=1;
          int averbose=verbose;
                  std::vector<cperm>      ccX = generateSingleConferenceExtensions ( al, ct, -1, averbose, 1, filterj2, ct.j3zero, filtersymminline );
                  if (verbose>=2) {
                      
showCandidates(ccX);
printf("\n-----------\n");
                  }
                  myassert(ccX.size() == 2, "number of candidnates generated");
                  {
          CandidateGenerator cgenerator ( array_link(), ct );
          int kz = maxz ( al ) +1;
          cgenerator.verbose=verbose;
          std::vector<cperm> ee = cgenerator.generateCandidatesZero ( al, kz );
          printf("ee.size() %d\n", (int)ee.size() );
          myassert(ee.size()==1, "number of candidnates generated");
          if (verbose>=2) {
          cgenerator.showCandidates(2);
          printf ( "generateCandidatesZero: %d\n-------------\n", ( int ) ee.size() );
          }
          }
          }
}

/// check transformation inverse. return 0 if test is good
int checkTransformationInverse ( const array_link &al ) {
    arraydata_t adataX = arraylink2arraydata ( al );
    array_transformation_t T1 ( &adataX );
    T1.randomize();

    array_transformation_t T1i = T1.inverse();
    array_transformation_t II = T1i * T1;

    myassert ( II.isIdentity(), "unittest error: inverse of array transformation\n" );

    return 0;
}
/// check composition operator. returns 0 if test id good
int checkTransformationComposition ( const array_link &al, int verbose=0 ) {
    arraydata_t adataX = arraylink2arraydata ( al );

    array_transformation_t T1 ( &adataX );
    T1.randomize();
    array_transformation_t T2 ( &adataX );
    T2.randomize();
    array_transformation_t T3;
    T3= T1 * T2;

    //T2.reset(); T2.randomize();
    array_link al2 = T2.apply ( al );

    //T1.show(); al2.showarray(); exit(0);

    array_link al12 = T1.apply ( al2 );
    array_link al3 = T3.apply ( al );

    if ( verbose ) {
        printfd ( "checkTransformationComposition\n" );
        T1.show();
        T2.show();
        T3.show();
        al.showarray();
        al2.showarray();
        al12.showarray();
        al3.showarray();
    }

    myassert ( al3==al12, "unittest error: composition of array transformations\n" );

    return 0;
}

enum {UGOOD, UERROR};

/** unittest for oapackage
 *
 * Returns UGOOD if all tests are ok.
 *
 */
int oaunittest ( int verbose, int writetests=0, int randval = 0 ) {
    double t0=get_time_ms();
    const char *bstr = "OA unittest";
    cprintf ( verbose, "%s: start\n", bstr );

    srand ( randval );

    int allgood=UGOOD;

    initncombscache ( 20 );

    /* constructors */
    {
        cprintf ( verbose,"%s: interaction matrices\n", bstr );

        array_link al = exampleArray ( 2 );
        Eigen::MatrixXd m1 = array2xfeigen ( al );
        Eigen::MatrixXd m2 = arraylink2eigen ( array2xf ( al ) );

        Eigen::MatrixXd dm = m1-m2;
        int sum = dm.sum();

        myassert ( sum==0, "unittest error: construction of interaction matrices\n" );


    }

#ifdef OADEV
    /* conference matrices */
    {
        cprintf ( verbose,"%s: conference matrices\n", bstr );

        int N=4;
        conference_t ctype ( N, N, 0 );

        arraylist_t kk;
        array_link al = ctype.create_root();
        kk.push_back ( al );

        for ( int extcol=2; extcol<N; extcol++ ) {
            kk = extend_conference ( kk, ctype,  0 );
        }
        myassert ( kk.size() ==1, "unittest error: conference matrices for N=4\n" );
    }
    
        {
        cprintf ( verbose,"%s: generators for conference matrix extensions\n", bstr );
        checkGenerators(verbose);
    }

    {
        cprintf ( verbose,"%s: conference matrix Fvalues\n", bstr );
        array_link  al = exampleArray ( 22, 0 );
        if ( verbose>=2 )
            al.show();
        std::vector<int> f3 = al.FvaluesConference ( 3 );
        if ( verbose>=2 ) {
            printf ( "F3: " );
            display_vector ( f3 );
            printf ( "\n" );
        }
        const int N = al.n_rows;
        jstructconference_t js ( N, 4 );
        std::vector<int> f4 = al.FvaluesConference ( 4 );
        std::vector<int> j4 = js.Jvalues();

        myassert ( j4[0]==28, "unittest error: conference matricex F values: j4[0]\n" );
        myassert ( f4[0]==0, "unittest error: conference matricex F values: f4[0] \n" );
        myassert ( f4[1]==0, "unittest error: conference matricex F values: j4[1]\n" );

        if ( verbose>=2 ) {
            printf ( "j4: " );
            display_vector ( j4 );
            printf ( "\n" );
            printf ( "F4: " );
            display_vector ( f4 );
            printf ( "\n" );
        }
    }

    {
        cprintf ( verbose,"%s: LMC0 check for arrays in C(4, 3)\n", bstr );

        array_link al= exampleArray ( 28,1 );
        if ( verbose>=2 )
            al.showarray();
        lmc_t r =  LMC0check ( al, verbose );
        if ( verbose>=2 )
            printf ( "LMC0check: result %d\n", r );
        myassert ( r>=LMC_EQUAL, "LMC0 check\n" )   ;

        al= exampleArray ( 29,1 );
        if ( verbose>=2 )
            al.showarray();
        r =  LMC0check ( al, verbose );
        if ( verbose>=2 )
            printf ( "LMC0check: result %d (LMC_LESS %d)\n", r, LMC_LESS );
        myassert ( r==LMC_LESS, "LMC0 check of example array 29\n" )   ;
    }

    {
        cprintf ( verbose,"%s: LMC0 check\n", bstr );

        array_link al= exampleArray ( 31,1 );
        if ( verbose>=2 )
            al.showarray();
        conference_transformation_t T ( al );

        for ( int i=0; i<80; i++ ) {
            T.randomize();
            array_link alx = T.apply ( al );

            lmc_t r =  LMC0check ( alx, verbose );

            if (verbose>=2) {
                printfd("%d: transformed array: r %d\n", i, r);
                alx.showarray();   
            }
            if ( alx==al )
                myassert ( r>=LMC_EQUAL, "result should be LMC_MORE\n" );
            else {
                myassert ( r==LMC_LESS, "result should be LMC_LESS\n" );
            }
        }
    }


    {
        cprintf ( verbose,"%s: random transformation for conference matrices\n", bstr );

        array_link al = exampleArray ( 19,1 );
        conference_transformation_t T ( al );
        //T.randomizerowflips();
        T.randomize();

        conference_transformation_t Ti = T.inverse();
        array_link alx = Ti.apply ( T.apply ( al ) );

        if ( 0 ) {
            printf ( "input array:\n" );
            al.showarray();
            T.show();
            printf ( "transformed array:\n" );
            T.apply ( al ).showarray();
            Ti.show();
            alx.showarray();
        }

        myassert ( alx==al, "transformation of conference matrix\n" )	;
    }
#endif

    /* constructors */
    {
        cprintf ( verbose,"%s: constructors\n", bstr );

        array_transformation_t t;
        conference_transformation_t ct;

    }

    /* J-characteristics */
    {
        cprintf ( verbose,"%s: J-characteristics\n", bstr );

        array_link al=exampleArray ( 8,1 );

        const int mm[]= {-1,-1,0, 0, 8,16,0,-1};

        for ( int jj=2; jj<7; jj++ ) {
            std::vector<int> jx = al.Jcharacteristics ( jj );
            int j5max = vectormax ( jx, 0 );
            if ( verbose>=2 ) {
                printf ( "oaunittest: jj %d: j5max %d\n", jj, j5max );
            }

            if ( j5max!=mm[jj] ) {
                printfd ( "j5max %d (should be %d)\n", j5max, mm[jj] );
                allgood=UERROR;
                return allgood;
            }
        }
    }
    {
        cprintf ( verbose, "%s: array transformations\n", bstr );

        const int N=9;
        const int t = 3;
        arraydata_t adataX ( 3, N, t, 4 );

        array_link al ( adataX.N,adataX.ncols,-1 );
        al.create_root ( adataX );

        if ( checkTransformationInverse ( al ) )
            allgood=UERROR;

        if ( checkTransformationComposition ( al, verbose>=2 ) )
            allgood=UERROR;

        al = exampleArray ( 5, 1 );
        if ( checkTransformationInverse ( al ) )
            allgood=UERROR;

        if ( checkTransformationComposition ( al ) )
            allgood=UERROR;


        for ( int i=0; i<15; i++ ) {
            al=exampleArray ( 18,0 );
            if ( checkConferenceComposition ( al ) )
                allgood=UERROR;
            if ( checkConferenceInverse ( al ) )
                allgood=UERROR;
            al=exampleArray ( 19,0 );
            if ( checkConferenceComposition ( al ) )
                allgood=UERROR;
            if ( checkConferenceInverse ( al ) )
                allgood=UERROR;

        }

    }

    {
        cprintf ( verbose, "%s: rank \n", bstr );

        const int idx[10]= {0, 1, 2, 3, 4, 6, 7, 8, 9};
        const int rr[10]= {4,11,13,18,16, 4, 4,29,29};
        for ( int ii=0; ii<9; ii++ ) {
            array_link al = exampleArray ( idx[ii], 0 );
            myassert ( al.is2level(), "unittest error: input array is not 2-level\n" );

            int r = arrayrankColPivQR ( array2xf ( al ) );

            int r3 = ( array2xf ( al ) ).rank();
            myassert ( r==r3, "unittest error: rank of array" );

            if ( verbose>=2 ) {
                al.showarray();
                printf ( "unittest: rank of array %d: %d\n", idx[ii], r );
            }

            //array2xf ( al ) .showarray();
            myassert ( rr[ii]==r, "unittest error: rank of example matrix\n" );
        }

    }

    {
        cprintf ( verbose, "%s: Doptimize \n", bstr );
        const int N=40;
        const int t = 0;
        arraydata_t arrayclass ( 2, N, t, 6 );
        std::vector<double> alpha ( 3 );
        alpha[0]=1;
        alpha[1]=1;
        alpha[2]=0;
        int niter=5000;
        double t00 =get_time_ms();
        DoptimReturn rr =  Doptimize ( arrayclass, 10, alpha, 0, DOPTIM_AUTOMATIC, niter );

        array_t ss[7] = {3,3,2,2,2,2,2};
        arraydata_t arrayclassmixed ( ss, 36, t, 5 );
        rr =  Doptimize ( arrayclassmixed, 10, alpha, 0, DOPTIM_AUTOMATIC, niter );

        cprintf ( verbose, "%s: Doptimize time %.3f [s] \n", bstr, get_time_ms() - t00 );
    }


    {
        cprintf ( verbose, "%s: J-characteristics for conference matrix\n", bstr );

        array_link al = exampleArray ( 19, 0 );
        std::vector<int> j2 = Jcharacteristics_conference ( al, 2 );
        std::vector<int> j3 = Jcharacteristics_conference ( al, 3 );

        myassert ( j2[0]==0, "j2 value incorrect" );
        myassert ( j2[1]==0, "j2 value incorrect" );
        myassert ( std::abs ( j3[0] ) ==1, "j3 value incorrect" );

        if ( verbose>=2 ) {
            al.showarray();
            printf ( "j2: " );
            display_vector ( j2 );
            printf ( "\n" );
            printf ( "j3: " );
            display_vector ( j3 );
            printf ( "\n" );
        }
    }

    {
// test PEC sequence
        cprintf ( verbose, "%s: PEC sequence\n", bstr );
        for ( int ii=0; ii<5; ii++ ) {
            array_link al  = exampleArray ( ii, 0 );
            std::vector<double> pec = PECsequence ( al );
            printf ( "oaunittest: PEC for array %d: ", ii );
            display_vector ( pec );
            printf ( " \n" );
        }
    }


    {
        cprintf ( verbose, "%s: D-efficiency test\n", bstr );
        //  D-efficiency near-zero test
        {
            array_link al = exampleArray ( 14 );
            double D= al.Defficiency();
            std::vector<double> dd = al.Defficiencies();
            printf ( "D %f, D (method 2) %f\n", D, dd[0] );
            assert ( fabs ( D-dd[0] ) <1e-4 );
        }
        {
            array_link al = exampleArray ( 15 );
            double D= al.Defficiency();
            std::vector<double> dd = al.Defficiencies();
            printf ( "D %f, D (method 2) %f\n", D, dd[0] );
            assert ( fabs ( D-dd[0] ) <1e-4 );
            assert ( fabs ( D-0.335063 ) < 1e-3 );
        }

    }

    arraydata_t adata ( 2, 20, 2, 6 );
    OAextend oaextendx;
    oaextendx.setAlgorithm ( ( algorithm_t ) MODE_ORIGINAL, &adata );

    std::vector<arraylist_t> aa ( adata.ncols+1 );
    printf ( "OA unittest: create root array\n" );
    create_root ( &adata, aa[adata.strength] );

    /** Test extend of arrays **/
    {
        cprintf ( verbose, "%s: extend arrays\n", bstr );

        setloglevel ( SYSTEM );

        for ( int kk=adata.strength; kk<adata.ncols; kk++ ) {
            aa[kk+1] = extend_arraylist ( aa[kk], adata, oaextendx );
            printf ( "  extend: column %d->%d: %ld->%ld arrays\n", kk, kk+1, aa[kk].size(), aa[kk+1].size() );

        }

        if ( aa[adata.ncols].size() != 75 ) {
            printf ( "extended ?? to %d arrays\n", ( int ) aa[adata.ncols].size() );
        }
        myassert ( aa[adata.ncols].size() ==75, "number of arrays is incorrect" );

        aa[adata.ncols].size();
        setloglevel ( QUIET );
    }



    {
        cprintf ( verbose, "%s: test LMC check\n", bstr );

        array_link al = exampleArray ( 1, 1 );

        lmc_t r = LMCcheckOriginal ( al );
        //printfd("r %d, strength %d\n", r, strength);

        myassert ( r!=LMC_LESS, "LMC check of array in normal form" );

        for ( int i=0; i<20; i++ ) {
            array_link alx = al.randomperm();
            if ( alx==al )
                continue;
            lmc_t r = LMCcheckOriginal ( alx );

            myassert ( r==LMC_LESS, "randomized array cannot be in minimal form" );
        }
    }

    {
        /** Test dof **/
        cprintf ( verbose, "%s: test delete-one-factor reduction\n", bstr );

        array_link al = exampleArray ( 4 );
        cprintf ( verbose>=2, "LMC: \n" );
        al.reduceLMC();
        cprintf ( verbose>=2, "DOP: \n" );
        al.reduceDOP();
    }

    arraylist_t lst;

    {
        /** Test different methods **/
        cprintf ( verbose, "%s: test 2 different methods\n", bstr );

        const int s = 2;
        arraydata_t adata ( s, 32, 3, 10 );
        arraydata_t adata2 ( s, 32, 3, 10 );
        OAextend oaextendx;
        oaextendx.setAlgorithm ( ( algorithm_t ) MODE_ORIGINAL, &adata );
        OAextend oaextendx2;
        oaextendx2.setAlgorithm ( ( algorithm_t ) MODE_LMC_2LEVEL, &adata2 );

        printf ( "OA unittest: test 2-level algorithm on %s\n", adata.showstr().c_str() );
        std::vector<arraylist_t> aa ( adata.ncols+1 );
        create_root ( &adata, aa[adata.strength] );
        std::vector<arraylist_t> aa2 ( adata.ncols+1 );
        create_root ( &adata, aa2[adata.strength] );

        setloglevel ( SYSTEM );

        for ( int kk=adata.strength; kk<adata.ncols; kk++ ) {
            aa[kk+1] = extend_arraylist ( aa[kk], adata, oaextendx );
            aa2[kk+1] = extend_arraylist ( aa2[kk], adata2, oaextendx2 );
            printf ( "  extend: column %d->%d: %ld->%ld arrays, 2-level method %ld->%ld arrays\n", kk, kk+1, aa[kk].size(), aa[kk+1].size(), aa2[kk].size(), aa2[kk+1].size() );

            if ( aa[kk+1]!=aa2[kk+1] ) {
                printf ( "oaunittest: error: 2-level algorithm unequal to original algorithm\n" );
                exit ( 1 );
            }
        }
        setloglevel ( QUIET );

        lst = aa[8];
    }

    {
        cprintf ( verbose, "%s: rank calculation using rankStructure\n", bstr );

        for ( int i=0; i<27; i++ ) {
            array_link al = exampleArray ( i,0 );
            if ( al.n_columns<5 )
                continue;
            al = exampleArray ( i,1 );

            rankStructure rs;
            rs.verbose=0;
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
    }
    {
        cprintf ( verbose, "%s: test dtable creation\n", bstr );
        double t0x=get_time_ms();

        for ( int i=0; i<4; i++ ) {
            array_link al = exampleArray ( 5 );
            array_link dtable = createJdtable ( al );
        }
    }

    {
        cprintf ( verbose, "%s: test Pareto calculation\n", bstr );
        double t0x=get_time_ms();

        int nn=lst.size();
        for ( int k=0; k<5; k++ ) {
            for ( int i=0; i<nn; i++ ) {
                lst.push_back ( lst[i] );
            }
        }
        Pareto<mvalue_t<long>,long> r = parsePareto ( lst, 1 );
        cprintf ( verbose,"%s: test Pareto %d/%d: %.3f [s]\n", bstr, r.number(), r.numberindices(), ( get_time_ms() - t0x ) );

    }

    {
        cprintf ( verbose,"%s: check reduction transformation\n", bstr );
        array_link al = exampleArray ( 6 ).reduceLMC();

        arraydata_t adata=arraylink2arraydata ( al );
        LMCreduction_t reduction ( &adata );
        reduction.mode=OA_REDUCE;

        reduction.init_state=COPY;
        OAextend oaextend;
        oaextend.setAlgorithm ( MODE_ORIGINAL, &adata );
        array_link alr = al.randomperm();

        array_link al2=reduction.transformation->apply ( al );

        //printf("input: \n"); al2.showarray();
        lmc_t tmp = LMCcheck ( alr,adata, oaextend, reduction );

        array_link alx = reduction.transformation->apply ( alr );

        bool c = alx==al;
        if ( !c ) {
            printf ( "oaunittest: error: reduction of randomized array failed!\n" );
            printf ( "-- al \n" );
            al.showarraycompact();
            printf ( "-- alr \n" );
            alr.showarraycompact();
            printf ( "-- alx \n" );
            alx.showarraycompact();
            allgood=UERROR;
        }
    }

    {
        cprintf ( verbose,"%s: reduce randomized array\n", bstr );
        array_link al = exampleArray ( 3 );

        arraydata_t adata=arraylink2arraydata ( al );
        LMCreduction_t reduction ( &adata );

        for ( int ii=0; ii<50; ii++ ) {
            reduction.transformation->randomize();
            //reduction.transformation->randomizecolperm();
            array_link al2=reduction.transformation->apply ( al );

            //printf("input: \n"); al2.showarray();

            array_link alr=al2.reduceLMC();
            if ( 0 ) {
                printf ( "\n reduction complete:\n" );
                al2.showarray();
                printf ( "	--->\n" );
                alr.showarray();
            }
            bool c= ( al==alr );
            if ( !c ) {
                printf ( "oaunittest: error: reduction of randomized array failed!\n" );
                allgood=UERROR;
            }

        }
    }



    /* Calculate symmetry group */
    {
        cprintf ( verbose,"%s: calculate symmetry group\n", bstr );

        const arraylist_t &arraylist = aa[5];

        array_link al=arraylist.at ( 0 );
        symmetry_group sg = al.row_symmetry_group();

        aa.resize ( 0 );
    }

    /* Test efficiencies */
    {
        cprintf ( verbose,"%s: efficiencies\n", bstr );

        std::vector<double> d;
        int vb=1;

        array_link al;
        if ( 1 ) {
            al = exampleArray ( 9, vb );
            al.showproperties();
            d = al.Defficiencies ( 0, 1 );
            //printf("verbose: %d\n", verbose);
            if ( verbose>=2 )
                printf ( "  efficiencies: D %f Ds %f D1 %f Ds0 %f\n", d[0], d[1], d[2], d[3] );
            if ( fabs ( d[0]-al.Defficiency() ) >1e-10 ) {
                printf ( "oaunittest: error: Defficiency not good!\n" );
                allgood=UERROR;

            }
        }
        al = exampleArray ( 8, vb );
        //al.showarray();
        al.showproperties();
        d = al.Defficiencies();
        if ( verbose>=2 )
            printf ( "  efficiencies: D %f Ds %f D1 %f\n", d[0], d[1], d[2] );
        if ( fabs ( d[0]-al.Defficiency() ) >1e-10 ) {
            printf ( "oaunittest: error: Defficiency of examlple array 8 not good!\n" );
        }

        al = exampleArray ( 13, vb );
        //al.showarray();
        //al.showproperties();

        d = al.Defficiencies ( 0,1 );
        if ( verbose>=2 )
            printf ( "  efficiencies: D %f Ds %f D1 %f\n", d[0], d[1], d[2] );

        if ( ( fabs ( d[0]-0.939014 ) > 1e-4 ) || ( fabs ( d[3]-0.896812 ) > 1e-4 )  || ( fabs ( d[2]-1 ) > 1e-4 ) ) {
            printf ( "ERROR: D-efficiencies of example array 13 incorrect! \n" );
            d = al.Defficiencies ( 2,1 );
            printf ( "  efficiencies: D %f Ds %f D1 %f Ds0 %f\n", d[0], d[1], d[2], d[3] );

            allgood=UERROR;
            exit ( 1 );
        }

        for ( int ii=11; ii<11; ii++ ) {
            printf ( "ii %d: ", ii );
            al = exampleArray ( ii, vb );
            al.showarray();
            al.showproperties();

            d = al.Defficiencies();
            //if ( verbose>=2 )
            printf ( "  efficiencies: D %f Ds %f D1 %f\n", d[0], d[1], d[2] );
        }

    }
    {
        cprintf ( verbose,"%s: test robustness\n", bstr );

        array_link A ( 0,8,0 );
        printf ( "should return an error\n  " );
        A.Defficiencies();

        A = array_link ( 1,8,0 );
        printf ( "should return an error\n  " );
        A.at ( 0,0 ) =-2;
        A.Defficiencies();
    }

    {
        cprintf ( verbose,"%s: test nauty\n", bstr );

        array_link alr = exampleArray ( 7, 0 );
        if ( unittest_nautynormalform ( alr, 1 ) ==0 ) {
            printf ( "oaunittest: error: unittest_nautynormalform returns an error!\n" );
        }

    }

#ifdef HAVE_BOOST
    if ( writetests ) {
        cprintf ( verbose,"OA unittest: reading and writing of files\n" );

        boost::filesystem::path tmpdir = boost::filesystem::temp_directory_path();
        boost::filesystem::path temp = boost::filesystem::unique_path ( "test-%%%%%%%.oa" );

        const std::string tempstr    = ( tmpdir / temp ).native(); // optional

        if ( verbose>=2 )
            printf ( "generate text OA file: %s\n", tempstr.c_str() );

        int nrows=16;
        int ncols=8;
        int narrays=10;
        arrayfile_t afile ( tempstr.c_str(), nrows, ncols, narrays, ATEXT );
        for ( int i=0; i<narrays; i++ ) {
            array_link al ( nrows, ncols, array_link::INDEX_DEFAULT );
            afile.append_array ( al );
        }
        afile.closefile();

        arrayfile_t af ( tempstr.c_str(), 0 );
        std::cout << "  "<<  af.showstr() << std::endl;
        af.closefile();

        // check read/write of binary file

        arraylist_t ll0;
        ll0.push_back ( exampleArray ( 22 ) );
        ll0.push_back ( exampleArray ( 22 ).randomperm() );
        writearrayfile ( tempstr.c_str(), ll0, ABINARY );
        arraylist_t ll = readarrayfile ( tempstr.c_str() );
        myassert ( ll0.size() ==ll.size(), "read and write of arrays: size of list" );
        for ( size_t i=0; i<ll0.size(); i++ ) {
            myassert ( ll0[i]==ll[i], "read and write of arrays: array unequal" );
        }

        ll0.resize ( 0 );
        ll0.push_back ( exampleArray ( 24 ) );
        writearrayfile ( tempstr.c_str(), ll0, ABINARY_DIFFZERO );
        ll = readarrayfile ( tempstr.c_str() );
        myassert ( ll0.size() ==ll.size(), "read and write of arrays: size of list" );
        for ( size_t i=0; i<ll0.size(); i++ ) {
            myassert ( ll0[i]==ll[i], "read and write of arrays: array unequal" );
        }
    }

#endif

    {
        cprintf ( verbose,"OA unittest: test nauty\n" );
        array_link al = exampleArray ( 5, 2 );
        arraydata_t arrayclass = arraylink2arraydata ( al );

        for ( int i=0; i<20; i++ ) {
            array_link alx = al;
            alx.randomperm();
            array_transformation_t t1 = reduceOAnauty ( al );
            array_link alr1 = t1.apply ( al );

            array_transformation_t t2 = reduceOAnauty ( alx );
            array_link alr2 = t2.apply ( alx );


            if ( alr1 != alr2 )
                printf ( "oaunittest: error: Nauty reductions unequal!\n" );
            allgood=UERROR;


        }
    }

    cprintf ( verbose,"OA unittest: complete %.3f [s]!\n", ( get_time_ms() - t0 ) );
    cprintf ( verbose,"OA unittest: also run ptest.py to perform checks!\n" );

    if (	allgood ) {
        printf ( "OA unittest: all tests ok\n" );
        return UGOOD;
    } else {
        printf ( "OA unittest: ERROR!\n" );
        return UERROR;
    }
}

/**
* @brief Read in files with arrays and join them into a single file
* @param argc
* @param argv[]
* @return
*/
int main ( int argc, char* argv[] ) {

    AnyOption opt;
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "verbose", 'v' );
    opt.setOption ( "random", 'r' );

    opt.addUsage ( "OA: unittest: Perform some checks on the code" );
    opt.addUsage ( "Usage: unittest [OPTIONS]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -v  --verbose  			Print documentation" );
    opt.addUsage ( " -r  --random  			Seed for random number generator" );

    opt.processCommandArgs ( argc, argv );
    int verbose = opt.getIntValue ( 'v', 1 );
    int random = opt.getIntValue ( 'r', 0 );

    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) ) {
        opt.printUsage();
        exit ( 0 );
    }


    if ( verbose ) {
        print_copyright();
    }
    if ( verbose>=2 ) {
        print_options ( std::cout );
    }

    oaunittest ( verbose, 1, random );


    return 0;
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
