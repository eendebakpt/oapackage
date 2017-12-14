
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
#include "graphtools.h"

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

void speedcheck_conf ( const char *input, int verbose, int nmax=2000 )
{
     double t0;
     arraylist_t ll= readarrayfile ( input );
     printf ( "### speedcheck: read from file (%d arrays)\n", ( int ) ll.size() );

     t0=get_time_ms();
     nmax = std::min ( ( int ) ( ll.size() ), nmax );
     for ( size_t i=0; i< ( size_t ) nmax; i++ ) {
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
          ll.resize ( nmax );
          double t0=get_time_ms();
          arraylist_t  out= selectConferenceIsomorpismClasses ( ll, 0,CONFERENCE_ISOMORPHISM );


          for ( size_t i=0; i< ( size_t ) nmax; i++ ) {
               array_link al = ll[i];
          }
          printf ( "dt nauty %.3f \n", get_time_ms()-t0 );
     }


     return ;

}


//std::vector<double> Defficiencies (const array_link &al, const arraydata_t & arrayclass, int verbose ) ;

array_link finalcheck ( const array_link &al,  const arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int optimmethod, int niterx, int nabortx )
{
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






array_link array2xf2 ( const array_link &al )
{
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


void paretoInfo ( const array_link & alx )
{
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
int arrayInPareto ( const  Pareto < mvalue_t < long >, array_link >  &pset, const array_link &al, int verbose=1 )
{

     std::vector<array_link> llx = pset.allindices();
     arraylist_t ll ( llx.begin(), llx.end() );
     int jx=arrayInList ( al, ll, 0 );
     if ( verbose )
          myprintf ( "arrayInPareto: index in pareto list %d\n", jx );
     return jx;
}

/// check composition operator. returns 0 if test id good
int checkConferenceComposition ( const array_link &al, int verbose=0 )
{
     conference_transformation_t T1 ( al );
     //T1.randomizecolperm();
     T1.randomizecolflips();
     //T1.randomizerowperm();
     T1.randomizerowflips();
     T1.randomize();

     conference_transformation_t T2 ( al );
     //T2.randomize();
     T2.randomizerowperm();
     //T2.randomizerowflips();
     //T2.randomizecolperm();
     //T2.randomizecolflips();

     conference_transformation_t T3 = T2 * T1;

     array_link al1 = T1.apply ( al );
     array_link al1t2 = T2.apply ( al1 );
     array_link al3 = T3.apply ( al );

     if ( verbose ) {
          printfd ( "checkTransformationComposition: transforms\n" );
          printf ( "-- T1 \n" );
          T1.show();
          printf ( "-- T2 \n" );
          T2.show();
          printf ( "-- T3 \n" );
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



int main ( int argc, char* argv[] )
{
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
     int r = opt.getIntValue ( 'r', 8 );
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

     //print_options(); exit(0);

     {
          arraylist_t ll = readarrayfile ( input );

          lmc_t r = LMC0checkDC ( ll[0], verbose ); 
          //printf ( "array %d: result %d (LMC_LESS %d, LMC_MORE %d)\n", ( int ) 0, ( int ) r, LMC_LESS, LMC_MORE ); exit(0);

          for ( size_t i=0; i<ll.size(); i++ ) {
               lmc_t r = LMC0checkDC ( ll[i], verbose>=2 );
               printf ( "array %d: result %d (LMC_LESS %d, LMC_MORE %d)\n", ( int ) i, ( int ) r, LMC_LESS, LMC_MORE );

               // randomize
               array_link al = ll[i];

               for ( int jj=0; jj<30; jj++ ) {
                    conference_transformation_t tr ( al );
                    tr.randomizecolperm();
                    tr.randomizerowperm();
                    tr.randomizecolflips();

                    array_link alx = tr.apply ( al );

                    if ( al!=alx ) {
                         lmc_t r = LMC0checkDC ( alx, verbose>=2 );

                         printf ( "  randomized array %d: result %d (LMC_LESS %d, LMC_MORE %d)\n", ( int ) i, ( int ) r, LMC_LESS, LMC_MORE );
                         if ( r!=LMC_LESS ) {
                              printf ( "error!?\n" );
                              exit ( 1 );
                         }
                    } else {
                         printf ( " randomized array %d: same array\n", ( int ) i );
                    }
               }
          }

          return 0;
     }
     {
          array_link al = exampleArray ( r );

          array_transformation_t tt = reduceOAnauty ( al );
          array_link alx = tt.apply ( al );
          exit ( 0 );
     }

     {
          bool addcol=true;

          int nx=3;
          array_link G ( nx,nx, array_link::INDEX_DEFAULT );
          G.setconstant ( 0 );
          G.at ( 0,1 ) =1;
          G.at ( 1,0 ) =1;
          if ( addcol )
               G.at ( nx-1,nx-1 ) =5;


          // arraydata_t ad = arraylink2arraydata(G);
          arraydata_t ad ( 6, nx, 2, nx );

          //ad.show();
          array_transformation_t T ( ad );

          T.reset();

          std::vector<int> perm ( nx );
          for ( size_t i=0; i<perm.size(); i++ ) perm[i]=i;
          random_perm ( perm );

          T.setrowperm ( perm );
          T.setcolperm ( perm );
          T.show();

          G.showarray();
          G = T.apply ( G );
          G.showarray();

          std::vector<int> colors ( G.n_rows );
          for ( size_t i=0; i<colors.size(); i++ ) colors[i]=0;
          if ( addcol )
               colors[nx-1]=1;
          std::vector<int> colorsx ( G.n_rows );
          colorsx = permute ( colors, perm );
          //colorsx = permuteback(colors, perm);

          printfd ( "colorsx.size %d\n", colorsx.size() );
          print_perm ( "colorsx", colorsx );

          std::vector<int> tr = nauty::reduceNauty ( G, colorsx, verbose );
          print_perm ( "tr        ", tr );
          print_perm ( "tr inverse", invert_permutation ( tr ) );
          tr=invert_permutation ( tr );

          myprintf ( "-- output:\n" );
          array_link Gx = transformGraph ( G, tr, verbose );
          Gx.showarray();

          exit ( 0 );
     }
     if ( 1 ) {
          array_link al = exampleArray ( r, 1 );
          conference_t ct ( al.n_rows, al.n_columns+4, 0 );
          ct.j3zero=0;

          if ( xx>0 )
               al=al.selectFirstColumns ( xx );

          if ( verbose>=1 )
               al.showarray();

          assert ( al.is_conference() );
          assert ( al.min() ==-1 );


          int filterj2=1;
          int filtersymminline=1;
          int averbose=verbose;
          std::vector<cperm>      ccX = generateSingleConferenceExtensions ( al, ct, -1, averbose, 1, filterj2, ct.j3zero, filtersymminline );
          showCandidates ( ccX );
          printf ( "\n-----------\n" );

          CandidateGenerator cgenerator ( array_link(), ct );
          int kz = maxz ( al ) +1;
          cgenerator.verbose=verbose;
          std::vector<cperm> ee = cgenerator.generateCandidatesZero ( al, kz );

          cgenerator.showCandidates ( 2 );
          printf ( "generateCandidatesZero: %d\n-------------\n", ( int ) ee.size() );


          exit ( 0 );
     }
     if ( 0 ) {
          const int N = 16;
          int j1zero=1;
          conference_t ct ( N, 4, j1zero );
          ct.ctype=conference_t::DCONFERENCE;
          ct.j3zero=1;

          if ( 0 ) {
               arraylist_t ll = ct.createDconferenceRootArrays();
               printfd ( "generated %d root arrays\n", ll.size() );
               //array_link al2 = ct.create_root();
               //        array_link al3 = ct.create_root_three();
               array_link al = ll[0];
          }
          array_link al = exampleArray ( r, 1 );
          al.showarray();
//exit(0);

          CandidateGeneratorDouble cgenerator ( array_link(), ct );
          cgenerator.verbose=2;
          for ( int i=0; i<2; i++ ) {
               {
                    printf ( "\n---------------------------------\n" );
                    const std::vector<cperm> &cl = cgenerator.generateCandidates ( al );
                    printfd ( "generated %d\n", cl.size() );
                    cgenerator.showCandidates ( 2 );
               }
          }
          exit ( 0 );
     }
     {
          int j1zero=0;
          const int N = r;


          conference_t ctype ( N, 4, j1zero );
          ctype.ctype=conference_t::CONFERENCE_DIAGONAL;
          ctype.ctype=conference_t::CONFERENCE_NORMAL;
          ctype.itype=CONFERENCE_ISOMORPHISM;
          ctype.j3zero=0;

          array_link al2 = ctype.create_root();
          array_link al3 = ctype.create_root_three();


          int ii=2;
          int ncstart=2;
          int ncmax=2;
          array_link root = al2;

          if ( 0 ) {
               int extcol=2;
               std::vector<cperm>        ee = generateConferenceExtensions ( al2, ctype, extcol, 0, 0, 1 );
               printfd ( "generated %d\n", ee.size() );
          }

          if ( 0 ) {
               int extcol=3;
               std::vector<cperm>        ee2 = generateConferenceExtensions ( al3, ctype, extcol, 0, 0, 1 );

               //    conf_candidates_t tmp = generateCandidateExtensions ( ctype, 2, ncstart, ncmax, root );

          }
          printf ( "------------------------------\n" );

          CandidateGenerator cgenerator ( array_link(), ctype );
          cgenerator.verbose=verbose;
          //cgenerator.generators[ii].verbose=2;
          cgenerator.showCandidates ( 2 );

          printf ( "------------------------------\n" );
          const std::vector<cperm> &cl = cgenerator.generateCandidatesZero ( al2, ii );
          cgenerator.showCandidates ( 2 );
          printfd ( " cache: generated %d\n", cl.size() );

          exit ( 0 );

          /*
          std::vector<cperm>  cltotal;
          for ( int i=0; i<ctype.N; i++ ) {
               CandidateGeneratorZero cgenerator2 ( array_link(), ctype, i );
               //cgenerator2.verbose=2;
               std::vector<cperm> clc;
               clc = cgenerator2.generateCandidates ( al2 );
               if ( verbose>=3 )
                    printfd ( "    cache %d: generated %d\n", ii, clc.size() );

               cltotal.insert ( cltotal.begin(), clc.begin(), clc.end() );
          }
          printfd ( " cache: generated %d\n", cltotal.size() );
          */
          printf ( "done\n" );
          exit ( 0 );
     }

     if ( 0 ) {
          array_link al=exampleArray ( 29,0 );
          checkConferenceComposition ( al, 1 ) ;
          exit ( 0 );

     }
     if ( 0 ) {
          array_link al= exampleArray ( 28,1 );
          al.showarray();
          lmc_t r =  LMC0check ( al, verbose );
          printf ( "result %d\n", r );
          al= exampleArray ( 29,1 );
          al.showarray();
          r =  LMC0check ( al, verbose );
          printf ( "result %d\n", r );
          al= exampleArray ( 30,1 );
          al.showarray();
          r =  LMC0check ( al, verbose );
          printf ( "result %d\n", r );
          exit ( 0 );
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
          conference_t ctype ( al.n_rows, al.n_rows , 1 );
          ctype.itype=CONFERENCE_RESTRICTED_ISOMORPHISM;
          ctype.ctype=conference_t::DCONFERENCE;

          CandidateGeneratorDouble cgenerator ( array_link() , ctype );
          cgenerator.verbose=verbose;

          for ( int i=0; i< ( int ) ll.size(); i++ ) {
               std::vector<cperm> cc = cgenerator.generateCandidates ( ll[i] );
               printfd ( "generated %d\n", cc.size() );
          }
          exit ( 0 );
     }


     return 0;

}

// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 
