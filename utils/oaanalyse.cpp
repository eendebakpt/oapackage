/** \file oajoin.cpp

 C++ program: oaanalyse

 oaanalyse: report on array data

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2008

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <numeric>

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"
#include "arrayproperties.h"

#ifdef _OPENMP
//#include "omp.h"
#endif

template <class Type>
/// Write a vector of vector elements to file
void vectorvector2file ( const std::string fname, const std::vector<std::vector<Type> > vals )
{
    FILE *fid = fopen ( fname.c_str(), "wb" );
    if ( fid==0 ) {
        cout << "error with file " << fname << std::endl;
        exit ( 1 );
    }
    for ( unsigned int i=0; i<vals.size(); i++ ) {
        std::vector<Type> x = vals[i];
        for ( unsigned int j=0; j<x.size(); j++ ) {
            fprintf ( fid, "%f ", x[j] );
        }
        fprintf ( fid, "\n" );
    }
    fclose ( fid );

}

int gma2jtype ( std::vector<double> gwlp )
{
    int jtype=0; // fold-over
    for ( size_t i=1; i<gwlp.size(); i+=2 ) {
        if ( fabs(gwlp[i])>1e-6 ) {
            jtype=i;
            break;
        }
    }
    return jtype;
}


template <class Type>
/// Write a vector of integer elements to file
void intvector2file ( std::string fname, std::vector<Type> vals )
{
    FILE *fid = fopen ( fname.c_str(), "wb" );
    if ( fid==0 ) {
        cout << "error with file " << fname << std::endl;
        exit ( 1 );
    }
    for ( unsigned int i=0; i<vals.size(); i++ ) {
        Type x = vals[i];
        fprintf ( fid, "%d\n", x );
    }
    fclose ( fid );

}




/// Helper function
bool testarray ( const jstruct_t &js, int maxj )
{
    //bool ret=1;
    for ( int x=0; x<js.nc; x++ ) {
        // printf("test array %d %d\n", maxj, js.vals[x]);
        if ( js.values[x]>maxj )
            return false;
    }
    return true;
}

/**
 * @brief Analyse arrays
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] )
{

    /* parse command line options */
    AnyOption opt;
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "output", 'o' );
    opt.setOption ( "analysis", 'a' );
    opt.setOption ( "maxj", 'm' );
    opt.setOption ( "jj", 'j' );
    opt.setFlag ( "rank", 'r' );
    opt.setFlag ( "gwp" );
    opt.setFlag ( "mixedlevel" );
    opt.setFlag ( "gma" );
    opt.setFlag ( "jtype" );
    opt.setFlag ( "cl2d" );
    opt.setFlag ( "Defficiency", 'A' );
    opt.setOption ( "verbose", 'v' );

    opt.addUsage ( "Orthonal Array Analyse: analyse arrays" );
    opt.addUsage ( "Usage: oaanalyse [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.addUsage ( " -r  					Calculate rank of matrix (only defined for 2-factor arrays (default: 0) " );
    opt.addUsage ( " -A  					Calculate D-efficiency of matrix (only defined for 2-factor arrays (default: 0) " );
    opt.addUsage ( " -v --verbose  			Verbose level (default: 1) " );
    opt.addUsage ( " --gwp  				Calculate GWLP (default: 0) " );
    opt.addUsage ( " --cl2d  				Calculate centered L2-discrepancy (default: 0) " );
    opt.addUsage ( " -f [NUMBER]  			Signed J-characteristic to compute (default: 4) " );
    opt.addUsage ( " --maxj [NUMBER] 			Maximum value for j-value (default: Inf) " );
    opt.addUsage ( " -o [FILE] --output [FILE]	Output file for filtered arrays (default: no output) " );
    opt.addUsage ( " -a [FILE] --analysis [FILE]	Output file for analysis data (default: no output) " );
    opt.addUsage ( " -j [NUMBER] 			Value (max) for J-characteristics (default: 4) " );
    opt.processCommandArgs ( argc, argv );

    print_copyright();

    /* parse options */
    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() ==0 ) {
        opt.printUsage();
        exit ( 0 );
    }

    const char *outputprefix = 0;
    if ( opt.getValue ( "output" ) !=0 )
        outputprefix = opt.getValue ( 'o' );

    const char *analysisfile = 0;
    if ( opt.getValue ( "analysis" ) !=0 )
        analysisfile = opt.getValue ( 'a' );

    int maxj =  opt.getIntValue ( "maxj", 1000000 );
    int dorank = opt.getFlag ( 'r' );
    int dojtype = opt.getFlag ( "jtype" );
    int doDefficiency = opt.getFlag ( 'A' );
    int dogma = opt.getFlag ( "gwp" );
    int mixedlevel = opt.getFlag ( "mixedlevel" );
    if ( !dogma ) {
        dogma = opt.getFlag ( "gma" );
    }
    int docl2d = opt.getFlag ( "cl2d" );
    int jj = opt.getIntValue ( 'j', 4 );
    int doj= ( jj>=0 );
    int verbose = opt.getIntValue ( 'v', 1 );

    bool sortarrays = 0;
    if ( opt.getFlag ( 's' ) !=0 )
        sortarrays = 1;

    /* read in the arrays */
    arraylist_t *arraylist = new arraylist_t;
    if ( verbose )
        std::cout << "oaanalyse: reading " << opt.getArgc() << " file(s)" << std::endl;

#ifdef _OPENMP
    if (verbose>=1)
    {
        printf("oaanalyse: openmp %d\n", omp_get_max_threads() );

    }
#endif
    int nr;
    int nc;
    int nn;
    for ( int i = 0 ; i < opt.getArgc() ; i++ ) {
        if ( verbose )
            std::cout << "file " <<  opt.getArgv ( i ) << std::endl ;
        int n = readarrayfile ( opt.getArgv ( i ), arraylist );

        arrayfileinfo ( opt.getArgv ( i ), nn, nr, nc );
        if ( verbose )
            std::cout << "   read " << n << printfstring ( " array(s) (size %d %d)", nr, nc ) << std::endl;
    }

    /* perform operations on arrays */
    if ( sortarrays ) {
        cout << "oaanalyse: sorting arrays" << endl;
        sort ( arraylist->begin(), arraylist->end() );
    }


    // analyse
    //jstruct_t *js ;

    vector<jstruct_t> results;

    if ( doj ) {
        results = analyseArrays ( *arraylist, verbose, jj );
    }
    vector<int> rankx;
    std::vector<std::vector<double> > rankvals;
    rankx.resize(arraylist->size());
    rankvals.resize(arraylist->size());
    std::vector<double> AA ( arraylist->size() );
    double maxA=0;

    /* filter */
    arraylist_t filtered;
    if ( doj ) {
        for ( unsigned int ii=0; ii<arraylist->size(); ii++ ) {
            if ( ii%5000==0 || verbose>=3 )
                printf ( "  calculation (j-filter) for array %d/%d\n", ii, ( int ) arraylist->size() );

            bool ret = testarray ( results[ii], maxj );
            if ( ret ) {
                filtered.push_back ( arraylist->at ( ii ) );
            }
        }
    }

    int loopcounter=0;
    
#ifdef _OPENMP
    if (verbose>=2) {
	// printing results should be ordered...
      omp_set_num_threads(1);
    }
#endif
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for ( int ii=0; ii<(int)arraylist->size(); ii++ ) {
#ifdef _OPENMP
        #pragma omp critical
        {
            if ( loopcounter%20000==0 || verbose>=3 ) {
                //  if (omp_get_thread_num()==0)
                printf ( "  calculation (rank, Defficiency) for array %d/%d\n", loopcounter, ( int ) arraylist->size() );

            }
            loopcounter++;
        }
#else
        if ( ii%5000==0 || verbose>=3 ) {
            printf ( "  calculation (rank, Defficiency) for array %d/%d\n", ii, ( int ) arraylist->size() );
        }
#endif
        if ( dorank ) {
            std::vector<double> rr;
            int r = array_rank_D_B ( arraylist->at ( ii ), &rr, verbose );
            if ( verbose>=2 ) {
                #pragma omp critical
                {
                    cout << "   rank values: ";
                    for ( std::vector<double>::iterator it = rr.begin(); it<rr.end(); it++ )
                        cout << *it << " ";
                    std::cout << std::endl;
                }
            }

            #pragma omp critical
            {
                //rankvals.push_back ( rr ); rankx.push_back ( r );
                rankvals[ii] = rr;
                rankx[ii] = r ;
            }
        }

        if ( doDefficiency ) {
	    
	   if (mixedlevel) {
            std::vector<double> dd = arraylist->at ( ii ).Defficiencies(0);
 AA[ii]=dd[0];
	   }
else {
            AA[ii]=Defficiency ( arraylist->at ( ii ), 0 );
}
	    if (0) {
            arraylist->at ( ii ).DsEfficiency(1);
            std::vector<double> dd = arraylist->at ( ii ).Defficiencies(1);
	  printf("  dd %f %f %f\n", dd[0], dd[1], dd[2] );
	    }
	    #pragma omp critical
            {
                if ( AA[ii]>maxA ) {
                    maxA=AA[ii];
                    if ( verbose>=2 ) {
                        printf ( "array: %d Defficiency %f\n", ii, AA[ii] );
                    }
                } else {
                    if ( verbose>=2 ) {
                        printf ( "array: %d Defficiency %f\n", ii, AA[ii] );
                    }
                }
            }
        }

    }



    //printf("dogma %d: GMA\n", dogma);
    int ncols=-1;


    std::vector<std::vector<double> > gmalist(arraylist->size());
    std::vector<double> c2list(arraylist->size());
    loopcounter=0;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (  int i=0; i<(int)arraylist->size(); i++ ) {
        #pragma omp critical
        {
            //printf("array %d: GMA\n", i);
            if ( loopcounter%4000==0 && verbose>=1 && (dogma || dojtype || docl2d) ) {
                printf ( "  statistics calculation (GWLP) for array %d/%d\n", loopcounter, ( int ) arraylist->size() );
            }
            loopcounter++;
        }
        array_link &al = arraylist->at ( i );
        if ( dogma || dojtype ) {
            int n = al.n_columns;
            std::vector<double> gwp = GWLP ( al, verbose>=3 );

            #pragma omp critical
            {
                ncols=n;
                if ( verbose>=3 )
                    al.showarray();
                if ( dojtype ) {
                    int jtype = gma2jtype ( gwp );
                    if ( verbose>=1 ) {
                        printf ( "array %d: jtype %d\n", i, jtype );
                    }

                }
                gmalist[i] =  gwp;

                if ( verbose>=2 ) {
                    printf ( "array %d: %dx%d GWP ", i, al.n_rows, al.n_columns );
                    for ( int x=0; x<=n; x++ ) {
                        printf ( "%.4f ", gwp[x] );
                    }
                    std::cout << std::endl;
                    //display_vector(gwp); std::cout << std::endl;
                    //double sum_of_elems =std::accumulate(gwp.begin(),gwp.end(),0.); printf("   sum %.3f\n", sum_of_elems);
                }
            }
        }
        if ( docl2d ) {
            double c2 = al.CL2discrepancy();

            #pragma omp critical
            {
                if ( verbose>=2 ) {
                    printf ( "array %d: centered L2 discrepancy %f\n ", i, c2 );
                }

                c2list[i]= ( c2 );
            }
        }
    }

    if ( analysisfile!=0 ) {
        if ( dogma ) {
            std::string abase = analysisfile;
            std::string gmafile = abase + printfstring ( "-gwlpfull.bin" );
            int ncb=ncols+1;	// number of columns
            vectorvector2binfile ( gmafile, gmalist, 1, ncb );
            if ( verbose>=2 )
                std::cout << "writing gma file " << gmafile << std::endl;
        }
        if ( docl2d ) {
            std::string c2file = analysisfile + printfstring ( "-cl2d.bin" );
            doublevector2binfile ( c2file, c2list );
            if ( verbose>=2 )
                std::cout << "writing centered L2 discrepancy file " << c2file << std::endl;
        }
    }

    /* write analysis data */
    if ( analysisfile!=0 ) {
        std::string abase = analysisfile;
        printf ( "oaanalyse: writing analysis file %s for %d array(s)\n", analysisfile, ( int ) arraylist->size() );

        std::string sj = abase + printfstring ( "-jvals%d.txt", jj );
        std::string sgma = abase + printfstring ( "-gma%d.txt", jj );

        if ( doj ) {
            FILE *fid = fopen ( sj.c_str(), "wb" );
            FILE *fidgma = fopen ( sgma.c_str(), "wb" );
            vector<jstruct_t>::iterator it;
            for ( unsigned int i=0; i<results.size(); i++ ) {
                if ( verbose>=2 ) {
                    printf ( " writing data for array %d\n", i );
                }
                //int *vv = results[i].vals[j];
                for ( int j=0; j<results[i].nc; j++ ) {
                    if ( j==results[i].nc-1 )
                        fprintf ( fid, "%d", results[i].values[j] );
                    else
                        fprintf ( fid, "%d,", results[i].values[j] );
                }
                fprintf ( fid, "\n" );

                fprintf ( fidgma, "%f\n", results[i].A );

            }
            fclose ( fid );
            fclose ( fidgma );
        }

        if ( dorank ) {
            std::string rankfile = abase + "-rank.txt";
            //intvector2file<int>(rankfile, rankx);
            //rankfile = abase + "-rankvals.txt";
            //vectorvector2file<double>(rankfile, rankvals);
            rankfile = abase + "-rankvalues.bin";
            vectorvector2binfile ( rankfile, rankvals, 1, 4 );
        }

        if ( doDefficiency ) {
            std::string DeffFile = abase + "-Dvalues.bin";
            if ( verbose>=2 )
                printf ( "writing: %s\n", DeffFile.c_str() );
            doublevector2binfile<double> ( DeffFile, AA, 1 );

        }
    }


    if ( 1 ) {
        /* write arrays to file */
        string outfile = "";

        if ( opt.getValue ( "output" ) !=0 ) {
            outfile += outputprefix;
            outfile += ".oa";
            cout << "Writing " << filtered.size() << " arrays (input " << arraylist->size() << ") to file " << outfile << endl;
            writearrayfile ( outfile.c_str(), &filtered );
        } else {
            outfile += "<standard output>";
        }

    }

    /* free allocated structures */
    delete arraylist;


    return 0;
}
