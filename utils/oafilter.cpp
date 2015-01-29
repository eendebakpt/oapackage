/** \file oafilter.cpp

 C++ program: oafilter

 oafilter: filter arrays

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "arrayproperties.h"

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"


template <class Type>
/// Write a vector of vector  elements to file
void vectorvector2file(const std::string fname, const std::vector<std::vector<Type> > vals)
{
    FILE *fid = fopen(fname.c_str(), "wb");
    if (fid==0)
    {
        cout << "error with file " << fname << std::endl;
        exit(1);
    }
    for (unsigned int i=0; i<vals.size(); i++) {
        std::vector<Type> x = vals[i];

        for (unsigned int j=0; j<x.size(); j++) {
            fprintf(fid, "%f ", x[j]);
        }
        fprintf(fid, "\n");
    }
    fclose(fid);

}

/// Write a vector of vector  elements to binary file
void vectorvector2binfile(const std::string fname, const std::vector<std::vector<double> > vals)
{
    FILE *fid = fopen(fname.c_str(), "wb");
    if (fid==0)
    {
        cout << "error with file " << fname << std::endl;
        exit(1);
    }
    for (unsigned int i=0; i<vals.size(); i++) {
        std::vector<double> x = vals[i];

        for (unsigned int j=0; j<x.size(); j++) {
            fwrite(&(x[j]), sizeof(double), 1, fid);
        }
    }
    fclose(fid);
}

template <class Type>
/// Write a vector of integer elements to file
void intvector2file(std::string fname, std::vector<Type> vals)
{
    FILE *fid = fopen(fname.c_str(), "wb");
    if (fid==0)
    {
        cout << "error with file " << fname << std::endl;
        exit(1);
    }
    for (unsigned int i=0; i<vals.size(); i++) {
        Type x = vals[i];
        fprintf(fid, "%d\n", x);
    }
    fclose(fid);

}


/**
 * @brief Read in files with arrays and join them into a single file
 * @param argc
 * @param argv[]
 * @return
 */
int main(int argc, char* argv[])
{
    AnyOption opt;

    /* parse command line options */
    opt.setFlag(  "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption("output", 'o');
    opt.setOption("verbose", 'v');
    opt.setOption("na", 'a');
    opt.setOption("index", 'i');
    opt.setOption("format", 'f');

    opt.addUsage( "Orthonal Array Filter: filter arrays" );
    opt.addUsage( "Usage: oaanalyse [OPTIONS] [INPUTFILE] [VALUESFILE] [THRESHOLD] [OUTPUTFILE]" );
    opt.addUsage( "  The VALUESFILE is a binary file with k*N double values. Here N is the number of arrays" );
    opt.addUsage( "  and k the number of analysis values." );

    opt.addUsage( "" );
    opt.addUsage( " -h --help  			Prints this help " );
    opt.addUsage( " -v --verbose  			Verbose level (default: 1) " );
    opt.addUsage( " -a 		 			Number of values in analysis file (default: 1) " );
    opt.addUsage( " --index 		 		Index of value to compare (default: 0) " );
    opt.addUsage( " -f [FORMAT]					Output format (default: TEXT, or BINARY; B) " );
    opt.processCommandArgs(argc, argv);

    print_copyright();

    /* parse options */
    if ( opt.getFlag( "help" ) || opt.getFlag( 'h' ) || opt.getArgc()<4 ) {
        opt.printUsage();
        exit(0);
    }
    int sortarrays = opt.getIntValue('s', 0);
    int verbose = opt.getIntValue('v', 1);
    int na = opt.getIntValue('a', 1);
    int index = opt.getIntValue("index", 0);

    std::string format = opt.getStringValue('f', "BINARY");
    arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString(format);

    //printf(" we have %d args\n", opt.getArgc());

    /* read in the arrays */
    std::string infile =  opt.getArgv( 0 );
    std::string anafile =  opt.getArgv( 1 );
    double threshold = atof(opt.getArgv( 2 ) );
    std::string outfile =  opt.getArgv( 3 );

    if (verbose)
        printf("oafilter: input %s, threshold %f (analysisfile %d values, index %d)\n", infile.c_str(), threshold, na, index);

    arraylist_t *arraylist = new arraylist_t;
    int n = readarrayfile( opt.getArgv( 0 ), arraylist);

    /* perform operations on arrays */
    if (sortarrays) {
        cout << "Sorting arrays" << endl;
        sort(arraylist->begin(), arraylist->end());
    }

    if (verbose)
        printf("oafilter: filtering %d arrays\n", n);

    // TODO: update to new binary format

    // read analysis file
    FILE *afid = fopen(anafile.c_str(),  "rb");
    if (afid==0)
    {
        printf("   could not read analysisfile %s\n", anafile.c_str() );
        exit(0);
    }

    double *a = new double[n*na];
    fread(a, sizeof(double), n*na, afid);
    fclose(afid);

    std::vector<int> gidx;;
    for(int jj=0; jj<n; jj++) {
        double val=a[jj*na+index];
        if (val>= threshold)
            gidx.push_back(jj);
        if (verbose>=2)
            printf("jj %d: val %f, threshold %f\n", val>=threshold, val, threshold);
    }

    // filter
    arraylist_t *filtered = new arraylist_t;
    selectArrays(*arraylist,  gidx, *filtered);


    /* write arrays to file */
    if (verbose)
        cout << "Writing " << filtered->size() << " arrays (input " << arraylist->size() << ") to file " << outfile << endl;
    writearrayfile(outfile.c_str(), filtered, mode);

    /* free allocated structures */
    delete [] a;
    delete arraylist;
    delete filtered;

    return 0;
}
