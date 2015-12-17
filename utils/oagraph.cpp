/** \file oagraph.cpp

 C++ program: oagraph

 oagraph: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "graphtools.h"


template<class Type>
std::vector<Type> sorthelper(std::vector<Type> &v)
{
std::vector<Type> w(v.size());
	std::copy(v.begin(), v.end(), w.begin() );
//			printf("sorthelper: "); display_vector(w); printf("\n");
	int rmin=w[0]; for(int j=0; j<w.size(); j++) {
//	printf("rmin %d, w[j] %d\n", rmin, w[j] );
		rmin = std::min(rmin, w[j]);
	}
//			printf("sorthelper: "); display_vector(w); printf("\n  rmin %d\n", rmin);

	for(int i=0; i<w.size(); i++) w[i] -= rmin;
//			printf("sorthelper: "); display_vector(w); printf("\n");
																return w;

}


/// From a relabelling of the graph return the corresponding array transformation 
array_transformation_t labelling2transformation(std::vector<int> &pp, arraydata_t &arrayclass, int verbose=1)
{
	///invert the labelling transformation
   std::vector<int> ppi(pp.size()); for(size_t i=0; i<pp.size(); i++) ppi[pp[i]]=i;

     // extract colperms and rowperm and levelperms from pp
   array_transformation_t ttr(arrayclass);

   if (verbose) {
	   printf("labelling2transformation: class %s\n", arrayclass.idstr().c_str() );
   }
   const int N= arrayclass.N;
	std::copy(pp.begin(), pp.begin()+N, ttr.rperm);
	int rmin=pp.size(); for(int j=0; j<N; j++) rmin = std::min(rmin,(int) ttr.rperm[j]);
	for(int i=0; i<N; i++) ttr.rperm[i] -= rmin;
ttr=ttr.inverse();

   if (verbose) {
	   printf("labelling2transformation: rowperm "); print_perm(ttr.rperm, N);
   }
   
int ncols=arrayclass.ncols;
   array_transformation_t ttc(arrayclass);
    std::vector<int> colperm(arrayclass.ncols); std::copy(pp.begin()+N, pp.begin()+N+ncols, colperm.begin() );
	colperm = sorthelper(colperm);
	ttc.setcolperm(colperm);

   if (verbose) {
	   printf("labelling2transformation: colperm "); display_vector(colperm); printf("\n");
   }
	
	std::vector<int> s = arrayclass.getS();
	
	int ns = std::accumulate(s.begin(), s.end(), 0);
	array_transformation_t ttl(arrayclass);
	ttl=ttl.inverse();
	
	std::vector<int> lvlperm(ns);
	std::copy(pp.begin()+N+ncols, pp.begin()+N+ncols+ns, lvlperm.begin() );
	lvlperm=sorthelper(lvlperm);

		std::vector<int>cs=cumsum0(s);

    //lp = []
    for(int ii=0; ii<ncols; ii++) {
        std::vector<int> ww(lvlperm.begin()+cs[ii], lvlperm.begin()+cs[ii+1] ); //  = lvlperm[cs[ii]:cs[ii + 1]]
		//printf("ii %d: ww ", ii); display_vector(ww); printf("\n");
		ww=sorthelper(ww);
		//printf("ww "); display_vector(ww); printf("\n");
        //ww = ww - ww.min();        ww = np.argsort(ww)
        //lp.append(ww)
		if (verbose)
		{
			printf("lvlperm %d: ",ii); display_vector(ww); printf("\n"); fflush(0);
		}
        ttl.setlevelperm(ii, ww);
}

 

	ttl=ttl.inverse();
	array_transformation_t tt = ttr*ttc*ttl;

	
	return tt;
}



int main ( int argc, char* argv[] )
{
	AnyOption opt;
	/* parse command line options */
	opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt.setOption ( "output", 'o' );
	opt.setOption ( "rand", 'r' );
	opt.setOption ( "verbose", 'v' );
	opt.setOption ( "ii", 'i' );
	opt.setOption ( "dverbose", 'd' );
	opt.setOption ( "rows" );
	opt.setOption ( "cols" );
	opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

	opt.addUsage ( "Orthogonal Array: oagraph: testing platform" );
	opt.addUsage ( "Usage: oagraph [OPTIONS] " );
	opt.addUsage ( "" );
	opt.addUsage ( " -h --help  			Prints this help " );
	opt.processCommandArgs ( argc, argv );


	int randvalseed = opt.getIntValue ( 'r', 1 );

	srand ( randvalseed );

	int verbose = opt.getIntValue ( 'v', 1 );
	int ix = opt.getIntValue ( 'i', 2 );

	print_copyright();
	setloglevel ( NORMAL );

	/* parse options */
	if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
		opt.printUsage();
		exit ( 0 );
	}


	setloglevel ( SYSTEM );

	// select a design
	array_link al = exampleArray ( ix, 2 );
	arraydata_t arrayclass = arraylink2arraydata(al);
	
	array_transformation_t trans(arrayclass);
trans.reset();
	
	trans.randomize();
	array_link alr = trans.apply(al);
	alr=al;
	
	std::vector<int> tr = nauty::reduceOAnauty ( alr, 1 );
	printf ( "canon: " );
	display_vector ( tr );
	printf ( "\n" );

array_transformation_t ttm = labelling2transformation(tr, arrayclass, verbose);
	ttm.show();
	
	array_link alm = ttm.apply(alr);
	
	al.showarraycompact();
	printf("---\n");
	alr.showarraycompact();
	printf("---\n");
alm.showarraycompact();



//return 0;
	
	array_link alx = al.randomperm();
	alx=alr;
	std::vector<int> trx = nauty::reduceOAnauty ( alx, 0 );
	array_transformation_t ttx = labelling2transformation(trx, arrayclass, 0);
	array_link alxm = ttx.apply(alx);

	printf("--- alxm \n");
	alxm.showarraycompact();
	
	if(alxm!=alm) {
	printf("error with nauty reduction...\n");	
	}
	// TODO: extract arraytransformation_t from canon
	// TODO: print minimal arrays and make unit test using 2 random transformations
	// TODO: clean up interfaces



	/*
	// convert design (with isomorphism type) to a colored graph
	nauty::NyGraph G = array2graph(al, 1);

	if (verbose)
	printf("creatd graph with %d vertices and ? edges\n", G.num_nodes );

	*/
	if ( verbose )
		printf ( "done\n" );


	return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
