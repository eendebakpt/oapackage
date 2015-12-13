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

template<class NumType>
std::vector<NumType> cumsum ( const std::vector<NumType> x )
{
	// initialize the result vector
	std::vector<NumType> res ( x.size() );
	std::partial_sum ( x.begin(), x.end(), res.begin() );
	return res;
}

template<class NumType>
std::vector<NumType> cumsum0 ( const std::vector<NumType> x )
{
	// initialize the result vector
	std::vector<NumType> res ( x.size() +1 );
	res[0]=0;
	std::partial_sum ( x.begin(), x.end(), res.begin() +1 );
	return res;
}

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

template<class Type, class InputType>
std::vector<Type> cumsum0(std::vector<InputType> s)
{
	std::vector<Type> c(s.size()+1);
	c[0]=0;
	for(int i=0; i<s.size();i++) c[i+1]=c[i]+s[i];
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

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair<array_link, std::vector<int> >  array2graph ( const array_link &al, int verbose=1 )
{
	arraydata_t arrayclass = arraylink2arraydata ( al );
	int nrows = al.n_rows;
	int ncols = al.n_columns;
	const std::vector<int> s = arrayclass.getS();

	int nRowVertices = nrows;
	int nColumnLevelVertices = std::accumulate ( s.begin(), s.end(), 0 );
	int nColVertices = ncols;

	int nVertices = nrows + ncols + nColumnLevelVertices; // number of vertices in the graph
	int colOffset = nrows;

	std::vector<int> vertexOffsets ( s.size() +1 );

	std::vector<int> cs = cumsum0 ( s );

	for ( size_t j=0; j<s.size(); j++ )
		vertexOffsets[j] = colOffset + ncols + cs[j];

	std::vector<int> colors ( nVertices );

	// row vertices: color 0
	for ( int i=0; i<nrows; i++ )
		colors[i]=0;
	// column vertices: color 1
	for ( int i=0; i<ncols; i++ )
		colors[i+nrows]=1;
	// other vertices: color 2
	for ( int i=0; i<nColumnLevelVertices; i++ )
		colors[i+nrows+ncols]=2;


	if (verbose)
		printf("array2graph: generating graph of size %d=%d+%d+%d\n", nVertices, nrows, ncols, nColumnLevelVertices);
	array_link G ( nVertices, nVertices, 0 ); // graph
G.setconstant(0);

//    im = np.zeros((nVertices, nVertices))  # incidence matrix


	for ( int r=0; r<nrows; r++ ) {
		for ( int c=0; c<ncols; c++ ) {
			int idx = al.at ( r, c ) + vertexOffsets[c];
			G.at ( r, idx )  = 1;
		}
	}


	if ( nColVertices > 0 ) {
		int colidx = 2;
		for ( int col =0; col<ncols; col++ ) {
			for ( int sx=0; sx<s[col]; sx++ ) {
				int sel = vertexOffsets[col] + sx;
				G.at ( colOffset+col, sel ) = colidx; // = colidx;
			}
		}
	}

	// The non-row vertices do not have any connections to other non-row vertices.

	int ss = std::accumulate ( s.begin(), s.end(), 0 ); // std::accumulate<int>(s.begin(), s.end(), 0);
	array_link xy ( 2, nrows+ ss, 0 ) ;


// calculate positions
	for ( int row=0; row<nrows; row++ ) {
		xy.at ( 0,row ) = 0 ;
		xy.at ( 1, row ) = row;
	}
	int pos = nrows;

	for ( int col=0; col<ncols; col++ ) {
		for ( int ss=0; ss<s[col]; ss++ ) {
			xy.at ( 0,pos ) = 2 + ss / s[col];
			xy.at ( 1,pos ) = col;
			pos = pos + 1;
		}
	}


	return std::pair<array_link, std::vector<int> > ( G, colors );
}


/* This program prints generators for the automorphism group of an
   n-vertex polygon, where n is a number supplied by the user.

   This version uses dynamic allocation.
*/

namespace nauty
{
#include "nauty.h"
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

template<class Type>
std::vector<Type> uniquevec(const std::vector<Type> &v)
{
std::vector<Type> w = v;
// fixme: better option for this?
typename std::vector< Type >::iterator last = std::unique(w.begin(), w.end());
                w.erase(last,  w.end());
return w;
	
}

/// set colors in nauty format
void setcolors ( std::vector<int> colors, int *lab, int *ptn )
{
	const int verbose=0;
	const int n = colors.size();
	std::vector<int> ucols =uniquevec( colors );
	if (verbose)
	printf ( "setcolors: found %d/%d unique colors\n" , (int)ucols.size(), n);
	if (verbose) {
display_vector(colors); printf("\n");
	}
	int x=-1;
	for ( size_t i=0; i<ucols.size(); i++ ) {
		if (verbose)
	printf ( "color %d: %d\n", (int)i, ucols[i] );
		for ( size_t j=0; j<colors.size(); j++ ) {
			if ( colors[j]==ucols[i] ) {
				x++;
				lab[x]=j;
				ptn[x]=1;
			}
		}
		ptn[x]=0;
	}
	if (verbose) {
	printf("lab and ptn:\n");
	print_perm(lab, n);
	print_perm(ptn, n);
	}	
}

	
/// testing function
std::vector<int> reduceNauty ( const array_link &G, std::vector<int> colors )
{
	int verbose=1;

	// FIXME: colors of vertices are ignored???
	
	//for(size_t j=0; j<colors.size(); j++) colors[j]=j;
	if (verbose) {
		printf("reduceNauty: %d vertices\n", G.n_rows);
		   printf("colors: "); display_vector(colors); printf("\n");

}
	int nvertices=G.n_rows;

	/* DYNALLSTAT declares a pointer variable (to hold an array when it
	   is allocated) and a size variable to remember how big the array is.
	   Nothing is allocated yet.  */

	DYNALLSTAT ( graph,g,g_sz );
	DYNALLSTAT ( int,lab,lab_sz );
	DYNALLSTAT ( int,ptn,ptn_sz );
	DYNALLSTAT ( int,orbits,orbits_sz );
	static DEFAULTOPTIONS_GRAPH ( options );
	statsblk stats;

	int m,v;
	set *gv;

	/* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
	   Here we change those options that we want to be different from the
	   defaults.  writeautoms=TRUE causes automorphisms to be written. */

	options.writeautoms = TRUE;
	options.writeautoms = FALSE;

	
	int n = nvertices;


	/* The nauty parameter m is a value such that an array of
	   m setwords is sufficient to hold n bits.  The type setword
	   is defined in nauty.h.  The number of bits in a setword is
	   WORDSIZE, which is 16, 32 or 64.  Here we calculate
	   m = ceiling(n/WORDSIZE). */

	m = SETWORDSNEEDED ( n );

	/* The following optional call verifies that we are linking
	   to compatible versions of the nauty routines. */

	nauty_check ( WORDSIZE,m,n,NAUTYVERSIONID );

	/* Now that we know how big the graph will be, we allocate
	 * space for the graph and the other arrays we need. */

	DYNALLOC2 ( graph,g,g_sz,m,n,"malloc" );
	DYNALLOC1 ( int,lab,lab_sz,n,"malloc" );
	DYNALLOC1 ( int,ptn,ptn_sz,n,"malloc" );
	DYNALLOC1 ( int,orbits,orbits_sz,n,"malloc" );

	EMPTYGRAPH ( g,m,n );

	for ( int ix=0; ix<nvertices; ix++ ) {
		for ( int iy=0; iy<nvertices; iy++ ) {
			if ( G.at ( ix,iy ) >0 ) {

				ADDONEEDGE ( g, ix, iy, m );
			}
		}
	}

	for ( v = 0; v < n; ++v )
		ADDONEEDGE ( g,v, ( v+1 ) %n,m );

	
	setcolors(colors, lab, ptn);
	options.defaultptn=false;
	printf("options.defaultptn: %d\n ", options.defaultptn);
	
	densenauty ( g,lab,ptn,orbits,&options,&stats,m,n,NULL );
	if ( verbose>=2 ) {
		printf ( "Generators for Aut(C[%d]):\n",n );

		printf ( "order = " );
		writegroupsize ( stdout,stats.grpsize1,stats.grpsize2 );
		printf ( "\n" );
	}




	std::vector<int> tr ( nvertices );
	std::copy ( lab, lab+nvertices, tr.begin() );

	// TODO: use colors in graph
	//  delete allocated variables
	DYNFREE ( g,g_sz );
	DYNFREE ( lab,lab_sz );
	DYNFREE ( ptn,ptn_sz );
	DYNFREE ( orbits,orbits_sz );

	return tr;
}


std::vector<int> reduceOAnauty(const array_link &al, int verbose=1)
{
	std::pair<array_link, std::vector<int> > Gc = array2graph ( al,  verbose );

	array_link &G = Gc.first;
	std::vector<int> &colors = Gc.second;

std::vector<int> r = reduceNauty ( G, colors );
	return r;
}
	
	
} // end of nauty namespace

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

	std::vector<int> tr = nauty::reduceOAnauty ( al );
	printf ( "canon: " );
	display_vector ( tr );
	printf ( "\n" );

array_transformation_t ttm = labelling2transformation(tr, arrayclass, verbose);
	ttm.show();
	
	array_link alm = ttm.apply(al);
	
	al.showarraycompact();
	printf("---\n");
alm.showarraycompact();

	
	array_link alx = al.randomperm();
	std::vector<int> trx = nauty::reduceOAnauty ( alx, 0 );
	array_transformation_t ttx = labelling2transformation(trx, arrayclass, 0);
	array_link alxm = ttx.apply(alx);
	
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
