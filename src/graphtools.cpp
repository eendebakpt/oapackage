#include <string>


#include "arraytools.h"
#include "mathtools.h"
#include "tools.h"

#include "graphtools.h"

#ifdef RPACKAGE
#define printf notallowed
#endif


/* Interface to Nauty code
 *
 */

namespace nauty
{
#include "nauty.h"
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

template<class Type>
std::vector<Type> uniquevec ( const std::vector<Type> &v )
{
	std::vector<Type> w = v;
// fixme: better option for this?
	typename std::vector< Type >::iterator last = std::unique ( w.begin(), w.end() );
	w.erase ( last,  w.end() );
	return w;

}

/// set colors in nauty format
void setcolors ( std::vector<int> colors, int *lab, int *ptn )
{
	const int verbose=0;
	const int n = colors.size();
	std::vector<int> ucols =uniquevec ( colors );
	if ( verbose )
		printf ( "setcolors: found %d/%d unique colors\n" , ( int ) ucols.size(), n );
	if ( verbose ) {
		display_vector ( colors );
		printf ( "\n" );
	}
	int x=-1;
	for ( size_t i=0; i<ucols.size(); i++ ) {
		if ( verbose )
			printf ( "color %d: %d\n", ( int ) i, ucols[i] );
		for ( size_t j=0; j<colors.size(); j++ ) {
			if ( colors[j]==ucols[i] ) {
				x++;
				lab[x]=j;
				ptn[x]=1;
			}
		}
		ptn[x]=0;
	}
	if ( verbose ) {
		printf ( "lab and ptn:\n" );
		print_perm ( lab, n );
		print_perm ( ptn, n );
	}
}


std::vector<int> reduceNauty ( const array_link &G, std::vector<int> colors )
{
	int verbose=1;


	//for(size_t j=0; j<colors.size(); j++) colors[j]=j;
	if ( verbose ) {
		printf ( "reduceNauty: %d vertices\n", G.n_rows );
		printf ( "colors: " );
		display_vector ( colors );
		printf ( "\n" );

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

	//for ( v = 0; v < n; ++v )
	//	ADDONEEDGE ( g,v, ( v+1 ) %n,m );


	setcolors ( colors, lab, ptn );
	options.defaultptn=false;
	
	if (verbose>=2) {
	printf ( "options.defaultptn: %d\n ", options.defaultptn );
	printf ( "lab: \n " );
	print_perm ( lab, n );
	printf ( "ptn: \n " );
	print_perm ( ptn, n );
	}
	
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


std::vector<int> reduceOAnauty ( const array_link &al, int verbose )
{
	std::pair<array_link, std::vector<int> > Gc = array2graph ( al,  verbose );

	array_link &G = Gc.first;
	std::vector<int> &colors = Gc.second;

	//for(int i=0; i<colors.size(); i++) colors[i]=i; // hack

	if (verbose>=2) {
	printf ( "graph:\n" );
	G.showarray();
	}

	std::vector<int> r = reduceNauty ( G, colors );
	return r;
}


} // end of nauty namespace

template<class IntType>
std::vector<int> indexvector ( const std::vector<IntType> s )
{
	//printf("indexvector: input "); print_perm(s);
	int n = std::accumulate ( s.begin(), s.end(), 0 );

	std::vector<int> v ( n );

	int x=0;
	for ( int i=0; i<s.size(); i++ ) {
		for ( int j=0; j<s[i]; j++ ) {
			v[x]=i;
			x++;
		}
	}
	//printf("indexvector: "); print_perm(v);
	return v;
}

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair<array_link, std::vector<int> >  array2graph ( const array_link &al, int verbose )
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


	// other vertices: color 2+...
	std::vector<int> v = indexvector ( s );
	for ( int i=0; i<nColumnLevelVertices; i++ )
		colors[i+nrows+ncols]=2+v[i];


	if ( verbose )
		printf ( "array2graph: generating graph of size %d=%d+%d+%d\n", nVertices, nrows, ncols, nColumnLevelVertices );
	array_link G ( nVertices, nVertices, 0 ); // graph
	G.setconstant ( 0 );

//    im = np.zeros((nVertices, nVertices))  # incidence matrix


	for ( int r=0; r<nrows; r++ ) {
		for ( int c=0; c<ncols; c++ ) {
			int idx = al.at ( r, c ) + vertexOffsets[c];
			G.at ( r, idx )  = 1;
			G.at ( idx, r )  = 1;
		}
	}


	if ( nColVertices > 0 ) {
		int colidx = 2;
		for ( int col =0; col<ncols; col++ ) {
			for ( int sx=0; sx<s[col]; sx++ ) {
				int sel = vertexOffsets[col] + sx;
				G.at ( colOffset+col, sel ) = colidx; // = colidx;
				G.at ( sel, colOffset+col ) = colidx; // = colidx;
			}
		}
	}

	// The non-row vertices do not have any connections to other non-row vertices.

	int ss = std::accumulate ( s.begin(), s.end(), 0 ); // std::accumulate<int>(s.begin(), s.end(), 0);

	if ( 0 ) {
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
	}

	return std::pair<array_link, std::vector<int> > ( G, colors );
}



// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
