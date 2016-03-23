/** \file conference.cpp

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>

#include "arraytools.h"
#include "arrayproperties.h"
#include "extend.h"

#include "oadevelop.h"
//#include "lmc.h"

#include "graphtools.h"

#include "conference.h"

/*

######## Extension numbers

When generating candidate extension columns many properties of the column are determined by N and the position of the zero in the column. (If the first 2 columns are in normal form).
Below we write down the conditions and solve the resulting equations.

#### Zero in middle block

Values in column k for k > 1 and k <= N/2 + 1 = q + 2:

First value 1, value at position k is 0.

From the remaining N-2 values we have q=(N-2)/2 positive. Inner product with column 0 is then satisfied.

Value v at position 1 is irrevant to innerproduct with column 1, this can be either +1 or -1

Top: 2 elements, Upper half: N/2-1 elements, bottom: N/2-1 elements

put q1 in upper half and q2 in lower half. we need:

[Case v=-1]

(inner prod with col 0 == 0):

q1+q2=q

(inner prod with col 1 == 0)

q1 - (N/2 - 1 - 1 -q1) - [ q2-(N/2-1-q2) ]  + 1= 0

--> 2q1 - 2q2  +2 = 0
--> q1 - q2  = -1

---> q1+(q1+1)=q --> 2q1 = q -1  --> q1 = q/2 - 1/2 = (q-1)/2

[Case v=1] (<-- q even/odd?)

(inner prod 0 == 0):

q1+q2+1=q

(inner prod 1 == 0)

q1 - (N/2 -1 -1 -q1) - [ q2-(N/2-1-q2) ]  + 1 = 0

--> 2q1 - 2q2 + 2 = 0
--> q1 - q2 + 1 =0

---> q1+(q1+1)+1=q --> 2q1 = q - 2 --> q1 = q/2 - 1

Example:

N=8
q=3
q1=1
q2=2
v=1

N=6
q=2
q1=0
q2=1
n1=1
v=1

#### Zero in last block

Values in column k for k > 1 and k > N/2 + 1 = 2 + q:
First value 1, value at position k is 0.

From the remaining N-2 values we have q=(N-2)/2 positive. Inner product with column 0 is then satisfied.

Value v at position 1 is irrevant to innerproduct with column 1, this can be either +1 or -1

Top: 2 elements, Upper half: N/2-1 elements, bottom: N/2-1 elements

The zero is in the bottom part. Let q1 and q2 be the number of +1 signs in the upper and lower half, respectively.
We need:

[Case v=-1]

(inner prod with col 0 == 0):

q1+q2=q

(inner prod with col 1 == 0)

1 + q1 + (q-1-q2) = q

[Case v=1] (<-- q even)

(inner prod 0 == 0):

q1+q2+1=q

(inner prod 1 == 0)

1 + q1 + (q-1-q2) = q

==> v=+1: q1=(q-1)/2, q2=q1
==> v=-1: q1=q/2, q2=q1

Examples:

N=8
q=3
q1=?
q2=?
v=?

N=6
q=2
q1=?0
q2=?
n1=?
v=1

 */

/// return true of the argument is even
inline bool iseven ( int q )
{
	return ( q%2 ) ==0;
}

/// return parameters of a conference design
void getConferenceNumbers ( int N,int k, int &q, int &q1, int &q2, int &v )
{
	q = ( N-2 ) /2;

	if ( k< 2+q ) {
		if ( iseven ( q ) ) {
			q1=q/2-1;
			v=1;
			q2=q1+1;
		} else {
			q1= ( q-1 ) /2;
			v=-1;
			q2=q1+1;
		}
	} else {

		if ( iseven ( q ) ) {
			q1=q/2;
			v=-1;
			q2=q1;
		} else {
			q1= ( q-1 ) /2;
			v=1;
			q2=q1;
		}
	}
	//printfd ( "getConferenceNumbers: k %d, q %d: q1 q2 %d, %d\n", k, q, q1, q2 );
}

conference_t::conference_t ( int N, int k )
{
	this->N = N;
	this->ncols = k;
	this->ctype = CONFERENCE_NORMAL;
}

array_link conference_t::create_root_three ( ) const
{
	array_link al ( this->N, 3, 0 ); // c.ncols
	
		al.at ( 0,0 ) =0;
	for ( int i=1; i<this->N; i++ ) {
		al.at ( i, 0 ) = 1;
	}
	for ( int i=0; i<this->N; i++ ) {
		if ( i==1 ) {
			al.at ( 1,1 ) =0;
			continue;
		}
		if ( i<=this->N/2 )
			al.at ( i, 1 ) = 1;
		else
			al.at ( i, 1 ) = -1;		
	}
	
	int q, q1, q2, v;
	const int k = 2;
	const int N = this->N;
	getConferenceNumbers ( this->N, k, q, q1, q2, v );

	for(int i=3;i<N; i++) al.at(i,2)=-1;
	al.at(0,2)=1;
	al.at(1,2)=v;
	al.at(2,2)=0;
	for(int i=0; i<q1; i++) al.at(2+1+i,2) = 1;
	for(int i=0; i<q2; i++) al.at(2+q+i,2) = 1;
	
	return al;
}
array_link conference_t::create_root ( ) const
{
	array_link al ( this->N, 2, 0 ); // c.ncols

	al.at ( 0,0 ) =0;
	for ( int i=1; i<this->N; i++ ) {
		al.at ( i, 0 ) = 1;
	}
	for ( int i=0; i<this->N; i++ ) {
		if ( i==1 ) {
			al.at ( 1,1 ) =0;
			continue;
		}
		if ( i<=this->N/2 )
			al.at ( i, 1 ) = 1;
		else
			al.at ( i, 1 ) = -1;
	}

	return al;
}

conference_transformation_t reduceConferenceTransformation ( const array_link &al, int verbose )
{
		const int nr = al.n_rows;
	const int nc=al.n_columns;
	const int nn = 2* ( nr+nc );
	/// create graph
	array_link G ( 2* ( nr+nc ), 2* ( nr+nc ), array_link::INDEX_DEFAULT );
	G.setconstant ( 0 );

	if ( verbose )
		printf ( "reduceConference: %d, %d\n", nr, nc );

	std::vector<int> colors ( 2* ( nr+nc ) );

	const int roffset0=0;
	const int roffset1=nr;
	const int coffset0=2*nr;
	const int coffset1=2*nr+nc;

	/*
	Add edges as follows:
	(1)  r[i]--r'[i] for i=0..nr-1;  c[j]--c'[j] for j=0..nc-1.
	(2)  r[i]--c[j] and r'[i]--c'[j] for all A[i,j] = +1
	(3)  r[i]--c'[j] and r'[i]--c[j] for all A[i,j] = -1.
	Zeros in A don't cause any edges.
	*/

	// set colors
	for ( int i=0; i<coffset0; i++ )
		colors[i]=0;
	for ( int i=coffset0; i<coffset0+2*nc; i++ )
		colors[i]=1;

	// (1)
	for ( int i=0; i<nr; i++ )
		G.atfast ( roffset0+i, i+roffset1 ) =1;
	for ( int i=0; i<nc; i++ )
		G.atfast ( coffset0+i, i+coffset1 ) =1;

	// (2), (3)
		for ( int c=0; c<nc; c++ ) {
	for ( int r=0; r<nr; r++ ) {
			if ( al.atfast ( r,c ) ==1 ) {
				G.atfast ( roffset0+r, coffset0+c ) =1;
				G.atfast ( roffset1+r, coffset1+c ) =1;
			} else {
			if ( al.atfast ( r,c ) ==-1 ) {
				G.atfast ( roffset0+r, coffset1+c ) =1;
				G.atfast ( roffset1+r, coffset0+c ) =1;
				//G.atfast ( coffset0+c, roffset1+r ) =1;
			}
			}
		}
	}

	// make symmetryic
	/*
	for ( int i=0; i<nn; i++ )
		for ( int j=i; j<nn; j++ )
			G.atfast ( j,i ) =G.atfast ( i, j );
*/
	const int nrg = G.n_rows;
	for ( int i=0; i<nn; i++ ) {
		
		array_t *x = G.array+i*nrg; // offset to column
		for ( int j=i; j<nn; j++ ) {
			x[j] =G.array[i+j*nrg];
		}
	}

	if ( verbose>=3 ) {
		printf ( "reduceConference: incidence graph:\n" );
		printf ( "   2x%d=%d row vertices and 2x%d=%d column vertices\n", nr, 2*nr, nc, 2*nc );
		G.showarray();
	}
	/// call nauty
	const std::vector<int> tr = nauty::reduceNauty ( G, colors, verbose>=2 );
	const std::vector<int> tri = invert_permutation ( tr );
	const std::vector<int> trx=tri;

	// extract transformation
	if ( verbose>=2 ) {
		if ( verbose>=3 ) {
			array_link Gx = transformGraph ( G, tri, 0 );
			printfd ( "transformed graph\n" );
			Gx.showarray();
		}
		std::vector<int> tr1 = std::vector<int> ( trx.begin () , trx.begin () + 2*nr );
		std::vector<int> tr2 = std::vector<int> ( trx.begin () +coffset0, trx.end() );
		printf ( "  row vertex transformations: " );
		display_vector ( tr1 );
		printf ( "\n" );
		printf ( "  col vertex transformations: " );
		display_vector ( tr2 );
		printf ( "\n" );
	}

	// ...

	// define conference matrix transformation object....
	conference_transformation_t t ( al );

	// extract transformation
	std::vector<int> rr ( nr );
	for ( int i=0; i<nr; i++ ) {
		rr[i] = std::min ( trx[i], trx[i+nr] );
	}

	if ( verbose>=2 ) {
		printf ( "rr: " );
		print_perm ( rr );
	}

	t.rperm = invert_permutation ( argsort ( rr ) );

	for ( int i=0; i<nr; i++ ) {
		t.rswitch[ t.rperm[i]] = 2* ( trx[i]<trx[i+nr] )-1;
	}

	std::vector<int> cc ( nc );
	for ( int i=0; i<nc; i++ ) {
		cc[i] = std::min ( trx[coffset0+i], trx[coffset0+i+nc] );
	}
	t.cperm = invert_permutation ( argsort ( cc ) );

	for ( int i=0; i<nc; i++ ) {
		t.cswitch[ t.cperm[i]] = 2* ( trx[coffset0+i]< trx[coffset0+i+nc] ) -1;
	}

	if ( verbose>=2 ) {
		printf ( "transform: \n" );
		t.show();
	}
	return t;
}

array_link reduceConference ( const array_link &al, int verbose )
{
	const int nr = al.n_rows;
	const int nc=al.n_columns;
	const int nn = 2* ( nr+nc );
	/// create graph
	array_link G ( 2* ( nr+nc ), 2* ( nr+nc ), array_link::INDEX_DEFAULT );
	G.setconstant ( 0 );

	if ( verbose )
		printf ( "reduceConference: %d, %d\n", nr, nc );

	std::vector<int> colors ( 2* ( nr+nc ) );

	const int roffset0=0;
	const int roffset1=nr;
	const int coffset0=2*nr;
	const int coffset1=2*nr+nc;

	/*
	Add edges as follows:
	(1)  r[i]--r'[i] for i=0..nr-1;  c[j]--c'[j] for j=0..nc-1.
	(2)  r[i]--c[j] and r'[i]--c'[j] for all A[i,j] = +1
	(3)  r[i]--c'[j] and r'[i]--c[j] for all A[i,j] = -1.
	Zeros in A don't cause any edges.
	*/

	// set colors
	for ( int i=0; i<coffset0; i++ )
		colors[i]=0;
	for ( int i=coffset0; i<coffset0+2*nc; i++ )
		colors[i]=1;

	// (1)
	for ( int i=0; i<nr; i++ )
		G.atfast ( roffset0+i, i+roffset1 ) =1;
	for ( int i=0; i<nc; i++ )
		G.atfast ( coffset0+i, i+coffset1 ) =1;

	// (2), (3)
	for ( int r=0; r<nr; r++ ) {
		for ( int c=0; c<nc; c++ ) {
			if ( al.atfast ( r,c ) ==1 ) {
				G.atfast ( roffset0+r, coffset0+c ) =1;
				G.atfast ( roffset1+r, coffset1+c ) =1;
			}
			if ( al.atfast ( r,c ) ==-1 ) {
				G.atfast ( roffset0+r, coffset1+c ) =1;
				G.atfast ( roffset1+r, coffset0+c ) =1;
				//G.atfast ( coffset0+c, roffset1+r ) =1;
			}
		}
	}

	// make symmetryic
	/*
	for ( int i=0; i<nn; i++ )
		for ( int j=i; j<nn; j++ )
			G.atfast ( j,i ) =G.atfast ( i, j );
*/
	const int nrg = G.n_rows;
	for ( int i=0; i<nn; i++ ) {
		
		array_t *x = G.array+i*nrg; // offset to column
		for ( int j=i; j<nn; j++ ) {
			x[j] =G.array[i+j*nrg];
		}
	}

	if ( verbose>=3 ) {
		printf ( "reduceConference: incidence graph:\n" );
		printf ( "   2x%d=%d row vertices and 2x%d=%d column vertices\n", nr, 2*nr, nc, 2*nc );
		G.showarray();
	}
	/// call nauty
	const std::vector<int> tr = nauty::reduceNauty ( G, colors, verbose>=2 );
	const std::vector<int> tri = invert_permutation ( tr );
	const std::vector<int> trx=tri;

	// extract transformation
	if ( verbose>=2 ) {
		if ( verbose>=3 ) {
			array_link Gx = transformGraph ( G, tri, 0 );
			printfd ( "transformed graph\n" );
			Gx.showarray();
		}
		std::vector<int> tr1 = std::vector<int> ( trx.begin () , trx.begin () + 2*nr );
		std::vector<int> tr2 = std::vector<int> ( trx.begin () +coffset0, trx.end() );
		printf ( "  row vertex transformations: " );
		display_vector ( tr1 );
		printf ( "\n" );
		printf ( "  col vertex transformations: " );
		display_vector ( tr2 );
		printf ( "\n" );
	}

	// ...

	// define conference matrix transformation object....
	conference_transformation_t t ( al );

	// extract transformation
	std::vector<int> rr ( nr );
	for ( int i=0; i<nr; i++ ) {
		rr[i] = std::min ( trx[i], trx[i+nr] );
	}

	if ( verbose>=2 ) {
		printf ( "rr: " );
		print_perm ( rr );
	}

	t.rperm = invert_permutation ( argsort ( rr ) );

	for ( int i=0; i<nr; i++ ) {
		t.rswitch[ t.rperm[i]] = 2* ( trx[i]<trx[i+nr] )-1;
	}

	std::vector<int> cc ( nc );
	for ( int i=0; i<nc; i++ ) {
		cc[i] = std::min ( trx[coffset0+i], trx[coffset0+i+nc] );
	}
	t.cperm = invert_permutation ( argsort ( cc ) );

	for ( int i=0; i<nc; i++ ) {
		t.cswitch[ t.cperm[i]] = 2* ( trx[coffset0+i]< trx[coffset0+i+nc] ) -1;
	}

	if ( verbose>=2 ) {
		printf ( "transform: \n" );
		t.show();
	}
	array_link alx = t.apply ( al );
	return alx;

}

// return vector of length n with specified positions set to one
cperm get_comb ( cperm p, int n, int zero=0, int one=1 )
{
	cperm c ( n, zero );
	//for ( int i=0; i<n ; i++ )
	//	c[i]=zero;
	for ( size_t i=0; i<p.size() ; i++ )
		c[p[i]]=one;
	return c;
}

// set vector of length n with specified positions set to one
inline void get_comb ( const cperm &p, int n, int zero, int one, cperm &c )
{
	//std::fill(c.begin(), c.end(), zero);
	for ( int i=0; i<n ; i++ )
		c[i]=zero;
	for ( size_t i=0; i<p.size() ; i++ )
		c[p[i]]=one;
}

/// return copy of vector with zero inserted at specified position
inline cperm insertzero ( const cperm &c, int pos, int value=0 )
{
	cperm cx(c.size()+1);
	std::copy(c.begin(), c.begin()+pos, cx.begin() );
	cx[pos]=value;
	std::copy(c.begin()+pos, c.end(), cx.begin() +pos+1);
	return cx;
}


/** Return all admissible columns (first part) for a conference array in normal form
 *
 *
 **/
std::vector<cperm> get_first ( int N, int extcol, int verbose=1 )
{
	int k1=-1;
	int n1=-1;
	int k = extcol;

	int q, q1, q2, v;
	getConferenceNumbers ( N, k, q, q1, q2, v );

	int haszero=extcol<q+2;
	if ( haszero ) {

		n1=q-1;
		//k1=n1/2;
	} else {
		//printf ( "conference array: extcol %d, N %d, n1 %d, not implemented...\n", extcol, N, n1 );
		//k1 = n1/2;
		n1=q;
	}

	cperm c ( q1 );
	for ( int i=0; i<q1 ; i++ )
		c[i]=i;


	int nc = ncombs<long> ( n1, q1 );
	if ( verbose>=2 )
		printf ( "get_first: conference array: extcol %d: N %d, n1 %d, q %d, v %d, q1 %d, q2 %d, nc %d\n", extcol, N, n1, q, v, q1, q2, nc );

	std::vector<cperm> ff;
	for ( long j=0; j<nc; j++ ) {
		cperm cc =get_comb ( c, n1, -1, 1 );

		cc=insertzero ( cc, 0, 1 );
		cc=insertzero ( cc, 1, v );

		if ( haszero )
			cc=insertzero ( cc, extcol );
		//printfd("get_first: add element of size %d =  2 + %d\n", cc.size(), q );
		ff.push_back ( cc );

		if ( j+1<nc ) {
			next_comb ( c, q1, n1 );
		}
	}

	return ff;
}

/** Return all admissible columns (block two) for a conference array in normal form */
std::vector<cperm> get_second ( int N, int extcol, int target, int verbose=0 )
{
	//verbose=2;
	if ( verbose )
		printfd ( "get_second: N %d, extcol %d, target %d\n" );
	int k = extcol;
	int q, q1, q2, v;
	getConferenceNumbers ( N, k, q, q1, q2, v );

	int n1=-1;
	int haszero=extcol>=q+2;
	if ( verbose )
		printf ( "get_second: extcol: %d, q %d, haszero %d\n", extcol, q, haszero );

	if ( haszero ) {
		//k1=n1/2;
		n1=q-1;
	} else {
		//k1 = n1- ( n1-target ) /2;
		n1=q;
	}
	int qx=q2;

	cperm c ( qx );
	for ( int i=0; i<qx ; i++ )
		c[i]=i;


	int nc = ncombs<long> ( n1, qx );
	if ( verbose )
		printf ( "get_second: N %d, n1 %d, qx %d, target %d, nc %d\n", N, n1, qx, target, nc );
	std::vector<cperm> ff;
	cperm ccx(n1 );
	cperm cc;
	for ( long j=0; j<nc; j++ ) {
		//printf("ccc: "); display_vector(cc); printf("\n");
		//printf("ccx: "); display_vector(ccx); printf("\n");
		
		if ( haszero ) {
		get_comb ( c, n1, -1, 1, ccx );
			cc=insertzero ( ccx, extcol- ( q+2 ) );
		} else {
		cc =get_comb ( c, n1, -1, 1 );			
		}
		if ( verbose>=2 ) {
			printfd ( "add element of size %d =   %d\n", cc.size(), q );
			display_vector ( cc );
			printf ( "\n" );
			printf ( "c: " );
			display_vector ( c );
			printf ( "\n" );
		}

		ff.push_back ( cc );
		if ( n1>0 ) // guard
			next_comb ( c, qx, n1 );
	}

	return ff;
}

/// calculate inner product between two permutations
int innerprod ( const cperm &a, const array_link &al, int col )
{
	int ip=0;
	size_t nn = a.size();
	const array_t *b = al.array+col*al.n_rows;

	for ( size_t i=0; i<nn; i++ ) {
		ip+= a[i] * b[i];
	}
	return ip;
}


/// calculate inner product between two permutations
int innerprod ( const cperm &a, const cperm &b )
{
	int ip=0;
	size_t nn = b.size();
	for ( size_t i=0; i<nn; i++ ) {
		//printf("innerprod %d: %d += %d + %d\n", (int)i, ip, a[i], b[i] );
		ip+= a[i] * b[i];
	}
	return ip;
}

/// helper function
int satisfy_symm ( const cperm &c, const symmdata & sd )
{
	const int verbose=0;

	if ( verbose>=2 ) {
		printf ( "satisfy_symm: sd: " );
		sd.show();
	}
	//return true; // hack
	int k = sd.rowvalue.n_columns-1;

	if ( verbose ) {
		printf ( "satisfy_symm: " );
		display_vector ( c );
		printf ( "\n" );
	}
	for ( size_t i=2; i<c.size()-1; i++ ) {
		if ( sd.rowvalue.atfast ( i, k ) ==sd.rowvalue.atfast ( i+1, k ) ) {
			if ( c[i]<c[i+1] && c[i]!=0 and c[i+1]!=0 ) {
				// discard

				if ( verbose ) {
					printf ( "satisfy_symm: perm: " );
					display_vector ( c );
					printf ( "\n" );
					printf ( "  discard i %d, k %d, c[i]=%d:   %d %d\n", ( int ) i, k, c[i], sd.rowvalue.atfast ( i, k ), sd.rowvalue.atfast ( i+1, k ) );
				}
				return false;
			}
		}
	}
	if ( verbose>=2 ) {
		printf ( "satisfy_symm: return true\n" );
	}
	return true;
}

/// return column of an array in cperm format
cperm getColumn ( const array_link &al, int c )
{
	cperm cx ( al.n_rows );
	std::copy ( al.array+c*al.n_rows, al.array+ ( c+1 ) *al.n_rows, cx.begin() );
	return cx;
}

// return true if the extension column satisfies the inner product check
int ipcheck ( const cperm &col, const array_link &al, int cstart=2, int verbose=0 )
{
	for ( int c=cstart; c<al.n_columns; c++ ) {
		if ( innerprod ( col, al, c ) !=0 ) {
			if ( verbose ) {
				printf ( "ipcheck: column %d to %d (inclusive), failed at col %d\n", c, cstart, al.n_columns+1 );
			}
			return false;
		}
	}
	return true;
}

int maxz ( const array_link &al, int k )
{
	int maxz=-1;
	const int nr=al.n_rows;
	if ( k==-1 ) {
		for ( int k=0; k<al.n_columns; k++ ) {
			for ( int r=nr-1; r>=maxz; r-- ) {
				//printf("r k %d %d\n", r, k); al.show();
				if ( al._at ( r, k ) ==0 ) {
					maxz = std::max ( maxz, r );
				}
			}
		}
		return maxz;
	} else {
		for ( int r=nr-1; r>=0; r-- ) {
			if ( al._at ( r, k ) ==0 ) {
				return r;
			}
		}
	}
	return maxz ;
}

std::vector<cperm> filterCandidates(const std::vector<cperm> &extensions, const array_link &als, int filtersymm, int filterip, int verbose )
{
	symmetry_group rs = als.row_symmetry_group();
	symmdata sd ( als );


	if ( verbose>=2 )
		sd.show ( 1 );

	std::vector<cperm> e2 ( 0 );
	for ( size_t i=0; i<extensions.size(); i++ ) {

		if ( filtersymm ) {
			if ( ! satisfy_symm ( extensions[i], sd ) ) {
				if ( verbose>=2 ) {
					printf ( "filterCandidates: reject due to row symm: " );
					display_vector ( extensions[i] );
					printf ( "\n" );
				}
				continue;
			}
		}
		if ( filterip ) {

			// perform inner product check for all columns
			if ( ! ipcheck ( extensions[i], als ) ) {
				if ( verbose>=2 ) {
					printf ( "   extension " );
					display_vector ( extensions[i] );
					printf ( "\n" );
					printf ( "filterCandidates: reject due to innerproduct (extension %d)\n", ( int ) i );
					ipcheck ( extensions[i], als, 2, 1 );
				}
				continue;
			}
		}
		e2.push_back ( extensions[i] );
	}
		return e2;
	}

std::vector<cperm> generateConferenceExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose , int filtersymm, int filterip)
{
	conference_extend_t ce;
	std::vector<cperm> extensions ( 0 );
	const int N = ct.N;
	const int extcol=al.n_columns;

	// loop over all possible first combinations
	std::vector<cperm> ff = get_first ( N, kz, verbose );

	if ( verbose>=2 ) {
		for ( size_t i=0; i<ff.size(); i++ ) {
			printf ( "extend1 %d: N %d: ", ( int ) i, N );
			display_vector ( ff[i] );
			printf ( "\n" );
		}
	}

	ce.first=ff;
	ce.second=ff;

	array_link als = al.selectFirstColumns ( extcol );

	cperm c0 = getColumn ( al, 0 );
	cperm c1 = getColumn ( al, 1 );
	for ( size_t i=0; i<ce.first.size(); i++ ) {
		int ip = innerprod ( c0, ce.first[i] );
		//int ip = innerprod(c1, ce.first[i]);
		//printf("extend1 %d: inner product %d\n", (int)i, ip);

		// TODO: cache this function call
		int target = -ip;

		std::vector<cperm> ff2 = get_second ( N, kz, target, verbose>=2 );
		ce.second=ff2;

		//printfd("ce.second[0] "); display_vector( ce.second[0]); printf("\n");

		for ( size_t j=0; j<ff2.size(); j++ ) {
			cperm c = ce.combine ( i, j );

#ifdef OADEBUG
#else
			if ( 1 ) {
				extensions.push_back ( c );
				continue;
			}
#endif
			int ip0 = innerprod ( c0, c );
			int ip1 = innerprod ( c1, c );
			//printf("extend %d: N %d ", (int)i, N); display_vector(c);	 printf("\n");

			if ( verbose>=2 ) {
				printf ( "extend_conference %d: ip %d %d\n", ( int ) i, ip0, ip1 );
			}

			if ( verbose>=2 ) {
				//alx.showarraycompact();
				array_link ecol=al.selectFirstColumns ( 1 );
				array_link alx = hstack ( al, ecol );
				alx.setColumn ( kz, c );
				alx.showarray();
			}
			// add array to good set if ip2 is zero
			if ( ip0==0 && ip1==0 ) {
				extensions.push_back ( c );
			} else {
				printfd ( "huh?" );

				printf ( " ip0 ip1 %d %d\n" , ip0, ip1 );
				al.show();
				printf ( "  c: %d\n", ( int ) c.size() );
				array_link alx = hstack ( al, c );
				alx.showarray();
			}
		}
	}

	ce.extensions=extensions;
	if ( verbose>=2 )
		printf ( "generateConferenceExtensions: after generation: found %d extensions\n", ( int ) extensions.size() );
	// perform row symmetry check


	std::vector<cperm> e2 = filterCandidates(extensions, al,  filtersymm,  filterip,  verbose );

	if ( verbose>=1 )
		printf ( "extend_conference: symmetry check %d + ip filter %d: %d->%d\n", filtersymm, filterip, ( int ) extensions.size(), ( int ) e2.size() );

	ce.extensions=e2;

	return e2;
}

int selectZmax(int maxzpos, const conference_t::conference_type &ctype, const array_link &al, int extcol)
{
	if (maxzpos<0) {
		if (ctype==conference_t::CONFERENCE_NORMAL)
			maxzpos = al.n_rows-1;
		else {
			if  (ctype==conference_t::conference_t::CONFERENCE_DIAGONAL) {
				maxzpos = extcol; // maxzval+2;		
				//printf("ct.ctype==conference_t::conference_t::CONFERENCE_DIAGONAL: maxzpos %d/%d, extcol %d\n", maxzpos, al.n_rows-1, extcol);
			}
				else {
				//
				printfd("not implemented...");
					maxzpos = al.n_rows-1;
			}
		}
	}
	return maxzpos;
}

conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t & ct, int extcol, int verbose, int maxzpos )
{
	conference_extend_t ce;
	ce.extensions.resize ( 0 );

	const int N = ct.N;
	const int k = extcol;
	const int maxzval = maxz ( al );

	if ( verbose )
		printf ( "--- extend_conference_matrix: extcol %d, maxz %d ---\n", extcol, maxzval );

	const int zstart=maxzval+1;
	
	maxzpos = selectZmax(maxzpos, ct.ctype, al, extcol);
		
	//for ( int ii=maxzval+1; ii<std::min<int>(al.n_rows, maxzval+2); ii++ ) {
	for ( int ii=zstart; ii<maxzpos+1; ii++ ) {
		if ( verbose>=2 )
			printf ( "array: kz %d: generate\n", ii );
		std::vector<cperm> extensionsX  = generateConferenceExtensions ( al, ct, ii, verbose, 1, 1 );

		if ( verbose>=2 )
			printf ( "array: kz %d: %d extensions\n", ii, ( int ) extensionsX.size() );
		ce.extensions.insert ( ce.extensions.end(), extensionsX.begin(), extensionsX.end() );
	}

	return ce;
}

conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t & ct, int extcol, int verbose, int maxzpos, const conf_candidates_t &cande )
{
	conference_extend_t ce;
	ce.extensions.resize ( 0 );

	const int N = ct.N;
	const int k = extcol;
	const int maxzval = maxz ( al );

	if ( verbose )
		printf ( "--- extend_conference_matrix: extcol %d, maxz %d ---\n", extcol, maxzval );

	const int zstart=maxzval+1;
	
	maxzpos = selectZmax(maxzpos, ct.ctype, al, extcol);
		
	//for ( int ii=maxzval+1; ii<std::min<int>(al.n_rows, maxzval+2); ii++ ) {
	for ( int ii=zstart; ii<maxzpos+1; ii++ ) {
		if ( verbose>=2 )
			printf ( "array: kz %d: generate\n", ii );
		//std::vector<cperm> extensionsX  = generateConferenceExtensions ( al, ct, ii, verbose, 1, 1 );
	std::vector<cperm> extensionsX  = filterCandidates ( cande.ce[ii],  al,1, 1, verbose);

		if ( verbose>=2 )
			printf ( "array: kz %d: %d extensions\n", ii, ( int ) extensionsX.size() );
		ce.extensions.insert ( ce.extensions.end(), extensionsX.begin(), extensionsX.end() );
	}

	return ce;
}


/// sort rows of an array based on the zero elements
array_link sortrows ( const array_link al )
{
	size_t nr=al.n_rows;
	// initialize original index locations
	std::vector<size_t> idx ( nr );
	for ( size_t i = 0; i != nr; ++i )
		idx[i] = i;

	//compfunc = ..;
	// sort indexes based on comparing values in v
// sort(idx.begin(), idx.end(), compfunc );

	printfd ( "not implemented...\n" );
	return al;
}



/**
void test_comb ( int n, int k )
{
	std::vector<int> c ( k );
	for ( int i=0; i<k ; i++ )
		c[i]=i;

	int nc = ncombs<long> ( n, k );

	for ( long j=0; j<nc; j++ ) {
		std::vector<int> cc =get_comb ( c, n );
		display_vector ( cc );
		printf ( "\n" );
		next_comb ( c, k, n );
	}
}
*/

template<typename T>
size_t vectorsizeof(const typename std::vector<T>& vec)
{
    return sizeof(T) * vec.size();
}


conf_candidates_t generateCandidateExtensions(const conference_t ctype, int verbose=1) {
	conf_candidates_t cande;

		cande.ce.resize(ctype.N);
		
	// TODO: use 3-column root array
	array_link al = ctype.create_root();
	array_link al3 = ctype.create_root_three();
	for(int i=3; i<ctype.N; i++) {
		std::vector<cperm> ee = generateConferenceExtensions ( al3, ctype, i, 0, 0, 1);
	
	//printf("al3:\n"); al3.showarray();
		
	if ( (long)vectorsizeof(ee)>(long(4)*1024*1024*1024)/(long)ctype.N) {
		printfd("generateCandidateExtensions: set of generated candidates too large, aborting");
			assert(0);
			exit(0);
	}
cande.ce[i] = ee;
	}
	
cande.ce[2] = generateConferenceExtensions ( al, ctype, 2, 0, 0, 1);
	
	cande.info(verbose);

	return cande;
}
arraylist_t extend_conference ( const arraylist_t &lst, const conference_t ctype, int verbose )
{
	arraylist_t outlist;

	if ( verbose>=2 ) {
		printfd ( "extend_conference: start %d\n", ( int ) lst.size() );

	}

	int vb=std::max ( 0, verbose-1 );

	conf_candidates_t cande = generateCandidateExtensions(ctype, verbose);

	for ( size_t i=0; i<lst.size(); i++ ) {
		const array_link &al = lst[i];
		int extcol=al.n_columns;
		conference_extend_t ce = extend_conference_matrix ( al, ctype, extcol, vb, -1, cande );
		//conference_extend_t ce = extend_conference_matrix ( al, ctype, extcol, vb, -1);


		arraylist_t ll = ce.getarrays ( al );
		const int nn = ll.size();

		arraylist_t::iterator it = outlist.end();

		outlist.insert ( it, ll.begin(), ll.end() );

		if ( verbose>=2 || (verbose>=1 && (i%50==0 || i==lst.size()-1)  ) ) {
			printf ( "extend_conference: extended array %d/%d to %d arrays\n", ( int ) i, ( int ) lst.size(), nn );
		}
	}
	return outlist;
}

std::pair<arraylist_t, std::vector<int> > selectConferenceIsomorpismHelper ( const arraylist_t lst, int verbose )
{
	const int nn = lst.size();

	arraylist_t lstr;
	arraylist_t lstgood;
	//printf ( "read %d arrays\n" , (int)lst.size());

	for ( int i=0; i< ( int ) lst.size(); i++ ) {
		if ( verbose>=1 && i%1500==0 )
			printf ( "selectConferenceIsomorpismClasses: reduce %d/%d\n", i, ( int ) lst.size() );
		array_link alx = reduceConference ( lst[i], verbose>=2 );
		lstr.push_back ( alx );
	}

	// perform stable sort
	indexsort sortidx ( lstr );

	const std::vector<int> &idx = sortidx.indices;

	std::vector<int> cidx ( nn );

	array_link prev;

	if ( lst.size() >0 )
		prev= lst[0];
	prev.setconstant ( -10 );

	int ci=-1;
	for ( size_t i=0; i<idx.size(); i++ ) {
		array_link al=lstr[idx[i]];
		if ( al!=prev ) {
			// new isomorphism class
			if ( verbose>=2 )
				printf ( "selectConferenceIsomorpismClasses: representative %d: index %d\n", ( int ) lstgood.size(), ( int ) idx[i] );

			lstgood.push_back (	lst[idx[i]] );
			prev=al;
			ci++;
		}
		cidx[i]=ci;
	}

	if ( verbose )
		myprintf ( "selectConferenceIsomorpismClasses: select classes %d->%d\n", ( int ) lst.size(), ( int ) lstgood.size() );

	return std::pair<arraylist_t, std::vector<int> > ( lstgood, cidx );
}

std::vector<int> selectConferenceIsomorpismIndices ( const arraylist_t lst, int verbose )
{

	std::pair<arraylist_t, std::vector<int> > pp = selectConferenceIsomorpismHelper ( lst, verbose ) ;
	return pp.second;
}

arraylist_t selectConferenceIsomorpismClasses ( const arraylist_t lst, int verbose )
{

	std::pair<arraylist_t, std::vector<int> > pp = selectConferenceIsomorpismHelper ( lst, verbose ) ;
	return pp.first;
}

/// testing function
bool compareLMC0x ( const array_link &alL, const array_link &alR )
{
	array_link L = alL;
	array_link R = alR;

	assert ( alL.n_rows==alR.n_rows );
	assert ( alL.n_columns==alR.n_columns );

	size_t nn = alL.n_columns*alL.n_rows;
	for ( size_t i=0; i<nn; i++ ) {
		if ( L.array[i]==0 )
			L.array[i]=-100;
		if ( R.array[i]==0 )
			R.array[i]=-100;
	}
	return L < R;
}

bool compareLMC0 ( const array_link &alL, const array_link &alR )
{
	assert ( alL.n_rows==alR.n_rows );
	assert ( alL.n_columns==alR.n_columns );

	for ( int c=0; c<alL.n_columns; c++ ) {
		// check position of zero in column c
		int zl = maxz ( alL, c );
		int zr = maxz ( alR, c );

		if ( zl<zr )
			return true;
		if ( zl>zr )
			return false;

		// zero is at same position in column, let LMC ordering decide

		const array_t *al = alL.array + c*alL.n_rows;
		const array_t *ar = alR.array + c*alR.n_rows;

		for ( int r=0; r<alL.n_rows; r++ ) {
			if ( al[r]> ar[r] )
				return true;	// note the reversed sign here
			if ( al[r]< ar[r] )
				return false;	// note the reversed sign here
		}
	}
	// the arrays are equal
	return false;
}


arraylist_t sortLMC0 ( const arraylist_t &lst )
{
	arraylist_t outlist = lst;
	sort ( outlist.begin(), outlist.end(), compareLMC0 );
	return outlist;
}

conference_options::conference_options(int maxpos) { maxzpos=-1; }

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
