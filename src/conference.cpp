/** \file conference.cpp

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "arraytools.h"
#include "arrayproperties.h"
#include "tools.h"
#include "extend.h"

#include "oadevelop.h"
#include "lmc.h"

#include "conference.h"

conference_t::conference_t ( int N, int k )
{
	this->N = N;
	this->ncols = k;
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

std::vector<int> get_comb ( std::vector<int> p, int n, int zero=0, int one=1 )
{
	std::vector<int> c ( n );
	for ( int i=0; i<n ; i++ )
		c[i]=zero;
	for ( size_t i=0; i<p.size() ; i++ )
		c[p[i]]=one;
	return c;
}
void next_combx ( std::vector<int> c, int n )
{
}

cperm insertzero ( cperm c, int pos, int value=0 )
{
	cperm cx =c;


	cx.insert ( cx.begin() +pos, value );
	return cx;

}


/*

 Values in column k for k > 1 and k <= N/2 + 1:

First value 1, value at position k is 0.

From the remaining N-2 values we have q=(N-2)/2 positive. Inner product with column 0 is then satisfied.

Value v at position 1 is irrevant to innerproduct with column 1, this can be either +1 or -1

Top: 2 elements, Upper half: N/2-1 elements, bottom: N/2-1 elements

put q1 in upper half and q2 in lower half. we need:

[Case v=-1]

(inner prod 0 == 0):

q1+q2=q

(inner prod 1 == 0)

q1 - (N/2 - 1 - 1 -q1) - [ q2-(N/2-1-q2) ]  + 1= 0

--> 2q1 - 2q2  +2 = 0
--> q1 - q2  = -1

---> q1+(q1+1)=q --> 2q1 = q -1  --> q1 = q/2 - 1/2

[Case v=1] (<-- q even)

(inner prod 0 == 0):

q1+q2+1=q

(inner prod 1 == 0)

q1 - (N/2 -1 -1 -q1) - [ q2-(N/2-1-q2) ]  + 1 = 0

--> 2q1 - 2q2 + 2 = 0
--> q1 - q2 + 1 =0

---> q1+(q1+1)+1=q --> 2q1 = q - 2 --> q1 = q/2 - 1

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

 */

std::vector<cperm> get_first ( int N, int extcol, int verbose=1 )
{
	int n1=N/2-1;
	int k1=-1;

	int haszero=extcol<n1+2;
	int q = -1;
	int q1=-1, q2=-1;
	int v=-100;
	if ( haszero ) {

		n1=n1-1;
		k1=n1/2;
		q  = ( N-2 ) /2;
	} else {
		printf ( "get_first: extcol %d, N %d, n1 %d, not implemented...\n", extcol, N, n1 );
		k1 = n1/2;
	}

	std::vector<int> c ( k1 );
	for ( int i=0; i<k1 ; i++ )
		c[i]=i;



	// if q is even
	if ( q%2==0 ) {
		q1 = q/2-1;
		q2 = q-q1-1;
		v=1;
	} else {
		q1 = ( q-1 ) /2;
		q2 = q1+1;

		v=-1;
	}
	int nc = ncombs<long> ( n1, q1 );
	if ( verbose )
		printf ( "get_first: N %d, n1 %d, q %d, v %d, q1 %d, q2 %d, nc %d\n", N, n1, q, v, q1, q2, nc );

	std::vector<cperm> ff;
	for ( long j=0; j<nc; j++ ) {
		cperm cc =get_comb ( c, n1, -1, 1 );

		cc=insertzero ( cc, 0, 1 );
		cc=insertzero ( cc, 1, v );

		if ( haszero )
			cc=insertzero ( cc, extcol );
		ff.push_back ( cc );

		if ( j+1<nc ) {
			next_comb ( c, q1, n1 );
		}
	}


	return ff;

}
std::vector<cperm> get_second ( int N, int extcol, int target, int verbose=0 )
{
	int n1=N/2-1;
	int k1=-1;

	int haszero=extcol>=n1+2;
	if ( haszero ) {
		printf ( "get_second: not  implemented...\n" );
		k1=n1/2;
	} else {
		k1 = n1- ( n1-target ) /2;
	}

	std::vector<int> c ( k1 );
	for ( int i=0; i<k1 ; i++ )
		c[i]=i;

	if ( verbose )
		printf ( "get_second: N %d, n1 %d, target %d, k1 %d\n", N, n1, target, k1 );

	int nc = ncombs<long> ( n1, k1 );
	std::vector<cperm> ff;
	for ( long j=0; j<nc; j++ ) {
		cperm cc =get_comb ( c, n1, -1, 1 );

		if ( haszero )
			cc=insertzero ( cc, extcol );
		ff.push_back ( cc );
		next_comb ( c, k1, n1 );
	}


	return ff;

}


/// calculate inner product between twee permutations
int innerprod ( const cperm a, const cperm b )
{
	int ip=0;
	for ( size_t i=0; i<b.size(); i++ ) {
		//printf("innerprod %d: %d += %d + %d\n", (int)i, ip, a[i], b[i] );
		ip+= a[i] * b[i];
	}
	return ip;
}

int satisfy_symm ( const cperm c, const symmdata sd )
{
	const int verbose=0;

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
					printf ( "  discard i %d, k %d, c[i]=%d:   %d %d\n", ( int ) i, k, c[i], sd.rowvalue.atfast ( i, k ), sd.rowvalue.atfast ( i+1, k ) );
				}
				return false;
			}
		}

	}
	return true;
}

cperm getColumn ( const array_link al, int c )
{
	cperm cx ( al.n_rows );
	std::copy ( al.array+c*al.n_rows, al.array+ ( c+1 ) *al.n_rows, cx.begin() );
	return cx;
}

// return true if the extension column satisfies the inner product check
int ipcheck ( const cperm col, const array_link &al, int cstart=2 )
{

	for ( int c=cstart; c<al.n_columns; c++ ) {
		cperm cx = getColumn ( al, c );
		if ( innerprod ( col, cx ) !=0 ) {
			return false;
		}
	}
	return true;
}

conference_extend_t extend_conference ( const array_link al, const conference_t ct, int extcol, int verbose )
{
	// first create two sets of partial extensions

	conference_extend_t ce;
	ce.extensions.resize ( 0 );
	std::vector<cperm> extensions ( 0 );

	int N = ct.N;

	int k = extcol;

	// loop over all possible first combinations
	std::vector<cperm> ff = get_first ( N, extcol );
	//std::vector<cperm> ff2 = get_second(N, extcol);

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

	printf ( "--- extend_conference: extcol %d ---\n", extcol );
	cperm c0 = getColumn ( al, 0 );
	cperm c1 = getColumn ( al, 1 );
	for ( size_t i=0; i<ce.first.size(); i++ ) {
		int ip = innerprod ( c0, ce.first[i] );
		//int ip = innerprod(c1, ce.first[i]);
		//printf("extend1 %d: inner product %d\n", (int)i, ip);

		int target = -ip;
		std::vector<cperm> ff2 = get_second ( N, extcol, target, verbose>=2 );
		ce.second=ff2;

		for ( size_t j=0; j<ff2.size(); j++ ) {
			cperm c = ce.combine ( i, j );
			int ip0 = innerprod ( c0, c );
			int ip1 = innerprod ( c1, c );
			//printf("extend %d: N %d ", (int)i, N); display_vector(c);	 printf("\n");
			if ( verbose>=2 ) {
				printf ( "extend %d: ip %d %d\n", ( int ) i, ip0, ip1 );
			}

			if ( verbose>=2 ) {
				//alx.showarraycompact();
				array_link ecol=al.selectFirstColumns ( 1 );
				array_link alx = hstack ( al, ecol );
				alx.setColumn ( extcol, c );
				alx.showarray();
			}
			// add array to good set if ip2 is zero
			if ( ip0==0 && ip1==0 )
				extensions.push_back ( c );
		}
	}

	ce.extensions=extensions;
	if ( verbose>=2 )
		printf ( "after generation: found %d extensions\n", ( int ) extensions.size() );
	// perform row symmetry check

	symmetry_group rs = al.row_symmetry_group();
	symmdata sd ( al.selectFirstColumns ( k ) );


	if ( verbose>=2 )
		sd.show ( 1 );

	std::vector<cperm> e2 ( 0 );
	for ( size_t i=0; i<extensions.size(); i++ ) {

		if ( satisfy_symm ( extensions[i], sd ) ) {

			// perform inner product check for all columns
			if ( ipcheck ( extensions[i], als ) ) {
				e2.push_back ( extensions[i] );
			} else {
				printf ( "  reject due to innerproduct\n" );
			}
		} else {
			if ( verbose>=2 ) {
				printf ( "  reject due to row symm: " );
				display_vector ( extensions[i] );
				printf ( "\n" );
			}
		}
	}

	if ( verbose>=1 )
		printf ( "extend_conference: symmetry check %d->%d\n", ( int ) extensions.size(), ( int ) e2.size() );


	ce.extensions=e2;

	return ce;
}

/// helper function
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


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
