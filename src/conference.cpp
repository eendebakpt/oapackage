/** \file conference.cpp

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stack>
#include <algorithm>

#include <vector>

#include "arraytools.h"
#include "arrayproperties.h"
#include "extend.h"

#include "graphtools.h"

#include "lmc.h"
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

conference_t::conference_t ( const conference_t  &rhs )
{
     this->N = rhs.N;
     this->ncols = rhs.ncols;
     this->ctype = rhs.ctype;
     this->itype = rhs.itype;

     this->j1zero = rhs.j1zero;
     this->j3zero = rhs.j3zero;
}

/// dummy constructor
conference_t::conference_t ( )
{
     this->N = -1;
     this->ncols = -1;
     this->ctype = CONFERENCE_NORMAL;
     this->itype = MATRIX_ISOMORPHISM;
     this->j1zero = 0;
     this->j3zero = 0;
}

conference_t::conference_t ( int N, int k, int _j1zero )
{

     this->N = N;
     this->ncols = k;
     this->ctype = CONFERENCE_NORMAL;
     this->itype = CONFERENCE_ISOMORPHISM;
     this->j1zero = _j1zero;
     this->j3zero = 0;
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

     for ( int i=3; i<N; i++ )
          al.at ( i,2 ) =-1;
     al.at ( 0,2 ) =1;
     al.at ( 1,2 ) =v;
     al.at ( 2,2 ) =0;
     for ( int i=0; i<q1; i++ )
          al.at ( 2+1+i,2 ) = 1;
     for ( int i=0; i<q2; i++ )
          al.at ( 2+q+i,2 ) = 1;

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

bool isConferenceFoldover ( const array_link &al, int verbose )
{
     // FIXME: implement reduce to conference matrix

     array_link alt = al.transposed();
     array_link alt2 = alt*-1;

     std::vector<int> ri ( al.n_rows );
     std::fill ( ri.begin(), ri.end(), -1 );

     for ( int i=0; i<al.n_rows; i++ ) {
          if ( ri[i]>-1 )
               continue;
          //array_link alx = alt.selectColumns ( i );
          int foundcol=0;
          for ( int j=i+1; j<al.n_rows; j++ ) {
               if ( ri[j]>-1 )
                    continue;
               //array_link alx2 = alt2.selectColumns ( j );
               //assert ( alt.columnEqual(i, alt2, j)==(alx==alx2) );
               if ( alt.columnEqual ( i, alt2, j ) ) {
                    //if ( alx==alx2 ) {
                    foundcol=1;
                    ri[i]=j;
                    ri[j]=i;
                    break;
               }
          }
          if ( !foundcol ) {
               if ( verbose ) {
                    printf ( "isConferenceFoldover: no foldover for row %d\n", i );
               }
               return false;
          }
     }
     return true;
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

     // make array symmetryic
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
     conference_transformation_t t = reduceConferenceTransformation ( al, verbose );
     array_link alx = t.apply ( al );
     return alx;

}

// return vector of length n with specified positions set to one
cperm get_comb ( const cperm p, int n, int zero=0, int one=1 )
{
     cperm c ( n, zero );
     //for ( int i=0; i<n ; i++ )
     //	c[i]=zero;
     for ( size_t i=0; i<p.size() ; i++ )
          c[p[i]]=one;
     return c;
}

// set vector of length n with specified positions set to one
inline void set_comb ( const cperm &p, cperm &c, int n, int zero=0, int one=1 )
{
     std::fill ( c.begin(), c.begin() +n, zero );
     for ( size_t i=0; i<p.size() ; i++ )
          c[p[i]]=one;
}

// set vector of length n with specified positions set to one
inline void get_comb ( const cperm &p, int n, int zero, int one, cperm &c )
{
     for ( int i=0; i<n ; i++ )
          c[i]=zero;
     for ( size_t i=0; i<p.size() ; i++ )
          c[p[i]]=one;
}

/// return copy of vector with zero inserted at specified position
inline cperm insertzero ( const cperm &c, int pos, int value=0 )
{
     cperm cx ( c.size() +1 );
     std::copy ( c.begin(), c.begin() +pos, cx.begin() );
     cx[pos]=value;
     std::copy ( c.begin() +pos, c.end(), cx.begin() +pos+1 );
     return cx;
}


/// set vector with zero inserted at specified position, no error checking
void insertzero ( cperm &target, const cperm &c, int pos, int value=0 )
{
     std::copy ( c.begin(), c.begin() +pos, target.begin() );
     target[pos]=value;
     std::copy ( c.begin() +pos, c.end(), target.begin() +pos+1 );
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
     } else {
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
          n1=q-1;
     } else {
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
     cperm ccx ( n1 );
     cperm cc ( n1+haszero );
     for ( long j=0; j<nc; j++ ) {
          //printf("ccc: "); display_vector(cc); printf("\n");

          if ( haszero ) {
               get_comb ( c, n1, -1, 1, ccx );
               insertzero ( cc, ccx, extcol- ( q+2 ) );
          } else {
               //cc =get_comb ( c, n1, -1, 1 );
               set_comb ( c, cc, n1, -1, 1 );

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

/// calculate inner product between partial two permutations
int partial_inner_product ( const cperm &a, const array_link &al, int col, int rmax )
{
     int ip=0;
     size_t nn = a.size();
     const array_t *b = al.array+col*al.n_rows;

     for ( int i=0; i<rmax; i++ ) {
          ip+= a[i] * b[i];
     }
     return ip;
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
inline int check_symm_zero ( const cperm &c, const std::vector<int>  & check_indices, int i )
{
     if ( check_indices[i] ) {
          if ( ( ( unsigned char ) c[i-1] ) > ( ( unsigned char ) c[i] ) ) {
               // discard
               return false;
          }
     }
     return true;
}

/// helper function, return true if a candidate extensions satisfies the symmetry test
int satisfy_symm ( const cperm &c, const std::vector<int>  & check_indices, int rowstart, int rowend )
{

     for ( int i=rowstart+1; i<rowend; i++ ) {
          if ( check_indices[i] ) {
               if ( ( ( unsigned char ) c[i-1] ) > ( ( unsigned char ) c[i] ) ) {
                    // discard
                    return false;
               }
          }
     }
     // accept
     return true;
}

/// helper function
int fix_symm ( cperm &c, const std::vector<int>  & check_indices, int rowstart, int rowend )
{
     for ( int i=rowstart+1; i<rowend; i++ ) {
          if ( check_indices[i] ) {
               if ( ( ( unsigned char ) c[i-1] ) > ( ( unsigned char ) c[i] ) ) {
                    // hack
                    if ( c[i-1]!=0 && c[i]!=0 ) {
                         // discard

                         cperm xx ( c.begin() +rowstart, c.begin() +rowend ) ;
                         printf ( "# fix_symm: input %d: ", i );
                         print_cperm ( xx );
                         printf ( "\n" );

                         std::swap ( c[i-1], c[i] );
                         cperm xx2 ( c.begin() +rowstart, c.begin() +rowend ) ;
                         printf ( "# fix_symm: output %d: ", i );
                         print_cperm ( xx2 );
                         printf ( "\n" );
                         return false;
                    }
               }
          }
     }
     // accept
     return true;
}


/// helper function, return true if a candidate extensions satisfies the symmetry test
int satisfy_symm ( const cperm &c, const std::vector<int>  & check_indices, int rowstart )
{
     //return true; // hack
//	int k = sd.rowvalue.n_columns-1;

     for ( size_t i=rowstart; i<c.size()-1; i++ ) {
          if ( check_indices[i+1] ) {
               if ( ( ( unsigned char ) c[i] ) > ( ( unsigned char ) c[i+1] ) ) {
                    // discard
                    return false;
               }
          }
     }
     // accept
     return true;
}

/// helper function, return true if a candidate extensions satisfies the symmetry test
int satisfy_symm ( const cperm &c, const symmdata & sd, int rowstart )
{
     const int verbose=0;

     if ( verbose>=2 ) {
          printf ( "satisfy_symm: sd: " );
          sd.show();
     }
     int k = sd.rowvalue.n_columns-1;

     if ( verbose ) {
          printf ( "satisfy_symm: " );
          display_vector ( c );
          printf ( "\n" );
     }
     for ( size_t i=rowstart; i<c.size()-1; i++ ) {
          // TODO: use the sd.checkIdx() for this
          if ( sd.rowvalue.atfast ( i, k ) ==sd.rowvalue.atfast ( i+1, k ) ) {
               //if ( c[i]<c[i+1] && c[i]!=0 && c[i+1]!=0 ) {
               //printf("c[i] %d, (char)c[i] %d\n", c[i], (unsigned char)c[i]);
               if ( ( ( unsigned char ) c[i] ) > ( ( unsigned char ) c[i+1] ) ) {
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
int ipcheck ( const cperm &col, const array_link &al, int cstart, int verbose )
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


int minz ( const array_link &al, int k )
{
     const int N = al.n_rows;
     int minzidx=-1;
     const int nr=al.n_rows;
     if ( k==-1 ) {
          for ( int k=0; k<al.n_columns; k++ ) {
               for ( int r=0; r<N; r++ ) {
                    if ( al._at ( r, k ) ==0 ) {
                         minzidx = std::min ( minzidx, r );
                    }
               }
          }
          return minzidx;
     } else {
          for ( int r=0; r<N; r++ ) {
               if ( al._at ( r, k ) ==0 ) {
                    return r;
               }
          }
     }
     return minzidx ;
}

int maxz ( const array_link &al, int k )
{
     int maxzidx=-1;
     const int nr=al.n_rows;
     if ( k==-1 ) {
          for ( int k=0; k<al.n_columns; k++ ) {
               for ( int r=nr-1; r>=maxzidx; r-- ) {
                    //printf("r k %d %d\n", r, k); al.show();
                    if ( al._at ( r, k ) ==0 ) {
                         maxzidx = std::max ( maxzidx, r );
                    }
               }
          }
          return maxzidx;
     } else {
          for ( int r=nr-1; r>=0; r-- ) {
               if ( al._at ( r, k ) ==0 ) {
                    return r;
               }
          }
     }
     return maxzidx ;
}

/// filter candidate extensions based on symmetry propery
std::vector<cperm> filterCandidatesSymm ( const std::vector<cperm> &extensions, const array_link &als, int verbose )
{
     const int N = als.n_rows;

     int mval=0;
     if ( N%4==2 ) {
          // S is symmetric
          mval=1;
     } else {
          // else S is anti-symmetric
          mval=-1;
     }
     if ( verbose>=2 )
          printf ( "N %d, mval %d\n", N, mval );

     const int k = als.n_columns;
     cperm tmp ( N );
     for ( int i=0; i<k; i++ )
          tmp[i] = als.at ( k, i );

     std::vector<cperm> e ( 0 );
     for ( size_t i=0; i<extensions.size(); i++ ) {
          const cperm &ex = extensions[i];

          int good=1;
          for ( int x=2; x<k; x++ ) {
               if ( tmp[x]!=mval*ex[x] ) {
                    good=0;
                    if ( k>3 && verbose>=3 ) {
                         printf ( "discard extension %d: N %d, k %d: x %d\n", ( int ) i, N, k, x );
                         printf ( "tmp: " );
                         printf_vector<signed char> ( tmp, "%d " );
                         printf ( "\n" );
                         printf ( "ex: " );
                         printf_vector<signed char> ( ex, "%d " );
                         printf ( "\n" );
                    }
                    break;
               }
          }
          if ( good )
               e.push_back ( extensions[i] );
     }
     if ( verbose>=1 )
          printf ( "filterCandidatesSymm: k %d, filter %d/%d\n", k, ( int ) e.size(), ( int ) extensions.size() );
     return e;

}

/// filter candidate extensions on J3 value (only pairs with extension candidate are checked)
std::vector<cperm> filterJ3 ( const std::vector<cperm> &extensions, const array_link &als, int verbose )
{

     const int N = als.n_rows;

     array_link dtable = createJ2tableConference ( als );

     if ( verbose>=2 ) {
          printf ( "array + dtable\n" );
          als.showarray();
          printfd ( "dtable:\n" );
          dtable.showarray();
     }

     int nc = dtable.n_columns;

     if ( verbose>=1 ) {
          printf ( "filterJ3: array %dx%d, nc %d\n", N, als.n_columns, nc );
     }

     std::vector<cperm> e2 ( 0 );
     for ( size_t i=0; i<extensions.size(); i++ ) {
          const cperm &c = extensions[i];

          //std::vector<int> cx(c.begin(), c.end() ); printf("i %d: ", (int)i); print_perm(cx);

          int jv=0;
          for ( int idx1=0; idx1<nc; idx1++ ) {
               jv=0;

               const array_t *o1 = dtable.array+dtable.n_rows*idx1;
               for ( int xr=0; xr<N; xr++ ) {
                    jv += ( o1[xr] ) * ( c[xr] );
               }

               if ( jv!=0 )
                    break;
          }

          if ( jv==0 )
               e2.push_back ( c );
     }

     if ( verbose>=1 ) {
          printf ( "filterJ3: %ld -> %ld extensions\n", ( long ) extensions.size(), ( long ) e2.size() );
     }
     return e2;
}



/** filter conferece matrix extension candidates
 *
 * Filtering is based in symmetry and ip
 */
std::vector<cperm> filterDconferenceCandidates ( const std::vector<cperm> &extensions, const array_link &als, int filtersymm, int filterip, int verbose )
{
//	symmetry_group rs = als.row_symmetry_group();
     symmdata sd ( als );
     DconferenceFilter dfilter ( als, filtersymm, filterip );
     dfilter.filterfirst=1;

     if ( verbose>=2 )
          sd.show ( 1 );

     std::vector<cperm> e2 ( 0 );
     for ( size_t i=0; i<extensions.size(); i++ ) {

          if ( dfilter.filter ( extensions[i] ) ) {
               e2.push_back ( extensions[i] );
          }
     }
     return e2;
}


/** filter conferece matrix extension candidates
 *
 * Filtering is based in symmetry and ip
 */
std::vector<cperm> filterCandidates ( const std::vector<cperm> &extensions, const array_link &als, int filtersymm, int filterip, int verbose )
{
     symmetry_group rs = als.row_symmetry_group();
     symmdata sd ( als );


     if ( verbose>=2 )
          sd.show ( 1 );

     std::vector<cperm> e2 ( 0 );
     for ( size_t i=0; i<extensions.size(); i++ ) {

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
          e2.push_back ( extensions[i] );
     }
     return e2;
}

std::vector<cperm> generateConferenceExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose , int filtersymm, int filterip )
{
     if ( ct.ctype==conference_t::DCONFERENCE ) {
          return generateDoubleConferenceExtensions ( al, ct, verbose,filtersymm,filterip );
     }

     if ( ct.itype==CONFERENCE_RESTRICTED_ISOMORPHISM ) {
          return generateConferenceRestrictedExtensions ( al, ct, kz, verbose,filtersymm,filterip );
     }


     int filterj2=1;
     int filterj3 = 0;
     int filtersymminline=1;
     std::vector<cperm> ee = generateSingleConferenceExtensions ( al, ct, kz, verbose, filtersymm,  filterj2,  filterj3,  filtersymminline );
     return ee;
}

std::vector<cperm> generateConferenceExtensionsOld ( const array_link &al, const conference_t & ct, int kz, int verbose , int filtersymm, int filterj2 )
{
     double t0=get_time_ms();

     /// legacy code

     assert ( kz>=0 ); // FIXME: refactor so we can accept kz==-1

     conference_extend_t ce;
     std::vector<cperm> extensions ( 0 );
     const int N = ct.N;
     const int extcol=al.n_columns;

     // loop over all possible first combinations
     std::vector<cperm> ff = get_first ( N, kz, verbose );

     // filter based on symmetry

     if ( filtersymm && 1 ) {
          array_link alx = al.selectFirstRows ( ( N+2 ) /2 );

          DconferenceFilter cfilter ( alx, 1, 0, 0 );
          ff= cfilter.filterList ( ff );
     }

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

     std::map<int, std::vector<cperm> > second_cache;
     std::vector<cperm> ff2;

     for ( size_t i=0; i<ce.first.size(); i++ ) {
          int ip = innerprod ( c0, ce.first[i] );
          //printfd("extend1 %d: inner product %d\n", (int)i, ip);
          int target = -ip;

          if ( second_cache.count ( target ) ==1 ) {
               std::map<int, std::vector<cperm> >::iterator it;
               it= second_cache.find ( target ); // use .find because .at is C++11
               ff2=it->second;
          }
          else {
               ff2 = get_second ( N, kz, target, verbose>=2 );
               second_cache.insert ( std::pair<int, std::vector<cperm> > ( target, ff2 ) );
          }
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

     //ce.extensions=extensions;
     if ( verbose>=2 )
          printf ( "generateConferenceExtensions: after generation: found %d extensions\n", ( int ) extensions.size() );
     // perform row symmetry check

     std::vector<cperm> e2 = filterCandidates ( extensions, al,  filtersymm,  filterj2,  verbose );

     if ( verbose>=1 )
          myprintf ( "%s: %.3f [s]: symmetry check %d + ip filter %d: %d->%d\n", __FUNCTION__, get_time_ms() - t0, filtersymm, filterj2, ( int ) extensions.size(), ( int ) e2.size() );

     ce.extensions=e2;

     return e2;
}

/// return number of +1 values in the first column of an array
int nOnesFirstColumn ( const array_link &al )
{
     int n=0;
     for ( int i=0; i<al.n_rows; i++ )
          n+= al.atfast ( i,0 ) ==1;
     return n;
}

/// helper function to sort rows of an array
indexsort rowsorter ( const array_link &al )
{
     std::vector<mvalue_t<int> > rr;
     for ( int i=0; i<al.n_rows; i++ ) {
          mvalue_t<int> m;
          for ( int k=0; k<al.n_columns; k++ )
               m.v.push_back ( al.at ( i, k ) );
          rr.push_back ( m );
     }
     indexsort is ( rr );
     return rr;
}



/// special structure for branch-and-bound generation of candidates
struct branch_t {
     int row; // current row
     int rval;
     int nvals[3];	/// number of 0, 1, and -1 remaining

     void show() const {
          myprintf ( "branch_t: row %d, val %d: %d %d %d\n", row, rval, nvals[0], nvals[1], nvals[2] );
     }
};

const int bvals[3] = {0,1,-1}; // these are ordered

template<class object_t>
/// stack object with minimal amount of memory allocations
class lightstack_t {
     object_t *stack;
private:
     int size;
     int n;
public:

     lightstack_t ( int sz ) : size ( sz ), n ( 0 ) {
          stack = new object_t[sz];
     }
     ~lightstack_t() {
          delete [] stack;
     }

     bool empty() const {
          return n==0;
     }
     void push ( const object_t &o ) {
#ifdef OADEBUG
          assert ( this->n<this->size );
#endif
          this->stack[n] = o;
          n++;
     }
     object_t & top() const {
          assert ( n>0 );
          return this->stack[n-1];
     }
     void pop() {
          n--;
     }
};

/// add new branches to a stack of branches
inline void branch ( lightstack_t<branch_t> &branches, branch_t &b, const int * const ( &bvals ), int istart = 0, int iend = 3 )
{
     for ( int i=istart; i<iend; i++ ) {

          if ( b.nvals[i]==0 )
               continue;

          branch_t bnew ( b );
          bnew.row++;
          bnew.rval = bvals[i];
          bnew.nvals[i]--;
          branches.push ( bnew );
     }
}

template <class NumType>
size_t countvalues ( const std::vector<NumType> &c, size_t start, size_t end, NumType value )
{
     size_t n=0;
     for ( size_t i=start; i<end; i++ ) {
          if ( c[i]==value )
               n++;
     }
     return n;
}

void debug_candidate ( const cperm &candidate, const std::vector<int> &check_indices, const char *str )
{
     printfd ( "%s:\n", str ) ;
     printf ( " : perm " );
     print_cperm ( candidate );
     printf ( "\n" );
     cperm xx ( check_indices.begin() , check_indices.end() ) ;
     printf ( " : chec " );
     print_cperm ( xx );
     printf ( "\n" );
}


std::vector<cperm> debug_branch ( cperm candidate, int gstart, int gend, int block, int blocksize, std::vector<int> check_indices, int showd )
{
     std::vector<cperm> cc;

     printf ( "#### debug_branch: range %d %d\n", gstart, gend );
     debug_candidate ( candidate, check_indices, "debug_branch" );
     std::sort ( candidate.begin() +gstart, candidate.begin() +gend );
     cperm candidatetmp ( candidate.begin(), candidate.end() );

     const int N = gend-gstart;
     // special construction for larger blocksizes
     lightstack_t<branch_t> branches ( 3* ( gend-gstart+2 ) );

     if ( showd )
          printfd ( "inflation of large block: block %d, blocksize %d\n", block, blocksize );

     // count items;
     // push initial branches
     branch_t b1 = {-1, -1, {-1,-1,-1} };
     b1.nvals[0] = countvalues<signed char> ( candidate, gstart, gend, bvals[0] );
     b1.nvals[1] = countvalues<signed char> ( candidate, gstart, gend, bvals[1] );
     b1.nvals[2] = countvalues<signed char> ( candidate, gstart, gend, bvals[2] );

     printf ( "  counts 1: %d 0: %d -1: %d\n", b1.nvals[0], b1.nvals[1], b1.nvals[2] );
     branches.push ( b1 );

     for ( int x=gstart; x<gend; x++ ) {
          candidatetmp[x]=-9;
     }

     double t0=get_time_ms();

     long n=0;
     do {

          branch_t b = branches.top();

          branches.pop();
          if ( b.row>=0 ) { // special check for dummy first branch element...
               candidatetmp[b.row+gstart]=b.rval;
          }

          if ( b.row>=0 ) { // special check for dummy first branch element...
               if ( check_symm_zero ( candidatetmp, check_indices, b.row+gstart ) ) {
                    // all good
                    debug_candidate ( candidatetmp, check_indices, printfstring ( "good    based on symmetry: block %d, row %d, range %d %d", block, b.row+gstart, gstart, gend ).c_str() );
               } else {
                    // discard branch

                    debug_candidate ( candidatetmp, check_indices, printfstring ( "discard based on symmetry: block %d, row %d", block, b.row+gstart ).c_str() );
                    continue;
               }
          }
          if ( b.row==N-1 ) {
               n++;
               // call the inflate function
               //inflateCandidateExtensionHelper ( list, basecandidate, candidatetmp, block+1, al, alsg, check_indices, ct, verbose, filter,ntotal );
               cc.push_back ( candidatetmp );
               continue;
          }

          // perform branching
          branch ( branches, b, bvals, 0, 3 );

     } while ( ! branches.empty() );
     if ( showd ) {
          printfd ( "inflation of large block: blocksize %d: %d branch calls\n", blocksize, n );
     }
     return cc;
}




void inflateCandidateExtensionHelper ( std::vector<cperm> &list, const cperm &basecandidate,  cperm &candidate, int block, const array_link &al,
                                       const symmetry_group &alsg, const std::vector<int> & check_indices, const conference_t & ct, int verbose , const DconferenceFilter &filter, long &ntotal )
{
     const symmdata &sd = filter.sd;
     // TODO: inline symmetry checks
     int nblocks = alsg.ngroups;

     if ( block==nblocks ) {
          // TODO: make this loop in n-1 case?

          ntotal++;
          const int j2start=std::max ( 0, al.n_columns-2 );
          // TODO: this can probably be a restricted filter (e.g. inner product check only last col and no symm check)
          if ( filter.filterJlast ( candidate, j2start ) ) {
               list.push_back ( candidate );
          }
          return;
     }

     const int blocksize = alsg.gsize[block];

     // FIXME: test this feature
     if ( block>nblocks-8 && blocksize>1 && 1 ) {
          int r =   alsg.gstart[block]-1;

          bool check = filter.filterJpartial ( candidate, r );
          if ( verbose>=2 ) {
               int j = partial_inner_product ( candidate, filter.als, filter.als.n_columns-1, r );
               printfd ( "partial check: block %d/%d:  r %d, j %d, check %d\n", block,nblocks, r, j, check );
          }
          if ( !check )
               return;
     }
     //printfd("sg.gstart.size() %d, block %d: blocksize %d\n", sg.gstart.size(), block, sg.gsize[block]);
     if ( verbose>=2 )
          printfd ( "inflateCandidateExtensionHelper: block %d/%d: blocksize %d\n", block, alsg.gsize.size(), blocksize );

     if ( blocksize==1 ) {
          // easy case
          inflateCandidateExtensionHelper ( list, basecandidate, candidate, block+1, al, alsg, check_indices, ct, verbose, filter,ntotal );
          return;
     }
     int gstart = alsg.gstart[block];
     int gend = alsg.gstart[block+1];

     if ( block<=-1 ) {
          printfd ( "row: %d to %d\n", gstart, gend );
          cperm tmp ( candidate.begin() +gstart, candidate.begin() +gend );
          printf ( "   current perm: " );
          print_cperm ( tmp );
          printf ( "\n" );
     }
     /*
         if ( blocksize>8 && 0) {

             double tx0=get_time_ms();
             std::vector<cperm> cc = debug_branch0 ( candidate,  gstart, gend, block, blocksize,  check_indices, 0 );
             double tx1=get_time_ms();
             std::vector<cperm> ccd = debug_branch ( candidate,  gstart, gend, block, blocksize,  check_indices, 0 );
             double tx2=get_time_ms();


             printfd ( "## debug branching %ld -> %ld (%.3f %.3f)\n",(long) cc.size(), ( long ) ccd.size(), tx1-tx0, tx2-tx1 );

             printf ( " ---> orig list\n" );
             showCandidates ( cc );
             printf ( " ---> debug list\n" );
             showCandidates ( ccd );
             exit(0);
         }
     */
     // FIXME: enable the other branch!!!!!!!
     if ( verbose>=2 )
          printfd ( "  split\n" );
     if ( blocksize<3 || ( block >1 && 0 ) ) {
          unsigned long iter=0;
          std::sort ( candidate.begin() +gstart, candidate.begin() +gend );
          unsigned long nbc=0;
          do {
               const int showd=0;

               // NOTE: for larger block sizes do not use naive generation
               if ( block<=4 && blocksize>1 && showd ) {
                    cperm xx ( candidate.begin() +gstart, candidate.begin() +gend ) ;
                    printf ( "  block %d, blocksize %d, iter %ld (k? %d): perm ", block, blocksize, iter, al.n_columns );
                    print_cperm ( xx );
                    printf ( "\n" );
                    cperm xxc ( check_indices.begin() +gstart, check_indices.begin() +gend );
                    cperm tmp ( check_indices.begin() +gstart, check_indices.begin() +gend );

                    printf ( "    : check: " );
                    print_cperm ( xxc );
                    printf ( " ---> %d\n",  satisfy_symm ( candidate, check_indices, gstart, gend ) );
                    //if (blocksize>20)
                    //	exit(0);
               }
               //cout << s1 << endl;
               iter++;

               // TODO: smart symmetry generation
               if ( satisfy_symm ( candidate, check_indices, gstart, gend ) ) {
                    nbc++;
                    inflateCandidateExtensionHelper ( list, basecandidate, candidate, block+1, al, alsg, check_indices, ct, verbose, filter,ntotal );
               } else {

               }
               // TODO: run inline filter
          } while ( std::next_permutation ( candidate.begin() +gstart, candidate.begin() +gend ) );
          if ( verbose>=2 )
               printfd ( "nbc block %d: %d/%ld\n", block, nbc, iter );

          if ( blocksize>10 && 0 ) {
               printfd ( "block %d: nbc %ld\n", block, ( long ) nbc ) ;
               debug_candidate ( candidate, check_indices, "..." );
          }

     } else {
          const int showd = 0; // show debugging output
          std::sort ( candidate.begin() +gstart, candidate.begin() +gend );
          cperm candidatetmp ( candidate.begin(), candidate.end() );

          for ( int x=gstart; x<gend; x++ ) {
               candidatetmp[x]=-9;
          }
          const int N = gend-gstart;
          // special construction for larger blocksizes
          lightstack_t<branch_t> branches ( 3* ( gend-gstart+2 ) );

          if ( showd )
               printfd ( "inflation of large block: block %d, blocksize %d\n", block, blocksize );

          // count items;
          // push initial branches
          branch_t b1 = {-1, -1, {-1,-1,-1} };
          b1.nvals[0] = countvalues<signed char> ( candidate, gstart, gend, bvals[0] );
          b1.nvals[1] = countvalues<signed char> ( candidate, gstart, gend, bvals[1] );
          b1.nvals[2] = countvalues<signed char> ( candidate, gstart, gend, bvals[2] );

          branches.push ( b1 );

          double t0=get_time_ms();

          long nbc=0;
          do {
               branch_t b = branches.top();

               branches.pop();
               if ( b.row>=0 ) { // special check for dummy first branch element...
                    candidatetmp[b.row+gstart]=b.rval;
               }

               if ( b.row>=0 ) { // special check for dummy first branch element...
                    if ( check_symm_zero ( candidatetmp, check_indices, b.row+gstart ) ) {
                         // all good
                         //debug_candidate ( candidatetmp, check_indices, printfstring ( "good    based on symmetry: block %d, row %d, range %d %d", block, b.row+gstart, gstart, gend ).c_str() );
                    } else {
                         // discard branch

                         //debug_candidate ( candidatetmp, check_indices, printfstring ( "discard based on symmetry: block %d, row %d", block, b.row+gstart ).c_str() );
                         continue;
                    }
               }
               if ( b.row==N-1 ) {
                    nbc++;
                    // call the inflate function
                    inflateCandidateExtensionHelper ( list, basecandidate, candidatetmp, block+1, al, alsg, check_indices, ct, verbose, filter,ntotal );
                    continue;
               }

               // perform branching
               branch ( branches, b, bvals, 0, 3 );

          } while ( ! branches.empty() );
          if ( verbose>=2 )
               printfd ( "nbc block %d: %d/%ld\n", block, nbc, nbc );

          if ( showd ) {
               printfd ( "inflation of large block: blocksize %d: %d branch calls\n", blocksize, nbc );
          }
     }
}

std::vector<cperm> inflateCandidateExtension ( const cperm &basecandidate,  const array_link &als,  const symmetry_group &alsg, const std::vector<int> &check_indices, const conference_t & ct, int verbose , const DconferenceFilter &filter )
{
     long ntotal=0;

     cperm candidate = basecandidate;
     int block=0;
     std::vector<cperm> cc;
     inflateCandidateExtensionHelper ( cc, basecandidate, candidate, block, als, alsg, check_indices, ct, verbose, filter, ntotal );

     if ( verbose>=2 || 0 ) {
          printfd ( "inflateCandidateExtension: generated %ld/%ld candidates (k %d)\n", ( long ) cc.size(), ntotal, als.n_columns );
     }
     return cc;
}

std::vector<cperm> generateSingleConferenceExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose , int filtersymm, int filterj2, int filterj3, int filtersymminline )
{
     double t0=get_time_ms();
     if ( verbose ) {
          myprintf ( "%s: filters: symmetry %d, symmetry inline %d, j2 %d, j3 %d\n", __FUNCTION__, filtersymm, filtersymminline, filterj2, filterj3 );
     }
     assert ( al.n_columns>1 );

     const int N = al.n_rows;
     DconferenceFilter dfilter ( al, filtersymm, filterj2, filterj3 );
     if ( verbose>=2 ) {
          dfilter.show()      ;
     }

     std::vector<int> sidx = dfilter.sd.checkIdx();
     if ( verbose>=2 ) {
          print_perm ( "sidx: ", sidx );
     }
     std::vector<cperm> cc;
     cperm c ( al.n_rows );

     lightstack_t<branch_t> branches ( 3*N );

     c[0]=1;

     // push initial branches
     branch_t b1 = {0, 1, {1,N/2-1,N/2-1} };    // branches starting with a 1
     branches.push ( b1 );
     //branch_t b0 = {0, 0, {1,N/2-1,N/2-1} };
     //branches.push ( b0 );

#ifdef OADEBUG
     std::vector<long> nb ( N+1 );
#endif

     // TODO: inline kz filtering, TODO: partial J2 filter with second column

     long n=0;
     do {

          branch_t b = branches.top();

          if ( verbose>=3 ) {
               b.show();
          }
#ifdef OADEBUG
          nb[b.row]++;
#endif
          branches.pop();
          c[b.row]=b.rval;      // update column vector

          if ( b.row==dfilter.inline_row && filterj3 ) {
               // TODO: inline_row can be one earlier?
               if ( ! dfilter.filterJ3inline ( c ) )
                    // discard branch
                    continue;
          }
          if ( b.row==N-1 ) {
               n++;
               //verbose=2;
               //dfilter.filterReason(c);
               if ( verbose>=3 ) {
                    myprintf ( "n %d: filter %d: ", ( int ) n, dfilter.filter ( c ) );
                    print_cperm ( c );
                    printf ( "\n" );

               }

               if ( dfilter.filter ( c ) ) {
                    if ( kz>0 ) {
                         if ( ! checkZeroPosition ( c, kz ) ) {

                              continue    ;
                         }
                    }

                    if ( verbose>=2 ) {
                         printfd ( "## push candidate   : " );
                         print_cperm ( c );
                         myprintf ( "\n" );
                    }
                    cc.push_back ( c );
               } else {
                    // printfd ( "## discard candindate: " );
               }
               continue;
          }
          int istart=0;
          const int iend = 3;
          if ( sidx[b.row+1] && filtersymminline ) {
               if ( c[b.row]==-1 ) {
                    istart=2;

               }
               if ( c[b.row]==1 && 1 ) {
                    istart=1;
               }
          }
          for ( int i=istart; i<iend; i++ ) {

               if ( b.nvals[i]==0 )
                    continue;

               // checks

               branch_t bnew ( b );
               bnew.row++;
               bnew.rval = bvals[i];
               bnew.nvals[i]--;
               //printf("push new branch: i %d\n", i);
               // NOTE: make direct push possible
               // NOTE: can we eliminate the branches object altogether?
               branches.push ( bnew );
          }

     } while ( ! branches.empty() )  ;

#ifdef OADEBUG
     if ( 1 ) {
          printf ( "branch count:\n" );
          for ( int i=0; i<=N; i++ ) {
               printf ( "  %d: %ld\n", i, nb[i] );
          }
     }
#endif
     if ( verbose || 0 ) {
          printfd ( "%s: %.3f [s]: generated %ld/%ld/%ld perms (len %ld)\n", __FUNCTION__, get_time_ms() - t0, ( long ) cc.size(), n, factorial<long> ( c.size() ), ( long ) c.size() );
          //al.show();
          //al.transposed().showarray(); showCandidates ( cc );
     }
     return cc;
}


std::vector<cperm> generateDoubleConferenceExtensions ( const array_link &al, const conference_t & ct, int verbose , int filtersymm, int filterj2, int filterj3, int filtersymminline )
{
     if ( verbose )
          myprintf ( "generateDoubleConferenceExtensions: filters: symmetry %d, symmetry inline %d, j2 %d, j3 %d\n", filtersymm, filtersymminline, filterj2, filterj3 );

     assert ( ct.j1zero==1 );

     const int N = al.n_rows;
     DconferenceFilter dfilter ( al, filtersymm, filterj2 );

     std::vector<int> sidx = dfilter.sd.checkIdx();
     if ( verbose>=2 ) {
          print_perm ( "sidx: ", sidx );
     }
     std::vector<cperm> cc;
     cperm c ( al.n_rows );

     lightstack_t<branch_t> branches ( 3*N );

     c[0]=1;

     // push initial branches
     branch_t b1 = {0, 1, {2,N/2-2,N/2-1} };
     branches.push ( b1 );
     branch_t b0 = {0, 0, {1,N/2-1,N/2-1} };
     branches.push ( b0 );
     double t0=get_time_ms();

#ifdef OADEBUG
     std::vector<long> nb ( N+1 );
#endif

     // TODO: do faster inline checks (e.g. abort with partial symmetry, take combined J2 check with many zeros)
     long n=0;
     do {

          branch_t b = branches.top();
#ifdef OADEBUG
          nb[b.row]++;
#endif
          branches.pop(); // TODO: use reference and pop later
          c[b.row]=b.rval;

          if ( b.row==dfilter.inline_row && filterj3 ) {
               // TODO: inline_row can be one earlier?
               if ( ! dfilter.filterJ3inline ( c ) )
                    // discard branch
                    continue;
          }
          if ( b.row==N-1 ) {
               n++;
               //verbose=2;
               //dfilter.filterReason(c);
               if ( dfilter.filter ( c ) ) {
                    if ( verbose>=2 ) {
                         printfd ( "## push candindate   : " );
                         print_cperm ( c );
                         myprintf ( "\n" );
                    }
                    cc.push_back ( c );
               } else {
                    // printfd ( "## discard candindate: " );
               }
               continue;
          }
          int istart=0;
          const int iend = 3;
          if ( sidx[b.row+1] && filtersymminline ) {
               if ( c[b.row]==-1 ) {
                    istart=2;

               }
               if ( c[b.row]==1 && 1 ) {
                    istart=1;
               }
          }
          for ( int i=istart; i<iend; i++ ) {

               if ( b.nvals[i]==0 )
                    continue;

               // checks

               branch_t bnew ( b );
               bnew.row++;
               bnew.rval = bvals[i];
               bnew.nvals[i]--;
               //printf("push new branch: i %d\n", i);
               // NOTE: make direct push possible
               // NOTE: can we eliminate the branches object altogether?
               branches.push ( bnew );
          }

     } while ( ! branches.empty() )  ;

#ifdef OADEBUG
     if ( 1 ) {
          printf ( "branch count:\n" );
          for ( int i=0; i<=N; i++ ) {
               printf ( "  %d: %ld\n", i, nb[i] );
          }
     }
#endif
     if ( verbose || 1 ) {
          printfd ( "generateDoubleConferenceExtensions: generated %ld/%ld/%ld perms (len %ld)\n", ( long ) cc.size(), n, factorial<long> ( c.size() ), ( long ) c.size() );
          //al.show();
          //al.transposed().showarray(); showCandidates ( cc );
     }
     return cc;
}

/** generate double conference matrices
 *
 * Old vesion that still can handle the j1zero=0 case
 *
 **/
std::vector<cperm> generateDoubleConferenceExtensions2 ( const array_link &al, const conference_t & ct, int verbose , int filtersymm, int filterip )
{
     assert ( ct.itype==CONFERENCE_RESTRICTED_ISOMORPHISM  || ct.itype==CONFERENCE_ISOMORPHISM );

     int j1zero = ct.j1zero;

     std::vector<cperm> cc;

     const int N = ct.N;
     cperm c ( N );


     DconferenceFilter dfilter ( al, filtersymm, filterip );
     dfilter.filterfirst=1;
     dfilter.filterj3=ct.j3zero;
     unsigned long n=0;
     for ( int i=0; i<N-2; i++ ) {
          if ( j1zero && i!= ( N-2 ) /2 )
               continue;

          // fill initial permutation
          std::fill ( c.begin(), c.end(), -1 );
          c[0]=0;
          c[1]=0;
          for ( int k=2; k<i+2; k++ )
               c[k]=1;

          std::sort ( c.begin(), c.end() );


          do {
               //cout << s1 << endl;
               n++;

               if ( dfilter.filter ( c ) ) {
                    cc.push_back ( c );
               }
          } while ( std::next_permutation ( c.begin(), c.end() ) );
     }

     //printfd ( "generateDoubleConferenceExtensions: before filter generated %d/%ld perms (len %ld)\n", n, factorial<long> ( c.size() ), ( long ) c.size() );
     //cc= filterDconferenceCandidates ( cc, al, filtersymm,  filterip, verbose );
     if ( verbose || 0 ) {
          printfd ( "generateDoubleConferenceExtensions: generated %ld/%ld/%ld perms (len %ld)\n", ( long ) cc.size(), n, factorial<long> ( c.size() ), ( long ) c.size() );
          //al.show();
          //al.transposed().showarray(); showCandidates ( cc );
     }
     return cc;
}

std::vector<cperm> generateConferenceRestrictedExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose , int filtersymm, int filterip )
{

     const int extcol=al.n_columns;
//const int N = al.n_rows;
     const int N = ct.N;

     // special case
     if ( extcol==1 ) {
          std::vector<cperm> ee;

          // we have no values +1 in the first column and
          int no = nOnesFirstColumn ( al );
          int n1=no-1;
          int n2 = N-n1-2;

          // in the second column we start with [1,0]^T and then in block1: k1 values +1, block2: k2 values +1
          // we must have k2 = k1+(n2-n1)/2 = k1 +N -2n1-2

          cperm cc ( N );
          cc[0]=1;
          cc[1]=0;
          for ( int k1=0; k1<=n1; k1++ ) {
               int k2 = k1+ ( n2-n1 ) /2;
               if ( k2>n2 )
                    continue;
               if ( k2<0 )
                    continue;
               printf ( "generateConferenceRestrictedExtensions: column 1: n1 %d, n2 %d, k1 %d, k2 %d\n", n1, n2, k1, k2 );

               std::fill ( cc.begin() +2, cc.end(), -1 );

               for ( int i=2; i<2+k1; i++ )
                    cc[i]=1;
               for ( int i=2+n1; i<2+n1+k2; i++ )
                    cc[i]=1;

               ee.push_back ( cc );
          }
          return ee;
     }

     std::vector<int> moi; // indices of -1 in first column
     for ( int i=0; i<N; i++ ) {
          if ( al.atfast ( i,0 ) ==-1 )
               moi.push_back ( i );
     }
     array_link alx = al.clone();

     // multiply
     for ( size_t i =0; i<moi.size(); i++ ) {
          alx.negateRow ( moi[i] );
     }

     // sort rows of array
     indexsort is = rowsorter ( alx );

     // now get candidate columns for the normal case, afterwards convert then using the rowsorter and row negations

     // loop over all possible first combinations
     std::vector<cperm> ff = get_first ( N, kz, verbose );

     if ( verbose>=2 ) {
          for ( size_t i=0; i<ff.size(); i++ ) {
               printf ( "extend1 %d: N %d: ", ( int ) i, N );
               display_vector ( ff[i] );
               printf ( "\n" );
          }
     }

     conference_extend_t ce;
     std::vector<cperm> extensions ( 0 );

     ce.first=ff;
     ce.second=ff;

     array_link als = al.selectFirstColumns ( extcol );

     cperm c0 = getColumn ( al, 0 );
     cperm c1 = getColumn ( al, 1 );
     for ( size_t i=0; i<ce.first.size(); i++ ) {
          int ip = innerprod ( c0, ce.first[i] );

          // TODO: cache this function call
          int target = -ip;

          std::vector<cperm> ff2 = get_second ( N, kz, target, verbose>=2 );
          ce.second=ff2;

          //printfd("ce.second[0] "); display_vector( ce.second[0]); printf("\n");

          for ( size_t j=0; j<ff2.size(); j++ ) {
               cperm c = ce.combine ( i, j );

               extensions.push_back ( c );
               continue;

          }
     }

     ce.extensions=extensions;
     if ( verbose>=2 )
          printf ( "generateConferenceExtensions: after generation: found %d extensions\n", ( int ) extensions.size() );

     // perform row symmetry check
     std::vector<cperm> e2 = filterCandidates ( extensions, al,  filtersymm,  filterip,  verbose );

     if ( verbose>=1 )
          printf ( "extend_conference: symmetry check %d + ip filter %d: %d->%d\n", filtersymm, filterip, ( int ) extensions.size(), ( int ) e2.size() );

     ce.extensions=e2;

     return e2;
}

int selectZmax ( int maxzpos, const conference_t::conference_type &ctype, const array_link &al, int extcol )
{
     if ( maxzpos<0 ) {
          switch ( ctype ) {
          case conference_t::CONFERENCE_NORMAL:
               maxzpos = al.n_rows-1;
               break;
          case conference_t::CONFERENCE_DIAGONAL:

               maxzpos = extcol;
               //printf("ct.ctype==conference_t::conference_t::CONFERENCE_DIAGONAL: maxzpos %d/%d, extcol %d\n", maxzpos, al.n_rows-1, extcol);
               break;
          case conference_t::DCONFERENCE:
               maxzpos = al.n_rows-1;

               break;
          default
                    :
               printfd ( "not implemented...\n" );
               maxzpos = al.n_rows-1;
          }
     }
     return maxzpos;
}

conference_extend_t extend_double_conference_matrix ( const array_link &al, const conference_t & ct,  CandidateGeneratorDouble &cgenerator, int extcol, int verbose, int maxzpos )
{
     conference_extend_t ce;
     ce.extensions.resize ( 0 );

     const int N = ct.N;
     const int k = extcol;
     const int maxzval = maxz ( al );

     if ( verbose )
          printf ( "--- extend_double_conference_matrix: extcol %d, maxz %d, itype %d ---\n", extcol, maxzval, ct.itype );

     //ct.j1zero

     int filterip=1;
     int filtersymm=1;
     std::vector<cperm> cc;

     //cgenerator.generateDoubleConferenceExtensions(al, ct, verbose, filterip, 1);

     if ( k>=3 && filtersymm && filterip && ct.j1zero==1 && 1 ) {
          //cgenerator.last_valid=0;
          // FIXME: for large symmetry blocks start with k higher!
          cc = cgenerator.generateCandidates ( al );
     } else {
          if ( k>3 && ct.j1zero==1 ) {
               cc = generateDoubleConferenceExtensionsInflate ( al, ct, verbose, filterip, 1 );
          } else {
               if ( ct.j1zero==1 )
                    cc= generateDoubleConferenceExtensions ( al, ct, verbose, filtersymm, filterip );
               else
                    cc= generateDoubleConferenceExtensions2 ( al, ct, verbose, filtersymm, filterip );
          }
     }

     if ( ct.j3zero ) {
          //printfd("filter on j3 values\n");
          cc = filterJ3 ( cc, al, verbose );
     }

     ce.extensions = cc;

     return ce;
}

//conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t &ct, int extcol, int verbose=1, int maxzpos=-1 );
conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t &ct, int extcol, int verbose, int maxzpos )
{
     conference_extend_t ce;
     ce.extensions.resize ( 0 );

     const int N = ct.N;
     const int k = extcol;
     const int maxzval = maxz ( al );

     if ( verbose )
          printf ( "--- extend_conference_matrix: extcol %d, maxz %d, itype %d ---\n", extcol, maxzval, ct.itype );

     const int zstart=maxzval+1;

     maxzpos = selectZmax ( maxzpos, ct.ctype, al, extcol );

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

conference_extend_t extend_conference_matrix_generator ( const array_link &al, const conference_t & ct, int extcol, int verbose, int maxzpos, const CandidateGenerator &cgenerator )
{
     conference_extend_t ce;
     ce.extensions.resize ( 0 );

     const int N = ct.N;
     const int k = extcol;
     const int maxzval = maxz ( al );

     const int zstart=maxzval+1;

     maxzpos = selectZmax ( maxzpos, ct.ctype, al, extcol );
     if ( verbose ) {
          printfd ( "--- extend_conference_matrix: extcol %d, zstart %d, maxz %d, maxzpos %d ---\n", extcol, zstart, maxzval, maxzpos );

          if (verbose>=2)
               al.showarray();
     }
     for ( int ii=zstart; ii<maxzpos+1; ii++ ) {
          if ( verbose>=2 )
               printfd ( "array: generate with zero at %d\n", ii );
          std::vector<cperm> extensionsX;

          {
               //const std::vector<cperm> &cl = cgenerator.cande.ce[ii];
               //FIXME: only if valid
               const std::vector<cperm> &cl = cgenerator.generateCandidatesZero ( al, ii );
               if ( verbose>2 ) {
                    printfd ( "-- ii %d, %d candidates \n", ii, cl.size() );
                    if ( verbose>=2 ) {
                         cgenerator.showCandidates ( 2 );
                    }
               }
               //printfd ( " cande %d --> cache %d\n", ( int ) cgenerator.cande.ce[ii].size(), ( int ) cl.size() );
               if ( ct.ctype==conference_t::CONFERENCE_DIAGONAL ) {
                    extensionsX = filterCandidatesSymm ( cl,  al, verbose ); // FIXME: is this needed?
                    extensionsX = filterCandidates ( extensionsX,  al,1, 1, verbose );
               } else {
                    // FIXME: for cached generated candidates we could remove the filterj2 check
                    extensionsX  = filterCandidates ( cl,  al,1, 1, verbose );
               }

               //extensionsX  = filterCandidates ( cande,  al,1, 1, verbose );
          }

          if ( verbose>=2 ) {
               printf ( "array: kz %d: %d extensions\n", ii, ( int ) extensionsX.size() );
          }
          if ( 0 ) {
               printfd ( "### generated candidates at kz %d...\n", ii );
               showCandidates ( extensionsX );

          }
          ce.extensions.insert ( ce.extensions.end(), extensionsX.begin(), extensionsX.end() );
     }


     return ce;
}

// FIXME: implement a local_symmetry check
// FIXME: implement partial j2 checking for cached conference matrices

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



template<typename T>
/// return size of std::vector in memory (in bytes)
size_t vectorsizeof ( const typename std::vector<T>& vec )
{
     return sizeof ( T ) * vec.size();
}


conf_candidates_t generateCandidateExtensions ( const conference_t ctype, int verbose=1, int ncstart=3, int ncmax=-1, int root = -1 )
{

     conf_candidates_t cande;

     cande.ce.resize ( ctype.N );

     array_link al2 = ctype.create_root();
     array_link al3 = ctype.create_root_three();

     if ( ncmax==-1 ) {
          ncmax=ctype.N;
          if ( ctype.ctype==conference_t::CONFERENCE_DIAGONAL ) {
               ncmax=ncstart;
          }
     }

     for ( int extcol=ncstart-1; extcol<ncmax; extcol++ ) {
          std::vector<cperm> ee;
          {
               switch ( root ) {
               case -1: {
                    if ( extcol==2 )
                         ee = generateConferenceExtensions ( al2, ctype, extcol, 0, 0, 1 );
                    else
                         ee = generateConferenceExtensions ( al3, ctype, extcol, 0, 0, 1 );
               }
               break;
               case 2:
                    ee = generateConferenceExtensions ( al2, ctype, extcol, 0, 0, 1 );
                    break;
                    break;
               case 3:
                    ee = generateConferenceExtensions ( al3, ctype, extcol, 0, 0, 1 );
                    break;
               default
                         :
                    printfd ( "not implemented!\n" );
                    break;
               }


          }
          //printf("al3:\n"); al3.showarray();

          if ( ( long ) vectorsizeof ( ee ) > ( long ( 1 ) *1024*1024*1024 ) / ( long ) ctype.N ) {
               printfd ( "generateCandidateExtensions: set of generated candidates too large, aborting (root %d, extcol %d, ee.size() %d)", root, extcol, ee.size() );
               assert ( 0 );
               exit ( 0 );
          }
          cande.ce[extcol] = ee;
     }

     cande.info ( verbose );
     return cande;
}



arraylist_t extend_double_conference ( const arraylist_t &lst, const conference_t ctype, int verbose )
{
     // TODO: cache candidate extensions
     arraylist_t outlist;
     if ( verbose>=2 ) {
          printfd ( "extend_double_conference: start with %d arrays\n", ( int ) lst.size() );
     }
     double t0=get_time_ms();

     int vb=std::max ( 0, verbose-1 );

     int ncstart=3;
     if ( lst.size() >0 )
          ncstart=lst[0].n_columns+1;

     CandidateGeneratorDouble cgenerator ( array_link() , ctype );

     for ( size_t i=0; i<lst.size(); i++ ) {
          const array_link &al = lst[i];
          int extcol=al.n_columns;
          conference_extend_t ce = extend_double_conference_matrix ( al, ctype, cgenerator, extcol, vb, -1 );


          arraylist_t ll = ce.getarrays ( al );
          const int nn = ll.size();

          outlist.insert ( outlist.end(), ll.begin(), ll.end() );

          if ( verbose>=2 || ( verbose>=1 && ( i%400==0 || i==lst.size()-1 ) ) ) {
               double dt = get_time_ms() - t0;
               printf ( "extend_conference: extended array %d/%d to %d/%d arrays (%.1f [s])\n", ( int ) i, ( int ) lst.size(), nn, ( int ) outlist.size(), dt );
               fflush ( 0 );
          }
     }
     return outlist;
}

arraylist_t extend_conference_restricted ( const arraylist_t &lst, const conference_t ctype, int verbose )
{
     arraylist_t outlist;

     if ( verbose>=2 ) {
          printfd ( "extend_conference: start %d\n", ( int ) lst.size() );
     }

     int vb=std::max ( 0, verbose-1 );

     int ncstart=3;
     if ( lst.size() >0 )
          ncstart=lst[0].n_columns+1;


     for ( size_t i=0; i<lst.size(); i++ ) {
          const array_link &al = lst[i];
          int extcol=al.n_columns;
          conference_extend_t ce = extend_conference_matrix ( al, ctype, extcol, vb, -1 );


          arraylist_t ll = ce.getarrays ( al );
          const int nn = ll.size();

          outlist.insert ( outlist.end(), ll.begin(), ll.end() );

          if ( verbose>=2 || ( verbose>=1 && ( i%200==0 || i==lst.size()-1 ) ) ) {
               printf ( "extend_conference: extended array %d/%d to %d arrays\n", ( int ) i, ( int ) lst.size(), nn );
               fflush ( 0 );
          }
     }
     return outlist;
}

/// select the unique arrays from a list of arrays. the indices of the unique arrays are returned
std::vector<int> selectUniqueArrayIndices ( const arraylist_t &lstr, int verbose )
{
     // perform stable sort
     indexsort sortidx ( lstr );

     const std::vector<int> &idx = sortidx.indices;

     const size_t nn = lstr.size();

     std::vector<int> cidx ( nn );
     std::vector<int> ridx;

     array_link prev;

     if ( lstr.size() >0 )
          prev= lstr[0];
     prev.setconstant ( -10 );

     int ci=-1;
     for ( size_t i=0; i<idx.size(); i++ ) {
          array_link al=lstr[idx[i]];
          if ( al!=prev ) {
               // new isomorphism class
               ci++;
               if ( verbose>=3 )
                    printf ( "selectConferenceIsomorpismClasses: representative %d: index %d\n", ( int ) ci, ( int ) idx[i] );
               ridx.push_back ( idx[i] );
               prev=al;
          }
          cidx[i]=ci;
     }
     return ridx;
}


/// FIXME: name
/// FIXME: reduce to more compact format
array_link reduceMatrix ( const array_link &al, matrix_isomorphism_t itype, int verbose )
{
     array_link alx;
     switch ( itype ) {
     case CONFERENCE_ISOMORPHISM: {
          alx= reduceConference ( al, verbose>=2 );
     }
     break;
     case CONFERENCE_RESTRICTED_ISOMORPHISM: {
          arraydata_t arrayclass ( 3, al.n_rows, 1, al.n_columns );
          //printfd("run %d ", i); arrayclass.show();
          array_transformation_t t = reduceOAnauty ( al+1, verbose>=2, arrayclass );
          alx=t.apply ( al+1 ) + ( - 1 );
          break;
     }
     default
               :
          printfd ( "error: isomorphism type not implemented\n" );
          break;
     }
     return alx;
}



class ConferenceIsomorphismSelector {
// OPTIONS: special array_link with reduced size --> factor 4 for char type
// online adding of single element or batches with isomorphism reduction
public:
     matrix_isomorphism_t itype;
     int verbose;
     int select_isomorphism_classes; /// if true then select only a single representative for each isomorphism class

     arraylist_t candidates;
     arraylist_t reductions;

private:
     int nadd ;
public:
     ConferenceIsomorphismSelector ( matrix_isomorphism_t itype, int verbose, int select_isomorphism_classes ) {
          this->itype=itype;
          this->verbose=verbose;
          this->select_isomorphism_classes=select_isomorphism_classes;
          this->nadd=0;

     }

     size_t size() const {
          return candidates.size();
     }

     void add ( const arraylist_t &lst ) {

          long nstart = candidates.size();
          long nextra = lst.size();
          candidates.insert ( candidates.end(), lst.begin(), lst.end() );

          nadd++;
          if ( select_isomorphism_classes ) {
               for ( size_t i =0; i<lst.size(); i++ ) {
                    array_link alr = reduceMatrix ( lst[i], itype, verbose );
                    reductions.push_back ( alr );
               }

               if ( nadd%1000==0 ) {
                    // reduce...
                    std::vector<int> ridx = selectUniqueArrayIndices ( reductions, verbose );

                    //FIXME: remove entries
                    // see http://stackoverflow.com/questions/596162/can-you-remove-elements-from-a-stdlist-while-iterating-through-it

                    arraylist_t tmp;
                    arraylist_t tmpr;
                    for ( size_t i=0; i<ridx.size(); i++ ) {
                         size_t ix = ridx[i];
                         tmp.push_back ( candidates[ix] );
                         tmpr.push_back ( reductions[ix] );
                    }
                    candidates=tmp;
                    reductions=tmpr;
                    if ( verbose )
                         printf ( "  reduce ... %ld -> %ld \n", nstart, ( long ) candidates.size() );

               }
          }

          long nfinal = candidates.size();

          if ( verbose ) {
               printf ( "ConferenceIsomorphismSelector: add %ld -> %ld -> %ld\n", ( long ) nstart, nstart + nextra, nfinal );
          }
     }
};


arraylist_t extend_conference_plain ( const arraylist_t &lst, const conference_t ctype, int verbose, int select_isomorphism_classes )
{
     double t0=get_time_ms();
     arraylist_t outlist;

     if ( verbose>=2 ) {
          printfd ( "extend_conference: start %d\n", ( int ) lst.size() );
     }

     int vb=std::max ( 0, verbose-1 );

     //int ncstart=3;
     //if ( lst.size() >0 )
     //     ncstart=lst[0].n_columns+1;

     ConferenceIsomorphismSelector selector ( ctype.itype, verbose>=2, select_isomorphism_classes );

     for ( size_t i=0; i<lst.size(); i++ ) {
          const array_link &al = lst[i];
          int extcol=al.n_columns;
          conference_extend_t ce = extend_conference_matrix ( al, ctype, extcol, vb, -1 );


          arraylist_t ll = ce.getarrays ( al );
          const int nn = ll.size();

          selector.add ( ll );

          if ( verbose>=2 || ( verbose>=1 && ( i%1000==0 || i==lst.size()-1 ) ) ) {
               printf ( "extend_conference: extended array %d/%d to %d arrays (total %ld, %.1f [s])\n", ( int ) i, ( int ) lst.size(), nn, ( long ) selector.size(), get_time_ms()-t0 );
               fflush ( 0 );
          }
     }

     return selector.candidates;
}

arraylist_t extend_conference ( const arraylist_t &lst, const conference_t ctype, int verbose, int select_isomorphism_classes )
{
     double t0=get_time_ms();
     arraylist_t outlist;

     if ( verbose>=2 ) {
          printfd ( "extend_conference: start with %d arrays\n", ( int ) lst.size() );
     }

     int vb=std::max ( 0, verbose-1 );

     //int ncstart=3;
     //if ( lst.size() >0 )
     //     ncstart=lst[0].n_columns+1;

     /// FIXME: move higher up in hierarchy
     CandidateGenerator cgenerator ( array_link(), ctype );
     //printfd ( "CandidateGenerator: constructed\n" );
     //cgenerator.verbose=2;

     if ( 0 ) {
          // debugging
          cgenerator.generateCandidatesZero ( lst[0],5 );
          cgenerator.generateCandidatesZero ( lst[0],7 );
          cgenerator.showCandidates();
     }

     // legacy code, fixed set of candidates, remove later
     //cgenerator.cande  = generateCandidateExtensions ( ctype, verbose>=2, ncstart );

     ConferenceIsomorphismSelector selector ( ctype.itype, verbose>=2, select_isomorphism_classes );

     for ( size_t i=0; i<lst.size(); i++ ) {
          if ( verbose>=2 )
               printfd ( "extend_conference: extend array %d\n", i );

          const array_link &al = lst[i];
          int extcol=al.n_columns;
          conference_extend_t ce = extend_conference_matrix_generator ( al, ctype, extcol, vb, -1, cgenerator );

          arraylist_t ll = ce.getarrays ( al );
          const int nn = ll.size();

          selector.add ( ll );

          if ( verbose>=2 || ( verbose>=1 && ( i%400==0 || i==lst.size()-1 ) ) ) {
               printf ( "extend_conference: extended array %d/%d to %d arrays (total %ld, %.1f [s])\n", ( int ) i, ( int ) lst.size(), nn, ( long ) selector.size(), get_time_ms()-t0 );
               fflush ( 0 );
          }
     }

     return selector.candidates;
}

std::pair<arraylist_t, std::vector<int> > selectConferenceIsomorpismHelper ( const arraylist_t &lst, int verbose, matrix_isomorphism_t itype )
{
     const int nn = lst.size();

     arraylist_t lstr;
     double t0=get_time_ms();

     // safety check
     if ( lst.size() >0 ) {
          if ( lst[0].min() <-1 ) {
               printfd ( "error: arrays should have positive integer values\n" );
               arraylist_t lstgood;
               std::vector<int> cidx;
               return std::pair<arraylist_t, std::vector<int> > ( lstgood, cidx );

          }
     }
     for ( int i=0; i< ( int ) lst.size(); i++ ) {
          if ( verbose>=1 && ( i%20000==0 || i== ( int ) lst.size()-1 ) )
               printf ( "selectConferenceIsomorpismClasses: reduce %d/%d\n", i, ( int ) lst.size() );
          //array_link alx;

          array_link alx = reduceMatrix ( lst[i], itype, 2* ( verbose>=3 ) );


          lstr.push_back ( alx );
     }

     // perform stable sort
     arraylist_t lstgood;
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
               if ( verbose>=3 )
                    printf ( "selectConferenceIsomorpismClasses: representative %d: index %d\n", ( int ) lstgood.size(), ( int ) idx[i] );

               lstgood.push_back (	lst[idx[i]] );
               prev=al;
               ci++;
          }
          cidx[i]=ci;
     }

     if ( verbose ) {
          double dt=get_time_ms()-t0;
          myprintf ( "selectConferenceIsomorpismClasses: select classes %d->%d (%.1f kArrays/h)\n", ( int ) lst.size(), ( int ) lstgood.size(), 3600*1e-3*double ( lst.size() ) /dt );
     }
     return std::pair<arraylist_t, std::vector<int> > ( lstgood, cidx );
}



std::vector<int> selectConferenceIsomorpismIndices ( const arraylist_t &lst, int verbose,  matrix_isomorphism_t itype )
{
     std::pair<arraylist_t, std::vector<int> > pp = selectConferenceIsomorpismHelper ( lst, verbose, itype ) ;
     return pp.second;
}

arraylist_t selectConferenceIsomorpismClasses ( const arraylist_t &lst, int verbose, matrix_isomorphism_t itype )
{
     std::pair<arraylist_t, std::vector<int> > pp = selectConferenceIsomorpismHelper ( lst, verbose , itype ) ;
     return pp.first;
}

arraylist_t  selectLMC0 ( const arraylist_t &list, int verbose,  const conference_t &ctype )
{

     arraylist_t out ;
     for ( size_t i=0; i<list.size(); i++ ) {
          lmc_t r=LMC0check ( list[i] );


          if ( verbose>=2 || ( verbose && i%10000==0 ) )
               printfd ( "selectLMC0: i %d/%d, r %d, total %d\n", i, list.size(), r, out.size() );
          if ( r==LMC_LESS ) {
               // pass, array is not in LMC0 format
          } else  {
               if ( verbose>=2 )
                    list[i].showarray();
               out.push_back ( list[i] );
          }
     }
     return out;

}


/// return true of alL is smaller than alR in LMC-0 ordering
bool compareLMC0 ( const array_link &alL, const array_link &alR )
{
     assert ( alL.n_rows==alR.n_rows );
     assert ( alL.n_columns==alR.n_columns );

     for ( int c=0; c<alL.n_columns; c++ ) {
          const array_t *al = alL.array + c*alL.n_rows;
          const array_t *ar = alR.array + c*alR.n_rows;

          // check position of zero(s) in column c
          for ( int r=0; r<alL.n_rows; r++ ) {
               if ( al[r]==0 &&  ar[r]!=0 )
                    return true;
               if ( al[r]!=0 && ar[r]==0 )
                    return false;
          }

          int zl = maxz ( alL, c );
          int zr = maxz ( alR, c );

          if ( zl<zr )
               return true;
          if ( zl>zr )
               return false;

          // zero is at same position(s) in column, let LMC ordering decide
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

bool compareLMC0_1 ( const array_link &alL, const array_link &alR )
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
     std::sort ( outlist.begin(), outlist.end(), compareLMC0 );
     return outlist;
}

conference_options::conference_options ( int maxpos )
{
     maxzpos=-1;
}

/// inflate a list of extensions
std::vector<cperm> conferenceReduce ( const std::vector<cperm> &ccX, const array_link &als, const array_link &alfull, const DconferenceFilter &filter, const conference_t &ct, int verbose )
{
     // FIXME: j2 check only for last pair
     //std::vector<cperm> cci = filter.filterList ( ccX );
     std::vector<cperm> cci = filter.filterListJ2last ( ccX );

     //printfd ( "inflate: %d: %d -> %d (als %d, alx %d, filter size %d)\n", als.n_columns, ( int ) ccX.size(), ( int ) cci.size() , als.n_columns, alfull.n_columns, filter.als.n_columns );
     return cci;
}


/// inflate a list of extensions
std::vector<cperm> extensionInflate ( const std::vector<cperm> &ccX, const array_link &als, const array_link &alfull, const DconferenceFilter &filter, const conference_t &ct, int verbose )
{
     std::vector<cperm> cci;
     std::vector<cperm> cc;

     symmetry_group alfullsg = alfull.row_symmetry_group();
     const std::vector<int> check_indices = alfullsg.checkIndices();
     symmetry_group alsg = als.row_symmetry_group();

     // loop over all candidinates with k columns and inflate to (k+1)-column candidates
     for ( size_t i=0; i<ccX.size(); i++ ) {
          const cperm &basecandidate = ccX[i];

          if ( verbose>2 )
               myprintf ( "### inflate candidate:" );
          //printf(" "); print_cperm( basecandidate); printf("\n");
          cc=  inflateCandidateExtension ( basecandidate, als, alsg, check_indices, ct, verbose, filter );

          if ( verbose>=2 ) {
               myprintf ( "inflate: array %d/%d: generated %ld candidates\n", ( int ) i, ( int ) ccX.size(), ( long ) cc.size() );
          }
          cci.insert ( cci.begin(), cc.begin(), cc.end() );
     }
     return cci;
}

std::vector<cperm> generateDoubleConferenceExtensionsInflate ( const array_link &al, const conference_t &ct, int verbose, int filterj2, int filterj3, int kstart )
{
     //kstart=1;
     double t00=get_time_ms();
     int kfinal=al.n_columns;
     if ( kstart<0 )
          kstart = al.n_columns-1;
     assert ( kstart>=1 );

     std::vector<cperm> cci;

     array_link als = al.selectFirstColumns ( kstart );

     double t0=get_time_ms();
     std::vector<cperm> ccX = generateDoubleConferenceExtensions ( als, ct, verbose, 1, filterj2, filterj3 );
     if ( verbose )
          printf ( "generateDoubleConferenceExtensionsInflate: extend array with %d columns (kfinal %d, initial array %d columns)\n", al.n_columns, kfinal, als.n_columns );
     if ( verbose )
          printf ( "   initial generation: dt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );

     for ( int kx=kstart; kx<kfinal; kx++ ) {
          als = al.selectFirstColumns ( kx );
          array_link alx = al.selectFirstColumns ( kx+1 );
          DconferenceFilter filter ( alx, 1, filterj2, filterj3 );

          if ( verbose )
               printf ( "## generateDoubleConferenceExtensionsInflate: at %d columns: start with %d extensions\n", kx+1, ( int ) ccX.size() );

          cci = extensionInflate ( ccX, als, alx, filter, ct, verbose );

          if ( verbose ) {
               printf ( "## generateDoubleConferenceExtensionsInflate: at %d columns: total inflated: %ld candidates for column %d\n", kx+1, cci.size(), kx+1 );
               printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t00 ) );
          }

          ccX=cci;
     }

     return cci;
}

CandidateGeneratorX::CandidateGeneratorX ( const array_link &al, const conference_t &ct_ ) : ct ( ct_ )   // , filter ( DconferenceFilter ( al, 1, 1, 1 ) )
{
     this->verbose=1;
     this->candidate_list_conf.resize ( ct.N+1 ); // set a safe max
     this->last_valid_conf.resize ( ct.N+1 );
     this->alz.resize ( ct.N+1 );

     for ( size_t i=0; i< ( size_t ) ct.N; i++ ) {
          this->candidate_list_conf[i].resize ( ct.N );
          this->alz[i] = al;
          this->last_valid_conf[i]=0;
     }
}

CandidateGeneratorBase::CandidateGeneratorBase ( const array_link &al, const conference_t &ct_ )
{
     this->ct = conference_t ( ct_ );
     this->al = al;
     this->verbose=1;
     this->last_valid=0;
     this->candidate_list.clear();
     this->candidate_list.resize ( ct.N+1 ); // set a safe max
}

CandidateGeneratorConference::CandidateGeneratorConference ( const array_link &al, const conference_t &ct_, int zero_position ) : CandidateGeneratorBase ( al, ct_ )
{
     this->zero_position = zero_position;
}

CandidateGeneratorZero::CandidateGeneratorZero ( const array_link &al, const conference_t &ct_, int zero_position ) : CandidateGeneratorBase ( al, ct_ )
{
     this->zero_position = zero_position;
}

CandidateGeneratorInflate::CandidateGeneratorInflate ( const array_link &al, const conference_t &ct_ ) : ct ( ct_ )
{
     for ( int i=0; i<ct.N; i++ ) {
          this->generators.push_back ( CandidateGeneratorZero ( al, ct, i ) );
     }
}

CandidateGeneratorDouble::CandidateGeneratorDouble ( const array_link &al, const conference_t &ct_ ) : CandidateGeneratorBase ( al, ct_ )   // , filter ( DconferenceFilter ( al, 1, 1, 1 ) )
{
}


const std::vector<cperm> & CandidateGeneratorX::generateConfCandidatesZero ( const array_link &al, int kz ) const
{
     std::vector<cperm> cx;

     assert ( ct.itype==CONFERENCE_ISOMORPHISM );
     assert ( ct.ctype==conference_t::CONFERENCE_DIAGONAL || ct.ctype==conference_t::CONFERENCE_NORMAL );
     assert ( ct.j1zero==0 );
     assert ( ct.j3zero==0 );

     const char *tag = "generateConfCandidates (cache)";
     const int filterj2=1;
     const int filterj3=0;
     double t00=get_time_ms();

     int startcol = this->startColumn ( al, kz );

     int kfinal=al.n_columns;
     int ncfinal = al.n_columns+1;
     int finalcol=kfinal;

     if ( verbose>=2 ) {
          printf ( "\n" );
          printf ( "## %s: kz %d, startcol %d, ncfinal %d, kfinal %d\n", tag, kz, startcol, ncfinal, kfinal );
     }

     std::vector<cperm> ccX, cci;
     int kstart=-1;

     if ( startcol==-1 ) {
          int targetcol =al.n_columns+1;
          int averbose=1;

          int root=3;
          if ( kfinal==2 )
               root=2;
          startcol=root;

          double t0=get_time_ms();
          const int ncmax=kz+1;
          const int ncstart=kz+1;
          conf_candidates_t tmp = generateCandidateExtensions ( ct, 0, ncstart, ncmax, root );
          ccX = tmp.ce[kz];

          //::showCandidates(tmp.ce[kz]);

          this->candidate_list_conf[kz][startcol] = ccX;
          this->last_valid_conf[kz]= startcol;
          kstart=startcol; // ??

          if ( verbose>=2 ) {
               printf ( "   ## %s: initial generate: %d candidates (kstart %d, kfinal %d)\n", tag, ( int ) ccX.size(), kstart, kfinal );
               this->startColumn ( al, kz , 1 );
          }

     } else {
          // FIXME: check this bound is sharp
          ccX = this->candidate_list_conf[kz][startcol];
          this->last_valid_conf[kz]=startcol;
          kstart=startcol; // ??

          if ( verbose>=2 || 0 ) {
               printf ( "   ## %s: initial cache: %d (startcol %d, kstart %d, kfinal %d)\n", tag, ( int ) ccX.size(), startcol, kstart, kfinal );

               conf_candidates_t tmp = generateCandidateExtensions ( ct, 0, kz+1, kz+1 );

               this->startColumn ( al, kz , 1 );

               if ( 0 ) {
                    printf ( "from cache:\n" );
                    ::showCandidates ( ccX );
                    printf ( "direct:\n" );
                    ::showCandidates ( ccX );
                    //exit(0);
               }
          }
     }

     for ( int kx=kstart; kx<kfinal; kx++ ) {
          // generates candidates for column index kx+1?

          array_link alx = al.selectFirstColumns ( kx+1 );
          DconferenceFilter filter ( alx, 0, filterj2, filterj3 );

          if ( verbose >=2 || 0 ) {
               printf ( "## %s: at %d columns: start with %d extensions, to generate extensions for column %d (%d column array)\n", tag, kx+1, ( int ) ccX.size(), kx+1 , kx+2 );

               printfd ( " --> inflate generates candindates for column %d\n", kx+1 );
          }
          size_t nprev=ccX.size();

          ccX = conferenceReduce ( ccX, alx, alx, filter, ct, verbose>=2 );

          if ( verbose >=2 ) {
               printf ( "## %s: at %d columns: total inflated: %ld->%ld\n",tag,  kx+1, ccX.size(), ( long ) nprev );
               printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t00 ) );
          }

          this->candidate_list_conf[kz][kx+1] = ccX;
          this->last_valid_conf[kz]=kx+1;
     }

     if ( verbose>=2 )
          printf ( "## %s: generated %d candidates for column index %d (with %d columns)\n",	tag, ( int ) this->candidate_list_conf[kz][ncfinal-1].size(), ncfinal-1, ncfinal );

     this->alz[kz] = al;

     // NOTE: consistency check
     if ( 0 ) {
          conf_candidates_t tmp = generateCandidateExtensions ( ct, 0, kz+1, kz+1 );
          cperm_list cctmp=tmp.ce[kz];
          DconferenceFilter filter ( al, 0, 1, 0 );

          cctmp=filter.filterList ( cctmp );
          if ( cctmp.size() != this->candidate_list_conf[kz][ncfinal].size()  || 1 ) {
               printf ( "DIFF: direct %d -> cached %d\n", ( int )	cctmp.size(), ( int ) this->candidate_list_conf[kz][ncfinal].size() );
          }
     }

     return 	this->candidate_list_conf[kz][kfinal];
}

const std::vector<cperm> & CandidateGeneratorConference::generateCandidates ( const array_link &al ) const
{
     // assert we have the right settings
     const char *tag = "generateCandidates (conference, zero fixed, cache)";
     const int filterj2=1;
     assert ( ct.j1zero==0 );
     const int filterj3=ct.j3zero;
     const int filtersymminline=1;
     double t00=get_time_ms();

     if ( verbose>=2 )
          myprintf ( "CandidateGenerator::%s: start\n",tag );

     int startcol = this->startColumn ( al, 0 );
     int kfinal=al.n_columns;
     int ncfinal = al.n_columns+1;
     int finalcol=kfinal;

     if ( verbose>=2 ) {
          myprintf ( "\n" );
          myprintf ( "## %s: startcol %d, ncfinal %d\n", tag, startcol, ncfinal );
     }

     std::vector<cperm> ccX, cci;
     int kstart=-1;

     /* select initial set of columns */
     if ( startcol==-1 ) {
          array_link als = al.selectFirstColumns ( START_COL );
          startcol=START_COL+1;
          int averbose=0;

          double t0=get_time_ms();
          ccX = generateSingleConferenceExtensions ( als, ct, -1, averbose, 1, filterj2, filterj3, filtersymminline );
          this->candidate_list[START_COL+1] = ccX;
          last_valid= START_COL+1;
          kstart=startcol-1;
     } else {
          // TODO: check this bound is sharp
          ccX = this->candidate_list[startcol];
          last_valid=startcol;
          kstart=startcol-1;
     }

     /* inflate the columns untill we have reached the target */
     array_link als;
     for ( int kx=kstart; kx<kfinal; kx++ ) {
          als = al.selectFirstColumns ( kx );
          array_link alx = al.selectFirstColumns ( kx+1 );
          DconferenceFilter filter ( alx, 1, filterj2, filterj3 );

          if ( verbose >=2 )
               myprintf ( "## %s: at %d columns: start with %d extensions, to generate extensions for column %d (%d column array)\n", tag, kx+1, ( int ) ccX.size(), kx+1 , kx+2 );

          cci = extensionInflate ( ccX, als, alx, filter, ct, ( verbose>=2 ) * ( verbose-1 ) );
          // FIXME: use the following
          // cci = filter.filterListZero ( cci );

          if ( verbose >=2 ) {
               myprintf ( "## %s: at %d columns: total inflated: %ld\n",tag,  kx+1, cci.size() );
               myprintf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t00 ) );
          }

          ccX=cci;

          this->candidate_list[kx+2] = ccX;
          this->last_valid=kx+2;
     }

     if ( verbose>=2 )
          printf ( "CandidateGenerator::%s: generated %d candidates with %d columns\n",tag, ( int ) this->candidate_list[ncfinal].size(), ncfinal );

     this->al = al;
     return this->candidate_list[ncfinal];
}


const std::vector<cperm> & CandidateGeneratorZero::generateCandidates ( const array_link &al ) const
{
     // assert we have the right settings
     const char *tag = "generateCandidates (conference, zero fixed, cache)";
     const int filterj2=1;
     assert ( ct.j1zero==0 );
     const int filterj3=ct.j3zero;
     double t00=get_time_ms();

     if ( verbose>=2 )
          printf ( "CandidateGenerator::%s: start\n",tag );


     int startcol = this->startColumn ( al, 1 );
     int kfinal=al.n_columns;
     int ncfinal = al.n_columns+1;
     int finalcol=kfinal;

     if ( verbose>=2 ) {
          printf ( "\n" );
          printf ( "## %s: startcol %d, ncfinal %d\n", tag, startcol, ncfinal );
     }

     std::vector<cperm> ccX, cci;
     int kstart=-1;

     /* select initial set of columns */
     if ( startcol==-1 ) {
          array_link als = al.selectFirstColumns ( START_COL );
          startcol=START_COL+1;
          int averbose=1;

          double t0=get_time_ms();
          ccX = generateConferenceExtensions ( al, ct, this->zero_position, averbose, 1, filterj2 );
          this->candidate_list[START_COL+1] = ccX;
          last_valid= START_COL+1;
          kstart=startcol-1;
     } else {
          // TODO: check this bound is sharp
          ccX = this->candidate_list[startcol];
          last_valid=startcol;
          kstart=startcol-1;
     }

     /* inflate the columns untill we have reached the target */
     array_link als;
     for ( int kx=kstart; kx<kfinal; kx++ ) {
          als = al.selectFirstColumns ( kx );
          array_link alx = al.selectFirstColumns ( kx+1 );
          DconferenceFilter filter ( alx, 1, filterj2, filterj3 );

          if ( verbose >=2 )
               printf ( "## %s: at %d columns: start with %d extensions, to generate extensions for column %d (%d column array)\n", tag, kx+1, ( int ) ccX.size(), kx+1 , kx+2 );

          cci = extensionInflate ( ccX, als, alx, filter, ct, ( verbose>=2 ) * ( verbose-1 ) );
          cci = filterZeroPosition ( cci, this->zero_position );

          if ( verbose >=2 ) {
               printf ( "## %s: at %d columns: total inflated: %ld\n",tag,  kx+1, cci.size() );
               printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t00 ) );
          }

          ccX=cci;

          this->candidate_list[kx+2] = ccX;
          this->last_valid=kx+2;
     }

     if ( verbose>=2 )
          printf ( "CandidateGenerator::%s: generated %d candidates with %d columns\n",tag, ( int ) this->candidate_list[ncfinal].size(), ncfinal );

     this->al = al;
     return this->candidate_list[ncfinal];
}

const std::vector<cperm> & CandidateGeneratorDouble::generateCandidates ( const array_link &al ) const
{
     // assert we have the right settings

     const char *tag = "generateCandidates (double conf matrices, cache)";
     const int filterj2=1;
     assert ( ct.j1zero==1 );
     const int filterj3=ct.j3zero;
     double t00=get_time_ms();

     int startcol = this->startColumn ( al, 0 );
     int kfinal=al.n_columns;
     int ncfinal = al.n_columns+1;
     int finalcol=kfinal;

     if ( verbose>=2 ) {
          printf ( "\n" );
          printf ( "## %s: startcol %d, ncfinal %d\n", tag, startcol, ncfinal );
     }

     std::vector<cperm> ccX, cci;
     int kstart=-1;

     /* select initial set of columns */
     if ( startcol==-1 ) {
          array_link als = al.selectFirstColumns ( START_COL );
          startcol=START_COL+1;
          int averbose=1;

          double t0=get_time_ms();
          ccX = generateDoubleConferenceExtensions ( als, ct, averbose, 1, filterj2, filterj3 );
          this->candidate_list[START_COL+1] = ccX;
          last_valid= START_COL+1;
          kstart=startcol-1;
     } else {
          // NOTE: check this bound is sharp?
          ccX = this->candidate_list[startcol];
          last_valid=startcol;
          kstart=startcol-1;
     }

     /* inflate the columns untill we have reached the target */
     array_link als;
     for ( int kx=kstart; kx<kfinal; kx++ ) {
          als = al.selectFirstColumns ( kx );
          array_link alx = al.selectFirstColumns ( kx+1 );
          DconferenceFilter filter ( alx, 1, filterj2, filterj3 );

          if ( verbose >=2 )
               printf ( "## %s: at %d columns: start with %d extensions, to generate extensions for column %d (%d column array)\n", tag, kx+1, ( int ) ccX.size(), kx+1 , kx+2 );

          cci = extensionInflate ( ccX, als, alx, filter, ct, ( verbose>=2 ) * ( verbose-1 ) );
          cci = filter.filterListZero ( cci );

          if ( verbose >=2 ) {
               printf ( "## %s: at %d columns: total inflated: %ld\n",tag,  kx+1, cci.size() );
               printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t00 ) );
          }
          ccX=cci;

          this->candidate_list[kx+2] = ccX;
          this->last_valid=kx+2;
     }

     if ( verbose>=2 )
          printf ( "CandidateGenerator::%s: generated %d candidates with %d columns\n",tag,	( int ) this->candidate_list[ncfinal].size(), ncfinal );


     this->al = al;
     return this->candidate_list[ncfinal];
}

void rowlevel_permutation ( const array_link &al, rowsort_t *rowperm, const std::vector<int> &colperm, std::vector<int> &rowsignperm, const rowindex_t n_rows, int column )
{

     for ( rowindex_t r=0; r < n_rows; r++ ) {
          if ( al.atfast ( r, colperm[column] ) < 0 )
               rowsignperm[ r ] = -1;
     }

}

void init_lmc0_rowsort ( const array_link &al, int sutk_col, rowsort_t *rowperm, std::vector<int> &rowsignperm, const rowindex_t n_rows, const colindex_t n_cols )
{

     for ( int i=0; i < n_rows; i++ ) {

          int rx = rowperm[i].r;
          int trans_val = ( rowsignperm[ rx ] ) *al.at ( rx, sutk_col );
          rowperm[ i ].val = ( ( trans_val+3 ) % 3 );

     }
     std::stable_sort ( rowperm, rowperm+al.n_rows );
}

int get_zero_pos_block ( const array_link &al, const int x1, const int x2, rowsort_t *rowperm, const std::vector<int> &colperm, int column, const int nrows )
{
     int position_zero = nrows+1;
     for ( int i = x1; i < x2; i++ ) {
          int current_val = al.atfast ( rowperm[i].r, colperm[column] );
          if ( current_val == 0 ) {
               position_zero = i;
          }
     }

     return position_zero;

}

lmc_t lmc0_compare_zeropos_block ( const array_link &al, const int x1, const int x2, rowsort_t *rowperm, const std::vector<int> &colperm, int column, const std::vector<int> &rowsignperm, const std::vector<int> &colsignperm, const int nrows )
{
     /* Get zero position in the original array*/
     int al_position_zero = nrows+1;

     for ( int i = x1; i < x2; i++ ) {
          if ( al.atfast ( i, column ) == 0 ) {
               al_position_zero = i; // changed from rowperm[r].r
          }
     }

     /* Check position of zeros */
     int position_zero = get_zero_pos_block ( al, x1, x2, rowperm, colperm, column, nrows );

     if ( position_zero > al_position_zero ) {
          return LMC_MORE;
     }
     if ( position_zero < al_position_zero ) {
          return LMC_LESS;
     }

     return LMC_NONSENSE;

}


/* Compare two columns with the zero element in the same position */
lmc_t compare_conf_columns ( const array_link &al, rowsort_t *rowperm, const std::vector<int> &colperm, int column, const std::vector<int> &rowsignperm, const std::vector<int> &colsignperm, const int nrows )
{

     for ( int i=0; i<nrows; i++ ) {
          int cp = colperm[column];
          int value_cdesign_trans = ( colsignperm[cp]*rowsignperm[rowperm[i].r] ) * al.atfast ( rowperm[i].r, cp );
          if ( ( ( al.atfast ( i, column ) +3 ) % 3 ) < ( ( value_cdesign_trans+3 ) % 3 ) ) // Transform the elements from (0, 1, -1) to (0, 1, 2)
               return LMC_MORE;
          if ( ( ( al.atfast ( i, column ) +3 ) % 3 ) > ( ( value_cdesign_trans+3 ) % 3 ) )
               return LMC_LESS;
     }
     return LMC_EQUAL;
}


lmc_t init_lmc0_sort_comp ( const array_link &al, int column, int sel_col, rowsort_t *rowperm, std::vector<int> &rowsignperm, const std::vector<int> &colperm, const std::vector<int> &colsignperm, const int n_rows )
{
     lmc_t r = LMC_NONSENSE;
     for ( int i=0; i < n_rows; i++ ) {

          int rx = rowperm[i].r; // play with current rowperm structure
          int posit_al = al.at ( rx, sel_col );
          int trans_val = ( rowsignperm[ rx ] ) *posit_al;
          int m = ( ( trans_val+3 ) % 3 );
          rowperm[ i ].val = m;

     }
     std::stable_sort ( rowperm, rowperm+n_rows ); // and now play with rowperm.rr

     // Compare zero position
     r = lmc0_compare_zeropos_block ( al, 0, n_rows, rowperm, colperm, column, rowsignperm, colsignperm, n_rows );
     if ( r==LMC_LESS || r==LMC_MORE ) {
          return r;
     }
     /* Compare whole Columns */
     r = compare_conf_columns ( al, rowperm, colperm, column, rowsignperm, colsignperm, n_rows );
     return r;

}

lmc_t LMC0_sortrows_compare ( const array_link &al, int column, rowsort_t *rowperm, const std::vector<int> &colperm, const std::vector<int> &rowsignperm, const std::vector<int> &colsignperm, const rowindex_t n_rows, const symmdata &sd )
{

     lmc_t r = LMC_NONSENSE;
     const int cp = colperm[column];
     for ( int i=0; i < n_rows; i++ ) {
          int rx = rowperm[i].r;
          int current_val = ( ( colsignperm[cp]*rowsignperm[ rx ] ) *al.atfast ( rx, cp ) );
          rowperm[ i ].val = ( ( current_val+3 ) % 3 );
     }

     /* Sort rows of the array in blocks*/
     int scol = column - 1;
     int nb = sd.ft.atfast ( sd.ft.n_rows - 1, scol ); // number of blocks
     /* we check in blocks determined by the ft */
     for ( int j = 0; j < nb; j++ ) {
          int x1 = sd.ft.atfast ( 2*j, scol );
          int x2 = sd.ft.atfast ( 2*j+1, scol );
          //std::stable_sort( rowperm+x1, rowperm+x2);
          shellSort ( rowperm, x1, x2-1 );

          // Compare blocks wrt zero position
          r = lmc0_compare_zeropos_block ( al, x1, x2, rowperm, colperm, column, rowsignperm, colsignperm, n_rows );
          if ( r==LMC_LESS || r==LMC_MORE ) {
               return r;
          }
     }

     /* Compare all elements */
     r = compare_conf_columns ( al, rowperm, colperm, column, rowsignperm, colsignperm, n_rows );
     return r;
}

lmc_t LMC0_columns ( const array_link &al, rowsort_t *rowperm, std::vector<int> colperm, int column, std::vector<int> &rowsignperm, std::vector<int> colsignperm, const int ncols, const int nrows, const symmdata &sd )
{

     lmc_t r = LMC_NONSENSE;

     for ( int c=column; c<ncols ; c++ ) {

          int col = colperm[c];
          colperm[c]=colperm[column];
          colperm[column]=col;

          /* i. Apply the correct column level permutation to make the element X(1,k) equal to 1*/
          int current_sign_col = colsignperm[colperm[column]];
          int current_val_firstrow = ( rowsignperm[rowperm[0].r]*current_sign_col ) * ( al.atfast ( rowperm[0].r, colperm[column] ) );
          colsignperm[ colperm[column] ] = colsignperm[ colperm[column] ] * current_val_firstrow;

          /* ii. Sort rows using the ordering 0, 1, -1 and compare */
          r = LMC0_sortrows_compare ( al, column, rowperm, colperm, rowsignperm, colsignperm, nrows, sd );

          if ( r==LMC_EQUAL ) {
               r = LMC0_columns ( al, rowperm, colperm, column+1, rowsignperm, colsignperm, ncols, nrows, sd );
          }
          if ( r==LMC_LESS ) {
               break; // array is not minimal form
          }

          colsignperm[ colperm[column] ] = current_sign_col;
          colperm[ column ] = colperm[ c ];
          colperm[ c ] = col;

     }
     return r;
}

lmc_t LMC0check ( const array_link &al, int verbose )
{
     /*0. Initialize data */
     lmc_t result = LMC_MORE;

     if ( ! al.is_conference() ) {
          printfd ( "error: input array is not a conference design" );
          return LMC_NONSENSE;
     }
     const int ncols = al.n_columns;
     const int nrows = al.n_rows;

     std::vector<int> colperm ( ncols );
     init_perm ( colperm );

     std::vector<int> colsignperm ( ncols );
     init_signperm ( colsignperm );

     std::vector<int> rowsignperm ( nrows );
     init_signperm ( rowsignperm );

     dyndata_t rowperm_data = dyndata_t ( nrows );
     rowsort_t *rowsort = rowperm_data.rowsort;

     for ( rowindex_t i = 0; i < nrows; i++ ) {
          rowsort[i].val = i;
     }

     symmdata sd ( al );

     for ( int sel_col = 0; sel_col < ncols; sel_col++ ) {

          /*1. Select the first (sel_col) column */
          int first_col = colperm[ 0 ];
          colperm[ 0 ] = colperm[ sel_col ];
          colperm[ sel_col ] = first_col;

          /*2. Find row-level permutation such that the first column only contains ones */
          rowlevel_permutation ( al, rowsort, colperm, rowsignperm, nrows, 0 );//

          /* 3. Find permutation to sort the array and compare column*/
          result = init_lmc0_sort_comp ( al, 0, sel_col, rowsort, rowsignperm, colperm, colsignperm, nrows );
          if ( result==LMC_LESS ) {
               return result;
          }

          //printf("--- sel_col %d\n",sel_col);
          //print_rowsort(rowsort, al.n_rows);

          /* 4. Select one of two possible sign permutations for the first row */
          for ( int r_sign = 0; r_sign < 2; r_sign++ ) {
               rowsignperm[ rowsort[0].r ] = 2*r_sign - 1;

               /* 5. Select the next column */
               result = LMC0_columns ( al, rowsort, colperm, 1, rowsignperm, colsignperm, ncols, nrows, sd );
               if ( result==LMC_LESS ) {
                    return result;

               }

          }

          colperm[ sel_col ] = colperm[ 0 ];
          colperm[ 0 ] = first_col;
          init_signperm ( rowsignperm );

     }

     return result;
}



// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 

