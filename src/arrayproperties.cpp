#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Core>

#ifdef _WIN32
#else
#include <stdbool.h>
#include <unistd.h>
#endif

#include "printfheader.h"	// with stdio
#include "mathtools.h"
#include "tools.h"
#include "arraytools.h"
#include "arrayproperties.h"


using namespace Eigen;


template<class Type>
/// helper class: n-dimensional array
class ndarray
{

public:
    Type *data;
    std::vector<int> dims; /// dimensions of the array
    int k; // dimension of the array
    int n; // total number of elements in the array
    std::vector<int> cumdims;
    std::vector<int> cumprod;

private:
    //std::vector<int> tmpidx;
public:
    ndarray ( std::vector<int> dimsx, int verbose=0 )
    {
        //verbose=2;
        k = dimsx.size();
        dims = dimsx;
        cumdims.resize ( k+1 );
        cumprod.resize ( k+1 );
        //tmpidx.resize ( k );
        n=1;
        for ( int i=0; i<k; i++ )
            n*= dims[i];
        cumdims[0]=0;
        for ( int i=0; i<k; i++ )
            cumdims[i+1] = cumdims[i]+dims[i];
        cumprod[0]=1;
        for ( int i=0; i<k; i++ )
            cumprod[i+1] = cumprod[i]*dims[i];

        if ( verbose ) {
            myprintf ( "ndarray: dimension %d, total %d\n", k, n );
            myprintf ( "  cumprod: " );
            printf_vector ( cumprod, "%d " );
            myprintf ( "\n" );
        }
        data = new Type[n];
        std::fill ( data, data+n, 0 );
    }


    std::string idxstring ( int i ) const
    {
        std::string s = "";
        std::vector<int> tmpidx ( this->k );
        linear2idx ( i, tmpidx );

        for ( int i=0; i<k; i++ ) {
            s+= printfstring ( "[%d]", tmpidx[i] ) ;
        }
        return s;
    }


    long totalsize() const
    {
        return n;
//	return this->cumprod[
    }
    void show() const
    {
        for ( int i=0; i<n; i++ ) {
            std::string idxstrx=idxstring ( i );
            myprintf ( "B[%d] = B%s = %f\n", i, idxstrx.c_str(), ( double ) data[i] );
        }
    }

    /// convert a linear index to normal indices
    inline void linear2idx ( int ndx, int *nidx = 0 ) const
    {

        if ( nidx==0 )
            return;

        for ( int i=k-1; i>=0; i-- ) {
            div_t xx = div ( ndx, cumprod[i] );
            int vi = xx.rem;
            int vj=xx.quot;
            nidx[i] = vj;
            //std::vector<int> *p = ( std::vector<int> * ) &tmpidx; p->at ( i ) = vj;
            ndx = vi;
        }
    }
    /// convert a linear index to normal indices
    inline void linear2idx ( int ndx, std::vector<int> & nidx ) const
    {

        assert ( (int) nidx.size() ==this->k );
        for ( int i=k-1; i>=0; i-- ) {
            div_t xx = div ( ndx, cumprod[i] );
            int vi = xx.rem;
            int vj=xx.quot;
            nidx[i] = vj;
            ndx = vi;
        }
    }

    inline int getlinearidx ( int *idx ) const
    {
        int lidx=0;
        for ( int i=0; i<k; i++ ) {
            //lidx = dims[i]*lidx+idx[i];
            lidx += idx[i]*cumprod[i];
        }
        return lidx;
    }

    /// set all values of the array to specified value
    void setconstant ( Type val )
    {
        std::fill ( this->data, this->data+this->n, val );
    }

    /// set value at position
    void set ( int *idx, Type val )
    {
        int lidx=getlinearidx ( idx );
        data[lidx]=val;
    }

    /// set value using linear index
    inline void setlinear ( int idx, Type val )
    {
        data[idx] = val;
    }
    /// get value using linear index
    inline void getlinear ( int idx, Type val ) const
    {
        return data[idx] ;
    }

    /// get value using n-dimensional index
    Type get ( int *idx ) const
    {
        int lidx=getlinearidx ( idx );
        return data[lidx];
    }

    ~ndarray()
    {
        delete [] data;
    }
};


/// Hamming distance
inline int dH ( const int N, int k, const array_link &al, int r1, int r2 )
{
    int dh=0;
    carray_t *D = al.array;
    carray_t *d1 = D+r1;
    carray_t *d2 = D+r2;

    for ( int c=0; c<k; c++ ) {
        dh += d1[c*N]!=d2[c*N];
    }
    return dh;
}

/// Hamming distance
inline void dHmixed ( const int N, int k, const array_link &al, int r1, int r2, int *dh, int ncolgroups, const std::vector<int> colgroupindex )
{
    for ( int i=0; i<ncolgroups; i++ )
        dh[i]=0;

    carray_t *D = al.array;
    carray_t *d1 = D+r1;
    carray_t *d2 = D+r2;

    for ( int c=0; c<k; c++ ) {
        int ci = colgroupindex[c];
        dh[ci] += d1[c*N]!=d2[c*N];
    }
}

/// Hamming distance (transposed array)
inline int dHx ( const int nr, int k, carray_t *data, int c1, int c2 )
{
    // nr is the OA column variable

    //printf("nr %d, k %d, c1 c2 %d %d\n", nr, k, c1, c2);
    int dh=0;
    carray_t *d1 = data+c1*nr;
    carray_t *d2 = data+c2*nr;

    for ( int c=0; c<nr; c++ ) {
        //printf("  c %d\n", c);
        dh += d1[c]!=d2[c];
    }
    return dh;
}

/// compare 2 GWPL sequences
int GWPcompare ( const std::vector<double> &a, const std::vector<double> &b )
{
    for ( size_t x=0; x<a.size(); x++ ) {
        if ( a[x]!=b[x] )
            return a[x]<b[x];
    }
    return 0;
}

/// calculate distance distrubution (array is transposed for speed)
std::vector<double> distance_distributionT ( const array_link &al, int norm=1 )
{
    int N = al.n_rows;
    int n = al.n_columns;

    // transpose array
    array_t *x = new array_t[N*n];
    array_t *xx = x;
    for ( int i=0; i<N; i++ )
        for ( int j=0; j<n; j++ ) {
            ( *xx ) = al.array[i+j*N];
            xx++;
        }

    // calculate distance distribution
    //printf("distance_distribution\n");
    std::vector<double> dd ( n+1 );

    for ( int r1=0; r1<N; r1++ ) {
        for ( int r2=0; r2<r1; r2++ ) {
            int dh = dHx ( n, N, x, r1, r2 );
            dd[dh]+=2; 	// factor 2: dH is symmetric
        }
    }
    // along diagonal
    dd[0] += N;

    if ( norm ) {
        for ( int x=0; x<=n; x++ ) {
            dd[x] /= N;
        }
    }

    delete [] x;
    return dd;
}

/// calculate distance distribution for mixed array
void distance_distribution_mixed ( const array_link &al, ndarray<double> &B, int verbose=1 )
{
    int N = al.n_rows;
    int n = al.n_columns;

    // calculate distance distribution

    //printf("distance_distribution\n");
    std::vector<double> dd ( n+1 );

    arraydata_t ad = arraylink2arraydata ( al );

    symmetry_group sg ( ad.getS(), false );


    int *dh = new int[sg.ngroups];
//std::fill(dh, dh+sg.ngroups, 0);
//printfd("dh: %d, %d\n", sg.ngroups, ad.ncolgroups);

    std::vector<int> dims ( ad.ncolgroups );
    for ( size_t i=0; i<dims.size(); i++ )
        dims[i] = ad.colgroupsize[i];
    if ( verbose>=3 ) {
        myprintf ( "distance_distribution_mixed before: \n" ); //display_vector(dd); std::cout << std::endl;
        B.show();
    }

    for ( int r1=0; r1<N; r1++ ) {
        for ( int r2=0; r2<r1; r2++ ) {
            dHmixed ( N, n, al, r1, r2, dh, sg.ngroups, sg.gidx );

            if ( verbose>=4 ) {
                myprintf ( "distance_distribution_mixed: rows %d %d: ", r1, r2 );
                print_perm ( dh, sg.ngroups );
            }
            int v = B.get ( dh );
            B.set ( dh, v+2 );

            if ( verbose>=3 ) {
                int w=B.getlinearidx ( dh );
                if ( w==0 ) {
                    myprintf ( " r1 %d, r2 %d\n", r1, r2 );
                }
            }
        }
    }
    if ( verbose>=3 ) {
        myprintf ( "distance_distribution_mixed low: \n" ); //display_vector(dd); std::cout << std::endl;
        B.show();
    }

    // along diagonal
    for ( unsigned int i=0; i<dims.size(); i++ )
        dh[i]=0;
    int v = B.get ( dh );
    B.set ( dh, v+N );

    if ( verbose>=3 ) {
        myprintf ( "distance_distribution_mixed integer: \n" );
        B.show();
    }

    for ( int x=0; x<B.n; x++ ) {
        // myprintf("xxxx %d/ %d\n", x, B.n);
        B.data[x] /= N;
    }

    if ( verbose ) {
        myprintf ( "distance_distribution_mixed: \n" );
        B.show();
    }

    delete [] dh;

    /*
    	if ( 0 ) {
    		myprintf ( "test: " );
    		int tmp[100];
    		int li0= 1;
    		B.linear2idx ( li0, tmp );
    		int li = B.getlinearidx ( tmp );
    		print_perm ( tmp, sg.ngroups );
    		myprintf ( " linear index: %d -> %d\n", li0, li );
    	}
    */
}

template <class Type>
/// calculate MacWilliams transform
std::vector<double> macwilliams_transform ( std::vector<Type> B, int N, int s )
{
    int n = B.size()-1;
    // printfd("macwilliams_transform: %d\n", n);
    std::vector<double> Bp ( n+1 );

    if ( s==2 ) {
        if ( n <= ncombscacheNumber() ) {
            // use cached version of krawtchouks
            //printfd("macwilliams_transform: using krawtchouksCache\n");
            for ( int j=0; j<=n; j++ ) {
                Bp[j]=0;
                for ( int i=0; i<=n; i++ ) {
                    Bp[j] +=  B[i] * krawtchouksCache<long> ( j, i, n ); //  calculate krawtchouk with dynamic programming
                }
                Bp[j] /= N;
            }
        } else {

            for ( int j=0; j<=n; j++ ) {
                Bp[j]=0;
                for ( int i=0; i<=n; i++ ) {
                    //myprintf("  B[i] %.1f krawtchouk(%d, %d, %d, %d) %ld \n", B[i], j,i,n,s, krawtchouk<long>(j, i, n, s));
                    Bp[j] +=  B[i] * krawtchouks<long> ( j, i, n ); // pow(s, -n)
                }
                Bp[j] /= N;
            }
        }

    } else {
        for ( int j=0; j<=n; j++ ) {
            Bp[j]=0;
            for ( int i=0; i<=n; i++ ) {
                Bp[j] +=  B[i] * krawtchouk<long> ( j, i, n, s ); // pow(s, -n)
            }
            Bp[j] /= N;
        }
    }

    return Bp;
}


std::vector<double> distance_distribution ( const array_link &al )
{
    int N = al.n_rows;
    int n = al.n_columns;

    // calculate distance distribution

    std::vector<double> dd ( n+1 );

    for ( int r1=0; r1<N; r1++ ) {
        for ( int r2=0; r2<r1; r2++ ) {
            int dh = dH ( N, n, al, r1, r2 );
            dd[dh]+=2;
        }
    }
    // along diagonal
    dd[0] += N;

    for ( int x=0; x<=n; x++ ) {
        dd[x] /= N;
    }

    return dd;
}



std::vector<double> macwilliams_transform_mixed ( const ndarray<double> &B, const symmetry_group &sg, std::vector<int> sx, int N, ndarray<double> &Bout, int verbose=0 )
{
// TODO: calculation is not very efficient
    if ( verbose ) {
        myprintf ( "macwilliams_transform_mixed:\n" );
#ifdef FULLPACKAGE
        myprintf ( "sx: " );
        display_vector ( sx );
#endif
        myprintf ( "\n" );
    }

    int jprod = B.n; // std::accumulate(sg.gsize.begin(), sg.gsize.begin()+sg.ngroups, 1, multiplies<int>());

    int *bi = new int[sg.ngroups];
    int *iout = new int[sg.ngroups];

    for ( int i=0; i<sg.ngroups; i++ )
        bi[i]=0;
    for ( int i=0; i<sg.ngroups; i++ )
        iout[i]=0;

//--------------

    for ( int j=0; j<Bout.n; j++ ) {
        //if (verbose)
//	myprintf("macwilliams_transform_mixed: Bout[%d]= %d\n", j, 0 );

        Bout.linear2idx ( j, iout );
        //myprintf("  j %d: iout[0] %d\n", j, iout[0]);
        Bout.setlinear ( j, 0 ); // [j]=0;

        for ( int i=0; i<B.n; i++ ) {
            B.linear2idx ( i, bi );

            long fac = 1 ;
            for ( int f=0; f<B.k; f++ ) {
                long ji=iout[f];
                long ii = bi[f];
                long ni = B.dims[f]-1;
                long si = sx[f];
                //myprintf("  f %d ji %d, ni %ld, ni %ld, si %d\n", f, (int)ji, ii, ni, si);
                fac *= krawtchouk ( ji,ii, ni, si );
                if ( verbose>=4 )
                    myprintf ( "  Bout[%d] += %.1f * %ld (%d %d %d)\n", j , ( double ) B.data[i] , fac, ( int ) ji, ( int ) ii, ( int ) ni );
            }
            Bout.data[j] += B.data[i] * fac;

        }
        Bout.data[j] /= N;
        if ( verbose>=2 )
            myprintf ( "macwilliams_transform_mixed: Bout[%d]=Bout%s= %f\n", j, Bout.idxstring ( j ).c_str(), Bout.data[j] );
    }

    if ( verbose>=1 ) {
        myprintf ( "Bout: \n" );
        Bout.show();
    }

	// use formula from page 555 in Xu and Wu  (Theorem 4.i)
    std::vector<double> A ( sg.n+1, 0 );


    for ( int i=0; i<jprod; i++ ) {
        Bout.linear2idx ( i, bi );
        int jsum=0;
        for ( int j=0; j<Bout.k; j++ )
            jsum += bi[j];
        if ( verbose>=2 )
            myprintf ( "   jsum %d/%d, i %d\n", jsum, ( int ) A.size(), i );
        A[jsum]+=Bout.data[i];
    }

    delete [] iout;
    delete [] bi;

    return A;
}

#ifdef FULLPACKAGE
/** Calculate D-efficiencies for all projection designs */
std::vector<double> projDeff ( const array_link &al, int kp, int verbose=0 )
{

    int kk = al.n_columns;
    //colperm_t cp = new_comb_init<int>(kp);
    std::vector<int> cp ( kp );
    for ( int i=0; i<kp; i++ )
        cp[i]=i;
    int64_t ncomb = ncombsm<int64_t> ( kk, kp );

    std::vector<double> dd ( ncomb );

    int m = 1 + kp + kp* ( kp-1 ) /2;
    int N = al.n_rows;

    if ( verbose )
        myprintf ( "projDeff: k %d, kp %d: start with %ld combinations \n", kk, kp, ( long ) ncomb );

    for ( int64_t i=0; i<ncomb; i++ ) {

        array_link alsub = al.selectColumns ( cp );
        if ( m>N )
            dd[i]=0;
        else
            dd[i] = alsub.Defficiency();

        if ( verbose>=2 )
            myprintf ( "projDeff: k %d, kp %d: i %ld, D %f\n", kk, kp, ( long ) i, dd[i] );
        next_comb ( cp ,kp, kk );

    }

    if ( verbose )
        myprintf ( "projDeff: k %d, kp %d: done\n", kk, kp );

    return dd;
}

/**
 *
 * Calculate the projection estimation capacity sequence for a design.
 *
 * See "Ranking Non-regular Designs", J.L. Loeppky
 *
 */
std::vector<double> PECsequence ( const array_link &al, int verbose )
{

    int N = al.n_rows;

    int kk = al.n_columns;
    std::vector<double> pec ( kk );

    if ( kk>=20 ) {
        myprintf ( "PECsequence: error: not defined for 20 or more columns\n" );
        pec[0]=-1;
        return pec;
    }

#ifdef DOOPENMP
    #pragma omp parallel for
#endif
    for ( int i=0; i<kk; i++ ) {
        int kp = i+1;
        int m = 1 + kp + kp* ( kp-1 ) /2;

        if ( m>N ) {
            pec[i]=0;
        } else {
            std::vector<double> dd = projDeff ( al, kp, verbose>=2 );

            double ec=0;
            for ( unsigned long j=0; j<dd.size(); j++ )
                ec += dd[j]>0;

            ec=ec/dd.size();

            pec[i]=ec;
        }
    }

    return pec;
}
#endif



/** calculate GWLP (generalized wordlength pattern)
 *
 * Based on: "GENERALIZED MINIMUM ABERRATION FOR ASYMMETRICAL
FRACTIONAL FACTORIAL DESIGNS", Xu and Wu, 2001
 *
 */
std::vector<double> GWLPmixed ( const array_link &al, int verbose, int truncate )
{
    arraydata_t adata = arraylink2arraydata ( al );
    symmetry_group sg ( adata.getS(), false );
    //myprintf(" adata.ncolgroups: %d/%d\n", adata.ncolgroups, sg.ngroups);


    std::vector<int> dims ( adata.ncolgroups );
    for ( unsigned int i=0; i<dims.size(); i++ )
        dims[i] = adata.colgroupsize[i]+1;

    ndarray<double> B ( dims );
    ndarray<double> Bout ( dims );

    distance_distribution_mixed ( al, B, verbose );
    if ( verbose>=3 ) {
        myprintf ( "GWLPmixed: distance distrubution\n" );
        B.show();
    }

    int N = adata.N;
	// calculate GWLP
    std::vector<int> ss = adata.getS();

    std::vector<int> sx;
    for ( int i=0; i<sg.ngroups; i++ )
        sx.push_back ( ss[sg.gstart[i]] );

    std::vector<double> gma = macwilliams_transform_mixed ( B, sg, sx, N, Bout, verbose );

    if ( truncate ) {
        for ( size_t i=0; i<gma.size(); i++ ) {
            gma[i]=round ( N*N*gma[i] ) / ( N*N );
            if ( gma[i]==0 )
                gma[i]=0;	 // fix minus zero
        }
    }
    return gma;


}

/// calculate GWLP (generalized wordlength pattern)
std::vector<double> GWLP ( const array_link &al, int verbose, int truncate )
{
    int N = al.n_rows;
    int n = al.n_columns;

    array_t *me = std::max_element ( al.array, al.array+ ( N*n ) );
    int s = *me+1;
    int k ;
    for ( k=0; k<al.n_columns; k++ ) {
        array_t *mine = std::max_element ( al.array+N*k, al.array+ ( N* ( k+1 ) ) );
        // myprintf( "s %d, %d\n", s, *mine+1);
        if ( *mine+1 != s )
            break;
    }
    int domixed= ( k<al.n_columns );
    //domixed=1;
    if ( verbose )
        myprintf ( "GWLP: N %d, s %d, domixed %d\n", N, s, domixed );


    if ( domixed ) {
        // see: "Algorithmic Construction of Efficient Fractional Factorial Designs With Large Run Sizes", Xu
        std::vector<double> gma = GWLPmixed ( al, verbose, truncate );
        return gma;
    } else {
        // calculate distance distribution
        std::vector<double> B = distance_distributionT ( al );
#ifdef FULLPACKAGE
        if ( verbose ) {
            myprintf ( "distance_distributionT: " );
            display_vector ( B );
            std::cout << std::endl;
        }
#endif
		// calculate GWP
        std::vector<double> gma = macwilliams_transform ( B, N, s );

        if ( truncate ) {
            for ( size_t i=0; i<gma.size(); i++ ) {
                gma[i]=round ( N*N*gma[i] ) / ( N*N );
                if ( gma[i]==0 )
                    gma[i]=0;	 // fix minus zero

            }
        }

        return gma;
    }
}

/// convert GWLP sequence to unique value
inline double GWPL2val ( GWLPvalue x )
{
    double r=0;
    for ( int i=x.size()-1; i>0; i-- )
        r = r/10 + x.v[i];
    if ( 0 ) {
        r=0;
        for ( int i=x.size()-1; i>=0; i-- ) {
            //double fac=choose(x.size(),i);
            // r = r/fac + x[i];
        }
    }
    return r;
}

/// convert GWLP sequence to unique value
inline double GWPL2val ( std::vector<double> x )
{
    double r=0;
    //for ( size_t i=0; i<x.size(); i++ ) r = r/10 + x[i];
    for ( int i=x.size()-1; i>0; i-- )
        r = r/10 + x[i];
    if ( 0 ) {
        r=0;
        for ( int i=x.size()-1; i>=0; i-- ) {
            //double fac=choose(x.size(),i);
            // r = r/fac + x[i];
        }
    }
    return r;
}

std::vector< GWLPvalue > sortGWLP ( const std::vector< GWLPvalue > in )
{
    std::vector< GWLPvalue > v = in;

    std::sort ( v.begin(), v.end() );
    return v;
}


std::vector< GWLPvalue > projectionGWLPs ( const array_link &al )
{
    int ncols=al.n_columns;

    std::vector< GWLPvalue > v ( ncols );
    for ( int i=0; i<ncols; i++ ) {
        array_link d = al.deleteColumn ( i );
        std::vector<double> gma = GWLP ( d );
        v[i]= gma;
    }
    return v;
}

std::vector<double> projectionGWLPvalues ( const array_link &al )
{
    int ncols=al.n_columns;

    std::vector<double> v ( ncols );
    for ( int i=0; i<ncols; i++ ) {
        array_link d = al.deleteColumn ( i );
        std::vector<double> gma = GWLP ( d );
        v[i]= GWPL2val ( gma );
    }
    return v;
}



/// convert array to Eigen matrix structure
Eigen::MatrixXd arraylink2eigen ( const array_link &al )
{
    int k = al.n_columns;
    int n = al.n_rows;
	assert(n>=0); assert(k>=0);

    Eigen::MatrixXd mymatrix = Eigen::MatrixXd::Zero ( n,k );

    for ( int c=0; c<k; ++c ) {
        //int ci = c*n;
        array_t *p = al.array+c*n;
        for ( int r=0; r<n; ++r ) {
            mymatrix ( r, c ) = p[r];
        }
    }
    return mymatrix;
}

/// return rank of an array based on Eigen::ColPivHouseholderQR
int arrayrankFullPivQR ( const array_link &al, double threshold )
{
    Eigen::MatrixXd mymatrix = arraylink2eigen ( al );
	FullPivHouseholderQR<Eigen::MatrixXd> decomp(mymatrix.rows(), mymatrix.cols());
	decomp.compute(mymatrix);
	if (threshold>0) {
		decomp.setThreshold(threshold);	
	}
    int rank = decomp.rank();
    return rank;
}

/// return rank of an array based on Eigen::ColPivHouseholderQR
int arrayrankColPivQR ( const array_link &al, double threshold  )
{
    Eigen::MatrixXd mymatrix = arraylink2eigen ( al );
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp ( mymatrix );
	if (threshold>0) {
		decomp.setThreshold(threshold);	
	}
    int rank = decomp.rank();
    return rank;
}

/// return rank of an array based on Eigen::FullPivLU
int arrayrankFullPivLU ( const array_link &al, double threshold  )
{
    Eigen::MatrixXd mymatrix = arraylink2eigen ( al );
    Eigen::FullPivLU<Eigen::MatrixXd> decomp ( mymatrix );
	if (threshold>0) {
		decomp.setThreshold(threshold);	
	}
    int rank = decomp.rank();
    return rank;
}

/// return rank of an array based on Eigen::JacobiSVD
int arrayrankSVD ( const array_link &al, double threshold  )
{
    Eigen::MatrixXd mymatrix = arraylink2eigen ( al );
    Eigen::JacobiSVD<Eigen::MatrixXd> decomp ( mymatrix );
	if (threshold>0) {
		decomp.setThreshold(threshold);	
	}
    int rank = decomp.rank();
    return rank;
}

/// return rank of an array based on Eigen::FullPivLU
int arrayrank ( const array_link &al )
{
    Eigen::MatrixXd mymatrix = arraylink2eigen ( al );
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp ( mymatrix );
	//printfd("threshold %e\n", lu_decomp.threshold() );
    int rank = lu_decomp.rank();
    return rank;
}


/// return rank, maxPivot etc.
int arrayrankInfo ( const Eigen::MatrixXd &mymatrix, int verbose )
{
    if (verbose) {
	printfd("arrayrankInfo\n");
	}
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp ( mymatrix );
    int rank = lu_decomp.rank();
    if (verbose) {
		double p = lu_decomp.maxPivot();
        printfd("arrayrankInfo: FullPivLU: rank %d, threshold %e, max pivot %e\n", rank, lu_decomp.threshold(), p );
    }
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr_decomp ( mymatrix );
    int rank2 = qr_decomp.rank();
    if (verbose) {
		Eigen::MatrixXd qr = qr_decomp.matrixQR();
		qr_decomp.colsPermutation();
		
		Eigen::FullPivHouseholderQR<Eigen::MatrixXd>::PermutationType P = qr_decomp.colsPermutation();
		 
		Eigen::MatrixXd R = qr_decomp.matrixQ().inverse() * mymatrix * P;
		if (verbose>=3) {
		eigenInfo(R, "arrayrankInfo: R ");
		std::cout << R << std::endl;
		std::cout << R.diagonal() << std::endl;
		}
		Eigen::VectorXd d = R.diagonal();
		
		double dmin = qr_decomp.maxPivot();
		double dfalse = 0;
		for(int i=0; i<d.size(); i++) 
		{
				double q = std::fabs((double)d(i));
				//printf("i %d, q %e\n", i, q);
				if (q<dmin && q > qr_decomp.threshold())
					dmin=q;
				if (q>dfalse && q < qr_decomp.threshold())
					dfalse=q;
		}
		double p = qr_decomp.maxPivot();
        printfd("arrayrankInfo: FullPivHouseholderQR: rank %d, threshold %e, max pivot %e, min non-zero pivot %e, false pivot %e\n", rank, qr_decomp.threshold(), p, dmin, dfalse );
    }

    return rank;
}

/// return rank of an array based on Eigen::FullPivLU
int arrayrankInfo ( const array_link &al, int verbose )
{
	printfd("xxx\n"); al.show(); fflush(stdout);
    Eigen::MatrixXd mymatrix = arraylink2eigen ( al );
    int rank = arrayrankInfo( mymatrix, verbose );
    return rank;
}

/* Helper functions for rankStructure */

#include <Eigen/Core>
#include <Eigen/Dense>

/// helper function
std::vector<int> subIndices ( int ks, int k )
{
    const int m = 1+k+k* ( k-1 ) /2;
    const int msub = 1+ks+ks* ( ks-1 ) /2;
    std::vector<int> idxsub ( msub );
    for ( int i=0; i<ks+1; i++ )
        idxsub[i]=i;
    int n = ks* ( ks-1 ) /2;
    for ( int i=0; i<n; i++ )
        idxsub[i+1+ks]=i+1+k;

    return idxsub;
}
/// helper function
std::vector<int> subIndicesRemainder ( int ks, int k )
{
    const int m = 1+k+k* ( k-1 ) /2;
    const int msub = 1+ks+ks* ( ks-1 ) /2;

    const int t2 = k* ( k-1 ) /2;
    const int t2s = ks* ( ks-1 ) /2;

    std::vector<int> idxsub ( m-msub );

    for ( int i=0; i< ( k-ks ); i++ )
        idxsub[i]=1+ks+i;
    for ( int i=0; i< ( t2-t2s ); i++ )
        idxsub[ ( k-ks ) +i]=1+k+t2s+i;

    return idxsub;
}
/// helper function
Eigen::MatrixXi permM ( int ks, int k, const Eigen::MatrixXi subperm, int verbose=1 )
{
    std::vector<int> idxsub = subIndices ( ks, k );
    std::vector<int> idxrem = subIndicesRemainder ( ks, k );

    if ( verbose ) {
        myprintf ( "ks: %d, k %d, idxsub: ",ks, k );
        print_perm ( idxsub ); //printf("\n");
        myprintf ( "ks: %d, k %d, idxrem: ",ks, k );
        print_perm ( idxrem ); //printf("\n");
    }

    const int m = 1+k+k* ( k-1 ) /2;
    const int msub = 1+ks+ks* ( ks-1 ) /2;

    std::vector<int> ww ( idxsub.size() );

    for ( size_t i=0; i<ww.size(); i++ )
        ww[i] = idxsub[subperm ( i )];

    Eigen::MatrixXi pm ( m,1 );
    for ( size_t i=0; i<ww.size(); i++ )
        pm ( i ) =ww[i];
    for ( int i=0; i< ( m-msub ); i++ )
        pm ( msub+i ) =idxrem[i];

    return pm;
}

/// return the condition number of a matrix
double conditionNumber ( const array_link &M ) {
    Eigen::MatrixXd A = arraylink2eigen ( M );
    Eigen::JacobiSVD<Eigen::Matrix<double,-1,-1> > svd ( A );
    double cond = svd.singularValues() ( 0 ) / svd.singularValues() ( svd.singularValues().size()-1 );
    return cond;
}

/// calculate the rank of the second order interaction matrix of an array using the cache system
int rankStructure::rankxf ( const array_link &al )
{
    this->ncalc++;

    ks = this->ks;
    int k = al.n_columns;
    const int m = 1+k+k* ( k-1 ) /2;
    const int msub = 1+ks+ks* ( ks-1 ) /2;
    const int N=al.n_rows;

    if (verbose) {
        if(al.n_columns!=alsub.n_columns+nsub)
            printf ( "rankStructure: rankxf: alsub %d, al %d (nsub %d, ks %d, id %d)\n", alsub.n_columns, al.n_columns, nsub, ks, this->id );
    }
    if ( verbose )
        printf ( "rankStructure: rankxf: alsub %d, al %d (ks %d)\n", alsub.n_columns, al.n_columns, ks );
    if ( al.selectFirstColumns ( ks ) == this->alsub ) {

    } else {
        // update structure
        if (verbose>=2)
            printf ( "rankStructure: update structure (current ks %d, al.n_columns %d)\n", ks, al.n_columns );
        updateStructure ( al.selectFirstColumns ( al.n_columns-nsub ) );

    }
    int rank0 = decomp.rank();
    if ( verbose>=2 )
        printfd ( "rankStructure: rank0 %d\n", rank0 );

	// special case: the same matrix!
	if (ks==al.n_columns) {
		if (verbose)
			printfd("special case: k==al.n_columns\n");
		return decomp.rank();
	}
    Eigen::MatrixXd A = array2xfeigen ( al );

    // caculate permutation
	
	
    EigenDecomp::PermutationType subperm = decomp.colsPermutation ();
    MatrixXi perm = permM ( ks, k, subperm.indices(), verbose );
    EigenDecomp::PermutationType ptmp ( perm );

    // transform second order interaction matrix into suitable format

    // method 2
    //Eigen::MatrixXd Zxp = (A * ptmp).block(0, rank0, N, m-rank0);
    //Eigen::MatrixXd ZxSub = this->Qi.block(rank0, 0, N-rank0, N)*(A * ptmp).block(0, rank0, N, m-rank0);
    // method 1
    Eigen::MatrixXd Zxp = A * ptmp;
    Eigen::MatrixXd ZxSub = this->Qi.block ( rank0, 0, N-rank0, N ) *Zxp.block ( 0, rank0, N, m-rank0 );

    // direct version
    //Eigen::MatrixXd Zx = this->Qi*A * ptmp;
    //Eigen::MatrixXd ZxSub = Zx.block ( rank0, rank0, Zx.rows()-rank0, m-rank0 );


    if ( verbose>=2 ) {
        printf ( "  rankStructure: k %d, m %d\n", k, m );
        //eigenInfo ( Zx, "Zx" );
        printf ( "  rankStructure: msub %d, m %d\n", msub, m );
    }

    if ( verbose>=2 ) {
        printfd ( "rankStructure: ZxSub\n" );
        eigenInfo ( ZxSub );
		
		arrayrankInfo(ZxSub, 1);
		
		if (verbose>=3) {
			fflush(stdout);
			printf("ZxSub\n");
		std::cout << ZxSub; printf("\n");	
		}
    }

    int rankx=rankdirect ( ZxSub );
    int rank = rank0 + rankx;

    if ( verbose ) {
        printf ( "rankStructure: rank %d + %d = %d (%d)\n", rank0, rankx, rank, rankxfdirect ( al ) );
    }
    return rank;
}


Eigen::MatrixXd array2xfeigen ( const array_link &al )
{
    const int k = al.n_columns;
    const int n = al.n_rows;
    const int m = 1 + k + k* ( k-1 ) /2;
    Eigen::MatrixXd mymatrix = Eigen::MatrixXd::Zero ( n,m );

    // init first column
    int ww=0;
    //mymatrix.col(0).setConstant(1);
    for ( int r=0; r<n; ++r ) {
        mymatrix ( r, 0 ) =1;
    }
    // init array
    ww=1;
    for ( int c=0; c<k; ++c ) {
        int ci = c*n;
        for ( int r=0; r<n; ++r ) {
            mymatrix(r, ww+c) = 2*al.array[r+ci]-1;
        }
    }

    //	mymatrix.block(0,0, n, k+1).array() -= 1;

    // init interactions
    ww=k+1;
    for ( int c=0; c<k; ++c ) {
        int ci = c+1;
        for ( int c2=0; c2<c; ++c2 ) {
            int ci2 = c2+1;

            for ( int r=0; r<n; ++r ) {
                mymatrix(r,ww) = -mymatrix(r, ci)*mymatrix(r,ci2);
            }
            ww++;
        }
    }

    return mymatrix;
}



/// convert 2-level design to second order interaction matrix
inline void array2eigenxf ( const array_link &al, Eigen::MatrixXd &mymatrix )
{
    int k = al.n_columns;
    int n = al.n_rows;
    int m = 1 + k + k* ( k-1 ) /2;

    mymatrix = Eigen::MatrixXd::Zero ( n,m );

    // init first column
    int ww=0;
    for ( int r=0; r<n; ++r ) {
        mymatrix ( r, ww ) = 1;
    }

    // init array
    ww=1;
    for ( int c=0; c<k; ++c ) {
        int ci = c*n;
        for ( int r=0; r<n; ++r ) {
            mymatrix ( r, ww+c ) = al.array[r+ci];
        }
    }

    // init interactions
    ww=k+1;
    for ( int c=0; c<k; ++c ) {
        int ci = c*n;
        for ( int c2=0; c2<c; ++c2 ) {
            int ci2 = c2*n;

            const array_t * p1 = al.array+ci;
            const array_t * p2 = al.array+ci2;
            for ( int r=0; r<n; ++r ) {
                mymatrix ( r, ww ) = ( *p1+*p2 ) %2;
                p1++;
                p2++;
                //mymatrix ( r, ww ) = ( p1[r]+p2[r] ) %2;
            }
            ww++;
        }
    }

    mymatrix.array() *= 2;
    mymatrix.array() -= 1;
}


array_link array2secondorder ( const array_link &al )
{
    int k = al.n_columns;
    int n = al.n_rows;
    int m = 1 + k + k* ( k-1 ) /2;
    int m2 = k* ( k-1 ) /2;
    array_link out ( n, m2, array_link::INDEX_DEFAULT );

    // init interactions
    int ww=0;
    for ( int c=0; c<k; ++c ) {
        int ci = c*n;
        for ( int c2=0; c2<c; ++c2 ) {
            int ci2 = c2*n;

            const array_t * p1 = al.array+ci;
            const array_t * p2 = al.array+ci2;
            array_t *pout = out.array+ww*out.n_rows;

            for ( int r=0; r<n; ++r ) {
                pout[r] = ( p1[r]+p2[r] ) %2 ;
            }
            ww++;
        }
    }

    out *= 2;
    out -= 1;

    return out;

}

array_link array2xf ( const array_link &al )
{
    const int k = al.n_columns;
    const int n = al.n_rows;
    const int m = 1 + k + k* ( k-1 ) /2;
	//assert(k>=0);
	//assert(n>=0);

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
                pout[r] = -( p1[r]*p2[r] )  ;
            }
            ww++;
        }
    }
    return out;
}

using namespace Eigen;

#include <Eigen/LU>


void DAEefficiecyWithSVD ( const Eigen::MatrixXd &x, double &Deff, double &vif, double &Eeff, int &rank, int verbose )
{
    // printfd("start\n");

    Eigen::FullPivLU<MatrixXd> lu_decomp ( x );
    rank = lu_decomp.rank();

    JacobiSVD<Eigen::MatrixXd> svd ( x );

    const Eigen::VectorXd S = svd.singularValues();
    int rank2 = svd.nonzeroSingularValues();
    if ( rank2!=rank ) {
        if ( verbose>=3 ) {
            myprintf ( "ABwithSVD: rank calculations differ, unstable matrix: ranklu %d, ranksvd: %d\n", rank, rank2 );

#ifdef FULLPACKAGE
            if ( verbose>=4 ) {
                std::cout << "Its singular values are:" << std::endl << S << std::endl;
            }
#endif
        }
    }
    int m = x.cols();
    int N= x.rows();

    if ( m>N ) {
        Deff=0;
        vif=0;
        Eeff=0;
#ifdef FULLPACKAGE
        if ( verbose>=3 ) {
            Eigen::MatrixXd Smat ( S );
            Eigen::ArrayXd Sa=Smat.array();
            double Deff = exp ( 2*Sa.log().sum() /m ) /N;

            std::cout << printfstring ( "  singular matrix: m (%d) > N (%d): rank ", m, N ) << rank << std::endl;
            int mx = std::min ( m, N );
            myprintf ( "   Deff %e, smallest eigenvalue %e, rank lu %d\n", Deff, S[mx-1], rank );

        }
#endif
        return;
    }
    if ( verbose>=3 )
        myprintf ( "N %d, m %d\n", N, m );

    if ( S[m-1]<1e-15 || rank < m ) {
        if ( verbose>=2 ) {
            myprintf ( "   array is singular, setting D-efficiency to zero\n" );

            if ( 0 ) {
                Eigen::MatrixXd mymatrix = x;
                Eigen::MatrixXd mm = mymatrix.transpose() * mymatrix;
                SelfAdjointEigenSolver<Eigen::MatrixXd> es;
                es.compute ( mm );
                const Eigen::VectorXd evs =  es.eigenvalues();
                Eigen::VectorXd S = evs; //sqrt(S);
                for ( int j=0; j<m; j++ ) {
                    if ( S[j]<1e-14 ) {
                        if ( verbose>=3 )
                            myprintf ( "  singular!\n" );
                        S[j]=0;
                    } else {

                        S[j]=sqrt ( S[j] );
                    }
                }

                Eigen::MatrixXd Smat ( S );
                Eigen::ArrayXd Sa=Smat.array();
                double DeffDirect = exp ( 2*Sa.log().sum() /m ) /N;

                myprintf ( "  DeffDirect %f\n", DeffDirect );

            }
        }
        Deff=0;
        vif=0;
        Eeff=0;

        //Eigen::FullPivLU<MatrixXd> lu(x);
        //rank=lu.rank();
#ifdef FULLPACKAGE
        int rankold=rank2;
        if ( verbose>=3 ) {
            Eigen::MatrixXd Smat ( S );
            Eigen::ArrayXd Sa=Smat.array();
            double Deff = exp ( 2*Sa.log().sum() /m ) /N;

            std::cout << "  singular matrix: the rank of A is " << rank << std::endl;
            myprintf ( "   Deff %e, smallest eigenvalue %e, rankold %d, rank lu %d\n", Deff, S[m-1], rankold, rank );

        }
        if ( verbose>=4 )
            std::cout << "Its singular values are:" << std::endl << S << std::endl;
#endif
        return;
    }
#ifdef FULLPACKAGE
    if ( verbose>=3 )
        std::cout << "Its singular values are:" << std::endl << S << std::endl;
#endif
//cout << "x:\n"<< x << std::endl;

    Eeff = S[m-1]*S[m-1]/N;
    //cout << "Its singular values are:" << endl << S/sqrt(N) << endl;

    vif=0;
    for ( int i=0; i<m; i++ )
        vif+= 1/ ( S[i]*S[i] );
    vif = N*vif/m;

//Eigen::VectorXf Si= S.inverse();
//cout << "Its inverse singular values are:" << endl << Si << endl;

    Eigen::MatrixXd Smat ( S );
    Eigen::ArrayXd Sa=Smat.array();
    Deff = exp ( 2*Sa.log().sum() /m ) /N;

    if ( verbose>=2 ) {
        myprintf ( "ABwithSVD: Defficiency %.3f, Aefficiency %.3f (%.3f), Eefficiency %.3f\n", Deff, vif, vif*m, Eeff );

        Eigen::FullPivLU<MatrixXd> lu ( x );
        int ranklu=lu.rank();

        myprintf ( "   Defficiency %e, smallest eigenvalue %e, rank %d, rank lu %d\n", Deff, S[m-1], rank, ranklu );

    }
}


/** Calculate the rank of an orthogonal array (rank of [I X X_2] )
 *
 * The vector ret is filled with the rank, Defficiency, VIF efficiency and Eefficiency
 */
int array_rank_D_B ( const array_link &al, std::vector<double> *ret  , int verbose )
{
    //printfd("start\n");
    int k = al.n_columns;
    int n = al.n_rows;
    int m = 1 + k + k* ( k-1 ) /2;

    Eigen::MatrixXd mymatrix ( n, m );
    array2eigenxf ( al, mymatrix );

    double Deff;
    double B;
    double Eeff;
    int rank;

    //ABold(mymatrix, A, B, rank, verbose);
    DAEefficiecyWithSVD ( mymatrix, Deff, B, Eeff, rank, verbose );

    if ( ret!=0 ) {
        ret->push_back ( rank );
        ret->push_back ( Deff );
        ret->push_back ( B );
        ret->push_back ( Eeff );
    }
    return rank;

//	   cout << "Here is the matrix A:\n" << A << endl;
    // FullPivLU<Matrix3f> lu_decomp(A);
    //cout << "The rank of A is " << lu_decomp.rank() << endl;
}

std::vector<double> Aefficiencies ( const array_link &al, int verbose )
{
    myassert ( al.is2level() );
    int N = al.n_rows;
    int k = al.n_columns;
    int m = 1 + k + ( k* ( k-1 ) ) /2;
    if ( verbose )
        printf ( "Aefficiencies: array %d, %d\n", N, k );
    MatrixFloat X = array2eigenModelMatrix ( al );

    Eigen::FullPivLU<MatrixXd> lu_decomp ( X );
    int rank = lu_decomp.rank();
    if ( rank<m ) {
        if ( verbose )
            printf ( "Aefficiencies: rank %d/%d\n", rank, m );
        std::vector<double> aa ( 3 );
        aa[0] = 0;
        aa[1] = 0;
        aa[2] = 0;
        return aa;
    }

    if ( verbose )
        printf ( "Aefficiencies: calculate information matrix\n" );

    MatrixFloat matXtX = ( X.transpose() * ( X ) ) /N;

    if ( verbose )
        printf ( "Aefficiencies: invert information matrix\n" );
    MatrixFloat M = matXtX.inverse();

    std::vector<double> aa ( 3 );
    double Ax=0;
    for ( int i=0; i<m; i++ ) {
        Ax += M ( i,i );
    }
    aa[0] = 1./ ( Ax/m );
    Ax=0;
    for ( int i=1; i<k+1; i++ ) {
        Ax += M ( i,i );
    }

    aa[1] = 1./ ( Ax/k );
    Ax=0;
    for ( int i=k+1; i<m; i++ ) {
        Ax += M ( i,i );
    }
    aa[2] = 1./ ( Ax/ ( k* ( k-1 ) /2 ) );
    return aa;
}

double VIFefficiency ( const array_link &al, int verbose )
{
    std::vector<double> ret;
    array_rank_D_B ( al, &ret, verbose );
    return ret[2];
}

double Aefficiency ( const array_link &al, int verbose )
{
    std::vector<double> ret;
    array_rank_D_B ( al, &ret, verbose );
    if ( ret[2]==0 )
        return 0;
    else
        return 1./ret[2];
}

double Eefficiency ( const array_link &al, int verbose )
{
    std::vector<double> ret;
    int r = array_rank_D_B ( al, &ret, verbose );
    return ret[3];
}

std::vector<int> Jcharacteristics ( const array_link &al, int jj, int verbose )
{
    jstruct_t js ( al, jj );
    return js.values;
}

/// calculate determinant of X^T X by using the SVD
double detXtX ( const Eigen::MatrixXd &mymatrix, int verbose )
{
    double dd=-1;
    double dd2=-1;
    //int n = mymatrix.rows();
    int m = mymatrix.cols();
    //int N = n;

    //    Eigen::FullPivLU< Eigen::MatrixXd> lu_decomp(mymatrix);
    //  int rank = lu_decomp.rank();

    Eigen::MatrixXd mm = mymatrix.transpose() * mymatrix;
    SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute ( mm );
    const Eigen::VectorXd evs =  es.eigenvalues();
    Eigen::VectorXd S = evs; //sqrt(S);

    if ( S[m-1]<1e-15 ) {
        if ( verbose>=2 ) {
            myprintf ( "   array is singular, setting det to zero\n" );

        }
        dd=0;
        return dd;
    }

    for ( int j=0; j<m; j++ ) {
        if ( S[j]<1e-14 ) {
            if ( verbose>=3 )
                myprintf ( "  singular!\n" );
            S[j]=0;
        } else {
            S[j]=sqrt ( S[j] );
        }
    }

    Eigen::MatrixXd Smat ( S );
    Eigen::ArrayXd Sa=Smat.array();
    dd = exp ( 2*Sa.log().sum() );

    //

    if ( S[0]<1e-15 ) {
        if ( verbose>=2 )
            myprintf ( "Avalue: singular matrix\n" );
        dd=0;
        return dd;
    }
    return dd;
}

typedef Eigen::MatrixXf MyMatrixf;
typedef Eigen::ArrayXf MyArrayf;
typedef Eigen::VectorXf MyVectorf;


double detXtXfloat ( const MyMatrixf &mymatrix, int verbose )
{
    double dd=-1;
    //double dd2=-1;
    //int n = mymatrix.rows();
    int m = mymatrix.cols();
    //int N = n;

    //    Eigen::FullPivLU< Eigen::MatrixXd> lu_decomp(mymatrix);
    //  int rank = lu_decomp.rank();

    MyMatrixf mm = mymatrix.transpose() * mymatrix;
    SelfAdjointEigenSolver<MyMatrixf> es;
    es.compute ( mm );
    const MyVectorf evs =  es.eigenvalues();
    MyVectorf S = evs; //sqrt(S);

    if ( S[m-1]<1e-15 ) {
        if ( verbose>=2 ) {
            myprintf ( "   array is singular, setting det to zero\n" );

        }
        dd=0;
        return dd;
    }

    for ( int j=0; j<m; j++ ) {
        if ( S[j]<1e-14 ) {
            if ( verbose>=3 )
                myprintf ( "  singular!\n" );
            S[j]=0;
        } else {
            S[j]=sqrt ( S[j] );
        }
    }

    MyMatrixf Smat ( S );
    MyArrayf Sa=Smat.array();
    dd = exp ( 2*Sa.log().sum() );

    //

    if ( S[0]<1e-15 ) {
        if ( verbose>=2 )
            myprintf ( "Avalue: singular matrix\n" );
        dd=0;
        return dd;
    }
    return dd;
}

//typedef Eigen::MatrixXd MyMatrix;
typedef MatrixFloat MyMatrix;

std::vector<double> Defficiencies ( const array_link &al, const arraydata_t & arrayclass, int verbose, int addDs0 )
{
    if ( ( al.n_rows>500 ) || ( al.n_columns>500 ) ) {
        myprintf ( "Defficiencies: array size not supported\n" );
        return std::vector<double> ( 3 );
    }

    //if ( al.min()<0 ) {
    //	myprintf ( "Defficiencies: array elements should range from 0 to s-1\n" );
    //	return std::vector<double> ( 3 );
    //}

    int k = al.n_columns;
    int k1 = al.n_columns+1;
    int n = al.n_rows;
    int N=n;
    int m = 1 + k + k* ( k-1 ) /2;

    MyMatrix X;

    int n2fi = -1; /// number of 2-factor interactions in contrast matrix
    int nme = -1;	/// number of main effects in contrast matrix

    if ( arrayclass.is2level() ) {
        //Eigen::MatrixXi Xint =  array2eigenModelMatrixInt ( al );
        //matXtX = ( Xint.transpose() *(Xint) ).cast<eigenFloat>() /n;

        X = array2eigenModelMatrix ( al );



        //X1i = array2eigenX1 ( al, 1 );
        //	X2 = array2eigenX2 ( al );
        //X2 = X.block ( 0, 1+k, N, k* ( k-1 ) /2 );
        //std::cout << "X2\n" << X2 << std::endl; std::cout << "X2a\n" << X2a << std::endl;
        //if (X2a==X2) myprintf("equal!"); exit(0);

        n2fi =  k* ( k-1 ) /2;
        nme = k;
    } else {
        if ( verbose>=2 )
            myprintf ( "Defficiencies: mixed design!\n" );
        std::pair<MyMatrix,MyMatrix> mm = array2eigenModelMatrixMixed ( al, 0 );
        const MyMatrix &X1=mm.first;
        const MyMatrix &X2=mm.second;
        //eigenInfo(X1i, "X1i");
        //X1i.resize ( N, 1+X1.cols() );
        //X1i << MyMatrix::Constant ( N, 1, 1 ),  X1;
        //X2=mm.second;
        X.resize ( N, 1+X1.cols() +X2.cols() );
        X << MyMatrix::Constant ( N, 1, 1 ),  X1, X2;

        n2fi = X2.cols();
        nme = X1.cols();


    }
    MyMatrix matXtX = ( X.transpose() * ( X ) ) /n;

    //Matrix X1i =  X.block(0,0,N, 1+nme);

    /*
    // https://forum.kde.org/viewtopic.php?f=74&t=85616&p=146569&hilit=multiply+transpose#p146569
    MyMatrix matXtX(X.cols(), X.cols());
    //		matXtX.setZero();
    //		myprintf("assign triangularView\n");
    matXtX.triangularView<Eigen::Upper>() = X.transpose() *X/n;

    //		matXtX.sefladjointView<Upper>().rankUpdate(X.transpose() );

    //		myprintf("assign triangularView (lower)\n");

    matXtX.triangularView<Eigen::StrictlyLower>() = matXtX.transpose();
    */

    if ( 0 ) {
        Eigen::MatrixXf dummy = matXtX.cast<float>();
        double dum1 = matXtX.determinant();
        double dum2 = dummy.determinant();
        myprintf ( "%f %f\n", dum1, dum2 );
    }
    double f1 = matXtX.determinant();
    //double f2 = ( matXtX.block(1+nme, 1+nme, n2fi, n2fi) ).determinant();


    int nm = 1+nme+n2fi;

    //if (1) {
    MyMatrix tmp ( nm,1+n2fi ); //tmp.resize ( nm, 1+n2fi);
    tmp << matXtX.block ( 0, 0, nm, 1 ), matXtX.block ( 0, 1+nme, nm, n2fi );
    MyMatrix mX2i ( 1+n2fi, 1+n2fi ); //X2i.resize ( 1+n2fi, 1+n2fi);
    mX2i << tmp.block ( 0, 0, 1, 1+n2fi ), tmp.block ( 1+nme, 0, n2fi, 1+n2fi );

    double f2i = ( mX2i ).determinant();
    //double f2i = ((X2i.transpose()*X2i)/n ).determinant();

    //if (fabs(f2i)<1e-15) { std::cout << mX2i << std::endl; exit(1); }


    //double f2 = ( X2.transpose() *X2/n ).determinant();
    double t = ( matXtX.block ( 0,0,1+nme, 1+nme ) ).determinant();
    //double t = ( X1i.transpose() *X1i/n ).determinant();

    double D=0, Ds=0, D1=0;
    int rank = m;
    if ( fabs ( f1 ) <1e-15 ) {
        Eigen::FullPivLU<MyMatrix> lu_decomp ( X );
        rank = lu_decomp.rank();

        if ( verbose>=1 ) {
            myprintf ( "Defficiencies: rank of model matrix %d/%d, f1 %e, f2i %e\n", rank, m, f1, f2i );
        }
    }
    if ( rank<m ) {
        if ( verbose>=1 ) {
            myprintf ( "Defficiencies: model matrix does not have max rank, setting D-efficiency to zero (f1 %e)\n", f1 );
            myprintf ( "   rank lu_decomp %d/%d\n", rank, m ) ;
            myprintf ( "   calculated D %f\n", pow ( f1, 1./m ) );
        }

    } else {
        if ( verbose>=2 ) {
            myprintf ( "Defficiencies: f1 %f, f2i %f, t %f\n", f1, f2i, t );
        }


        Ds = pow ( ( f1/f2i ), 1./k );
        D = pow ( f1, 1./m );
    }
    D1 = pow ( t, 1./k1 );

    if ( verbose>=2 ) {
        myprintf ( "Defficiencies: D %f, Ds %f, D1 %f\n", D, Ds, D1 );
    }

    /*
    if ( 0 ) {
    	double dd = detXtX( X/sqrt ( double ( n ) ) );
    	double Dnew = pow ( dd, 1./m );
    	myprintf ( "D %.15f -> %.15f\n", D, Dnew );
    } */

    std::vector<double> d ( 3 );
    d[0]=D;
    d[1]=Ds;
    d[2]=D1;

    if ( addDs0 ) {
        double f2 = ( matXtX.block ( 1+nme, 1+nme, n2fi, n2fi ) ).determinant();
        double Ds0=0;
        if ( fabs ( f1 ) >=1e-15 ) {
            Ds0 = pow ( ( f1/f2 ), 1./k1 );
        }
        d.push_back ( Ds0 );
    }
    return d;

}

typedef MatrixFloat DMatrix;
typedef VectorFloat DVector;
typedef ArrayFloat DArray;

//typedef Eigen::MatrixXd DMatrix; typedef Eigen::VectorXd DVector; typedef  Eigen::ArrayXd DArray;

double Defficiency ( const array_link &al, int verbose )
{
    int k = al.n_columns;
    int n = al.n_rows;
    int m = 1 + k + k* ( k-1 ) /2;
    int N = n;
    double Deff = -1;

    //DMatrix mymatrix(n, m); array2eigenxf(al, mymatrix);
    DMatrix mymatrix = array2eigenModelMatrix ( al );

    Eigen::FullPivLU<DMatrix> lu_decomp ( mymatrix );
    int rank = lu_decomp.rank();

    DMatrix mm = mymatrix.transpose() * mymatrix;
    SelfAdjointEigenSolver<DMatrix> es;
    es.compute ( mm );
    const DVector evs =  es.eigenvalues();
    DVector S = evs; //sqrt(S);

    if ( S[m-1]<1e-15 || rank < m ) {
        if ( verbose>=2 ) {

            Eigen::FullPivLU<DMatrix> lu_decomp2 ( mm );
            int rank2 = lu_decomp2.rank();

            myprintf ( "Defficiency: array is singular (rank %d/%d/%d), setting D-efficiency to zero (S[m-1] %e)\n", rank, rank2, m, S[m-1] );

        }
        Deff=0;
        return Deff;
    }

    for ( int j=0; j<m; j++ ) {
        if ( S[j]<1e-14 ) {
            if ( verbose>=3 )
                myprintf ( "  singular!\n" );
            S[j]=0;
        } else {
            S[j]=sqrt ( S[j] );
        }
    }

    if ( verbose>=2 ) {
        JacobiSVD<DMatrix> svd ( mymatrix );
        const DVector S2 = svd.singularValues();

        if ( ! ( S2[m-1]<1e-15 ) ) {
            for ( int ii=0; ii<m-3; ii++ ) {
                myprintf ( "ii %d: singular values sqrt(SelfAdjointEigenSolver) SelfAdjointEigenSolver svd %f %f %f\n", ii, ( double ) S[m-ii-1], ( double ) evs[m-ii-1], ( double ) S2[ii] );
            }
            for ( int ii=m-3; ii<m-1; ii++ ) {
                myprintf ( "ii %d: %f %f %f\n", ii, ( double ) S[m-ii-1], ( double ) evs[m-ii-1], ( double ) S2[ii] );
            }

            DMatrix Smat ( S );
            DArray Sa=Smat.array();
            Deff = exp ( 2*Sa.log().sum() /m ) /N;
            myprintf ( "  Aold: %.6f\n", Deff );
        }
    }

    //   int rank = svd.nonzeroSingularValues();

    //myprintf("Avalue: S size %d %d\n", (int) S.cols(), (int) S.rows());

    if ( S[0]<1e-15 ) {
        if ( verbose>=2 )
            myprintf ( "Avalue: singular matrix\n" );
        Deff=0;
        return Deff;
    }

    DMatrix Smat ( S );
    DArray Sa=Smat.array();
    Deff = exp ( 2*Sa.log().sum() /m ) /N;

    if ( verbose>=2 ) {
        myprintf ( "Avalue: A %.6f (S[0] %e)\n", Deff, S[0] );
    }

    Deff = std::min ( Deff,1. );	// for numerical stability
    return Deff;
}

double CL2discrepancy ( const array_link &al )
{
    const int m = al.n_columns;

    std::vector<double> gwp = al.GWLP();

    double v = 1;
    for ( int k=1; k<=m; k++ ) {
        v += gwp[k]/std::pow ( double ( 9 ), k );
    }
    double w = std::pow ( double ( 13 ) /12, m ) - 2*pow ( 35./32, m ) + std::pow ( 9./8, m ) * v;

    return w;
}

#ifdef FULLPACKAGE

typedef Pareto<mvalue_t<long>, array_link >::pValue ( *pareto_cb ) ( const array_link &, int ) ;

/** Calculate the Pareto optimal arrays from a list of array files

    Pareto optimality is calculated according to (rank; A3,A4; F4)
*/
void calculateParetoEvenOdd ( const std::vector<std::string> infiles, const char *outfile, int verbose, arrayfilemode_t afmode, int nrows, int ncols, paretomethod_t paretomethod )
{
    pareto_cb paretofunction = calculateArrayParetoRankFA<array_link>;
    switch ( paretomethod ) {
    case PARETOFUNCTION_J5:
        paretofunction = calculateArrayParetoJ5<array_link>;
        break;
    default
            :
        break;
    }
    Pareto<mvalue_t<long>,array_link> pset;

    long ntotal=0;
    for ( size_t i=0; i<infiles.size(); i++ ) {
        // open arrayfile
        arrayfile_t af ( infiles[i].c_str() );

        if ( verbose ) {
            myprintf ( "calculateParetoEvenOdd: read file %s (%d arrays)\n", af.filename.c_str(), af.narrays );
        }
        int narrays = af.narrays;
        //#pragma omp parallel for num_threads(4) schedule(dynamic,1)
        #pragma omp parallel for
        for ( int k=0; k<narrays; k++ ) {
            array_link al;
            #pragma omp critical
            {
                al= af.readnext();
                ntotal++;
            }
            Pareto<mvalue_t<long>, array_link >::pValue p = paretofunction ( al, verbose>=3 );

            if ( verbose>=2 ) {
                printf ( "values: " );
                Pareto<mvalue_t<long>, array_link >::showvalue ( p );

            }
            #pragma omp critical
            {
                // add the new tuple to the Pareto set
                pset.addvalue ( p, al );
            }

            #pragma omp critical
            {
                if ( verbose>=2 || ( k%500000==0 && k>0 ) ) {
                    printf ( "calculateParetoEvenOdd: file %d/%d, array %d/%d\n", ( int ) i, ( int ) infiles.size(), k, narrays );
                }
            }
        }
    }

    if ( verbose )
        printf ( "calculateParetoEvenOdd: %ld arrays -> %d pareto values, %d pareto arrays \n", ntotal, pset.number(),pset.numberindices() );

    if ( verbose ) {
        pset.show ( verbose );
    }

    arraylist_t lst = pset.allindicesdeque();

    // write files to disk
    if ( verbose )
        printf ( "calculateParetoEvenOdd: writing arrays to file %s\n", outfile );
    if ( verbose>=3 ) {
        printf ( "calculateParetoEvenOdd: afmode %d (TEXT %d)\n", afmode, ATEXT );

    }

    if ( outfile!=0 ) {
        writearrayfile ( outfile, &lst, afmode, nrows, ncols );
    }
    return;
}


Pareto<mvalue_t<long>,long> parsePareto ( const arraylist_t &arraylist, int verbose, paretomethod_t paretomethod )
{
    pareto_cb paretofunction = calculateArrayParetoRankFA<array_link>;
    switch ( paretomethod ) {
    case PARETOFUNCTION_J5:
        paretofunction = calculateArrayParetoJ5<array_link>;
        break;
    default
            :
        break;
    }


    Pareto<mvalue_t<long>,long> pset;
    pset.verbose=verbose;

    #pragma omp parallel for num_threads(4) schedule(dynamic,1)
    for ( size_t i=0; i<arraylist.size(); i++ ) {
        if ( verbose>=2 || ( ( i%2000==0 ) && verbose>=1 ) ) {
            myprintf ( "parsePareto: array %ld/%ld\n", ( long ) i, ( long ) arraylist.size() );
        }
        if ( ( ( i%10000==0 ) && verbose>=1 ) ) {
            pset.show ( 1 );
        }
        const array_link &al = arraylist.at ( i );

        Pareto<mvalue_t<long>, long>::pValue p = paretofunction ( al, verbose );
        #pragma omp critical
        {
            // add the new tuple to the Pareto set
            pset.addvalue ( p, i );
        }
    }
    return pset;
}

void addArray(Pareto < mvalue_t < long >, long > &pset, const array_link &al, int idx, int verbose, paretomethod_t paretomethod )
{
    pareto_cb paretofunction = calculateArrayParetoRankFA<array_link>;
    switch ( paretomethod ) {
    case PARETOFUNCTION_J5:
        paretofunction = calculateArrayParetoJ5<array_link>;
        break;
    default
            :
        break;
    }
        Pareto<mvalue_t<long>, long>::pValue p = paretofunction ( al, verbose );
    pset.addvalue(p, idx);
}

template  Pareto < mvalue_t < long >, array_link >::pValue calculateArrayParetoJ5<array_link> (const array_link & al, int verbose);
template  Pareto < mvalue_t < long >, int >::pValue calculateArrayParetoJ5<int> (const array_link & al, int verbose);
template  Pareto < mvalue_t < long >, long >::pValue calculateArrayParetoJ5<long> (const array_link & al, int verbose);


#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
