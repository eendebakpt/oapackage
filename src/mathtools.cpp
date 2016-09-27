#include "mathtools.h"

#include <numeric>


#ifdef RPACKAGE
#else
void set_srand ( unsigned int s )
{
    srand ( s );
}
#endif

template<class Type>
void symmetry_group::init ( const std::vector<Type> vals, bool ascendingx, int verbose )
{
    //verbose=2;
#ifdef RPACKAGE
#else
    if ( verbose>=2 ) {
        myprintf ( "symmetry_group::init: %ld elements: ", vals.size() );
        for ( size_t i=0; i<vals.size(); i++ )
            std::cout << vals[i] << ", ";
        myprintf ( "\n" );
    }
#endif
    n=vals.size();
    ascending= ascendingx;

    // check we are sorted
    if ( verbose>=2 )
        myprintf ( "symmetry_group::init: check sorting\n" );

    indexsort is ( vals.size() );

    if ( ascending )
        is.sort ( vals );
    else
        is.sortdescending ( vals );

    //      myprintf("symmetry_group::init: xxx\n");

    if ( verbose>=2 || 0 ) {
        if ( ascending ) {
            if ( ! is.issorted() ) {
                myprintf ( "symmetry_group: input group was not sorted!\n" );
                is.show();
                myprintf ( "\n" );
            }
        } else {
            if ( ! is.issorted() ) {
                myprintf ( "symmetry_group: input group was not sorted!\n" );
                is.show();
                myprintf ( "\n" );
            }
        }
    }
    // calc group
    int nsg=0;
    Type prev; //=vals[0];

    prev = std::numeric_limits<Type>::quiet_NaN();
    /* count number of symmetry groups */
    for ( int i=0; i<n; i++ ) {
        if ( vals[i]!=prev || i==0 ) {
            nsg++;
#ifdef FULLPACKAGE
            if ( verbose ) {
                myprintf ( "  symmetry_group: %d: add ", i );
                std::cout << prev << "->";
                std::cout << vals[i] << std::endl;
            }
#endif
        }
        prev=vals[i];
    }

    if ( verbose ) {
        myprintf ( "symmetry_group: %d elements, %d groups\n", n, nsg );

    }

    this->ngroups=nsg;

// initialize structures
    gidx.resize ( n );
    gstart.resize ( nsg+1 );
    gsize.resize ( nsg+1 );


    /* find starting positions of symmetry group */
    if ( verbose>=2 )
        myprintf ( "symmetry_group::init:  find starting positions\n" );

    prev =  std::numeric_limits<Type>::quiet_NaN();
    nsg = 0;
    int i;
    for ( i=0; i<n; i++ ) {
        //printf("symm_group_index: nsg: %d i %d\n", nsg, i);
        if ( vals[i]!=prev || i==0 ) {
            gstart[nsg]=i;
            nsg++;
        }
        gidx[i]=nsg-1;
        prev=vals[i];
    }
    gstart[nsg]=n;	/* add dummy element */

    for ( int i=0; i<nsg; i++ ) {
        gsize[i] = gstart[i+1]-gstart[i];
    }
    gsize[nsg]=0;
}

symmetry_group::symmetry_group ( const std::vector<float> &vals, bool ascending, int verbose )
{
    this->init<float> ( vals, ascending, verbose );
}

symmetry_group::symmetry_group ( const std::vector<double> &vals, bool ascending, int verbose )
{
    this->init<double> ( vals, ascending, verbose );
}

symmetry_group::symmetry_group ( const std::vector<short int> &vals, bool ascending, int verbose )
{
    this->init ( vals, ascending, verbose );
}

symmetry_group::symmetry_group ( const std::vector<unsigned int> &vals, bool ascending, int verbose )
{
    this->init ( vals, ascending, verbose );
}

symmetry_group::symmetry_group ( const std::vector<int> &vals, bool ascending, int verbose )
{
    this->init<int> ( vals, ascending, verbose );
}

symmetry_group::symmetry_group ( const std::vector< mvalue_t<double> > &vals, bool ascending, int verbose )
{
    this->init ( vals, ascending, verbose );
}

symmetry_group::symmetry_group ( const std::vector<mvalue_t<int> > &vals, bool ascending, int verbose )
{
    //printf("symmetry_group::symmetry_group: type <mvalue_t<int>: %zu, %zu\n", vals.size(), vals[0].size() );
    this->init ( vals, ascending, verbose );
}

symmetry_group& symmetry_group::operator= ( const symmetry_group &sgx )
{
    gidx = sgx.gidx;
    gstart = sgx.gstart;
    gsize = sgx.gsize;
    ngroups = sgx.ngroups;
    n = sgx.n;
    ascending = sgx.ascending;

    return *this;
}

symmetry_group::symmetry_group ( const symmetry_group &sgx )
{
    gidx = sgx.gidx;
    gstart = sgx.gstart;
    gsize = sgx.gsize;
    ngroups = sgx.ngroups;
    n = sgx.n;
    ascending = sgx.ascending;
}
symmetry_group::symmetry_group ( )
{
    ngroups=0;
    n=0;
    ascending=0;
}


void symmetry_group_walker::show ( int verbose ) const
{
    myprintf ( "symmetry_group_walker: " );
    if ( verbose>=2 )
        myprintf ( "\n" );
    for ( size_t i=0; i< ( size_t ) perms.size(); i++ ) {
        if ( verbose>=2 ) {
            myprintf ( "  block %ld: ", i );
            print_perm ( perms[i] );
        } else {
            print_perm ( perms[i], 100, false );
            myprintf ( " " );
        }
    }
    if ( verbose==1 )
        myprintf ( "\n" );
}

bool symmetry_group_walker::nextsub ( int g )
{
//myprintf("symmetry_group_walker::nextsub: %d\n", g);
//print_perm(perms[g]);
    bool of = next_perm ( perms[g] );
    if ( of && g>0 )
        of = nextsub ( g-1 );

    return of;

}

std::vector<int> symmetry_group_walker::fullperm() const
{
    std::vector<int> w ( sg.n );
    for ( int i=0; i<sg.n; i++ )
        w[i]=i;

    std::vector<int> ww ( sg.n );
    for ( size_t j=0; j<perms.size(); j++ ) {
        for ( size_t k=0; k<perms[j].size(); k++ ) {
            int offset=sg.gstart[j];
            ww[offset+k] = w[offset+perms[j][k]];
        }
    }
    return ww;
}

// instantiate classes
//template symmetry_group::symmetry_group<int>(std::vector<int>, int);

/* Random number generators */

int g_seed = 123;
void seedfastrand ( int s )
{
    g_seed=s;
}
int fastrand()
{
    g_seed = ( 214013*g_seed+2531011 );
    return ( g_seed>>16 ) &0x7FFF;
}
int fastrandK ( int K )
{
    return fastrand() % K;
}

// cached data

long **ncombsdata = 0;
int ncombscachemax=0;

int ncombscacheNumber()
{
    return 	ncombscachemax;
}
void initncombscache(int N)
{


    if(N<=ncombscacheNumber() )
        return;
#ifdef OADEBUG
    myprintf("initncombscache: value %d\n", N);
#endif

    #pragma omp critical
    {
        const int rowsize=N+1;
        const int nrows=N+1;
        if(ncombsdata!=0) {
            delete [] ncombsdata[0];
            delete [] ncombsdata;
        }

        ncombsdata = new long* [nrows];
        // if debugging, check for memory allocation
        assert(ncombsdata);

        //myprintf("nrows*rowsize: %d * %d = %d\n", nrows, rowsize, nrows*rowsize);
        ncombsdata[0] = new long [nrows*rowsize];

        int offset = 0;
        for ( int i = 0; i < nrows; i++ ) {
            ncombsdata[i] = ncombsdata[0]+offset;
            offset += rowsize;
        }

        for(int i=0; i<nrows; i++) {
            for(int j=0; j<rowsize; j++) {
                //myprintf("i j %d %d, %d\n", i, j, N);
                ncombsdata[i][j]=ncombs(i,j);
            }

            ncombscachemax=N;
        }
    }

}
long ncombscache(int n, int k)
{
#ifdef OADEBUG
    assert(ncombsdata!=0);
    //assert(n<=ncombscacheValue() );
    //assert(k<=ncombscacheValue() );
#endif
    //printf("ncombscache: %d %d\n", n, k);
    return ncombsdata[n][k];
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
