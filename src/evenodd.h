/*! \file evenodd.h
 *  \brief Contains functions to generate even-odd designs
 *
 */

#pragma once

#ifdef OADEBUG
#else
//#define DOOPENMP
#endif

#include <algorithm>
#include <map>

#include "arraytools.h"
#include "lmc.h"
#include "extend.h"


/// structure containing current position in search tree
struct depth_path_t {
    std::vector < int >ncurr;	// vector with current position
    std::vector < int >nmax;	// vector with target
    std::vector < int >necols;	// number of extension columns
    std::vector < int >ngecols;	// number of good extension columns
    int depthstart;

    depth_path_t () {
    }

    void updatePositionGEC ( int k, int goodextensioncols ) {
        ngecols[k] = goodextensioncols;

    }
    void updatePosition ( int k, int c, int m, int extensioncols,
                          int goodextensioncols ) {
        ncurr[k] = c;
        nmax[k] = m;
        necols[k] = extensioncols;
        ngecols[k] = goodextensioncols;

    }
    void show ( int depth, int maxentries = 8 ) const {
        for ( int i = depthstart; i <= depth; i++ ) {

            if ( ( i - depthstart ) == maxentries ) {
                printf ( "... " );
                break;
            }

            printf ( "%d: %d/%d (%d->%d) ", i, ncurr[i], nmax[i], necols[i],
                     ngecols[i] );
        }
        printf ( "\n" );

    }
    void init ( int ncols, int _depthstart = 9 ) {
        ncurr.resize ( ncols + 1 );
        nmax.resize ( ncols + 1 );
        necols.resize ( ncols + 1 );
        ngecols.resize ( ncols + 1 );
        depthstart = _depthstart;
    }
};

/// structure to count and show number of arrays generated, the structure is thread safe
struct counter_t {

    std::vector < int >nfound;	// vector with number of arrays found

    counter_t ( int n ) {
        nfound.resize ( n + 1 );
    }

    void addNfound ( int col, int num ) {
#ifdef DOOPENMP
        #pragma omp atomic
#endif
        this->nfound[col] += num;
    }

    long nArrays () const {
        long na =
            std::accumulate ( this->nfound.begin (), this->nfound.end (), 0 );;
        return na;
    }
    void addNumberFound ( int n, int k ) {
#ifdef DOOMP
        #pragma omp critical (DEXTEND_NFOUND)
#endif
        {
            this->nfound[k] += n;
        }
    }

    void clearNumberFound () {
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            for ( size_t k = 0; k < this->nfound.size (); k++ ) {
                this->nfound[k] = 0;
            }
        }
    }

    void addNumberFound ( const counter_t & de ) {
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            for ( size_t k = 0; k < this->nfound.size (); k++ ) {
                this->nfound[k] += de.nfound[k];
            }
        }
    }

    /// show information about the number of arrays found
    inline void showcountscompact () const {
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            printf ( "depth_extend: counts " );
            display_vector ( this->nfound );
            printf ( "\n" );
        }
    }

    /// show information about the number of arrays found
    inline void showcounts ( const arraydata_t & ad ) const {
        printf ( "--results--\n" );
        for ( size_t i = ad.strength; i <= ( size_t ) ad.ncols; i++ ) {
            printf ( "depth_extend: column %ld: found %d\n", i, this->nfound[i] );
        }
    }

    /// show information about the number of arrays found
    inline void showcounts ( const char *str, int first, int last ) const {
        printf ( "--results--\n" );
        for ( size_t i = first; i <= ( size_t ) last; i++ ) {
            printf ( "%s: column %ld: found %d\n", str, i, this->nfound[i] );
        }
    }
};



/** Helper structure for dynamic extension
 *
 * In this structure we keep track of pointers to valid column extensions
 *
 */
struct depth_extend_sub_t {
public:
    std::vector < int >lmctype;
    /// last column changed in lmc check
    std::vector < int >lastcol;

    std::vector < double >strengthcheck;
    std::vector < int >valididx;
    // param
    int verbose;


    depth_extend_sub_t ( int nn = 0 ) :verbose ( 0 ) {
        resize ( nn );
    };

    void resize ( int nn ) {
        this->lmctype.resize ( nn );
        this->lastcol.resize ( nn );
        this->strengthcheck.resize ( nn );
        // this->tmp.resize(nn);

        std::fill ( this->lmctype.begin (), this->lmctype.begin () + nn, LMC_MORE );
        std::fill ( this->lastcol.begin (), this->lastcol.begin () + nn, -1 );
    }


    inline size_t n () const {
        return lmctype.size ();
    }

    std::vector < int > updateExtensionPointers ( int extcol ) {
        if ( verbose >= 3 )
            printf
            ( "updateExtensionPointers: determine extensions that can be used at the next stage\n" );

        std::vector < int >pointers;
        //
        for ( size_t i = 0; i < lmctype.size (); i++ ) {
            // we need proper strength
            if ( strengthcheck[i] ) {
                //printf("##  depth_extend_sub_t.updateExtensionPointers: i %zu, lastcol %d extcol %d \n", i, lastcol[i], extcol);
                if ( lastcol[i] >= extcol || lastcol[i] == -1 || extcol < 5 ) {
                    // NOTE: extcol < 5 condition --> make generic
                    // good candidate
                    //      printf("  %d %d \n", lastcol[i], extcol);
                    pointers.push_back ( valididx[i] );
                }
            }
        }
        if ( verbose >= 2 )
            printf ( "updateExtensionPointers: extcol %d, kept %ld/%ld pointers\n",
                     extcol, pointers.size (), lmctype.size () );
        return pointers;
    }

    /// initialize the new list of extension columns
    arraylist_t initialize ( const arraylist_t & alist, const arraydata_t & adf,
                             const OAextend & oaextend );

    /// select the arrays with are LMC and hence need to be written to disk
    inline arraylist_t selectArraysZ ( const arraylist_t & alist ) const {
        if ( verbose >= 2 )
            printf
            ( "depth_extend_sub_t: selectArrays: alist.size() %ld, lmctype %ld\n",
              alist.size (), lmctype.size () );
        arraylist_t ga;
        for ( size_t i = 0; i < lmctype.size (); i++ ) {
            //size_t ii = valididx[i];
            if ( verbose >= 3 )
                printf
                ( "  depth_extend_sub_t.selectArraysZ: array %ld: lmctype %d\n", i,
                  lmctype[i] );
            if ( lmctype[i] >= LMC_EQUAL ) {
                array_link ee = alist[i];

                ga.push_back ( ee );
                // printf ( "  selectArraysZ: selected array %zu/%zu\n", i, lmctype.size() );
            }
        }
        if ( verbose )
            printf ( "dextend_sub_t: selected %d/%d arrays\n", ( int ) ga.size (),
                     ( int ) alist.size () );
        return ga;
    }

    inline arraylist_t selectArraysXX ( const array_link & al,
                                        const arraylist_t & elist ) const {
        if ( verbose >= 2 )
            printf
            ( "depth_extend_sub_t: selectArraysXX: alist.size() %ld, lmctype %ld\n",
              elist.size (), lmctype.size () );
        arraylist_t ga;
        for ( size_t i = 0; i < n (); i++ ) {
            if ( verbose >= 3 ) {
                printf ( "  selectArraysXX lmctype %d\n", lmctype[i] );
            }
            if ( lmctype[i] >= LMC_EQUAL ) {
                array_link ee = hstacklastcol ( al, elist[valididx[i]] );

                ga.push_back ( ee );
            }
        }
        if ( verbose >= 1 )
            printf ( "dextend_sub_t: selected %d/%d arrays\n", ( int ) ga.size (),
                     ( int ) elist.size () );
        return ga;
    }

    inline void info () const {
        size_t nl = 0;
        for ( size_t t = 0; t < lmctype.size (); t++ ) {
            nl += ( lmctype[t] > +LMC_EQUAL );
        }
        //int ngood =std::accumulate(lmctype.begin(),lmctype.end(),0);//#include <numeric>
        if ( verbose ) {
            printf ( "lmc %ld/%d\n", nl, ( int ) lmctype.size () );
            printf ( "valididx size %ld\n", valididx.size () );
        }
    }
};

/** @brief Helper structure for dynamic extension
 *
 * This structure allows for writing the generated arrays to disk.
 * It also contains functions to print progress of the extension.
 *
 * Multiple copies of this class are made, but they all share the same counter_t and arraywriter_t object. Also t0 and tp are shared
 *
 */
struct depth_extend_t {
public:
    int verbose;
    OAextend oaextend;
    const arraydata_t *ad;

    int loglevelcol;
    /// print progress every x seconds
    double logtime;

    arraylist_t extension_column_list;	// list of possible extensions


    /// if set to true write arrays to disk
    int writearrays;

    int discardJ5;		/// if true, then we discard the designs which have J5 maximal
    long discardJ5number;

    // shared by mutiple instances of dextend_t (could be made static)
    arraywriter_t *arraywriter;
    counter_t *counter;

    static double t0;		// time since start of calculation
    static double tp;		// time since last progress report

private:
    long narraysmax;
    depth_path_t searchpath;


public:

    // constructure function
    depth_extend_t ( const arraydata_t * ad_, double _logtime = 10000000, int _discardJ5 = -1 ) :verbose ( 1 ), ad ( ad_ ),
        discardJ5
        ( _discardJ5 ) {
        loglevelcol = -1;
        t0 = get_time_ms ();
        tp = get_time_ms ();

        logtime = _logtime;
        if ( ad == 0 ) {
            printf ( "depth_extend_t: pointer to arraydata_t is zero!" );
        }

        writearrays = 1;
        narraysmax = LONG_MAX - 1;

        arraywriter = 0;
        counter = 0;

        discardJ5number = 0;

        searchpath.init ( ad->ncols );
    };

    depth_extend_t ( const depth_extend_t & de ) {
        //printf("depth_extend_t: copy constructor\n"); printf(" searchpath: "); de.searchpath.show(16);
        verbose = de.verbose;
        oaextend = de.oaextend;
        ad = de.ad;
        loglevelcol = de.loglevelcol;
        logtime = de.logtime;
        //cache_extensions=de.cache_extensions;
        extension_column_list = de.extension_column_list;
        arraywriter = de.arraywriter;
        writearrays = de.writearrays;
        narraysmax = de.narraysmax;
        searchpath = de.searchpath;
        discardJ5 = de.discardJ5;
        discardJ5number = de.discardJ5number;

        counter = de.counter;
    }
    ~depth_extend_t () {
        //closeafiles();
    }


public:

    void show () {
        printf ( "depth_extend_t: logtime %.1f [s]\n", logtime );
    }

    void setNarraysMax ( long n ) {
        this->narraysmax = n;
    }

    // helper function, thread safe
    void maxArrayCheck () {
        if ( arraywriter == 0 ) {
            return;
        }
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            if ( arraywriter->nwritten > this->narraysmax ) {
                /// HACK
                printf ( "dextend_t: number of arrays written: %d, quitting\n",
                         arraywriter->nwritten );
                this->counter->showcounts ( *this->ad );
                this->arraywriter->closeafiles ();
                //printfd ( "symmetry time: %.1f [s]\n", getdtsymm() );
                //printfd ( "lmc time: %.3f [s]\n", getdtlmc() );

                exit ( 0 );
            }
        }
    }


public:

    void showsearchpath ( int depth ) const {
        searchpath.show ( depth );
    }

    /// show information about the progress of the loop
    bool showprogress ( int showtime = 1, int depth = 0, int forcelog = 0 ) {
        {
            double currenttime = get_time_ms ();
            double dt = currenttime - tp;
            if ( ( dt > logtime ) || forcelog ) {
#ifdef DOOPENMP
                #pragma omp critical
#endif
                tp = get_time_ms ();
                this->arraywriter->flush ();
                double dt0 = currenttime - t0;
                if ( showtime ) {
                    int na = this->counter->nArrays ();

#ifdef DOOPENMP
                    printf
                    ( "-- depth_extend: progress: %.1f [s], narrays %d (%.1f arrays/s), thread %d/%d\n",
                      dt0, na, na / dt0, omp_get_thread_num (),
                      omp_get_num_threads () );
#else
                    printf
                    ( "-- depth_extend: progress: %.1f [s], narrays %d (%.1f arrays/s)\n",
                      dt0, na, na / dt0 );
#endif

                    if ( depth > 0 ) {
                        printf ( "-- depth %2d: %.1f [s]: ", depth, dt0 );
                        searchpath.show ( depth );
                    }
                }
                return true;
            } else {
                return false;
            }
        }
    }
    inline void info () const {
        printf ( "depth_extend: " );
        ad->show ();
    }


    /// set the position in the dextend structure
    void setposition ( int k, int c, int m, int extensioncols =
                           -1, int goodextensioncols = -1 ) {
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            searchpath.updatePosition ( k, c, m, extensioncols, goodextensioncols );
        }
    }
    /// set the position in the dextend structure
    void setpositionGEC ( int k, int goodextensioncols ) {
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            searchpath.updatePositionGEC ( k, goodextensioncols );
        }
    }
};

enum depth_alg_t
{ DEPTH_DIRECT, DEPTH_EXTENSIONS };


/// Helper structure for the even-odd depth extension
struct depth_extensions_storage_t {

    void resize ( size_t s ) {
        columnextensionsList.resize ( s );
        goodarrayslist.resize ( s );
        depthalglist.resize ( s );
        dextendsubList.resize ( s );
    }

    void set ( int ai, const arraylist_t & goodarrays,
               const arraylist_t & extension_column_list, depth_alg_t depthalg,
               const depth_extend_sub_t & dextendsub ) {
        this->goodarrayslist[ai] = ( goodarrays );
        this->columnextensionsList[ai] = ( extension_column_list );
        this->depthalglist[ai] = ( depthalg );
        this->dextendsubList[ai] = ( dextendsub );

    }
    //std::vector<arraylist_t>  extensions0list;
    std::vector < arraylist_t > columnextensionsList;
    std::vector < arraylist_t > goodarrayslist;
    std::vector < depth_alg_t > depthalglist;
    std::vector < depth_extend_sub_t > dextendsubList;
};

/** Extend arrays using a depth-first or breadth-first approach
 *
 * @param goodarrays List of arrays to extend
 * @param depthalg Extend using depth-first or breadth-first
 * @param dextend Option structure for the extension
 * @param dextendsublight Data structure for the extensions
 * @param extensioncol Column to extend
 * @param verbose Verbosity level
 */
void processDepth ( const arraylist_t & goodarrays, depth_alg_t depthalg,
                    depth_extend_t & dextend,
                    depth_extend_sub_t & dextendsublight, int extensioncol,
                    int verbose = 0 );

/// depth-first extension of arrays. depending on the symmetry group of the array to be extended a direct method is used or a method with caching of candidate columns
void depth_extend_hybrid ( const arraylist_t & alist, depth_extend_t & dextend,
                           int extcol, const OAextend & oaextendx,
                           int verbose );

/// variation of depth_extend for arrays with large symmetry groups
void depth_extend_direct ( const arraylist_t & alist, depth_extend_t & dextend,
                           int extcol, const OAextend & oaextendx,
                           int verbose );




/** @brief perform depth-first extension
 *
 * The arrays generated are pruned by keeping a list of possible extension values
 *
 */
//void depth_extend ( const arraylist_t &alist,  depth_extend_t &dextend, const depth_extend_sub_t &dextendsub, int col, int verbose=1 );


/// depth extend a single array
void depth_extend_array ( const array_link & al, depth_extend_t & dextend,
                          const arraydata_t & adfull, int verbose,
                          depth_extensions_storage_t * ds = 0, int = 0 );


/// callback function for Pareto calculations
typedef Pareto < mvalue_t < long >,
        array_link >::pValue ( *pareto_cb ) ( const array_link &, int );
/// callback function for Pareto calculations with cache
typedef Pareto < mvalue_t < long >,
        array_link >::pValue ( *pareto_cb_cache ) ( const array_link &, int,
                rankStructure & rs );

template < class IndexType >
inline typename Pareto < mvalue_t < long >, IndexType >::pValue
calculateArrayParetoJ5Cache ( const array_link & al, int verbose,
                              rankStructure & rs )
{

    const int N = al.n_rows;
    int r = rs.rankxf ( al );
    mvalue_t < long >wm = A3A4 ( al );
    mvalue_t < long >f4 = F4 ( al );

    typename Pareto < mvalue_t < long >, IndexType >::pValue p;
    p.push_back ( r );		// rank of second order interaction matrix
    p.push_back ( wm );		// A4
    p.push_back ( f4 );		// F
    addJmax < IndexType > ( al, p, verbose );

    if (0) {
        printf ( "calculateArrayParetoJ5Cache: %d ; ", r );
        wm.show_integer();
        printf ( " ; " );
        f4.show_integer();
        printf ( " ; " );
        p[3].show_integer();
        printf ( " ; " );
        p[4].show_integer();
        printf ( "\n" );
    }
    return p;
}


/// add arrays to set of Pareto results
void addArraysToPareto ( Pareto<mvalue_t<long>,array_link> &pset, pareto_cb paretofunction, const arraylist_t & arraylist, int jj, int verbose );
/// add arrays to set of Pareto results
void addArraysToPareto ( Pareto < mvalue_t < long >, array_link > &pset, pareto_cb_cache paretofunction, const arraylist_t & arraylist, int jj, int verbose );

/** helper class for indexing statistics of designs
 *
 * The index consists of the number of columns and the value for the J-characteristic
 */
struct jindex_t {
    int k;			// number of columns
    int j;			// J-value

    jindex_t ( int colindex, int jvalue ) :k ( colindex ), j ( jvalue ) {
    }

public:

    bool operator< ( const jindex_t & rhs ) const {
        if ( this->k < rhs.k ) {
            return true;
        }
        if ( this->k > rhs.k ) {
            return false;
        }

        return ( this->j < rhs.j );
    }
    std::string toString () const {
        return printfstring ( "k%d-j%d", this->k, this->j );
    }
};


/// object to hold counts of maximum J_k-values
class Jcounter
{
public:
    int N;			/// number of rows
    int jj;
    std::vector < int >fvals;
    std::map < jindex_t, long >maxJcounts;
    double dt;			/// time needed for calculation

    Jcounter () :N ( -1 ), jj ( -1 ) {
    }

    Jcounter ( int N, int jj = 5, int k = -1 ) {
        this->init ( N, jj, k );

    }

    bool validData() {
     if (N==-1 && jj==-1) 
         return false;
     else return true;
    }
    
    /// return true if specified column is in the data
    bool hasColumn ( int col ) const {
        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {

            if ( it->first.k == col ) {
                return true;
            }
        }
        return false;
    }

    bool isOpen () const {
        return N > 0;
    }
    void showPerformance () const {
        myprintf ( "Jcounter: %.1f Marrays/h\n",
                   ( 1e-6 * 3600. ) * double ( narrays () ) / this->dt );
    }
    long narrays () const {

        long r = 0;
        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            r += it->second;
        }

        return r;
    }

    /// show statistics of the object
    void show () const {

        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            myprintf ( "k %d: max(J%d) %d: %ld\n", it->first.k, this->jj,
                       it->first.j, it->second );
        }
    }

    int maxCols() const {
        int kmax=-1;
        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            kmax=std::max ( kmax,  it->first.k );
        }

        return kmax;
    }
    
    long getCount(int k, int j) const {
        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            if( it->first.j==j && it->first.k==k) {
                return it->second;
            }
        }
            return -1;
    }
    
    std::vector<long> getTotalsJvalue(int jval) const {
        int nmax=maxCols();
        std::vector<long> k ( nmax+1 );

        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            if ( it->second<0 ) {
                printf ( "Jcounter::getTotals: value -1 for index %s\n", it->first.toString().c_str() );
            } else {
                if( it->first.j==jval)
                    k[it->first.k] += it->second;
            }
        }
        return k;
    }
    std::vector<long> getTotals() const {
        int nmax=maxCols();
        std::vector<long> k ( nmax+1 );

        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            if ( it->second<0 ) {
                printf ( "Jcounter::getTotals: value -1 for index %s\n", it->first.toString().c_str() );
            } else {
                k[it->first.k] += it->second;
            }
        }
        return k;
    }
    /// show statistics of the object
    void showcompact () const {

        int kprev = -1;
        long nt = 0;
        for ( std::map < jindex_t, long >::const_iterator it = maxJcounts.begin ();
                it != maxJcounts.end (); ++it ) {
            if ( it->first.k == kprev ) {
                nt += it->second;
                myprintf ( "; %d: %ld", it->first.j, it->second );
            } else {
                if ( kprev != -1 ) {
                    myprintf ( "; total %ld\n", nt );
                }
                nt = 0;
                myprintf ( "k %d: max(J%d) %d: %ld", it->first.k, this->jj,
                           it->first.j, it->second );
                kprev = it->first.k;
            }
        }
        myprintf ( "; total %ld\n", nt );
    }

    Jcounter & operator += ( Jcounter & jc );


    /// add list of arrays to object
    void addArrays ( const arraylist_t & arraylist, int verbose = 0 );

    /// add single array to statistics object
    void addArray ( const array_link & al, int verbose = 0 ) {
        //jstruct_t js ( al, this->jj );
        jstruct_t js ( al.selectFirstColumns ( 5 ), this->jj );

        int maxJ = js.maxJ ();

        int k = al.n_columns;

        if ( verbose ) {
            jstruct_t js ( al, this->jj );
            std::vector < int >FF = js.calculateF ();
            printf ( "addArray: maxJ %d: ", maxJ );
            display_vector ( FF );
            printf ( "\n" );
        }
        jindex_t ji = jindex_t ( k, maxJ );
        #pragma omp critical
        maxJcounts[ji]++;
    }

private:
    void init ( int N, int jj, int k = -1 ) {
        this->N = N;
        this->jj = jj;
        this->fvals = Fval ( N, 3 );
        this->dt = 0;

        maxJcounts.clear ();

        if ( k > 0 ) {
            for ( size_t j = 0; j < fvals.size (); j++ ) {
                jindex_t ji ( k, fvals[j] );
                //printf("adding val %s\n", ji.toString().c_str());
                maxJcounts[ji] = 0;
            }
        }
    }

};

/// read statistics object from disk
Jcounter readStatisticsFile ( const char *numbersfile, int verbose );
/// write statistics object to disk
void writeStatisticsFile ( const char *numbersfile, const Jcounter & jc, int verbose );


/// calculate J-value statistics
Jcounter calculateJstatistics ( const char *afile, int jj = 5, int verbose =
                                    1 );
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
