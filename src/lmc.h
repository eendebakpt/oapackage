/** \file lmc.h

\brief This file contains definitions and functions to perform LMC tests and reductions

 C++ Interface: lmc

 This file contains definitions are functions to perform LMC tests and reductions

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2008

 Copyright: See COPYING file that comes with this distribution
*/

/*! \mainpage	Orthogonal Arrays

This package contains code to calculate orthogonal arrays with a specified number of runs, factor levels and strength. The main program is oaextendsingle (or oaextendmpi for the multi-core version), which starts with the specifications of an OA and creates the root of th arrray. Then the array is extended column wise.
For a more complete description of the algoritm see the article
"Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays", E.D. Schoen, P.T. Eendebak, M.V.M. Nguyen.

The programs in this package are:

- oaextendsingle/oaextendmpi: Extend LMC arrays with additional columns
- oainfo: Return information about files containing arrays
- oasplit: Split array files into multiple files
- oajoin: Join multiple files
- oacheck: Check a file with arrays using the LMC check. Optionally the arrays are reduced to LMC form
- oafilter: Filter a file with arrays using a binary file with indices
- oaanalyse: Calculate statistics for each array in an array file

For more information please contact Pieter Eendebak, <pieter.eendebak@gmail.com>.

  @sa extend_array(), LMCreduce(), array_link, arraydata_t
  */

#ifndef LMC_H
#define LMC_H

#include <list>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <set>
#include <vector>

#include "mathtools.h"
#include "tools.h"
#include "arrayproperties.h"

#ifdef SWIG
#else
#define ORDER_J5_SMALLER >
#define ORDER_J5_GREATER <
#define ORDER_J45_SMALLER >
#define ORDER_J45_GREATER <
#endif

// forward declaration
class OAextend;

#define stringify( name ) # name

#ifdef SWIG
%ignore check_root_update;
%ignore LMCreduceFull;
%ignore dyndata_t::dyndata_t(dyndata_t const &);
#endif


#ifdef LMCSTATS
// function is not thread safe
void lmc_stats();
#endif

/* constants and structures */

/** Possible results for the LMC check
 *
 * LMC_LESS: Found a permutation which leads to a lexicographically smaller array
 * LMC_EQUAL: Found a permutation which leads to a lexicographically equal array
 * LMC_MORE: Found a permutation which leads to a lexicographically larger array
 */
enum lmc_t {LMC_LESS, LMC_EQUAL, LMC_MORE, LMC_NONSENSE};

/// different algorithms
enum algorithm_t { MODE_ORIGINAL, MODE_J4, MODE_J5ORDER, MODE_J5ORDERX, MODE_INVALID, MODE_AUTOSELECT, MODE_LMC_SYMMETRY, MODE_LMC_2LEVEL, MODE_LMC_DEBUG, MODE_J5ORDERXFAST};
#define MODE_LMC MODE_ORIGINAL

inline std::string algorithm_t_list()
{
    std::string ss = printfstring("%d (automatic), %d (original), %d (check j4), %d (j5 order), %d (j5 order dominant), %d (MODE_J5ORDERXFAST)", MODE_AUTOSELECT,
                                  MODE_ORIGINAL, MODE_J4, MODE_J5ORDER, MODE_J5ORDERX, MODE_J5ORDERXFAST);
    ss+= printfstring(", %d (MODE_LMC_SYMMETRY), %d (MODE_LMC_2LEVEL)", MODE_LMC_SYMMETRY, MODE_LMC_2LEVEL);
    return ss;
}

/// method used for initialization of columns
enum initcolumn_t {INITCOLUMN_ZERO, INITCOLUMN_PREVIOUS, INITCOLUMN_J5};

/// variations of the J45 structures
enum j5structure_t {J5_ORIGINAL, J5_45};

/// return name of the algorithm
std::string algnames ( algorithm_t m );



struct dyndata_t;

// NOTE: unsigned long is enough for 2-factor arrays up to 60 columns
//typedef unsigned char rowsort_value_t; /** type for value for sorting rows*/
//typedef unsigned long rowsort_value_t; /** type for value for sorting rows*/
typedef unsigned int rowsort_value_t; /** type for value for sorting rows*/
/*!
 * @brief structure to perform row sorting
 */
struct rowsort_t {
    //! index of row
    rowindex_t r;
    //! value of row
    rowsort_value_t val;
};

typedef rowindex_t rowsortlight_t;

/**
 * @brief Comparision operator for the rowsort_t structure
 */
static inline bool operator< ( const rowsort_t& a, const rowsort_t& b )
{
#ifdef OAOVERFLOW
    assert(a.val>=0);
    assert(b.val>=0);
#endif
    return a.val < b.val;
}

/**
 * @brief Comparision operator for the rowsort_t structure
 */
static inline bool operator> ( const rowsort_t& a, const rowsort_t& b )
{
#ifdef OAOVERFLOW
    assert(a.val>=0);
    assert(b.val>=0);
#endif

    return a.val > b.val;
}


/// Apply Hadamard transformation to orthogonal array
void apply_hadamard ( const arraydata_t *ad, array_t *array, colindex_t hcol );

/// Apply Hadamard transformation to orthogonal array
void apply_hadamard ( array_link &al, colindex_t hcol );



/**
 * @brief Contains initialization data for static allocations
 *
 *  Part of the allocations is for structures that are constant and are re-used each time an LMC calculation is performed.
 *  Some other structures are temporary buffers that are written to all the time.
 */
struct LMC_static_struct_t {
public:

    /* variable data */
    int LMC_non_root_init;
    int LMC_root_init;
    //int LMC_reduce_non_root_init;
    //int LMC_reduce_root_init;
    int LMC_reduce_root_rowperms_init;

    /* data structures used by functions */
    arraydata_t *ad;

    /* used at the root_level_perm stage */
    int LMC_root_rowperms_init;
    int nrootrowperms;
    rowperm_t* rootrowperms;

    int LMC_root_rowperms_init_full;
    int nrootrowperms_full;
    rowperm_t* rootrowperms_full;


    array_t* colbuffer;	/** buffer for a single column */

    dyndata_t ** dyndata_p;	/** dynamic data; row permutations */
    colindex_t ** colperm_p;	/** column permutations */
    colindex_t ** localcolperm_p;	/** local column permutation */

    array_transformation_t *current_trans; /* not used at the moment? */

#ifdef OADEBUG
    std::string ref;
    inline void setRef(const std::string s)
    {
        ref=s;
        printf("LMC_static_struct_t: set ref %s\n", ref.c_str() );
    }
    int id;
#endif

    /* methods */
    LMC_static_struct_t();
    ~LMC_static_struct_t();

    void show(int verbose=1) const
    {
        printf("LMC_static_struct_t: ad %ld, LMC_non_root_init %d\n", long(this->ad), LMC_non_root_init );
    }

    void init ( const arraydata_t *adp );
    void freeall ( );

    /// update structure with new design specification
    int update ( const arraydata_t *adp );
    int needUpdate ( const arraydata_t *adp ) const;

    void init_root_stage ( levelperm_t * &lperm_p, colperm_t * &colperm_p, const arraydata_t *adp );
    void init_nonroot_stage ( levelperm_t * &lperm_p, colperm_t * &colperm_p, colperm_t * &localcolperm_p, dyndata_t ** &dynd_p, int &dynd_p_nelem,  array_t * &colbuffer, const arraydata_t *adp ) const;


/// Static initialization of root row permutations
    void init_rootrowperms ( int &totalperms, rowperm_t * &rootrowperms, levelperm_t * &lperm_p )
    {
        /* no static update, we assume this has been done already */

        totalperms = this->nrootrowperms;
        rootrowperms = this->rootrowperms;
        lperm_p =  this->current_trans->lperms;

        this->LMC_root_rowperms_init = 1; // FIXME: is this needed?
    }

    /** @brief Static initialization of root row permutations (full group)
     */
    void init_rootrowperms_full ( int &totalperms, rowperm_t * &rootrowperms, levelperm_t * &lperm_p )
    {
        /* no static update, we assume this has been done already */

        totalperms = this->nrootrowperms_full;
        rootrowperms = this->rootrowperms_full;
        lperm_p =  this->current_trans->lperms;

        this->LMC_root_rowperms_init_full = 1;
    }

};


// pool of available structures
LMC_static_struct_t * getGlobalStaticIndexed(int n);
void cleanGlobalStaticIndexed();

/// return static structure from dynamic global pool, return with releaseGlobalStatic
LMC_static_struct_t * getGlobalStatic();
void releaseGlobalStatic(LMC_static_struct_t *p);
void cleanGlobalStatic();
#ifdef OADEBUG
int getGlobalStaticNumber(LMC_static_struct_t * p);
#endif
LMC_static_struct_t &getGlobalStaticOne();


/// variable indicating the state of the reduction process
enum REDUCTION_STATE {REDUCTION_INITIAL, REDUCTION_CHANGED};
//! main mode for the LMC routine: test, reduce or reduce with initialization
enum OA_MODE {OA_TEST, OA_REDUCE, OA_REDUCE_PARTIAL};

typedef larray<rowindex_t> rowpermtypelight;
typedef larray<colindex_t> colpermtypelight;

typedef std::vector<int> colpermtype;
typedef std::vector< colpermtype > colpermset;
//typedef std::set<std::vector<int> > colpermset;

#define USE_ROWPERM_POOL

/// class representing array symmetries
struct arraysymmetry {
    larray<rowindex_t> *rowperm;

    larray<colindex_t> colperm;

public:
    static object_pool< larray<rowindex_t> > rowpermpool;

public:
    ~arraysymmetry();

    arraysymmetry(const dyndata_t *dyndata);

/// copy constructor
    arraysymmetry(const arraysymmetry &rhs);

    //Copy assignment operator
    arraysymmetry &operator= ( const arraysymmetry &rhs );

    inline int operator== (const arraysymmetry &as2)
    {
        return (( (*rowperm) == *(as2.rowperm) ) && (this->colperm==as2.colperm) );
    }


    int N() const
    {
        return rowperm->size();
    }
    void show() const
    {
        printf("arraysymmetry: rowperm ");
        print_perm<rowindex_t>(*rowperm, 24);
        printf("             : colperm ");
        print_perm(colperm);
    }
};





/// helper function
inline array_link rootPlus(const arraydata_t &ad)
{
    array_link al(ad.N, ad.ncols, -1); //al.setvalue(100);
    al.create_root(ad);
    for(int i=ad.strength; i<ad.ncols; i++) {
        for(int r=0; r<ad.N; r++) {
            al.at(r, i)=1000;
        }
    }
    //reduction.setArray(al);
    //reduction.init_state=INIT;
    //printf("LMCcheck: updated reduction.array to\n");
    //reduction.getArray().showarray();
    return al;
}


#ifdef _WIN32
// on windows do not use smart pointers, it is a mess
//#if _MSC_VER >= 1600
#else
#define SDSMART
#endif

#ifdef SDSMART
#ifdef WIN32

#include <memory>
typedef std::shared_ptr<symmdata> symmdataPointer;
//typedef std::shared_ptr<symmdata> symmdataPointer;
#else
#include <tr1/memory>
typedef std::tr1::shared_ptr<symmdata> symmdataPointer;
//typedef std::shared_ptr<symmdata> symmdataPointer;
#endif

#else
typedef symmdata * symmdataPointer;
#endif

//typedef std::vector<std::vector<int> > arraysymmetry;
typedef std::vector< arraysymmetry > symmetryset;

// hack: add documentation here!
enum INIT_STATE {INIT_STATE_INVALID, COPY, INIT, SETROOT};

/// helper function
template <class Type>
void insertUnique(std::vector<Type> &cp, const Type &cpv)
{
    if (cp.size()>0) {
        if(! (cp.back()==cpv) )
            //cp.emplace_back(cpv);
            cp.push_back(cpv);
    } else {
        cp.push_back(cpv);
    }
}

/** @brief Class to describe an LMC reduction
 *
 * The most important variable is the transformation itself, contained in transformation.
 * The state contains information about how the reduction was performed.
 */
struct LMCreduction_t {
    array_t *array;
    /// pointer to transformation_t structure
    array_transformation_t *transformation;
    OA_MODE mode;
    REDUCTION_STATE state;

    INIT_STATE init_state;	// initalization mode: INIT: initialized by user, COPY: copy from array argument

    //! maximum depth for search tree
    int maxdepth;

    /// last column visited in algorithm
    int lastcol;

    //! counter for number of reductions made
    long nred;

    int targetcol;	// NOTE: document this
    int mincol;  // used in debugging

    int nrows, ncols;

    LMC_static_struct_t *staticdata;

    //! store column permutations from array symmetry group
    //std::vector< std::vector< std::vector<int> > > colperms;
    class symm_t
    {
    public:
        int store;
        int ncols;
        std::vector< colpermset > colperms;
        std::vector< colpermset > colcombs;
        std::vector< symmetryset > symmetries;

        symm_t()
        {
            ncols = -1;
        }
        void show(int verbose=1) const
        {
            long ns=0, ncp=0, ncc=0;
            for(size_t i=0; i<symmetries.size(); i++) {
                if (verbose>=2) {
                    printf("  symm_t: k %ld: symms %ld\n", (long)i, (long)symmetries[i].size() );
                }
                ns+=symmetries[i].size();
            }
            for(size_t i=0; i<colperms.size(); i++)
                ncp+=colperms[i].size();
            for(size_t i=0; i<colcombs.size(); i++)
                ncc+=colcombs[i].size();

            printf("symm_t: store %d, symms %ld, colperms %ld, colcombs %ld, ncols %d\n", store, ns, ncp, ncc, ncols);
        }
        bool valid() const
        {
            return ncols>0;
        }
        void makeColpermsUnique(int dverbose=0)
        {
            symm_t &xx = *this;
            // create unique elements
            for(size_t i=0; i<xx.colperms.size(); i++) {
                std::sort(xx.colperms[i].begin(), xx.colperms[i].end());

                std::vector< colpermtype >::iterator last = std::unique(xx.colperms[i].begin(), xx.colperms[i].end());
                if (dverbose>=2) {
                    printf("makeColpermsUnique: i %ld vals %ld, unique %d\n", (long)i, (long) xx.colperms[i].size(), (int)( last - xx.colperms[i].begin()) );
                }
                xx.colperms[i].erase(last,  xx.colperms[i].end());


                for(size_t j=0; j<xx.colperms[i].size(); j++) {
                }
            }

        }

        void storeColumnCombination(colpermtype cpv)
        {
            if (store<2) return;
            int n=cpv.size();
            insertUnique(colcombs[n], cpv);
        }
        void storeColumnCombination(const colperm_t cp, int n)
        {
            if (store<2) return;
            colpermtype cpv = array2vector<int, colindex_t>(cp, n);
            insertUnique(colcombs[n], cpv);
        }
        void storeColumnPermutation(const colperm_t cp, int n)
        {
            if (store<2) return;
            colpermtype cpv = array2vector<int, colindex_t>(cp, n);
            insertUnique(colperms[n], cpv);
        }
        void showColperms(int verbose=1) const
        {
            for(size_t i=0; i<=(size_t)ncols; i++) {
                printf("LMCreduction: column permutations with %d cols: %ld/%ld\n", (int)i, (long)colperms[i].size(), ncombs<long>(ncols, i) );
                if (verbose>=2) {
                    for( colpermset::const_iterator it = colperms[i].begin(); it != colperms[i].end(); it++) {
                        print_perm( *it );
                    }
                }
            }
        }

        void showColcombs(int verbose=1) const
        {
            for(size_t i=0; i<=(size_t)ncols; i++) {
                printf("LMCreduction: column combinations with %d cols: %d/%ld\n", (int)i, (int)colcombs[i].size(), ncombs<long>(ncols, i) );
                if (verbose>=2) {
                    for( colpermset::const_iterator it = colcombs[i].begin(); it != colcombs[i].end(); it++) {
                        print_perm( *it );
                    }
                }
            }
        }
        void showSymmetries(int verbose=1) const
        {
            for(size_t i=0; i<(size_t)symmetries.size() ; i++) {
                printf("LMCreduction: symmetries with %ld cols: %ld/%ld\n", (long)i, (long)symmetries[i].size(), ncombs<long>(ncols, i) );
                if (verbose>=2 || (i==60 )) {
                    for( symmetryset::const_iterator it = symmetries[i].begin(); it != symmetries[i].end(); it++) {
                        it->show();
                    }
                }
            }
        }
        void storeSymmetryPermutation(const dyndata_t *dyndata);


    } symms;

    //symmdata *sd;
    symmdataPointer sd;

public:
    LMCreduction_t ( const LMCreduction_t  &at );	/// copy constructor
    LMCreduction_t ( const arraydata_t *ad );
    ~LMCreduction_t();

    LMCreduction_t & operator= ( const LMCreduction_t &at );	/// Assignment operator

    array_link getArray() const
    {
        array_link al(this->nrows, this->ncols, 0, this->array);
        return al;

    }

    void setArray(const array_link al)
    {
        init_state=INIT;
        std::copy(al.array, al.array+this->ncols*this->nrows, this->array);
    }
    void setArray(const array_t *array, int nrows, int ncols)
    {
        // FIXE: check on sizes
        init_state=INIT;
        std::copy(array, array+ncols*nrows, this->array);

    }


    void updateSDpointer(const array_link al, bool cache=false)
    {
        //reduction.sd = symmdataPointer(new symmdata(al) );
#ifdef SDSMART
        symmdata *sdp = this->sd.get();
#else
        symmdata *sdp = sd;
#endif
        if (sdp!=0 && cache ) { //&& (reduction.sd->orig == al) ) {
            //do nothing
            //  printf("using cached symmdata!\n");
        } else {
            // update symmetry data
            this->sd = symmdataPointer(new symmdata(al, 1) );
        }
    }

    void clearSymmetries();

    void releaseStatic()
    {
        if(this->staticdata!=0) {
            releaseGlobalStatic(this->staticdata);
            this->staticdata=0;
        }
    }

    /// acquire a reference to a LMC_static_struct_t object
    void initStatic()
    {
        if(this->staticdata==0) {
            // printfd("staticdata==0, allocating new structure\n");
            this->staticdata=getGlobalStatic();
        }
    }

    /// return a reference to a LMC_static_struct_t object
    LMC_static_struct_t & getStaticReference()
    {
        if(this->staticdata==0) {
            // printfd("problem! getStaticReference() calls getGlobalStaticOne!\n");
            return getGlobalStaticOne();
        } else {
            //printfd("return internal pointer...\n");
            return *(this->staticdata);
        }
    }

/// reset the reduction: clears the symmetries and sets the transformation to zero
    void reset();

    void show ( int verbose=2 ) const
    {
        printf ( "LMCreduction_t: mode %d, state %d (REDUCTION_INITIAL %d, REDUCTION_CHANGED %d), init_state %d, lastcol %d\n", this->mode, this->state, REDUCTION_INITIAL, REDUCTION_CHANGED, this->init_state, this->lastcol );
        if ( verbose>=1 ) {
            printf ( "LMCreduction_t: nred %ld\n", nred );
            print_array ( "array:\n", this->array, this->transformation->ad->N, this->transformation->ad->ncols );
        }
        if ( verbose>=2 )
            this->transformation->show();
    }

    std::string __repr__() const
    {
        std::string ss = printfstring ( "LMCreduction_t: mode %d, state %d (REDUCTION_INITIAL %d, REDUCTION_CHANGED %d), init_state %d, lastcol %d\n", this->mode, this->state, REDUCTION_INITIAL, REDUCTION_CHANGED, this->init_state, this->lastcol );
        return ss;
    }

    /// called whenever we find a reduction
    void updateFromLoop ( const arraydata_t &ad, const dyndata_t &dynd, levelperm_t *lperms, const array_t *original );
    void updateTransformation ( const arraydata_t &ad, const dyndata_t &dynd, levelperm_t *lperms, const array_t *original );

    inline bool doBreak(lmc_t ret)
    {
        if ( this->mode>=OA_REDUCE )
            return false;
        else {
            return ret==LMC_LESS;
        }
    }
    inline void updateLastCol(int col)
    {
        this->lastcol=col ;
    }

    // non-public

private:
    void free();



};


/** @brief Contains dynamic data of an array
 *
 * The dynamic data are used in the inner loops of the LMC algorithm. In particular
 * they keep track of the current row ordering and column permutation.
 * By not applying these transformations to the array we can save calculation time.
 *
 * We try to prevent copying the object, so it is re-used at different levels in the algorithm.
 *  - N: static
 * 	- col: changes at each column level
 *  - rowsort: changes at each column level, used mainly in non-root stage
 *  - colperm: changes at all levels
 * @sa arraydata_t
 */
struct dyndata_t {
    //! active column
    colindex_t col;
    //! number of rows
    rowindex_t N;
//private:
    //! ordering of rows
    rowsort_t *rowsort;
    rowsortlight_t *rowsortl;
public:
    //! current column permutation
    colperm_t colperm;

    dyndata_t ( int N, int col=0 );
    dyndata_t ( const dyndata_t *dd );
    dyndata_t ( const dyndata_t & );
    ~dyndata_t();


    dyndata_t & operator= ( const dyndata_t & );
    void show() const;

    void reset();
    void setColperm(const colperm_t perm, int n)
    {
        copy_perm ( perm, this->colperm, n );
    }
    void setColperm(const larray<colindex_t> &perm )
    {
        // todo: make check on size
        std::copy(perm.d, perm.d+perm.n, this->colperm);
    }

    void setColperm(const std::vector<colindex_t> &perm )
    {
        // todo: make check on size
        std::copy(perm.begin(), perm.end(), this->colperm);
    }

    /// initialize the rowsort structure from an arraysymmetry object
    void initsymmetry(const arraysymmetry &arraysymm, const symmdata &sd, int ncols)
    {

        const array_t *w = sd.rowvalue.array+(ncols-1)*N;
        for(int i=0; i<N; i++) {
            rowsort[i].r=arraysymm.rowperm->at(i);
            //rowsort[i].val=sd.rowvalue.at(i, ncols-1);
            rowsort[i].val=w[i];

            // TODO: where to set the level permutations? these are needed because we store them finally in the array_transformation_t!
        }

    }


    /// set lightweight row permutation
    void getRowperm(rowpermtypelight &rp) const
    {
        rp.resize(this->N);
        if(this->rowsortl==0) {
            for(int i=0; i<this->N; i++)
                rp[i]=this->rowsort[i].r;
        } else {
            for(int i=0; i<this->N; i++)
                rp[i]=this->rowsortl[i];
        }
    }

    /// get row permutation
    void getRowperm(rowperm_t &rperm) const
    {
        if(this->rowsortl==0) {
            for ( rowindex_t x=0; x<this->N; x++ )
                rperm[x] = this->rowsort[x].r;
        } else {
            for ( rowindex_t x=0; x<this->N; x++ )
                rperm[x] = this->rowsortl[x];
        }
        // printfd("getRowperm: after "); print_perm(rperm, this->N);
    }

    /// return lightweight row permutation
    rowpermtypelight getRowperm() const
    {
        rowpermtypelight rp(this->N);
        this->getRowperm(rp);
        return rp;
    }


    /// return column permutation
    colpermtypelight getColperm() const
    {
        colpermtypelight cp(this->colperm, this->col+1);
        //  printf("created colpermtypelight: cp.n %d\n", cp.n);
        return cp;
    }
    /// set column permutation
    void getColperm(colpermtypelight &cp) const
    {
        cp.resize(this->col+1);
        std::copy(this->colperm, this->colperm+this->col+1, cp.d);

    }

    /// allocate lightweight rowsort structure
    void allocrowsortl()
    {
        if(this->rowsortl==0) {
            this->rowsortl=new_perm<rowindex_t>(this->N);
        }
    }

    void deleterowsortl()
    {
        if(this->rowsortl!=0) {
            delete_perm(this->rowsortl);
            //delete [] this->rowsortl;
            this->rowsortl=0;
        }

    }
    void initrowsortl()
    {
        if(this->rowsortl!=0) {
            //delete [] this->rowsortl;
            for(int i=0; i<this->N; i++) {
                this->rowsortl[i]=this->rowsort[i].r;
            }
        } else {
            this->rowsortl=new_perm<rowindex_t>(this->N);
            for(int i=0; i<this->N; i++) {
                this->rowsortl[i]=this->rowsort[i].r;
            }
        }
    }

    /// helper function
    void rowsortl2rowsort()
    {
        for(int i=0; i<this->N; i++) {
            this->rowsort[i].r=this->rowsortl[i];
        }
    }

    void copydata ( const dyndata_t &dd );
private:
    void initdata ( const dyndata_t &dd );

};


/// return 0 if target is equal to original, otherwise return 1 and copy root initialization + 1
inline int check_root_update ( array_t *array, const arraydata_t &ad )
{
    int changed=0;
    array_t *root = create_array ( ad.N, ad.strength );
    create_root ( root, &ad );
    if ( ! std::equal ( array, array+ad.N*ad.strength, root ) ) {
        //printf("check_root_update: copying!\n");
        copy_array ( root, array, ad.N, ad.strength );
        for ( int j=0; j<ad.N; j++ )
            array[ad.N*ad.strength+j]=ad.s[ad.strength];
        changed=1;

    }
    destroy_array ( root );

    return changed;
}

/// return true if target is in root form, otherwise return false
inline bool check_root_form ( const array_t *array, const arraydata_t &ad )
{
    array_t *root = create_array ( ad.N, ad.strength );
    create_root ( root, &ad );
    if ( std::equal ( array, array+ad.N*ad.strength, root ) ) {
        destroy_array ( root );
        return true;
    } else {
        destroy_array ( root );
        return false;
    }

}

/// return 0 if target is equal to original, otherwise return 1 and copy root initialization + 1
inline int check_root_update ( carray_t *original, const arraydata_t &ad, array_t *target )
{
    int changed=0;
    array_t *root = create_array ( ad.N, ad.strength );
    create_root ( root, &ad );
    if ( ! std::equal ( original, original+ad.N*ad.strength, root ) ) {
        //printf("check_root_update: copying!\n");
        copy_array ( root, target, ad.N, ad.strength );
        for ( int j=0; j<ad.N; j++ )
            target[ad.N*ad.strength+j]=ad.s[ad.strength]+100;
        changed=1;

    }
    destroy_array ( root );

    return changed;
}


typedef double jj45_t ;

/// return value based on J4-J5 ordering
jj45_t jj45val ( carray_t *array, rowindex_t N, int jj, const colperm_t comb, int j5val=-1, int dosort=1 );


/** Apply a random transformation to an array **/
void random_transformation ( array_t *array, const arraydata_t *adp );

/* main functions for LMC reduction */
lmc_t LMCreduction_train ( const array_link &al, const arraydata_t* ad, LMCreduction_t *reduction,const OAextend &oaextend ) ;
lmc_t LMCreduction_train ( const array_t* original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction,const OAextend &oaextend ) ;

/// helper function
lmc_t LMCreduce ( array_t const *original,  array_t const *array, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction,const OAextend &oaextend );
/// Perform reduction or LMC check without root trick
lmc_t LMCreduceFull ( carray_t* original, const array_t *array, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction,const OAextend &oaextend, LMC_static_struct_t &tmpStatic );


/// generic LMCcheck function
lmc_t LMCcheck ( const array_t *array, const arraydata_t &ad,  const OAextend &oaextend,LMCreduction_t &reduction );

//lmc_t LMCcheck ( const array_link &al, const arraydata_t &ad,  const OAextend &oaextend,LMCreduction_t &reduction );

/// generic LMCcheck function
inline lmc_t LMCcheck ( const array_link &al, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction )
{
    myassert(ad.N==al.n_rows, "LMCcheck: wrong number of rows");
    mycheck(ad.ncols<=al.n_columns, "LMCcheck: wrong number of columns al %d, adata %d", al.n_columns, ad.ncols);

    return LMCcheck ( al.array, ad,  oaextend, reduction );
}

/// direct LMC check using the original LMC check
lmc_t LMCcheckOriginal ( const array_link &al );


/// helper function
LMCreduction_t calculateSymmetryGroups(const array_link &al, const arraydata_t &adata, const OAextend &oaextend, int verbose=1, int hack=0);
lmc_t LMCcheckSymmetryMethod(const array_link &al, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction,  LMCreduction_t &reductionsub, int dverbose);

//#include <array>

template<class numtype>
/// Convert selection of elements to extended permutation
void combadd2perm(const larray<numtype> &comb, int newidx, int n, larray<numtype> &target, larray<numtype> &wtmp)
{
//    larray<numtype> w(n);
//    larray<numtype> wtmp(n);	// FIXME: allocation is expensive!
    for(int i=0; i<n; i++) {
        wtmp[i]=i;
    }
    int j=comb.size();
    for(int i=0; i<j; i++) {
        target[i]=comb[i];
        wtmp[comb[i]]=-1;
    }
    if (j<0) printf("combadd2perm: j<0! j %d\n", j);
    target[j]=newidx;
    wtmp[newidx]=-1;
    j++;

    for(int i=0; i<n; i++) {
        if (wtmp[i]>=0) {
            target[j]=wtmp[i];
            j++;
        }
    }
}

template<class numtype>
/// Convert selection of elements to extended permutation
larray<numtype> comb2perm(const larray<numtype> comb, int n)
{
    larray<numtype> w(n);
    larray<numtype> wtmp(n);
    for(int i=0; i<n; i++) {
        wtmp[i]=i;
    }
    int j=comb.size();
    for(int i=0; i<j; i++) {
        w[i]=comb[i];
        wtmp[comb[i]]=-1;
    }
    for(int i=0; i<n; i++) {
        if (wtmp[i]>=0) {
            w[j]=wtmp[i];
            j++;
        }
    }
    // std::array<numtype, n> a;
    return w;
}

template<class numtype>
/// Convert selection of elements to extended permutation
std::vector<numtype> comb2perm(const std::vector<numtype> comb, int n)
{
    std::vector<numtype> w(n);
    std::copy(comb.begin(), comb.end(), w.begin());
    numtype j = comb.size();
    for(int i=0; i<n; i++) {
        bool c = std::find(comb.begin(), comb.end(), i) != comb.end();
        if (!c) {
            w[j] = i;
            j++;
        }
    }
    // std::array<numtype, n> a;
    return w;
}


/// reduce arrays to canonical form using delete-1-factor ordering
void reduceArraysGWLP ( const arraylist_t *input_arrays, arraylist_t &reduced_arrays, int verbose, int dopruning = 1, int strength = 2, int dolmc=1 );

array_transformation_t reductionDOP(const array_link &al, int verbose=0);


void selectUniqueArrays(arraylist_t &xlist, arraylist_t &earrays, int verbose=1);

/* Public interface */

/// reduce an array to canonical form using LMC ordering
array_link reduceLMCform(const array_link &al);

/// reduce an array to canonical form using delete-1-factor ordering
array_link reduceDOPform(const array_link &al, int verbose=0);

/* Helper functions */

/** Apply LMC check (original mode) to a list of arrays */
std::vector<int> LMCcheckLex ( arraylist_t const &list, arraydata_t const &ad, int verbose=0 );

/// Perform LMC check lexicographically
lmc_t LMCcheckLex ( array_link const &al, arraydata_t const &ad );

lmc_t LMCcheckj4 ( array_link const &al, arraydata_t const &ad, LMCreduction_t &reduction, const OAextend &oaextend, int jj=4);

/// Perform LMC check with J5 ordering
lmc_t LMCcheckj5 ( array_link const &al, arraydata_t const &ad, LMCreduction_t &reduction, const OAextend &oaextend, int hack=0 );

/* internal functions LMC reduction */
#ifdef OAMEM
lmc_t LMCreduce_root_level_perm_ME ( carray_t const *original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction );
#endif

void root_row_permutation_from_index ( int permindex,  const arraydata_t *ad, levelperm_t *lperms );
//void root_row_permutation(int *perm, arraydata_t *ad, int num);

/* helper functions */

rowperm_t* create_root_permutations_index ( const arraydata_t *ad, int &totalpermsr );
void create_root_permutations_index_helper ( rowperm_t *rperms, levelperm_t *lperms, const arraydata_t *ad, int level, int &permcounter );

void print_rowsort ( rowsort_t *rowsort, int N );
void print_column_rowsort ( const array_t *arraycol, rowsort_t *rowsort, int N );
//void print_column_rowsort2 ( const array_t *arraycol, const array_t *arraycol2, rowsort_t *rowsort, levelperm_t lperm, int N );

void print_fracs ( int logl = NORMAL );
void clear_fracs();

#ifdef SWIG
%ignore print_rowsort;
%ignore root_row_permutation_from_index;
%ignore create_root_permutations_index;
%ignore create_root_permutations_index_helper;
#endif


#endif
