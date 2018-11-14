/** \file lmc.h

\brief This file contains definitions and functions to perform LMC tests and reductions


Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2008

Copyright: See LICENSE file that comes with this distribution
*/

/*! \mainpage	Orthogonal Arrays

This package contains code to calculate orthogonal arrays with a specified number of runs, factor levels and strength.
The main program is oaextendsingle (or oaextendmpi for the multi-core version), which starts with the specifications of
an OA and creates the root of th arrray. Then the array is extended column wise.
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

#pragma once

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "arrayproperties.h"
#include "mathtools.h"
#include "tools.h"

#ifdef SWIG
#else
#define ORDER_J5_SMALLER >
#define ORDER_J5_GREATER <
#define ORDER_J45_SMALLER >
#define ORDER_J45_GREATER <
#endif

// forward declaration
class OAextend;

#define stringify(name) #name

#ifdef SWIG
%ignore check_root_update;
%ignore LMCreduceFull;
%ignore dyndata_t::dyndata_t (dyndata_t const &);
#endif

#ifdef LMCSTATS
// function is not thread safe
void lmc_stats ();
#endif

/* constants and structures */

/** Possible results for the LMC check
 *
 * LMC_LESS: Found a permutation which leads to a lexicographically smaller array
 * LMC_EQUAL: Found a permutation which leads to a lexicographically equal array
 * LMC_MORE: Found a permutation which leads to a lexicographically larger array
 */
enum lmc_t { LMC_LESS, LMC_EQUAL, LMC_MORE, LMC_NONSENSE };

/// different algorithms
enum algorithm_t {
        MODE_ORIGINAL,
        MODE_J4,
        MODE_J5ORDER,
        MODE_J5ORDERX,
        MODE_INVALID,
        MODE_AUTOSELECT,
        MODE_LMC_SYMMETRY,
        MODE_LMC_2LEVEL,
        MODE_LMC_DEBUG,
        MODE_J5ORDERXFAST
};
#define MODE_LMC MODE_ORIGINAL

inline std::string algorithm_t_list () {
        std::string ss =
            printfstring ("%d (automatic), %d (original), %d (check j4), %d (j5 order), %d (j5 order dominant), %d "
                          "(MODE_J5ORDERXFAST)",
                          MODE_AUTOSELECT, MODE_ORIGINAL, MODE_J4, MODE_J5ORDER, MODE_J5ORDERX, MODE_J5ORDERXFAST);
        ss += printfstring (", %d (MODE_LMC_SYMMETRY), %d (MODE_LMC_2LEVEL)", MODE_LMC_SYMMETRY, MODE_LMC_2LEVEL);
        return ss;
}

/// method used for initialization of columns
enum initcolumn_t { INITCOLUMN_ZERO, INITCOLUMN_PREVIOUS, INITCOLUMN_J5 };

/// variations of the J45 structures
enum j5structure_t { J5_ORIGINAL, J5_45 };

/// return name of the algorithm
std::string algnames (algorithm_t m);

struct dyndata_t;

// NOTE: unsigned long is enough for 2-factor arrays up to 60 columns
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
static inline bool operator< (const rowsort_t &a, const rowsort_t &b) {
#ifdef OAOVERFLOW
        assert (a.val >= 0);
        assert (b.val >= 0);
#endif
        return a.val < b.val;
}

/**
 * @brief Comparision operator for the rowsort_t structure
 */
static inline bool operator> (const rowsort_t &a, const rowsort_t &b) {
#ifdef OAOVERFLOW
        assert (a.val >= 0);
        assert (b.val >= 0);
#endif

        return a.val > b.val;
}

/// Apply Hadamard transformation to orthogonal array
void apply_hadamard (const arraydata_t *ad, array_t *array, colindex_t hcol);

/// Apply Hadamard transformation to orthogonal array
void apply_hadamard (array_link &al, colindex_t hcol);

/**
 * @brief Contains structures used by the LMC reduction or LMC check
 *
 *  Part of the allocations is for structures that are constant and are re-used each time an LMC calculation is
 * performed. Some other structures are temporary buffers that are written to all the time.
 */
struct LMC_static_struct_t {
      public:
        /* variable data */
        int LMC_non_root_init;
        int LMC_root_init;
        int LMC_reduce_root_rowperms_init;

        /* data structures used by functions */
        arraydata_t *ad;

        /* used at the root_level_perm stage */
        int LMC_root_rowperms_init;
        int nrootrowperms;
        rowperm_t *rootrowperms;

        int LMC_root_rowperms_init_full;
        int nrootrowperms_full;
        rowperm_t *rootrowperms_full;

        array_t *colbuffer; /** buffer for a single column */

        dyndata_t **dyndata_p;       /** dynamic data; row permutations */
        colindex_t **colperm_p;      /** column permutations */
        colindex_t **localcolperm_p; /** local column permutation */

        array_transformation_t *current_trans; /* not used at the moment? */

        LMC_static_struct_t ();
        ~LMC_static_struct_t ();

        void show (int verbose = 1) const {
                myprintf ("LMC_static_struct_t: ad %p, LMC_non_root_init %d\n", (void *)(this->ad),
                          LMC_non_root_init);
        }

        void init (const arraydata_t *adp);
        void freeall ();

        /// update structure with new design specification
        int update (const arraydata_t *adp);
        int needUpdate (const arraydata_t *adp) const;

        void init_root_stage (levelperm_t *&lperm_p, colperm_t *&colperm_p, const arraydata_t *adp);
        void init_nonroot_stage (levelperm_t *&lperm_p, colperm_t *&colperm_p, colperm_t *&localcolperm_p,
                                 dyndata_t **&dynd_p, int &dynd_p_nelem, array_t *&colbuffer,
                                 const arraydata_t *adp) const;

        /// Static initialization of root row permutations
        void init_rootrowperms (int &totalperms, rowperm_t *&rootrowperms, levelperm_t *&lperm_p) {
                /* no static update, we assume this has been done already */

                totalperms = this->nrootrowperms;
                rootrowperms = this->rootrowperms;
                lperm_p = this->current_trans->lperms;

                this->LMC_root_rowperms_init = 1; 
        }

        /** @brief Static initialization of root row permutations (full group)
         */
        void init_rootrowperms_full (int &totalperms, rowperm_t *&rootrowperms, levelperm_t *&lperm_p) {
                /* no static update, we assume this has been done already */

                totalperms = this->nrootrowperms_full;
                rootrowperms = this->rootrowperms_full;
                lperm_p = this->current_trans->lperms;

                this->LMC_root_rowperms_init_full = 1;
        }
};

//LMC_static_struct_t *getGlobalStaticIndexed (int n);
//void cleanGlobalStaticIndexed ();

/// return static structure from dynamic global pool, return with releaseGlobalStatic
LMC_static_struct_t *getGlobalStatic ();
void releaseGlobalStatic (LMC_static_struct_t *p);

/// release all objects in the pool
void cleanGlobalStatic ();

/// variable indicating the state of the reduction process
enum REDUCTION_STATE { REDUCTION_INITIAL, REDUCTION_CHANGED };
//! main mode for the LMC routine: test, reduce or reduce with initialization
enum OA_MODE { OA_TEST, OA_REDUCE, OA_REDUCE_PARTIAL };

typedef larray< rowindex_t > rowpermtypelight;
typedef larray< colindex_t > colpermtypelight;

typedef std::vector< int > colpermtype;
typedef std::vector< colpermtype > colpermset;

#ifdef _WIN32
// on windows do not use smart pointers, it is a mess
//#if _MSC_VER >= 1600
#elif __APPLE__
#define SDSMART
#else
#define SDSMART
#endif

#ifdef SDSMART
#ifdef WIN32

#include <memory>
typedef std::shared_ptr< symmdata > symmdataPointer;

#elif defined(__APPLE__)

#include <memory>
typedef std::shared_ptr< symmdata > symmdataPointer;

/** @brief An operator to test whether a pointer holds a null value or not
 *
 * It is recomended to not compare pointers to integers and use explicit
 * conversion `ptr to bool` instead or as another approach to compare pointers against
 * nullptr which has type of std::nullptr. But such methods would force migration to C++11. Hence,
 * the purpose of this method is to avoid ether rewritting the whole project or changing the code
 * that already works on other platforms (aka Windows) by adding a missing in <memory> header
 * comparison `bool operator!=(smart_ptr, int)`
 */
bool operator!= (symmdataPointer const &ptr, int x);

#else
#include <tr1/memory>
typedef std::tr1::shared_ptr< symmdata > symmdataPointer;
// typedef std::shared_ptr<symmdata> symmdataPointer;
#endif

#else
typedef symmdata *symmdataPointer;
#endif

enum INIT_STATE { INIT_STATE_INVALID, COPY, INIT, SETROOT };

/// Append element to vector if the element the element is not at the end of vector
template < class Type > void insert_if_not_at_end_of_vector (std::vector< Type > &cp, const Type &cpv) {
        if (cp.size () > 0) {
                if (!(cp.back () == cpv))
                        cp.push_back (cpv);
        } else {
                cp.push_back (cpv);
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

        INIT_STATE init_state; // initalization mode: INIT: initialized by user, COPY: copy from array argument

        //! maximum depth for search tree
        int maxdepth;

        /// last column visited in algorithm
        int lastcol;

        //! counter for number of reductions made
        long nred;

        int targetcol; 
        int mincol;    

        int nrows, ncols;

        LMC_static_struct_t *staticdata;

        symmdataPointer sd;

      public:
        LMCreduction_t (const LMCreduction_t &at); /// copy constructor
        LMCreduction_t (const arraydata_t *ad);
        ~LMCreduction_t ();

        LMCreduction_t &operator= (const LMCreduction_t &at); /// Assignment operator

        array_link getArray () const {
                array_link al (this->nrows, this->ncols, 0, this->array);
                return al;
        }

        void setArray (const array_link al) {
                init_state = INIT;
                std::copy (al.array, al.array + this->ncols * this->nrows, this->array);
        }
        void setArray (const array_t *array, int nrows, int ncols) {
                init_state = INIT;
                std::copy (array, array + ncols * nrows, this->array);
        }

        /// update the pointer to the symmetry data based on the specified array
        void updateSDpointer (const array_link al, bool cache = false) {
#ifdef SDSMART
                symmdata *sdp = this->sd.get ();
#else
                symmdata *sdp = sd;
#endif
                if (sdp != 0 && cache) { 
                } else {
                        // update symmetry data
                        this->sd = symmdataPointer (new symmdata (al, 1));
                }
        }

        void releaseStatic () {
                if (this->staticdata != 0) {
                        releaseGlobalStatic (this->staticdata);
                        this->staticdata = 0;
                }
        }

        /// acquire a reference to a LMC_static_struct_t object
        void initStatic () {
                if (this->staticdata == 0) {
                        this->staticdata = getGlobalStatic ();
                }
        }

        /// return a reference to a object with LMC reduction data
        LMC_static_struct_t &getStaticReference () {
                if (this->staticdata == 0) {
		        this->initStatic();
		}
		return *(this->staticdata);
        }

        /// reset the reduction: clears the symmetries and sets the transformation to zero
        void reset ();

		void show(int verbose = 2) const;

        std::string __repr__ () const {
                std::string ss = printfstring ("LMCreduction_t: mode %d, state %d (REDUCTION_INITIAL %d, "
                                               "REDUCTION_CHANGED %d), init_state %d, lastcol %d\n",
                                               this->mode, this->state, REDUCTION_INITIAL, REDUCTION_CHANGED,
                                               this->init_state, this->lastcol);
                return ss;
        }

        /// called whenever we find a reduction
        void updateFromLoop (const arraydata_t &ad, const dyndata_t &dynd, levelperm_t *lperms,
                             const array_t *original);
        void updateTransformation (const arraydata_t &ad, const dyndata_t &dynd, levelperm_t *lperms,
                                   const array_t *original);

        inline bool doBreak (lmc_t ret) {
                if (this->mode >= OA_REDUCE)
                        return false;
                else {
                        return ret == LMC_LESS;
                }
        }
        inline void updateLastCol (int col) { this->lastcol = col; }

      private:
        void free ();
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
        // private:
        //! ordering of rows
        rowsort_t *rowsort;
        rowsortlight_t *rowsortl;

      public:
        //! current column permutation
        colperm_t colperm;

        dyndata_t (int N, int col = 0);
        dyndata_t (const dyndata_t *dd);
        dyndata_t (const dyndata_t &);
        ~dyndata_t ();

        dyndata_t &operator= (const dyndata_t &);
        void show () const;

        void reset ();
        void setColperm (const colperm_t perm, int n) { copy_perm (perm, this->colperm, n); }
        void setColperm (const larray< colindex_t > &perm) { std::copy (perm.data_pointer, perm.data_pointer + perm.data_size, this->colperm); }

        void setColperm (const std::vector< colindex_t > &perm) {
                std::copy (perm.begin (), perm.end (), this->colperm);
        }

        /// get lightweight row permutation
		void getRowperm(rowpermtypelight &rp) const;

        /// get row permutation
		void getRowperm(rowperm_t &rperm) const;

        /// return lightweight row permutation
		rowpermtypelight getRowperm() const;

        /// return column permutation
		colpermtypelight getColperm() const;

        /// set column permutation
		void getColperm(colpermtypelight &cp) const;

        /// allocate lightweight rowsort structure
        void allocate_rowsortl () {
                if (this->rowsortl == 0) {
                        this->rowsortl = new_perm< rowindex_t > (this->N);
                }
        }

        void deleterowsortl () {
                if (this->rowsortl != 0) {
                        delete_perm (this->rowsortl);
                        this->rowsortl = 0;
                }
        }

		/// initialize rowsortl from rowsort
        void initrowsortl () {
                if (this->rowsortl != 0) {
                        for (int i = 0; i < this->N; i++) {
                                this->rowsortl[i] = this->rowsort[i].r;
                        }
                } else {
                        this->rowsortl = new_perm< rowindex_t > (this->N);
                        for (int i = 0; i < this->N; i++) {
                                this->rowsortl[i] = this->rowsort[i].r;
                        }
                }
        }

        /// copy rowsortl variable to rowsrt
        void rowsortl2rowsort () {
                for (int i = 0; i < this->N; i++) {
                        this->rowsort[i].r = this->rowsortl[i];
                }
        }

        void copydata (const dyndata_t &dd);

      private:
        void initdata (const dyndata_t &dd);
};

/// return true if target is in root form, otherwise return false
inline bool check_root_form (const array_t *array, const arraydata_t &ad) {
        array_t *root = create_array (ad.N, ad.strength);
        create_root (root, &ad);
        if (std::equal (array, array + ad.N * ad.strength, root)) {
                destroy_array (root);
                return true;
        } else {
                destroy_array (root);
                return false;
        }
}

/// return 0 if target is equal to original, otherwise return 1 and copy root initialization + 1
inline int check_root_update (carray_t *original, const arraydata_t &ad, array_t *target) {
        int changed = 0;

        array_t *root = create_array (ad.N, ad.strength);
        create_root (root, &ad);
        if (!std::equal (original, original + ad.N * ad.strength, root)) {
                copy_array (root, target, ad.N, ad.strength);
                for (int j = 0; j < ad.N; j++)
                        target[ad.N * ad.strength + j] = ad.s[ad.strength] + 100;
                changed = 1;
        }
        destroy_array (root);

        return changed;
}

typedef double jj45_t;


/* helper function for LMC reduction */
lmc_t LMCreduction_train (const array_link &al, const arraydata_t *ad, LMCreduction_t *reduction,
                          const OAextend &oaextend);
/* helper function for LMC reduction */
lmc_t LMCreduction_train (const array_t *original, const arraydata_t *ad, const dyndata_t *dyndata,
                          LMCreduction_t *reduction, const OAextend &oaextend);

/// helper function
lmc_t LMCreduce (array_t const *original, array_t const *array, const arraydata_t *ad, const dyndata_t *dyndata,
                LMCreduction_t *reduction, const OAextend &oaextend);
/// Perform reduction or LMC check without root trick
lmc_t LMCreduceFull (carray_t *original, const array_t *array, const arraydata_t *ad, const dyndata_t *dyndata,
                     LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic);

/// generic LMCcheck function
lmc_t LMCcheck (const array_t *array, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction);

/// generic LMCcheck function
lmc_t LMCcheck (const array_link &al, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction);

/// direct LMC check using the original LMC check
lmc_t LMCcheckOriginal (const array_link &al);

/// reduce arrays to canonical form using delete-1-factor ordering
void reduceArraysGWLP (const arraylist_t &input_arrays, arraylist_t &reduced_arrays, int verbose, int dopruning = 1,
                       int strength = 2, int dolmc = 1);

/// Caculate the transformation to delete-on-factor normal for from an array
array_transformation_t reductionDOP (const array_link &array, int verbose = 0);

/// select the unique arrays in a list, the original list is sorted in place. the unique arrays are append to the output list
void selectUniqueArrays (arraylist_t &input_arrays, arraylist_t &output_arrays, int verbose = 1);

/** Calculate projection values for delete-of-factor algorithm
 *
 */
std::vector< GWLPvalue > projectionDOFvalues (const array_link &array, int verbose = 0);

/// reduce an array to canonical form using LMC ordering
array_link reduceLMCform (const array_link &al);

/// reduce an array to canonical form using delete-1-factor ordering
array_link reduceDOPform (const array_link &al, int verbose = 0);

/** Apply LMC check (original mode) to a list of arrays */
std::vector< int > LMCcheckLex (arraylist_t const &list, arraydata_t const &ad, int verbose = 0);

/// Perform LMC check lexicographically
lmc_t LMCcheckLex (array_link const &al, arraydata_t const &ad);

lmc_t LMCcheckj4 (array_link const &al, arraydata_t const &ad, LMCreduction_t &reduction, const OAextend &oaextend,
                  int jj = 4);

/// Perform LMC check with J5 ordering
lmc_t LMCcheckj5 (array_link const &al, arraydata_t const &ad, LMCreduction_t &reduction, const OAextend &oaextend,
                  int hack = 0);

/* internal functions LMC reduction */
#ifdef OAMEM
lmc_t LMCreduce_root_level_perm_ME (carray_t const *original, const arraydata_t *ad, const dyndata_t *dyndata,
                                    LMCreduction_t *reduction);
#endif


/**
* @brief Print the contents of a rowsort structure
* @param rowsort Pointer to rowsort structure
* @param N Number of elements
*/
void print_rowsort (rowsort_t *rowsort, int N);
void print_column_rowsort (const array_t *arraycol, rowsort_t *rowsort, int N);


