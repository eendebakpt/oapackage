/*! \file strength.h
        \brief Contains code to perform a strength check on an OA

   Author: Pieter Eendebak <pieter.eendebak@gmail.com>
   Copyright: See LICENSE.txt file that comes with this distribution
*/

#ifndef STRENGTH_H
#define STRENGTH_H

#include "arrayproperties.h"
#include "arraytools.h"
#include "printfheader.h"
#include "tools.h"

#ifdef FULLPACKAGE
#include "extend.h"
#endif

struct rev_index;

typedef int freq_t; /* used for counting t-tuples in strength check */

/** Create a table with frequencies for t-tuples for a set of column combinations */
typedef freq_t **strength_freq_table;

// strength_freq_table new_strength_freq_table ( int ncolcombs, int *nvalues, int &nelements );

rev_index *create_reverse_colcombs_fixed (const int ncolcombs);
// rev_index *create_reverse_colcombs ( colindex_t **colcombs, const int ncols, const int strength );
// void free_colcombs_fixed ( colindex_t **colcombs, int *lambda, int *nvalues );

/*!
  Structure serves as an index, gives the combinations in which a certain column participates. Is used to for the
  strength check, the program needs to know which column participate in the combination
  \brief Reverse Index
 */
struct rev_index {
        /// Number of combination a column is involved in, equal to ncombs(n-1, k-1)
        int nr_elements;
        /// List of combination a column is involved in
        int *index;
        /// pointer to struct with more information+settings for combination
};

struct strength_check_t {
        int freqtablesize;
        strength_freq_table freqtable;
        vindex_t **indices;
        rev_index *r_index;

        colindex_t **colcombs;
        int *nvalues;
        int *lambda;
        int ncolcombs;
        const int strength;

        strength_check_t (int str)
            : freqtable (0), indices (0), r_index (0), colcombs (0), nvalues (0), lambda (0), strength (str){};
        ~strength_check_t ();

        void set_colcombs (const arraydata_t &ad) {
                int verbose = 0;

                {
                        int prod;
                        int n = ad.ncols;
                        int N = ad.N;
                        const carray_t *s = ad.s;

                        int k = strength; // we keep 1 column fixed, choose k columns
                        ncolcombs = ncombs (n, k);

                        if (verbose)
                                myprintf ("strength_check_t: set_colcombs: ncolcombs %d, strength %d\n", ncolcombs,
                                          strength);
                        this->colcombs = malloc2d< colindex_t > (ncolcombs, strength);
                        this->lambda = (int *)malloc (ncolcombs * sizeof (int));
                        this->nvalues = (int *)malloc (ncolcombs * sizeof (int));

                        // set initial combination
                        for (int i = 0; i < k; i++)
                                colcombs[0][i] = i;

                        for (int i = 1; i < ncolcombs; i++) {
                                memcpy (colcombs[i], colcombs[i - 1], k * sizeof (colindex_t));
                                next_combination< colindex_t > (colcombs[i], k, n);

                                prod = 1;
                                for (int j = 0; j < strength; j++) {
                                        if (verbose >= 2)
                                                myprintf ("i %d j %d: %d\n", i, j, colcombs[i][j]);
                                        prod *= s[colcombs[i][j]];
                                }
                                nvalues[i] = prod;
                                lambda[i] = N / prod;
                        }

                        prod = 1; // First lambda manually, because of copy-for-loop
                        for (int j = 0; j < strength; j++)
                                prod *= s[colcombs[0][j]];
                        nvalues[0] = prod;
                        lambda[0] = N / prod;
                }
        }

        void info () const {
                myprintf ("strength_check_t: %d column combintations: \n", ncolcombs);
                for (int i = 0; i < ncolcombs; i++) {
                        myprintf ("   ");
                        print_perm (colcombs[i], strength);
                }
        }

        void create_reverse_colcombs_fixed () { r_index = ::create_reverse_colcombs_fixed (ncolcombs); }

        void print_frequencies () const {
                register int i, j;

                for (i = 0; i < ncolcombs; i++) {
                        myprintf ("%i:\t", i);
                        for (j = 0; j < nvalues[i]; j++) {
                                myprintf ("%2i ", freqtable[i][j]);
                        }
                        myprintf ("\n");
                }
                myprintf ("\n");
        }
};

/** @brief Contains static data for the extend loop
 *
 */
struct extend_data_t {

        /* static data: not changed during calculations */
        const arraydata_t *adata;
        const colindex_t extcolumn;
        rowindex_t oaindextmin; /** index of t-1 columns */

        //! number of rows
        const rowindex_t N;

        //! column combinations used in strength check
        colindex_t **colcombs;

        int **indices; /* todo: type? */

        int ncolcombs;      /** number of relevant column combinations */
        rev_index *r_index; /** reverse pointer to column combinations */
        rev_index *r_index_total;

        //! index of each column
        int *lambda;
        int lambda2lvl; // specialized version for 2level arrays

        //!
        int *nvalues;

        //! for row symmetry calculations
        rowindex_t *gidx;
        rowindex_t *gstart;
        rowindex_t *gsize;

/* dynamic data: can be changed during extension */

#ifdef COUNTELEMENTCHECK
        //! elements count, for strength 1 check
        int *elements;
#endif

        int freqtablesize;
        //! frequency table for strength check. For each column combination this table contains the frequencies of
        //! tuples found so far
        strength_freq_table freqtable;

        //! strength check, cache
        strength_freq_table *freqtable_cache;

        //! used strength check, cache for each element
        int **freqtable_elem;

        //! used strength check, for each row+element combination and column combination give a pointer to the position
        //! in the tuple frequence table
        int **element2freqtable;
        // freq_t* **element2freqtable;

        //! used for setting the range, range_high is inclusive
        array_t range_low, range_high;
        // the range is are the values that can be taken at a certain position in the array
        // rangemax is equal max(values in current col)+1
        // rangemin only occurs if row symmetries are used (define: USE_ROW_SYMMETRY)

        /* functions */

        extend_data_t (const arraydata_t *ad, colindex_t extcol);
        ~extend_data_t ();
};

/** @brief Copy the frequency count table
 */
inline void copy_freq_table (strength_freq_table source, strength_freq_table target, int ftsize) {
        memcpy ((void *)target[0], (void *)source[0], ftsize * sizeof (freq_t));
}

/// check whether an array passes divisibility test
bool check_divisibility (const arraydata_t *);

void print_frequencies (int **frequencies, const int nelements, const int *lambda, const int N);

/**
 * Return all column combinations including a fixed column.
 * At the same time allocate space for the number of values these columns have
 * @param xlambda
 * @param nvalues
 * @param ncolcombs
 * @param s
 * @param strength
 * @param fixedcol
 * @param N
 * @return
 */
colindex_t **set_colcombs_fixed (int *&xlambda, int *&nvalues, int &ncolcombs, const array_t *s, const int strength,
                                 const int fixedcol, const int N);

/// Add row to frequency table using cache system
void add_element_freqtable (extend_data_t *es, rowindex_t activerow, carray_t *array, strength_freq_table freqtable);
/// fast version of add_element_freqtable
void add_element_freqtable_col (extend_data_t *es, rowindex_t activerow, carray_t *arraycol,
                                strength_freq_table freqtable);

/// Initialize the table of t-tuple frequencies
void init_frequencies (extend_data_t *es, array_t *array);

void recount_frequencies (int **frequencies, extend_data_t *es, colindex_t currentcol, rowindex_t rowstart,
                          rowindex_t rowlast, carray_t *array);

/** perform strength check on an array
 *
 * Special case for extension of an array with proper strength
 *
 */
bool strength_check (const arraydata_t &ad, const array_link &al, int verbose = 1);

/// perform strength check on an array
bool strength_check (const array_link &al, int strength, int verbose = 0);

#ifdef FULLPACKAGE
/**
 * @brief Determine whether an element passes the strength test
 * @param es
 * @param p
 * @param array
 * @return
 */
bool valid_element (const extend_data_t *es, const extendpos *p, carray_t *array);

/** Determine whether an element passes the strength test, specialized for 2-level array */
bool valid_element_2level (const extend_data_t *es, const extendpos *p);

#endif

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
