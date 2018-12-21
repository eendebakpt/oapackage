/*! \file strength.h
        \brief Contains code to perform a strength check on an OA

   Author: Pieter Eendebak <pieter.eendebak@gmail.com>
   Copyright: See LICENSE.txt file that comes with this distribution
*/

#pragma once

#include "arrayproperties.h"
#include "arraytools.h"
#include "printfheader.h"
#include "tools.h"

#ifdef FULLPACKAGE
#include "extend.h"
#endif

struct rev_index;

/** Create a table with frequencies for t-tuples for a set of column combinations */
typedef int **strength_freq_table;

rev_index *create_reverse_colcombs_fixed (const int ncolcombs);

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

		void set_colcombs(const arraydata_t &ad);
		
        void info () const {
                myprintf ("strength_check_t: %d column combintations: \n", ncolcombs);
                for (int i = 0; i < ncolcombs; i++) {
                        myprintf ("   ");
                        print_perm (colcombs[i], strength);
                }
        }

		void create_reverse_colcombs_fixed();
		
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
        //! frequency table for strength check. For each column combination this table contains the frequencies of tuples found so far
        strength_freq_table freqtable;

        //! strength check, cache
        strength_freq_table *freqtable_cache;

        //! used strength check, cache for each element
        int **freqtable_elem;

        //! used strength check, for each row+element combination and column combination give a pointer to the position
        //! in the tuple frequence table
        int **element2freqtable;

        //! used for setting the range, range_high is inclusive
        array_t range_low, range_high;
        // the range is are the values that can be taken at a certain position in the array
        // rangemax is equal max(values in current col)+1
        // rangemin only occurs if row symmetries are used (define: USE_ROW_SYMMETRY)

        extend_data_t (const arraydata_t *ad, colindex_t extcol);
        ~extend_data_t ();

		/// Initialize the table of t-tuple frequencies
		void init_frequencies(array_t *array);
};

/// Checks on the divisibility of the number of runs by the product of the levels in the factors for all t-tuples
bool check_divisibility(int N, int ncols, int strength, const array_t * s);

/// check whether an array passes divisibility test
bool check_divisibility (const arraydata_t *);


/// Add row to frequency table using cache system
void add_element_freqtable (extend_data_t *es, rowindex_t activerow, carray_t *array, strength_freq_table freqtable);
/// fast version of add_element_freqtable
void add_element_freqtable_col (extend_data_t *es, rowindex_t activerow, carray_t *arraycol,
                                strength_freq_table freqtable);

void recount_frequencies (int **frequencies, extend_data_t *es, colindex_t currentcol, rowindex_t rowstart,
                          rowindex_t rowlast, carray_t *array);

/** Perform strength check on an array
 *
 * Special case for extension of an array with proper strength
 *
 **/
bool strength_check (const arraydata_t &ad, const array_link &al, int verbose = 1);

/// perform strength check on an array
bool strength_check (const array_link &al, int strength, int verbose = 0);

/**
 * @brief Determine whether an element passes the strength test
 * @param es
 * @param position
 * @param array Pointer to array
 * @return
 */
bool valid_element (const extend_data_t *es, const extendpos *position, carray_t *array);

/** Determine whether an element passes the strength test, specialized for 2-level array */
bool valid_element_2level (const extend_data_t *es, const extendpos *p);

