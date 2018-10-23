/** \file conference.h

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2016

 Copyright: See LICENSE.txt file that comes with this distribution
*/
#pragma once

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "arrayproperties.h"
#include "arraytools.h"
#include "graphtools.h"

/// print a candidate extension
void print_cperm(const conference_column &c, const char *msg = 0);

/** Show a list of candidate extensions
 *
 * \param cc List of candidates to show
 */
void showCandidates (const std::vector< conference_column > &cc);

/// convert conference design to definitive screening design
array_link conference2DSD(const array_link &conf, bool add_zeros = 1);

/// structure to cache a list of candidate extensions
struct conf_candidates_t {
      public:
        /// list of candidate extentions for each number of columns
        std::vector< std::vector< conference_column > > ce;

        /// print information about the set of candidate extentions
        void info (int verbose = 1) const {
                for (int i = 2; i < (int)ce.size (); i++) {
                        if (verbose) {
                                myprintf ("generateCandidateExtensions: k %d: %d candinates\n", i, (int)ce[i].size ());
                        }
                }
        }
};

/// Structure representing the type of conference designs
class conference_t {
      public:
        rowindex_t N;     /// number of runs
        colindex_t ncols; /// total number of columns (factors) in the design

        /// Type of conference design
        enum conference_type { CONFERENCE_NORMAL, CONFERENCE_DIAGONAL, DCONFERENCE };
        conference_type ctype;      /// defines the type of designs
        matrix_isomorphism_t itype; /// defines the isomorphism type

        bool j1zero; /// if true then J1 values should be zero
        bool j3zero; /// if true then J3 values should be zero

      public:
        /// create new conference_t object
        conference_t ();
        conference_t (int N, int k, int j1zero);
        conference_t (const conference_t &rhs);

        // return short string describing the class
        std::string idstr () const;

        /// create the unique representative of the 2 column design (for conference matrices)
        array_link create_root () const;

        /// create the unique representative of the 3 column design
        array_link create_root_three () const;

        /// create the root arrays with 1 column for the double conference matrices
        arraylist_t createDconferenceRootArrays () const;

        /// add the root arrays to a list
        void addRootArrays (arraylist_t &lst) const;

        /// return string representation of the object
        std::string __repr__ () const {
                return printfstring ("conference type: N %d, ncols %d", this->N, this->ncols);
        }
};

/// reduce conference matrix to normal form
array_link reduceConference (const array_link &, int verbose = 0);

/// reduce conference matrix to normal form
conference_transformation_t reduceConferenceTransformation (const array_link &al, int verbose);

/** Helper structure containing extensions of conference designs
 */
struct conference_extend_t {
        std::vector< conference_column > first;      /// list of first block candidate extensions
        std::vector< conference_column > second;     /// list of first block candidate extensions
        std::vector< conference_column > extensions; /// list of candidate extensions

      public:
        // combine first and second section into a single column
        conference_column combine (int i, int j) const {
                conference_column c = vstack (this->first[i], this->second[j]);
                return c;
        }

        size_t nExtensions () const { return this->extensions.size (); }

        /// return the set of extension arrays
        arraylist_t getarrays (const array_link al) const {
                arraylist_t ll;

                for (size_t i = 0; i < this->extensions.size (); i++) {
                        array_link alx = hstack (al, extensions[i]);
                        ll.push_back (alx);
                }
                return ll;
        }
};

/// return true if zero is a specified position
inline bool checkZeroPosition (const conference_column &p, int zero_position) {
        if (zero_position <= 0)
                return false;

        if (p[zero_position] == 0) {
                return true;
        } else
                return false;
}

/** Class to generate candidate extensions with caching
 *
 * We assume that the designs to be extended are run ordered, so that the caching has maximal effect.
 *
 * The key idea used is that any valid extension of a design A with k columns is a permutation of a valid extension
 * of the design B obtained by taking the first l < k columns of A. The permutations that are allowed are called the
 * symmetry inflations. All the j2 checks performed for the extension of B do not have to be repeated for the
 * permutations of this extension.
 **/
class CandidateGeneratorBase {
      public:
        /// type of designs to generate
        conference_t ct;
        /// verbosity level
        int verbose;

        /// last array analyzed
        mutable array_link al;
        /// index of last valid column
        mutable int last_valid;

      protected:
        /// list of candidate extensions. the elements of candidate_list[k] correspond to columns with index k-1
        mutable std::vector< conference_column_list > candidate_list;

      public:
        CandidateGeneratorBase (const array_link &al, const conference_t &ct);

        /** Show the candidate extensions for each column
         *
         */
        void showCandidates (int verbose = 1) const {
                myprintf ("CandidateGenerator: N %d\n", this->ct.N);
                for (int i = 2; i <= last_valid; i++) {
                        myprintf ("CandidateGenerator: number of candidates for %dth column: %ld\n", i,
                                  (long)candidate_list[i].size ());
                        if (verbose >= 2) {
                                ::showCandidates (candidate_list[i]);
                        }
                }
        }

        /// return all candidates for the kth column
        conference_column_list candidates (int k);

      protected:
        static const int START_COL = 2;

        /** Find the starting column for the extension of a design
         *
         * The static variable START_COL is the number of columns for which is the starting point if the cache is
         * empty. Therefore
         * for a design with initial columns the same, START_COL+1 is the first number of columns with valid entries.
         *
         * For startcol k the elements in candidate_list[k] are valid, e.g. we can start with extensions valid for
         * index k-1
         */
        int startColumn (const array_link &alx, int verbose = 0) const {
                if (this->al.n_columns != alx.n_columns) {
                        int startcol = -1;
                        if (verbose) {
                                printfd ("startColumn: return -1\n");
                        }
                        return startcol;
                }
                int startcol = al.firstColumnDifference (alx) + 1;
                if (verbose) {
                        printfd ("startColumn: firstColumnDifference %d, last_valid %d\n", startcol, last_valid);
                }

                startcol = std::min (startcol, last_valid);
                if (startcol < this->START_COL) {
                        startcol = -1;
                }
                return startcol;
        }
};

/// Class to generate conference candidate extensions
class CandidateGeneratorConference : public CandidateGeneratorBase {

      public:
        CandidateGeneratorConference (const array_link &al, const conference_t &ct);

        /// Generate a list of candidate extensions for the specified design
        const std::vector< conference_column > &generateCandidates (const array_link &al) const;

        /// generate all candidate extensions with a zero at the specified position
        std::vector< conference_column > generateCandidatesZero (const array_link &al, int kz) const;
};

typedef CandidateGeneratorConference CandidateGenerator;

/// Class to generate double conference candidate extensions with caching
class CandidateGeneratorDouble : public CandidateGeneratorBase {

      public:
        CandidateGeneratorDouble (const array_link &al, const conference_t &ct);

        /** Generate a list of candidate extensions for the specified design
         *
         * This method uses symmetry inflation, assumes j1=0 and j2=0. Optimal performance is achieved when the arrays
         * to be extended have identical first columns.
         *
         */
        const std::vector< conference_column > &generateCandidates (const array_link &al) const;
};

/** Extend a single conference design with candidate columns */
conference_extend_t extend_conference_matrix (const array_link &al, const conference_t &ct, int extcol,
                                              int verbose = 1, int maxzpos = -1);

/** Extend a list of conference designs with a single column.
 *
 */
arraylist_t extend_conference (const arraylist_t &lst, const conference_t ctype, int verbose,
                               int select_isomorphism_classes = 0);

/// plain version without caching
arraylist_t extend_conference_plain (const arraylist_t &lst, const conference_t ctype, int verbose,
                                     int select_isomorphism_classes = 0);

/** Extend a list of conference designs with a single column */
arraylist_t extend_conference_restricted (const arraylist_t &lst, const conference_t ctype, int verbose);

// extend a list of double conference matrices
arraylist_t extend_double_conference (const arraylist_t &lst, const conference_t ctype, int verbose);

/// select representatives for the isomorphism classes of a list of conference arrays
arraylist_t selectConferenceIsomorpismClasses (const arraylist_t &list, int verbose,
                                               matrix_isomorphism_t itype = CONFERENCE_ISOMORPHISM);

/// select representatives for the isomorphism classes of a list of conference arrays, return indices of classes
std::vector< int > selectConferenceIsomorpismIndices (const arraylist_t &lst, int verbose,
                                                      matrix_isomorphism_t itype = CONFERENCE_ISOMORPHISM);

/// select double conference arrays in LMC0 form
arraylist_t selectLMC0doubleconference (const arraylist_t &list, int verbose, const conference_t &ctype);

/// select conference arrays in LMC0 form
arraylist_t selectLMC0 (const arraylist_t &list, int verbose, const conference_t &ctype);

/** Generate candidate extensions (wrapper function)
 *
 * \param al Design to be extended
 * \param conference_type Class of conference designs
 * \param zero_index index of zero in candidate column
 * \param verbose Verbosity level
 * \param filtersymm If True, filter based on symmetry
 * \param filterj2 If True, filter based on J2 values
 */
std::vector< conference_column > generateConferenceExtensions (const array_link &al, const conference_t &conference_type,
                                                   int zero_index, int verbose = 1, int filtersymm = 1,
                                                   int filterj2 = 1);

/** Generate candidate extensions for restricted isomorphism classes */
std::vector< conference_column > generateConferenceRestrictedExtensions (const array_link &al, const conference_t &conference_type,
                                                             int zero_index, int verbose = 1, int filtersymm = 1,
                                                             int filterip = 1);

// generate extensions for double conference matrices in LMC0 form
std::vector< conference_column > generateDoubleConferenceExtensions (const array_link &al, const conference_t &conference_type,
                                                         int verbose = 1, int filtersymm = 1, int filterip = 1,
                                                         int filterJ3 = 0, int filtersymminline = 1);

// generate extensions for conference matrices in LMC0 form
std::vector< conference_column > generateSingleConferenceExtensions (const array_link &al, const conference_t &conference_type,
                                                         int zero_index, int verbose, int filtersymm, int filterj2,
                                                         int filterj3, int filtersymminline = 0);

/** return max position of zero in array, returns -1 if no zero is found
 *
 * The parameter k specifies the column to search in. For k=-1 all columns are searched.
 */
int maxz (const array_link &al, int column_index = -1);

/** Return true of the array is smaller in LMC-0 ordering
 *
 */
bool compareLMC0 (const array_link &alL, const array_link &alR);

/// sort list of arrays according to LMC-0 ordering
arraylist_t sortLMC0 (const arraylist_t &lst);

/// LMC0 check for double conference matrix
lmc_t LMC0checkDC (const array_link &al, int verbose = 0);

/// check if array is in LM0 form
lmc_t LMC0check (const array_link &al, int verbose = 0);

/// return true if the design is a foldover array
bool isConferenceFoldover (const array_link &al, int verbose = 0);

// return true if the extension column satisfies the inner product check
int ipcheck (const conference_column &col, const array_link &al, int cstart = 2, int verbose = 0);

/// return minimal position of zero in design
int minz (const array_link &al, int column_index);

/// class to filter single or double conference designs
class DconferenceFilter {
      public:
        array_link als;

        /// filter based on symmetry
        int filtersymm;
        /// filter based on j2 value
        int filterj2;
        /// filter based on j3 value
        int filterj3;
        /// filter only columns with first value >=0
        int filterfirst;
        /// filter based on first occurence of zero in a column
        int filterzero;


        mutable long ngood;

      private:
        /// table of J2 vectors for J3 filter
        array_link dtable;
        /// table of J2 vectors for inline J3 filter
        array_link inline_dtable;

        /// indices to check for symmetry check
        std::vector< int > check_indices;
        /// used for filtering based on zero
        int minzvalue;

      public:
        int inline_row;
        symmdata sd;

      public:
        DconferenceFilter (const array_link &_als, int filtersymm_, int filterj2_, int filterj3_ = 1)
            : als (_als), filtersymm (filtersymm_), filterj2 (filterj2_), filterj3 (filterj3_), filterfirst (0),
              filterzero (0), ngood (0), sd (als) {

                check_indices = sd.checkIdx ();
                dtable = createJ2tableConference (als);

                if (als.n_columns >= 2) {
                        inline_dtable = als.selectColumns (0) - als.selectColumns (1);
                        inline_dtable = hstack (inline_dtable, als.selectColumns (0) + 1);
                        inline_dtable = hstack (inline_dtable, als.selectColumns (0) * als.selectColumns (0) - 1);
                        inline_dtable = hstack (inline_dtable, als.selectColumns (1) * als.selectColumns (1) - 1);

                        minzvalue = minz (als, als.n_columns - 1);

                        inline_row = als.n_rows;
                        int br = 0;
                        for (int i = als.n_rows - 1; i >= 0; i--) {
                                for (int c = 0; c < als.n_columns; c++) {
                                        if (inline_dtable.at (i, 0) != 0) {
                                                br = 1;
                                                break;
                                        }
                                }
                                if (br) {
                                        break;
                                }
                                inline_row = i;
                        }
                } else {
                        inline_row = -1;
                }
        }

        /// print object to stdout
        void show () const {
                myprintf ("DconferenceFilter: filterj1 -, filterj2 %d, filterj3 %d, filtersymm %d\n", filterj2,
                          filterj3, filtersymm);
        }

        /// filter a list of cperms using the filter method
        std::vector< conference_column > filterList (const std::vector< conference_column > &lst, int verbose = 0) const {
                std::vector< conference_column > out;
                for (size_t i = 0; i < lst.size (); i++) {
                        if (this->filter (lst[i])) {
                                out.push_back (lst[i]);
                        }
                }
                if (verbose) {
                        printfd ("filterList: %d -> %d\n", lst.size (), out.size ());
                }
                return out;
        }

        std::vector< conference_column > filterListJ2last (const std::vector< conference_column > &lst) const {
                std::vector< conference_column > out;
                for (size_t i = 0; i < lst.size (); i++) {
                        if (this->filterJ2last (lst[i])) {
                                out.push_back (lst[i]);
                        }
                }
                // printfd("filterListZero: minzvalue %d: %d -> %d\n", minzvalue, lst.size(), out.size() );
                return out;
        }

        /// filter a list of cperms using the filterZero method
        std::vector< conference_column > filterListZero (const std::vector< conference_column > &lst) const {
                std::vector< conference_column > out;
                for (size_t i = 0; i < lst.size (); i++) {
                        if (this->filterZero (lst[i])) {
                                out.push_back (lst[i]);
                        }
                }
                return out;
        }

        /// return True of the extension satisfies all checks
        bool filter (const conference_column &c) const;

        /** filter on partial column (only last col)
         *
         * r (int): the number of rows that are valid
         **/
        bool filterJpartial (const conference_column &c, int r) const;

        /// return True of the extension satisfies all J-characteristic checks
        bool filterJ (const conference_column &c, int j2start = 0) const {
                if (filterj2) {
                        // perform inner product check for all columns
                        if (!ipcheck (c, als, j2start)) {
                                return false;
                        }
                }
                if (filterj3) {
                        // perform inner product check for all columns
                        if (!this->filterJ3 (c)) {
                                return false;
                        }
                }
                ngood++;
                return true;
        }

        /// return True of the extension satisfies all J-characteristic checks for the last columns
        bool filterJlast (const conference_column &c, int j2start = 0) const {
                if (filterj2) {
                        // perform inner product check for all columns
                        if (!ipcheck (c, als, j2start)) {
                                return false;
                        }
                }
                int startidx = this->dtable.n_columns - this->als.n_columns;
                if (filterj3) {
                        // perform inner product check for all columns
                        if (!this->filterJ3s (c, startidx)) {
                                return false;
                        }
                }
                ngood++;
                return true;
        }
        /// return True of the extension satisfies all checks
        bool filterReason (const conference_column &c) const;

        /// return True of the candidate satisfies the J3 check
        bool filterJ3 (const conference_column &c) const {
                const int nc = dtable.n_columns;
                const int N = als.n_rows;
                int jv = 0;
                for (int idx1 = 0; idx1 < nc; idx1++) {
                        jv = 0;

                        const array_t *o1 = dtable.array + dtable.n_rows * idx1;
                        for (int xr = 0; xr < N; xr++) {

                                jv += (o1[xr]) * (c[xr]);
                        }

                        if (jv != 0) {
                                return false;
                        }
                }
                return true;
        }

        /// return True of the candidate satisfies the J3 check for specified pairs
        bool filterJ3s (const conference_column &c, int idxstart) const {
                const int nc = dtable.n_columns;
                const int N = als.n_rows;
                int jv = 0;
                for (int idx1 = nc - 1; idx1 >= idxstart; idx1--) {
                        jv = 0;

                        const array_t *o1 = dtable.array + dtable.n_rows * idx1;
                        for (int xr = 0; xr < N; xr++) {

                                jv += (o1[xr]) * (c[xr]);
                        }

                        if (jv != 0) {
                                return false;
                        }
                }
                return true;
        }
        /// return True of the candidate satisfies the J3 check
        bool filterJ3r (const conference_column &c) const {
                const int nc = dtable.n_columns;
                const int N = als.n_rows;
                int jv = 0;
                for (int idx1 = nc - 1; idx1 >= 0; idx1--) {
                        jv = 0;

                        const array_t *o1 = dtable.array + dtable.n_rows * idx1;
                        for (int xr = 0; xr < N; xr++) {

                                jv += (o1[xr]) * (c[xr]);
                        }

                        if (jv != 0) {
                                return false;
                        }
                }
                return true;
        }
        /// return True of the candidate satisfies the J3 check
        bool filterJ3inline (const conference_column &c) const {
                const int nc = inline_dtable.n_columns;
                const int N = als.n_rows;
                int jv = 0;
                for (int idx1 = 0; idx1 < nc; idx1++) {
                        jv = 0;

                        const array_t *o1 = inline_dtable.array + inline_dtable.n_rows * idx1;
                        for (int xr = 0; xr < N; xr++) {

                                jv += (o1[xr]) * (c[xr]);
                        }

                        if (jv != 0) {
                                return false;
                        }
                }
                return true;
        }

        /// return True of the candidate satisfies the symmetry check
        bool filterSymmetry (const conference_column &c) const;

        /// return True of the candidate extension satisfies the J2 check
        bool filterJ2 (const conference_column &c) const { return ipcheck (c, als, 0); }
        /// return True of the candidate extension satisfies the J2 check for the last column of the array checked
        /// against
        bool filterJ2last (const conference_column &c) const { return ipcheck (c, als, als.n_columns - 1); }
        /** return True of the candidate extension satisfies the zero check
         *
         * This means that the first entries of the extension do not contain a zero.
         */
        bool filterZero (const conference_column &c) const {
                for (int i = 0; i < minzvalue - 1; i++) {
                        if (c[i] == 0) {
                                return false;
                        }
                }
                return true;
        }

      private:
};

/** Inflate a candidate column
 *
 * The extensions are generated according to the symmertry specified by the symmetry group. Filtering is performed
 *using the filter object.
 *
 * From the filtering object only the J2 filtering is used.
 *
 * \return List of inflated extensions
 **/
std::vector< conference_column > inflateCandidateExtension (const conference_column &basecandidate, const array_link &als,
                                                const symmetry_group &alsg, const std::vector< int > &check_indices,
                                                const conference_t &ct, int verbose, const DconferenceFilter &filter);

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;
