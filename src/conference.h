/** \file conference.h

Contains functionality to generate and analyse conference designs.

For more information see:

* https://en.wikipedia.org/wiki/Conference_matrix
* "A Classification Criterion for Definitive Screening Designs", Schoen et al., The Annals of Statistics, 2018

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2018

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
void print_column(const conference_column &column, const char *msg = 0);

/** Show a list of candidate extensions
 *
 * \param column_candidates List of candidates to show
 */
void showCandidates (const std::vector< conference_column > &column_candidates);

/** Convert conference design to definitive screening design
 *
 * The DSD is created by appending the negated design to the conference design and then appending a row of zeros.
 *
 * \param conference_design Array with the conference design
 * \param add_zeros If True, then append a row of zeros
 */
array_link conference2DSD(const array_link &conference_design, bool add_zeros = 1);

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
		/// number of runs
        rowindex_t N;     
		/// total number of columns (factors) in the design
        colindex_t ncols; 

        /// Type of conference design
        enum conference_type { 
            /// normal conference design
            CONFERENCE_NORMAL,
            /// conference design with zeros only on diagonal
            CONFERENCE_DIAGONAL,
            /// double conference design
            DCONFERENCE };
		/// defines the type of designs
		conference_type ctype;      
		/// defines the isomorphism type
        matrix_isomorphism_t itype; 

		/// if true then J1 values should be zero
        bool j1zero; 
		/// if true then J3 values should be zero
        bool j3zero; 

      public:
        /// create new conference_t object
        conference_t ();
        conference_t (int N, int k, int j1zero);
        conference_t (const conference_t &rhs);

        // return short string describing the class
        std::string idstr () const;

        /// create the unique representative of the 2 column conference design 
        array_link create_root () const;

        /// create the unique representative of the 3 column conference design
        array_link create_root_three () const;

        /// create the root arrays with 1 column for the double conference matrices
        arraylist_t createDconferenceRootArrays () const;

        /// add the root arrays to a list
        void addRootArrays (arraylist_t &lst) const;

        /// return string representation of the object
        std::string __repr__ () const {
                return printfstring ("conference class: number of rows %d, number of columns %d", this->N, this->ncols);
        }
};

/** Reduce conference matrix to normal form using Nauty
 *
 * @see reduceConferenceTransformation
 */
array_link reduceConference (const array_link &, int verbose = 0);

/** Reduce conference matrix to normal form using Nauty
 *
 * The design is converted to a graph representation. The graph is then reduced using Nauty 
 * to normal form and the resulting graph translated back to a conference design.
 *
 * \param conference_design Design to be reduced to normal form
 * \param verbose Verbosity level
 * \returns A transformation that converts the input design to normal form
 *
 */
conference_transformation_t reduceConferenceTransformation (const array_link &conference_design, int verbose);

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
        void showCandidates (int verbose = 1) const;

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

/** Extend a list of conference designs with a single column.
 *
 */
arraylist_t extend_conference (const arraylist_t &lst, const conference_t conference_type, int verbose,
                               int select_isomorphism_classes = 0);

/// plain version without caching
arraylist_t extend_conference_plain (const arraylist_t &lst, const conference_t conference_type, int verbose,
                                     int select_isomorphism_classes = 0);

/** Extend a list of conference designs with a single column */
arraylist_t extend_conference_restricted (const arraylist_t &lst, const conference_t conference_type, int verbose);

// extend a list of double conference matrices with an additional column
arraylist_t extend_double_conference (const arraylist_t &lst, const conference_t conference_type, int verbose);

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

/// sort list of conference designs according to LMC-0 ordering
arraylist_t sortLMC0 (const arraylist_t &lst);

/// LMC0 check for double conference matrix
lmc_t LMC0checkDC (const array_link &al, int verbose = 0);

/// check if array is in LM0 form
lmc_t LMC0check (const array_link &al, int verbose = 0);

/// return true if the design is a foldover array
bool isConferenceFoldover (const array_link &al, int verbose = 0);

/** For a double conference design return a row permutation to a single conference design
 * 
 * If the design is not a foldover design then the first element of the returned permutation is -1.
 * 
 * \param double_conference A double conference design
 * \returns Permutation
 * 
 */
std::vector<int> double_conference_foldover_permutation(const array_link &double_conference);

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
		  DconferenceFilter(const array_link &_als, int filtersymm_, int filterj2_, int filterj3_ = 1);

        /// print object to stdout
		void show() const;

        /// filter a list of columns using the filter method
		std::vector< conference_column > filterList(const std::vector< conference_column > &lst, int verbose = 0) const;

		std::vector< conference_column > filterListJ2last(const std::vector< conference_column > &column_list) const;

        /// filter a list of cperms using the filterZero method
		std::vector< conference_column > filterListZero(const std::vector< conference_column > &lst) const;

        /// return True if the extension satisfies all checks
        bool filter (const conference_column &c) const;

        /** filter on partial column (only last col)
         *
         * maxrow (int): the number of rows that are valid
         **/
        bool filterJpartial (const conference_column &c, int maxrow) const;

        /// return True if the extension satisfies all J-characteristic checks
		bool filterJ(const conference_column &column, int j2start = 0) const;

        /// return True if the extension satisfies all J-characteristic checks for the last columns
		bool filterJlast(const conference_column &c, int j2start = 0) const;

        /// return True if the extension satisfies all checks. prints the reason for returning True or False to stdout
        bool filterReason (const conference_column &column) const;

        /// return True if the candidate satisfies the J3 check
		bool filterJ3(const conference_column &column) const;

        /// return True if the candidate satisfies the J3 check for specified pairs
		bool filterJ3s(const conference_column &column, int idxstart) const;
        /// return True if the candidate satisfies the J3 check
		bool filterJ3inline(const conference_column &column) const;

        /// return True of the candidate satisfies the symmetry check
        bool filterSymmetry (const conference_column &column) const;

        /// return True of the candidate extension satisfies the J2 check
        bool filterJ2 (const conference_column &c) const;
        /// return True of the candidate extension satisfies the J2 check for the last column of the array checked against
        bool filterJ2last (const conference_column &c) const;
        /** return True of the candidate extension satisfies the zero check
         *
         * This means that the first entries of the extension do not contain a zero.
         */
		bool filterZero(const conference_column &c) const;

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

