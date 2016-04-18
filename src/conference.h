/** \file oatest.cpp

 C++ program: oatest

 oatest: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "arraytools.h"
#include "graphtools.h"
#include "arrayproperties.h"
//#include "extend.h"


/// structure to cache a list of candidate extensions
struct conf_candidates_t {
public:
	std::vector<std::vector<cperm> > ce;

	void info ( int verbose=1 ) const {
		for ( int i=2; i< ( int ) ce.size(); i++ ) {
			if ( verbose ) {
				printf ( "generateCandidateExtensions: k %d: %d candinates\n", i, ( int ) ce[i].size() );
			}
		}
	}
};


/// Structure representing the type of conference designs
class conference_t
{
public:
	rowindex_t N;	/** number of runs */
	colindex_t ncols;	/** total number of columns (factors) in the design */

	enum conference_type {CONFERENCE_NORMAL, CONFERENCE_DIAGONAL, DCONFERENCE};
	conference_type ctype;
	matrix_isomorphism_t itype;

public:
	/// create new conference_t object
	conference_t ( int N, int k );

	/// create the unique representative of the 2 column design (for conference matrices)
	array_link create_root ( ) const;

	/// create the unique representative of the 3 column design
	array_link create_root_three ( ) const;

	arraylist_t createDconferenceRootArrays() const {
		arraylist_t lst;
		array_link al ( this->N, 1, array_link::INDEX_DEFAULT );
		for ( int i=N/2+1; i<=N; i++ ) {
			al.setconstant ( -1 );
			al.at ( 0,0 ) =0;
			al.at ( 1,0 ) =0;
			for ( int k=2; k<i; k++ )
				al.at ( k,0 ) =1;
			lst.push_back ( al );
		}

		return lst;
	}

	void addRootArrays ( arraylist_t &lst ) const {
		switch ( this->ctype ) {
		case CONFERENCE_NORMAL:
		case CONFERENCE_DIAGONAL:
			switch ( this->itype ) {
			case CONFERENCE_ISOMORPHISM:
				lst.push_back ( this->create_root() );
				break;
			case CONFERENCE_RESTRICTED_ISOMORPHISM: {
				array_link al ( this->N, 1, array_link::INDEX_DEFAULT );
				for ( int i=N/2+1; i<=N; i++ ) {
					al.setconstant ( -1 );
					al.at ( 0,0 ) =0;
					for ( int k=1; k<i; k++ )
						al.at ( k,0 ) =1;
					lst.push_back ( al );
				}

			}
			break;
			default
					:
				printfd ( "not implemented (itype %d)\n", this->itype );
			}
			break; //
		case DCONFERENCE: {
			switch ( this->itype ) {
			case CONFERENCE_RESTRICTED_ISOMORPHISM: {
				arraylist_t tmp = this->createDconferenceRootArrays();
				lst.insert ( lst.end(), tmp.begin(), tmp.end() );
			}
			break;
			case CONFERENCE_ISOMORPHISM:
			default
					:
				printfd ( "not implemented (itype %d)\n", this->itype );
			}
		}
		}
		if ( 0 ) {
			for ( size_t i=0; i<lst.size(); i++ ) {
				printf ( "root array %d:\n", int ( i ) );
				lst[i].showarray();
			}
		}
	}
/// return string representation of the object
	std::string __repr__() const {
		return printfstring ( "conference type: N %d, ncols %d", this->N, this->ncols );
	}

};

/// reduce conference matrix to normal form
array_link reduceConference ( const array_link &, int verbose = 0 );

/// reduce conference matrix to normal form
conference_transformation_t reduceConferenceTransformation ( const array_link &al, int verbose );



/// Helper structure
struct conference_extend_t {
	std::vector<cperm> first;
	std::vector<cperm> second;
	std::vector<cperm> extensions;

public:

	// combine first and second section into a single column
	cperm combine ( int i, int j ) const {
		cperm c =vstack ( this->first[i], this->second[j] );

		//printfd("c.size() %d = %d + %d\n", c.size(),  this->first[i].size(),  this->second[i].size() );
		return c;
	}

	size_t nExtensions() const {
		return this->extensions.size();
	}

	arraylist_t getarrays ( const array_link al ) const {
		arraylist_t ll;

		for ( size_t i=0; i<this->extensions.size(); i++ ) {
			array_link alx = hstack ( al, extensions[i] );
			ll.push_back ( alx );
		}
		return ll;
	}
};

struct conference_options {
	int maxzpos;

	//conference_options() { maxzpos=-1; };
	conference_options ( int maxpos = -1 ); // { maxzpos=-1; };
} ;

/** Extend a single conference design with candidate columns */
conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t &ct, int extcol, int verbose=1, int maxzpos=-1 );

/// helper function
conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t & ct, int extcol, int verbose, int maxzpos, const conf_candidates_t &cande );

/** Extend a list of conference designs with a single column.
 *
 */
arraylist_t extend_conference ( const arraylist_t &lst, const conference_t ctype, int verbose );

/** Extend a list of conference designs with a single column */
arraylist_t extend_conference_restricted ( const arraylist_t &lst, const conference_t ctype, int verbose );

// extend a list of double conference matrices
arraylist_t extend_double_conference ( const arraylist_t &lst, const conference_t ctype, int verbose ) ;


/// select representatives for the isomorphism classes of a list of conference arrays
arraylist_t  selectConferenceIsomorpismClasses ( const arraylist_t list, int verbose,  matrix_isomorphism_t itype = CONFERENCE_ISOMORPHISM );

/// select representatives for the isomorphism classes of a list of conference arrays, return indices of classes
std::vector<int> selectConferenceIsomorpismIndices ( const arraylist_t lst, int verbose,  matrix_isomorphism_t itype = CONFERENCE_ISOMORPHISM );

/** Generate candidate extensions
 *
 * \param al design to be extended
 * \param kz index of zero in candidate column
 */
std::vector<cperm> generateConferenceExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose = 1, int filtersymm= 1, int filterip=1 );

/** Generate candidate extensions for restricted isomorphism classes */
std::vector<cperm> generateConferenceRestrictedExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose=1 , int filtersymm=1, int filterip=1 );

std::vector<cperm> generateDoubleConferenceExtensions ( const array_link &al, const conference_t & ct, int verbose=1 , int filtersymm=1, int filterip=1 );

/** return max position of zero in array, returns -1 if no zero is found
 *
 * The parameter k specifies the column to search in. For k=-1 all columns are searched.
 */
int maxz ( const array_link &al, int k = -1 );

/** Return true of the array is smaller in LMC-0 ordering
 *
 */
bool compareLMC0 ( const array_link &alL, const array_link &alR );

/// sort list of arrays according to LMC-0 ordering
arraylist_t sortLMC0 ( const arraylist_t &lst );


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
