/** \file conference.h

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2016

 Copyright: See LICENSE.txt file that comes with this distribution
*/
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "arraytools.h"
#include "graphtools.h"
#include "lmc.h"
#include "arrayproperties.h"


/// print a candidate extension
inline void print_cperm ( const cperm &c )
{
	for ( size_t i=0; i<c.size(); i++ ) {
		printf ( "%3d", c[i] );
	}
}

/// show a list of candidate extensions
inline void showCandidates ( const std::vector<cperm> &cc )
{
	for ( size_t i=0; i<cc.size(); i++ ) {
		printf ( "%d: ", ( int ) i );
		print_cperm ( cc[i] );
		printf ( "\n" );
	}
}

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

	bool j3zero;
	bool j1zero; /// for the double conference type matrices

public:
	/// create new conference_t object
	conference_t ( int N, int k );
	conference_t ( const conference_t &rhs );

	/// create the unique representative of the 2 column design (for conference matrices)
	array_link create_root ( ) const;

	/// create the unique representative of the 3 column design
	array_link create_root_three ( ) const;

	arraylist_t createDconferenceRootArrays ( ) const {

		arraylist_t lst;
		array_link al ( this->N, 1, array_link::INDEX_DEFAULT );
		if ( j1zero ) {
			al.setconstant ( -1 );
			al.at ( 0,0 ) =0;
			al.at ( 1,0 ) =0;
			for ( int k=2; k<N/2+1; k++ )
				al.at ( k,0 ) =1;
			lst.push_back ( al );

		} else {
			for ( int i=N/2+1; i<=N; i++ ) {
				al.setconstant ( -1 );
				al.at ( 0,0 ) =0;
				al.at ( 1,0 ) =0;
				for ( int k=2; k<i; k++ )
					al.at ( k,0 ) =1;
				lst.push_back ( al );
			}
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
				//const int j1zero=1;
				arraylist_t tmp = this->createDconferenceRootArrays ( );
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

/// Class to generate candidate extensions with caching
class CandidateGenerator
{
public:
	conference_t ct;
	int verbose;
	
	// legacy structure
		conf_candidates_t cande;

	mutable std::vector<array_link> alz;

private:

	/// candidate_list_conf[kz][k] are the valid candidates with a zero at position kz and at column index k
	mutable std::vector< std::vector<cperm_list > > candidate_list_conf;
	/// index of last valid entry in cache
			mutable std::vector<int> last_valid_conf;

public:

	CandidateGenerator ( const array_link &al, const conference_t &ct );

	/** generate candidates with caching
	 * this method uses j2 filtering
	 */
	const std::vector<cperm> & generateConfCandidates ( const array_link &al, int kz ) const;

	void showCandidates() const {
		
		printf ( "CandidateGenerator: N %d\n", this->ct.N );
		for(int kz=0; kz<this->ct.N;kz++) {
		for ( int i =2; i<=this->last_valid_conf[kz]; i++ ) {
			printf ( "CandidateGenerator: %d columns: %ld elements\n", i, ( long ) candidate_list_conf[kz][i].size() );
		}
		}
	}

	void updateLastValid(int lv) const {
		// special case: if the data is invalid for one value of kz, we make all entries invalid
		// not efficient in general, but efficient for our orderly generation
		for(size_t i=0; i<this->last_valid_conf.size(); i++) {
		this->last_valid_conf[i]=std::min(this->last_valid_conf[i], lv);	
		}
	}
private:
	/** find the starting column for the extension
	 *
	 * Startcol is the index of the first column with valid elements, if corresponds to
	 * startcol=k the elements in candidate_list_conf[kz][k] are valid
	 */
	int startColumn ( const array_link &alx, int kz, int verbose=0 ) const {
		if ( this->alz[kz].n_columns!=alx.n_columns ) {
			int startcol=-1;
			//printfd(" ### reset because of column diff!\n");
			return startcol;
		}
		int startcol = alz[kz].firstColumnDifference ( alx );

		startcol = std::min ( startcol, last_valid_conf[kz] );
		if ( startcol<2 ) {
			startcol=-1;
		}
		
		if (verbose || startcol==-2) {
				printfd("startColumn: startcol %d\n", startcol);
				printf("-- cache --\n");
				this->alz[kz].transposed().showarray();
				printf("-- new --\n");
				alx.transposed().showarray();
		}
		
		return startcol;
	}

};


/// Class to generate candidate extensions with caching
class CandidateGeneratorDouble
{
public:
	conference_t ct;
	int verbose;
	

	mutable array_link al;
			mutable int last_valid;

private:
	mutable std::vector< cperm_list > candidate_list_double;


public:

	CandidateGeneratorDouble ( const array_link &al, const conference_t &ct );

	/** generate candidates with caching
	 * this method uses symmetry inflation
	 */
	const std::vector<cperm> & generateDoubleConfCandidates ( const array_link &al ) const;

	void showCandidates() const {
		printf ( "CandidateGenerator: N %d\n", this->ct.N );
		for ( int i =2; i<=last_valid; i++ ) {
			printf ( "CandidateGenerator: %d columns: %ld elements\n", i, ( long ) candidate_list_double[i].size() );
		}
	}
	static const  int START_COL = 2;
private:

	/** find the starting column for the extension
	 *
	 * For startcol k the elements in candidate_list[k] are valid
	 */
	int startColumn ( const array_link &alx ) const {
		if ( this->al.n_columns!=alx.n_columns ) {
			int startcol=-1;
			return startcol;
		}
		int startcol = al.firstColumnDifference ( alx );


		if ( 0 ) {
			int rx, ry;
			al.firstDiff ( alx, rx, ry, 1 );
			printfd ( " ---> startcol %d, last_valid %d\n", startcol, last_valid );
			printf ( " ---> cache array\n" );
			this->al.showarray();
			printf ( " --->alx\n" );
			alx.showarray();
		}
		startcol = std::min ( startcol, last_valid );
		if ( startcol<2 )
			startcol=-1;
		return startcol;
	}
};



/** Extend a single conference design with candidate columns */
conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t &ct, int extcol, int verbose=1, int maxzpos=-1 );

/// helper function
conference_extend_t extend_conference_matrix_generator ( const array_link &al, const conference_t & ct, int extcol, int verbose, int maxzpos, const CandidateGenerator &cgenerator );
//conference_extend_t extend_conference_matrix_generator ( const array_link &al, const conference_t & ct, int extcol, int verbose, int maxzpos, const CandidateGeneratorDouble &cgenerator );

/** Extend a list of conference designs with a single column.
 *
 */
arraylist_t extend_conference ( const arraylist_t &lst, const conference_t ctype, int verbose, int select_isomorphism_classes = 0 );

/** Extend a list of conference designs with a single column */
arraylist_t extend_conference_restricted ( const arraylist_t &lst, const conference_t ctype, int verbose );

// extend a list of double conference matrices
arraylist_t extend_double_conference ( const arraylist_t &lst, const conference_t ctype, int verbose ) ;


/// select representatives for the isomorphism classes of a list of conference arrays
arraylist_t  selectConferenceIsomorpismClasses ( const arraylist_t &list, int verbose,  matrix_isomorphism_t itype = CONFERENCE_ISOMORPHISM );

/// select representatives for the isomorphism classes of a list of conference arrays, return indices of classes
std::vector<int> selectConferenceIsomorpismIndices ( const arraylist_t &lst, int verbose,  matrix_isomorphism_t itype = CONFERENCE_ISOMORPHISM );

/** Generate candidate extensions
 *
 * \param al design to be extended
 * \param kz index of zero in candidate column
 */
std::vector<cperm> generateConferenceExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose = 1, int filtersymm= 1, int filterip=1 );

/** Generate candidate extensions for restricted isomorphism classes */
std::vector<cperm> generateConferenceRestrictedExtensions ( const array_link &al, const conference_t & ct, int kz, int verbose=1 , int filtersymm=1, int filterip=1 );

std::vector<cperm> generateDoubleConferenceExtensions ( const array_link &al, const conference_t & ct, int verbose=1 , int filtersymm=1, int filterip=1, int filterJ3=0, int filtersymminline = 1 );


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


/// return true if the design is a foldover array
bool isConferenceFoldover ( const array_link &al, int verbose = 0 );


/// helper function, return true if a candidate extensions satisfies the symmetry test
int satisfy_symm ( const cperm &c, const symmdata & sd, int rowstart=2 );

/// helper function, return true if a candidate extensions satisfies the symmetry test
int satisfy_symm ( const cperm &c, const std::vector<int>  & check_indices, int rowstart=2 );

/// helper function, return true if a candidate extensions satisfies the symmetry test
int satisfy_symm ( const cperm &c, const std::vector<int>  & check_indices, int rowstart, int rowend );

// return true if the extension column satisfies the inner product check
int ipcheck ( const cperm &col, const array_link &al, int cstart=2, int verbose=0 );


/// return minimal position of zero in design
int minz ( const array_link &al, int k );

/// class to filter designs
class DconferenceFilter
{
public:
	array_link als;
	int filtersymm; // filter based on symmetry
	int filterj2; /// filter based on j2 value
	int filterj3; /// filter based on j3 value
	int filterfirst; /// filter only columns with first value >=0
	int filterzero; /// filter based on first occurence of zero in a column

	mutable long ngood;
private:
	array_link dtable; /// table of J2 vectors for J3 filter
	array_link inline_dtable; /// table of J2 vectors for inline J3 filter

	/// indices to check for symmetry check
	std::vector<int> check_indices;
	/// used for filtering based on zero
	int minzvalue;
public:
	int inline_row;
	symmdata sd;

public:
	DconferenceFilter ( const array_link &_als, int filtersymm_, int filterj2_, int filterj3_ = 1 ) : als ( _als ), filtersymm ( filtersymm_ ), filterj2 ( filterj2_ ), filterj3 ( filterj3_ ), filterfirst ( 0 ), filterzero ( 0 ), ngood ( 0 ), sd ( als ) {
		//sd = symmdata( als );

		check_indices = sd.checkIdx();

		dtable = createJ2tableConference ( als );

		if ( als.n_columns>=2 ) {
			inline_dtable = als.selectColumns ( 0 )-als.selectColumns ( 1 ); // createJ2tableConference ( als.selectFirstColumns(2) );
			inline_dtable = hstack ( inline_dtable, als.selectColumns ( 0 ) +1 );
			inline_dtable = hstack ( inline_dtable, als.selectColumns ( 0 ) *als.selectColumns ( 0 )-1 );
			inline_dtable = hstack ( inline_dtable, als.selectColumns ( 1 ) *als.selectColumns ( 1 )-1 );

			//inline_dtable = als.selectColumns(0)*als.selectColumns(0)-1;
			//inline_dtable = als.selectColumns(0)+1; // createJ2tableConference ( als.selectFirstColumns(2) );


			minzvalue = minz ( als, als.n_columns-1 );

			inline_row = als.n_rows;
			int br=0;
			for ( int i=als.n_rows-1; i>=0; i-- ) {
				for ( int c=0; c<als.n_columns; c++ ) {
					if ( inline_dtable.at ( i,0 ) !=0 ) {
						br=1;
						break;
					}
				}
				if ( br ) {
					break;
				}
				inline_row=i;
			}
			//inline_dtable.showarray();
			//printfd("  inline J3 check: inline_row %d\n", inline_row);
		} else {
			inline_row=-1;
		}
	}

	/// filter a list of cperms using the filter method
	std::vector<cperm> filterList ( const std::vector<cperm> &lst, int verbose=0 ) const {
		std::vector<cperm> out;
		for ( size_t i=0; i<lst.size(); i++ ) {
			if ( this->filter ( lst[i] ) )
				out.push_back ( lst[i] );
		}
		if (verbose)
			printfd("filterList: %d -> %d\n", lst.size(), out.size() );
		return out;
	}

	std::vector<cperm> filterListJ2last ( const std::vector<cperm> &lst ) const {
		std::vector<cperm> out;
		for ( size_t i=0; i<lst.size(); i++ ) {
			if ( this->filterJ2last ( lst[i] ) )
				out.push_back ( lst[i] );
		}
		//printfd("filterListZero: minzvalue %d: %d -> %d\n", minzvalue, lst.size(), out.size() );
		return out;
	
	}
	
	/// filter a list of cperms using the filterZero method
	std::vector<cperm> filterListZero ( const std::vector<cperm> &lst ) const {
		std::vector<cperm> out;
		for ( size_t i=0; i<lst.size(); i++ ) {
			if ( this->filterZero ( lst[i] ) )
				out.push_back ( lst[i] );
		}
		//printfd("filterListZero: minzvalue %d: %d -> %d\n", minzvalue, lst.size(), out.size() );
		return out;
	}

	/// return True of the extension satisfies all checks
	bool filter ( const cperm &c ) const {
		if ( filterfirst ) {
			if ( c[0]<0 ) {
				return false;
			}
		}
		if ( filtersymm ) {
			if ( ! satisfy_symm ( c, check_indices, 0 ) ) {
				return false;
			}
		}
		if ( filterj2 ) {
			// perform inner product check for all columns
			if ( ! ipcheck ( c, als, 0 ) ) {
				return false;
			}
		}
		if ( filterj3 ) {
			// perform inner product check for all columns
			if ( ! this->filterJ3 ( c ) ) {
				return false;
			}
		}
		ngood++;
		return true;
	}
	/// return True of the extension satisfies all J-characteristic checks
	bool filterJ ( const cperm &c, int j2start=0 ) const {
		if ( filterj2 ) {
			// perform inner product check for all columns
			if ( ! ipcheck ( c, als, j2start ) ) {
				return false;
			}
		}
		if ( filterj3 ) {
			// perform inner product check for all columns
			if ( ! this->filterJ3 ( c ) ) {
				return false;
			}
		}
		ngood++;
		return true;
	}
	/// return True of the extension satisfies all checks
	bool filterReason ( const cperm &c ) const {
		if ( filterfirst ) {
			if ( c[0]<0 ) {
				printf ( "filterfirst\n" );
				return false;
			}
		}
		if ( filtersymm ) {
			if ( ! satisfy_symm ( c, sd, 0 ) ) {
				printf ( "symmetry\n" );
				return false;
			}
		}
		if ( filterj2 ) {
			// perform inner product check for all columns
			if ( ! ipcheck ( c, als, 0 ) ) {
				printf ( "j2\n" );
				return false;
			}
		}
		if ( filterj3 ) {
			// perform inner product check for all columns
			if ( ! this->filterJ3 ( c ) ) {
				printf ( "j3\n" );
				return false;
			}
		}
		ngood++;
		printf ( "filter check good\n" );

		return true;
	}

	/// return True of the candidate satisfies the J3 check
	bool filterJ3 ( const cperm &c ) const {
		const int nc = dtable.n_columns;
		const int N = als.n_rows;
		int jv=0;
		for ( int idx1=0; idx1<nc; idx1++ ) {
			jv=0;

			const array_t *o1 = dtable.array+dtable.n_rows*idx1;
			for ( int xr=0; xr<N; xr++ ) {

				jv += ( o1[xr] ) * ( c[xr] );
			}

			if ( jv!=0 )
				return false;
		}
		return true;
	}

	/// return True of the candidate satisfies the J3 check
	bool filterJ3inline ( const cperm &c ) const {
		const int nc = inline_dtable.n_columns;
		const int N = als.n_rows;
		int jv=0;
		for ( int idx1=0; idx1<nc; idx1++ ) {
			jv=0;

			const array_t *o1 = inline_dtable.array+inline_dtable.n_rows*idx1;
			for ( int xr=0; xr<N; xr++ ) {

				jv += ( o1[xr] ) * ( c[xr] );
			}

			if ( jv!=0 )
				return false;
		}
		return true;
	}


	/// return True of the candidate satisfies the symmetry check
	bool filterSymmetry ( const cperm &c ) const {
		return  satisfy_symm ( c, check_indices, 0 );
	}
	/// return True of the candidate extension satisfies the J2 check
	bool filterJ2 ( const cperm &c ) const {
		return ipcheck ( c, als, 0 );
	}
	/// return True of the candidate extension satisfies the J2 check for the last column of the array checked against
	bool filterJ2last ( const cperm &c ) const {
		return ipcheck ( c, als, als.n_columns-1 );
	}
	/// return True of the candidate extension satisfies the zero
	bool filterZero ( const cperm &c ) const {
		// TODO: minzvalue-1?
		for ( int i=0; i<minzvalue-1; i++ ) {
			if ( c[i]==0 )
				return false;
		}
		return true;

	}

private:

};

std::vector<cperm> generateDoubleConferenceExtensionsInflate ( const array_link &al, const conference_t &ct, int verbose, int filterj2, int filterj3, int kstart=2 );

// inflate candindate column: FIXME: documentation
std::vector<cperm> inflateCandidateExtension ( const cperm &basecandidate,  const array_link &als,  const symmetry_group &alsg, const std::vector<int> &check_indices, const conference_t & ct, int verbose , const DconferenceFilter &filter );
//std::vector<cperm> inflateCandidateExtension ( const cperm &basecandidate,  const array_link &als, const array_link &alx, const conference_t & ct, int verbose , const DconferenceFilter &filter );



// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
