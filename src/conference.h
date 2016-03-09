/** \file oatest.cpp

 C++ program: oatest

 oatest: tool for testing new algorithms

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "arraytools.h"
#include "arrayproperties.h"
#include "tools.h"
#include "extend.h"


/// Structure representing the type of conference designs
class conference_t
{
public:
	rowindex_t N;	/** number of runs */
	colindex_t ncols;	/** total number of columns (factors) in the design */

public:
	/// create new conference_t object
	conference_t ( int N, int k );

	array_link create_root ( ) const;

	std::string __repr__() const
	{
	return printfstring("conference type: N %d, ncols %d", this->N, this->ncols);
	}
	
};

/// reduce conference matrix to normal form
	array_link reduceConference(const array_link &, int verbose = 0);
	


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

	arraylist_t getarrays ( const array_link al ) {
		arraylist_t ll;

		for ( size_t i=0; i<this->extensions.size(); i++ ) {
			array_link alx = hstack ( al, extensions[i] );
			ll.push_back ( alx );
		}
		return ll;
	}
};


/** Extend a list of conference designs with a single column.
 *
 */
conference_extend_t extend_conference_matrix ( const array_link &al, const conference_t &ct, int extcol, int verbose=1 );

arraylist_t extend_conference ( const arraylist_t &lst, const conference_t ctype, int verbose );


/// select representatives for the isomorphism classes of a list of conference arrays
arraylist_t  selectConferenceIsomorpismClasses(const arraylist_t list, int verbose);

/// select representatives for the isomorphism classes of a list of conference arrays, return indices of classes
std::vector<int> selectConferenceIsomorpismIndices(const arraylist_t lst, int verbose);

/** return max position of zero in array, returns -1 if no zero is found
 * 
 * The parameter k specifies the column to search in. For k=-1 all columns are searched.
 */
int maxz(const array_link &al, int k = -1);

	
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
