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

};


/// Helper structure
struct conference_extend_t {
	std::vector<cperm> first;
	std::vector<cperm> second;
	std::vector<cperm> extensions;

public:

	cperm combine ( int i, int j ) const {
		cperm c =vstack ( this->first[i], this->second[j] );

		return c;
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
conference_extend_t extend_conference ( const array_link al, const conference_t ct, int extcol, int verbose=1 );


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
