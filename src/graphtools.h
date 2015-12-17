/** \file graphtools.h

\brief This file contains definitions and functions related to graphs and designs.


 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See COPYING file that comes with this distribution
*/


#pragma once


/* Interface to Nauty code
 *
 */

namespace nauty
{
#include "nauty.h"
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */


	
/// reduce graph to Nauty minimal form
std::vector<int> reduceNauty ( const array_link &G, std::vector<int> colors );


/// reduce design to Nauty minimal form
std::vector<int> reduceOAnauty(const array_link &al, int verbose=1);
	
	
} // end of nauty namespace

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair<array_link, std::vector<int> >  array2graph ( const array_link &al, int verbose=1 );
