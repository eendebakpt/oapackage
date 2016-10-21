/** \file graphtools.h

\brief This file contains definitions and functions related to graphs and designs.


 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2016

 Copyright: See COPYING file that comes with this distribution
*/


#pragma once

#include <vector>
#include "arraytools.h"


/** Isomorphism types for matrices
 *
 * Isotopy: permute rows, columns and symbols
 * Matrix isomorphism: permute rows and columns
 * Conference isomorphism: permute rows, columns and to row and column negations (values in 0, +1, -1)
 * Orthogonal array isomorphism: permutations of rows, columns and column symbol permutations
 */
enum matrix_isomorphism_t {ISOTOPY, MATRIX_ISOMORPHISM, CONFERENCE_ISOMORPHISM, OA_ISOMORPHISM};


/// isomorphism type for column and row permtations and column permutations
const matrix_isomorphism_t CONFERENCE_RESTRICTED_ISOMORPHISM = OA_ISOMORPHISM;


/* Interface to Nauty code
 *
 */

namespace nauty
{
#include "nauty.h"
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

/** reduce a colored graph to Nauty minimal form
 *
 * The transformation returned is from the normal form to the specified graph.
 *
 */
std::vector<int> reduceNauty ( const array_link &G, std::vector<int> colors, int verbose=0 );

} // end of nauty namespace

/// apply a vertex permutation to a graph
array_link transformGraph ( const array_link &G, const std::vector<int> tr, int verbose = 1 );

/// reduce an orthogonal array to Nauty minimal form. the array transformation is returned
array_transformation_t reduceOAnauty(const array_link &al, int verbose=0);

array_transformation_t reduceOAnauty(const array_link &al, int verbose, const arraydata_t &ad);

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair<array_link, std::vector<int> >  array2graph ( const array_link &al, int verbose=1 );

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair<array_link, std::vector<int> >  array2graph ( const array_link &al, int verbose, const arraydata_t &ad );

/// From a relabelling of the graph return the corresponding array transformation
array_transformation_t oagraph2transformation ( const std::vector<int> &pp, const arraydata_t &arrayclass, int verbose=1 );

int unittest_nautynormalform(const array_link &al, int verbose=1);


