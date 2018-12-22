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
enum matrix_isomorphism_t { 
  /// isotopy: permute rows, columns and symbols
  ISOTOPY,
  /// permute rows and columns
  MATRIX_ISOMORPHISM,
  /// permute rows, columns and to row and column negations (values in 0, +1, -1)
  CONFERENCE_ISOMORPHISM,
  /// permutations of rows, columns and column symbol permutations
  OA_ISOMORPHISM };

/// isomorphism type for column and row permtations and column permutations
const matrix_isomorphism_t CONFERENCE_RESTRICTED_ISOMORPHISM = OA_ISOMORPHISM;

namespace nauty {
#include "nauty.h"

/** Reduce a colored graph to Nauty minimal form
 *
 * The transformation returned is from the normal form to the specified graph.
 *
 * \param graph Graph in incidence matrix form
 * \param colors Colors of the graph nodes
 * \param verbose Verbosity level
 * \return Relabelling of the graph vertices
 *
 */
std::vector< int > reduceNauty (const array_link &graph, std::vector< int > colors, int verbose = 0);

} 

/// Apply a vertex permutation to a graph
array_link transformGraph (const array_link &graph, const std::vector< int > vertex_permutation, int verbose = 1);

/// Reduce an orthogonal array to Nauty minimal form. the array transformation is returned
array_transformation_t reduceOAnauty (const array_link &array, int verbose = 0);

/// Reduce an orthogonal array to Nauty minimal form. the array transformation is returned
array_transformation_t reduceOAnauty (const array_link &array, int verbose, const arraydata_t &arrayclass);

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair< array_link, std::vector< int > > array2graph (const array_link &array, int verbose = 1);

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
*/
std::pair< array_link, std::vector< int > > array2graph (const array_link &array, int verbose, const arraydata_t &arrayclass);

/// From a relabelling of the graph return the corresponding array transformation
array_transformation_t oagraph2transformation (const std::vector< int > &pp, const arraydata_t &arrayclass,
                                               int verbose = 1);

