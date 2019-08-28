/** \file conference.cpp

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <algorithm>
#include <iostream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "arrayproperties.h"
#include "arraytools.h"
#include "extend.h"
#include "graphtools.h"
#include "conference.h"


void print_column(const conference_column &c, const char *msg) {
	if (msg != 0)
		myprintf("%s: ", msg);
	for (size_t i = 0; i < c.size(); i++) {
		myprintf("%3d", c[i]);
	}
	myprintf("\n");
}

/// return true of the argument is even
inline bool is_even (int q) { return (q % 2) == 0; }

/** return structure parameters of a conference design
 *
 * When generating candidate extension columns many properties of the column are
 * determined by N and the position of the zero in the column. (If the first 2 columns are in normal form).
 *
 * This method returns the sizes q1 and q2 of the second and third block. For details see "A Classification Criterion for Definitive Screening Designs", Schoen et al.
 */

void getConferenceNumbers (int N, int k, int &q, int &q1, int &q2, int &v) {
        q = (N - 2) / 2;

        if (k < 2 + q) {
                if (is_even (q)) {
                        q1 = q / 2 - 1;
                        v = 1;
                        q2 = q1 + 1;
                } else {
                        q1 = (q - 1) / 2;
                        v = -1;
                        q2 = q1 + 1;
                }
        } else {

                if (is_even (q)) {
                        q1 = q / 2;
                        v = -1;
                        q2 = q1;
                } else {
                        q1 = (q - 1) / 2;
                        v = 1;
                        q2 = q1;
                }
        }
}

array_link conference2DSD(const array_link &conf, bool add_zeros)
{
	array_link dsd(2 * conf.n_rows + add_zeros, conf.n_columns, array_link::INDEX_DEFAULT);
	
	for(int row=0; row<conf.n_rows; row++) {
		for (int column = 0; column < conf.n_columns; column++) {
			dsd.atfast(row, column) = conf.atfast(row, column);
			dsd.atfast(conf.n_rows+row, column) = -conf.atfast(row, column);
		}
	}
	return dsd;
}

/// show a list of candidate extensions
void showCandidates (const std::vector< conference_column > &column_candidates) {
        for (size_t i = 0; i < column_candidates.size (); i++) {
                myprintf ("%d: ", (int)i);
                print_column (column_candidates[i]);
        }
}

conference_t::conference_t (const conference_t &rhs) {
        this->N = rhs.N;
        this->ncols = rhs.ncols;
        this->ctype = rhs.ctype;
        this->itype = rhs.itype;

        this->j1zero = rhs.j1zero;
        this->j3zero = rhs.j3zero;
}

conference_t::conference_t () {
        this->N = -1;
        this->ncols = -1;
        this->ctype = CONFERENCE_NORMAL;
        this->itype = MATRIX_ISOMORPHISM;
        this->j1zero = 0;
        this->j3zero = 0;
}

conference_t::conference_t (int N, int k, int _j1zero) {

        this->N = N;
        this->ncols = k;
        this->ctype = CONFERENCE_NORMAL;
        this->itype = CONFERENCE_ISOMORPHISM;
        this->j1zero = _j1zero;
        this->j3zero = 0;
}

array_link conference_t::create_root_three_columns () const {
        array_link array (this->N, 3, 0); 

        array.at (0, 0) = 0;
        for (int i = 1; i < this->N; i++) {
                array.at (i, 0) = 1;
        }
        for (int i = 0; i < this->N; i++) {
                if (i == 1) {
                        array.at (1, 1) = 0;
                        continue;
                }
                if (i <= this->N / 2)
                        array.at (i, 1) = 1;
                else
                        array.at (i, 1) = -1;
        }

        int q, q1, q2, v;
        const int k = 2;
        const int N = this->N;
        getConferenceNumbers (this->N, k, q, q1, q2, v);

        for (int i = 3; i < N; i++)
                array.at (i, 2) = -1;
        array.at (0, 2) = 1;
        array.at (1, 2) = v;
        array.at (2, 2) = 0;
        for (int i = 0; i < q1; i++)
                array.at (2 + 1 + i, 2) = 1;
        for (int i = 0; i < q2; i++)
                array.at (2 + q + i, 2) = 1;

        return array;
}
arraylist_t conference_t::createRootArrays () const {
	arraylist_t root_arrays;
	switch (this->ctype) {
        case CONFERENCE_NORMAL:
        case CONFERENCE_DIAGONAL:
                switch (this->itype) {
                case CONFERENCE_ISOMORPHISM:
                        root_arrays.push_back (this->create_root ());
                        break;
                case CONFERENCE_RESTRICTED_ISOMORPHISM: {
                        array_link al (this->N, 1, array_link::INDEX_DEFAULT);
                        for (int i = N / 2 + 1; i <= N; i++) {
                                al.setconstant (-1);
                                al.at (0, 0) = 0;
                                for (int k = 1; k < i; k++) {
                                        al.at (k, 0) = 1;
                                }
                                root_arrays.push_back (al);
                        }
                } break;
                default:
                        printfd ("not implemented (itype %d)\n", this->itype);
                }
                break; //
        case DCONFERENCE: {
                switch (this->itype) {
                case CONFERENCE_RESTRICTED_ISOMORPHISM: {
                        arraylist_t tmp = this->createDoubleConferenceRootArrays ();
                        root_arrays.insert (root_arrays.end (), tmp.begin (), tmp.end ());
                } break;
                case CONFERENCE_ISOMORPHISM: {
                        myassert (this->j1zero == 0, "j1zero should be zero for CONFERENCE_ISOMORPHISM type");
						myassert(this->j3zero == 0, "j3zero should be zero for CONFERENCE_ISOMORPHISM type");
                        arraylist_t tmp = this->createDoubleConferenceRootArrays ();
                        root_arrays.insert (root_arrays.end (), tmp.begin (), tmp.end ());
                } break;
                default:
			throw_runtime_exception(printfstring("root array geneated not implemented for class with itype=%d", this->itype));
                }
        }
        }
	return root_arrays;
}

arraylist_t conference_t::createDoubleConferenceRootArrays () const {
        arraylist_t lst;
        array_link al (this->N, 1, array_link::INDEX_DEFAULT);
        if (j1zero) {
                al.setconstant (-1);
                al.at (0, 0) = 0;
                al.at (1, 0) = 0;
                for (int row = 2; row < N / 2 + 1; row++) {
                        al.at (row, 0) = 1;
                }
                lst.push_back (al);
        } else {
                for (int i = N / 2 + 1; i <= N; i++) {
                        al.setconstant (-1);
                        al.at (0, 0) = 0;
                        al.at (1, 0) = 0;
                        for (int row = 2; row < i; row++) {
                                al.at (row, 0) = 1;
                        }
                        lst.push_back (al);
                }
        }
        return lst;
}

std::string conference_t::idstr () const {

        std::string s = printfstring ("conf-N%d-k%d-j-%d-%d-%d", this->N, this->ncols, this->j1zero, 1, this->j3zero);
        return s;
}

array_link conference_t::create_root () const {
        array_link al (this->N, 2, 0);

        al.at (0, 0) = 0;
        for (int i = 1; i < this->N; i++) {
                al.at (i, 0) = 1;
        }
        for (int i = 0; i < this->N; i++) {
                if (i == 1) {
                        al.at (1, 1) = 0;
                        continue;
                }
                if (i <= this->N / 2)
                        al.at (i, 1) = 1;
                else
                        al.at (i, 1) = -1;
        }

        return al;
}

/** Helper structure containing extensions of conference designs
*/
struct conference_extend_t {
	std::vector< conference_column > first;      /// list of first block candidate extensions
	std::vector< conference_column > second;     /// list of first block candidate extensions
	std::vector< conference_column > extensions; /// list of candidate extensions

public:
	// combine first and second section into a single column
	conference_column combine(int i, int j) const {
		conference_column c = vstack(this->first[i], this->second[j]);
		return c;
	}

	size_t nExtensions() const { return this->extensions.size(); }

	/// return the set of extension arrays
	arraylist_t getarrays(const array_link al) const {
		arraylist_t ll;

		for (size_t i = 0; i < this->extensions.size(); i++) {
			array_link alx = hstack(al, extensions[i]);
			ll.push_back(alx);
		}
		return ll;
	}
};

std::vector<int> double_conference_foldover_permutation(const array_link &double_conference) {
        int N = double_conference.n_rows/2;
        
        if (! double_conference.is_conference() ) {
         throw_runtime_exception("input array should be conference design");    
             
        }
        
        if (2*N != double_conference.n_rows) {
         throw_runtime_exception("double conference design should have even number of rows");    
        }
        
        array_link alt = double_conference.transposed ();
        array_link alt_minus = alt * -1;

        std::vector< int > ri (double_conference.n_rows);
        std::fill (ri.begin (), ri.end (), -1);
        std::vector< int > foldover_permutation (double_conference.n_rows);
        foldover_permutation[0]=-1;
        
        for (int row = 0; row < double_conference.n_rows; row++) {
                if (ri[row] > -1)
                        continue;
                int foundcol = 0;
                for (int j = row + 1; j < double_conference.n_rows; j++) {
                        if (ri[j] > -1)
                                continue;
                        if (alt.columnEqual (row, alt_minus, j)) {
                                foundcol = 1;
                                ri[row] = j;
                                ri[j] = row;
                                break;
                        }
                }
                if (!foundcol) {
                        return foldover_permutation;
                }
        }
        int row_index=0;
        for (int row = 0; row < double_conference.n_rows; row++) {
             if (row<ri[row]) {
               foldover_permutation[row_index] = row;
               foldover_permutation[row_index+N] = ri[row];
               row_index++;
             }
        }
        
        return foldover_permutation;
 
     
}
bool isConferenceFoldover (const array_link &array, int verbose) {
        array_link alt = array.transposed ();
        array_link alt_minus = alt * -1;

        std::vector< int > ri (array.n_rows);
        std::fill (ri.begin (), ri.end (), -1);

        for (int i = 0; i < array.n_rows; i++) {
                if (ri[i] > -1)
                        continue;
                int foundcol = 0;
                for (int j = i + 1; j < array.n_rows; j++) {
                        if (ri[j] > -1)
                                continue;
                        if (alt.columnEqual (i, alt_minus, j)) {
                                foundcol = 1;
                                ri[i] = j;
                                ri[j] = i;
                                break;
                        }
                }
                if (!foundcol) {
                        if (verbose) {
                                printf ("isConferenceFoldover: no foldover for row %d\n", i);
                        }
                        return false;
                }
        }
        return true;
}

/// reduce double conference matrix to normal form (work in progress)
conference_transformation_t reduceDoubleConferenceTransformation (const array_link &array, int verbose) {
        if (!array.is_conference (2)) {
                myprintf ("reduceConferenceTransformation: error: design is not a double conference design\n");
                conference_transformation_t t (array);
                return t;
        }

        arraydata_t arrayclass (3, array.n_rows, 1, array.n_columns);
        array_transformation_t at = reduceOAnauty (array + 1, verbose >= 2, arrayclass);

        conference_transformation_t t (array);
        t.rperm = at.rowperm ();
        t.cperm = at.colperm ();

        for (int c = 0; c < array.n_columns; c++) {
                std::vector< int > lp = at.lvlperm (c);
                myassert (lp[1] == 1, "error in reduction");                  // 0 should go to 0
                t.cswitch[c] = (lp[0] == 0) ? 1 : -1; 
        }

        return t;
}

/** Convert a conference design to a colored graph representation
 *
 * The conversion is such that conferrence design isomorphisms correspond to isomorphisms of the colored graph 
 * (e.g. permutations of the nodes respecting the edges and colors).
 *
 * The colored graph consists of N=(2*number_of_rows + 2*number_of_columns). The first section of size 2*number_of_rows encodes +1, -1 value
 * in the rows (the zeros are implicitly encoded) and the second section of size 2*number_of_columns corresponds to the columns.
 *
 *	Add edges as follows:
 *	(1)  r[i]--r'[i] for i=0..nr-1;  c[j]--c'[j] for j=0..nc-1.
 *	(2)  r[i]--c[j] and r'[i]--c'[j] for all A[i,j] = +1
 *	(3)  r[i]--c'[j] and r'[i]--c[j] for all A[i,j] = -1.
 *	Zeros in A don't cause any edges.
 * 
 * The first 2*number_of_rows are colored with 0, the remaining nodes with color 1.
 *
 * \param al Conference design
 * \param Verbose Vebosity level
 * \returns Pair of the graph in edge adjacency  matrix representation and the colors
 */
std::pair<array_link, std::vector<int> > conference_design2colored_graph(const array_link &conference_design, int verbose) {
	const int nr = conference_design.n_rows;
	const int nc = conference_design.n_columns;
	const int number_vertices = 2 * (nr + nc);
	/// create graph
	array_link G(number_vertices, number_vertices, array_link::INDEX_DEFAULT);
	G.setconstant(0);

	if (verbose)
		myprintf("reduceConferenceTransformation: reduce design with %d rows, %d columns\n", nr, nc);

	std::vector< int > colors(2 * (nr + nc));

	const int roffset0 = 0;
	const int roffset1 = nr;
	const int coffset0 = 2 * nr;
	const int coffset1 = 2 * nr + nc;

	// set colors
	for (int i = 0; i < coffset0; i++)
		colors[i] = 0;
	for (int i = coffset0; i < coffset0 + 2 * nc; i++)
		colors[i] = 1;

	// (1)
	for (int i = 0; i < nr; i++)
		G.atfast(roffset0 + i, i + roffset1) = 1;
	for (int i = 0; i < nc; i++)
		G.atfast(coffset0 + i, i + coffset1) = 1;

	// (2), (3)
	for (int c = 0; c < nc; c++) {
		for (int r = 0; r < nr; r++) {
			if (conference_design.atfast(r, c) == 1) {
				G.atfast(roffset0 + r, coffset0 + c) = 1;
				G.atfast(roffset1 + r, coffset1 + c) = 1;
			}
			else {
				if (conference_design.atfast(r, c) == -1) {
					G.atfast(roffset0 + r, coffset1 + c) = 1;
					G.atfast(roffset1 + r, coffset0 + c) = 1;
				}
			}
		}
	}

	// make array symmetryic
	const int nrg = G.n_rows;
	for (int i = 0; i < number_vertices; i++) {

		array_t *x = G.array + i * nrg; // offset to column
		for (int j = i; j < number_vertices; j++) {
			x[j] = G.array[i + j * nrg];
		}
	}

	if (verbose >= 3) {
		printf("reduceConference: incidence graph:\n");
		printf("   2x%d=%d row vertices and 2x%d=%d column vertices\n", nr, 2 * nr, nc, 2 * nc);
		G.showarray();
	}

	return std::pair<array_link, std::vector<int> >(G, colors);
}

/// From a graph tranformation calculate the corresponding tranformation of a conference design
conference_transformation_t graph_transformation2conference_transformation(int nrows, int ncolumns, const array_link &G, const std::vector< int > vertex_permutation, int verbose)
{
	const int number_vertices = 2 * (nrows + ncolumns);
	myassert(number_vertices == G.n_rows, "conference design specification does not match graph size");

	// extract transformation
	const int roffset0 = 0;
	const int roffset1 = nrows;
	const int coffset0 = 2 * nrows;
	const int coffset1 = 2 * nrows + ncolumns;

	if (verbose >= 2) {
		if (verbose >= 3) {
			array_link Gx = transformGraph(G, vertex_permutation, 0);
			printfd("transformed graph\n");
			Gx.showarray();
		}
		std::vector< int > tr1 = std::vector< int >(vertex_permutation.begin(), vertex_permutation.begin() + 2 * nrows);
		std::vector< int > tr2 = std::vector< int >(vertex_permutation.begin() + coffset0, vertex_permutation.end());
		printf("  row vertex transformations: ");
		display_vector(tr1);
		printf("\n");
		printf("  col vertex transformations: ");
		display_vector(tr2);
		printf("\n");
	}

	// define conference matrix transformation object....
	conference_transformation_t transformation(nrows, ncolumns);

	// extract transformation
	std::vector< int > rr(nrows);
	for (int i = 0; i < nrows; i++) {
		rr[i] = std::min(vertex_permutation[i], vertex_permutation[i + nrows]);
	}

	if (verbose >= 2) {
		printf("rr: ");
		print_perm(rr);
	}

	transformation.rperm = invert_permutation(argsort(rr));

	for (int i = 0; i < nrows; i++) {
		transformation.rswitch[transformation.rperm[i]] = 2 * (vertex_permutation[i] < vertex_permutation[i + nrows]) - 1;
	}

	std::vector< int > cc(ncolumns);
	for (int i = 0; i < ncolumns; i++) {
		cc[i] = std::min(vertex_permutation[coffset0 + i], vertex_permutation[coffset0 + i + ncolumns]);
	}
	transformation.cperm = invert_permutation(argsort(cc));

	for (int i = 0; i < ncolumns; i++) {
		transformation.cswitch[transformation.cperm[i]] = 2 * (vertex_permutation[coffset0 + i] < vertex_permutation[coffset0 + i + ncolumns]) - 1;
	}

	if (verbose >= 2) {
		printf("transform: \n");
		transformation.show();
	}
	return transformation;
}

conference_transformation_t reduceConferenceTransformation (const array_link &al, int verbose) {

        if (!al.is_conference (1)) {
                myprintf ("reduceConferenceTransformation: error: design is not a conference design\n");
                conference_transformation_t t (al);
                return t;
        }

		std::pair<array_link, std::vector<int> > colored_graph = conference_design2colored_graph(al, verbose);
		array_link G = colored_graph.first;
		std::vector<int> colors = colored_graph.second;

        const std::vector< int > vertex_permutation = nauty::reduceNauty (G, colors, verbose >= 2);
        const std::vector< int > inverse_vertex_permutation = invert_permutation (vertex_permutation);

		conference_transformation_t t = graph_transformation2conference_transformation(al.n_rows, al.n_columns, G, inverse_vertex_permutation, verbose);
        return t;
}

array_link reduceConference (const array_link &al, int verbose) {
        conference_transformation_t t = reduceConferenceTransformation (al, verbose);
        array_link alx = t.apply (al);
        return alx;
}

// return vector of length n with specified positions set to one
conference_column get_comb (const conference_column p, int n, int zero = 0, int one = 1) {
        conference_column c (n, zero);
        for (size_t i = 0; i < p.size (); i++)
                c[p[i]] = one;
        return c;
}

// set vector of length n with specified positions set to one
void set_comb (const conference_column &positions, conference_column &c, int n, int zero = 0, int one = 1) {
        std::fill (c.begin (), c.begin () + n, zero);
        for (size_t i = 0; i < positions.size (); i++)
                c[positions[i]] = one;
}

inline void get_comb (const conference_column &p, int n, int zero, int one, conference_column &c) {
        for (int i = 0; i < n; i++)
                c[i] = zero;
        for (size_t i = 0; i < p.size (); i++)
                c[p[i]] = one;
}

/// return copy of vector with zero inserted at specified position
inline conference_column insertzero (const conference_column &c, int pos, int value = 0) {
        conference_column cx (c.size () + 1);
        std::copy (c.begin (), c.begin () + pos, cx.begin ());
        cx[pos] = value;
        std::copy (c.begin () + pos, c.end (), cx.begin () + pos + 1);
        return cx;
}

/// set vector with zero inserted at specified position, no error checking
void insertzero (conference_column &target, const conference_column &c, int pos, int value = 0) {
        std::copy (c.begin (), c.begin () + pos, target.begin ());
        target[pos] = value;
        std::copy (c.begin () + pos, c.end (), target.begin () + pos + 1);
}

/** Return all admissible columns (first part) for a conference array in normal form
 *
 *
 **/
std::vector< conference_column > get_first (int N, int extcol, int verbose = 1) {
        int k1 = -1;
        int n1 = -1;
        int k = extcol;

        int q, q1, q2, v;
        getConferenceNumbers (N, k, q, q1, q2, v);

        int haszero = extcol < q + 2;
        if (haszero) {

                n1 = q - 1;
        } else {
                n1 = q;
        }

        conference_column c (q1);
        for (int i = 0; i < q1; i++)
                c[i] = i;

        int nc = ncombs< long > (n1, q1);
        if (verbose >= 2)
                printf ("get_first: conference array: extcol %d: N %d, n1 %d, q %d, v %d, q1 %d, q2 %d, nc %d\n",
                        extcol, N, n1, q, v, q1, q2, nc);

        std::vector< conference_column > ff;
        for (long j = 0; j < nc; j++) {
                conference_column cc = get_comb (c, n1, -1, 1);

                cc = insertzero (cc, 0, 1);
                cc = insertzero (cc, 1, v);

                if (haszero)
                        cc = insertzero (cc, extcol);
                ff.push_back (cc);

                if (j + 1 < nc) {
                        next_comb (c, q1, n1);
                }
        }

        return ff;
}

/** Return all admissible columns (block two) for a conference array in normal form */
std::vector< conference_column > get_second (int N, int extcol, int target, int verbose = 0) {
        // verbose=2;
        if (verbose)
                printfd ("get_second: N %d, extcol %d, target %d\n");
        int k = extcol;
        int q, q1, q2, v;
        getConferenceNumbers (N, k, q, q1, q2, v);

        int n1 = -1;
        int haszero = extcol >= q + 2;
        if (verbose)
                printf ("get_second: extcol: %d, q %d, haszero %d\n", extcol, q, haszero);

        if (haszero) {
                n1 = q - 1;
        } else {
                n1 = q;
        }
        int qx = q2;

        conference_column c (qx);
        for (int i = 0; i < qx; i++)
                c[i] = i;

        int nc = ncombs< long > (n1, qx);
        if (verbose)
                printf ("get_second: N %d, n1 %d, qx %d, target %d, nc %d\n", N, n1, qx, target, nc);
        std::vector< conference_column > ff;
        conference_column ccx (n1);
        conference_column cc (n1 + haszero);
        for (long j = 0; j < nc; j++) {
                if (haszero) {
                        get_comb (c, n1, -1, 1, ccx);
                        insertzero (cc, ccx, extcol - (q + 2));
                } else {
                        set_comb (c, cc, n1, -1, 1);
                }
                if (verbose >= 2) {
                        printfd ("add element of size %d =   %d\n", cc.size (), q);
                        display_vector (cc);
                        printf ("\n");
                        printf ("current_column: ");
                        display_vector (c);
                        printf ("\n");
                }

                ff.push_back (cc);
                if (n1 > 0) // guard
                        next_comb (c, qx, n1);
        }

        return ff;
}

/// calculate inner product between partial two permutations
int partial_inner_product (const conference_column &a, const array_link &al, int column_idx, int rmax) {
        int ip = 0;
        size_t nn = a.size ();
        const array_t *b = al.array + column_idx * al.n_rows;

        for (int row = 0; row < rmax; row++) {
                ip += a[row] * b[row];
        }
        return ip;
}

/// calculate inner product between two permutations
int innerprod (const conference_column &a, const array_link &al, int col) {
        int ip = 0;
        size_t nn = a.size ();
        const array_t *b = al.array + col * al.n_rows;

        for (size_t i = 0; i < nn; i++) {
                ip += a[i] * b[i];
        }
        return ip;
}

/// calculate inner product between two columns
int innerprod (const conference_column &a, const conference_column &b) {
        int ip = 0;
        size_t nn = b.size ();
        for (size_t i = 0; i < nn; i++) {
                ip += a[i] * b[i];
        }
        return ip;
}

int check_symm_zero (const conference_column &c, const std::vector< int > &check_indices, int i) {
        if (check_indices[i]) {
                if (((unsigned char)c[i - 1]) > ((unsigned char)c[i])) {
                        // discard
                        return false;
                }
        }
        return true;
}

/// filter list of columns, only return columns with zero at specified position
inline std::vector< conference_column > filterZeroPosition (const std::vector< conference_column > &lst, int zero_position) {
        std::vector< conference_column > out;
        for (size_t i = 0; i < lst.size (); i++) {
                if (lst[i][zero_position] == 0) {
                        out.push_back (lst[i]);
                }
        }
        return out;
}

/** Return true if a candidate extensions satisfies the symmetry test
 *
 * The method checks the candidate at positions given by check_indices and compares the element to the next position
 *
 */
int satisfy_lmc0_symmetry(const conference_column &c, const std::vector< int > &check_indices, int rowstart, int rowend) {

	for (int i = rowstart + 1; i < rowend; i++) {
		if (check_indices[i]) {
			if (((unsigned char)c[i - 1]) > ((unsigned char)c[i])) {
				// discard
				return false;
			}
		}
	}
	// accept
	return true;
}

/** Return true if a candidate extensions satisfies the symmetry test
 *
 * The method checks the candidate at positions given by check_indices and compares the element to the next position
 *
 */
int satisfy_lmc0_symmetry (const conference_column &column, const std::vector< int > &check_indices, int rowstart = 2) {
	return satisfy_lmc0_symmetry(column, check_indices, rowstart, column.size() - 1);
}

/** helper function, return true if a candidate extensions satisfies the symmetry test
 *
 * The candidate extension fails the symmetry test if a simple row permutation of two succesive rows leads to a smaller
 * design.
 *
 */
int satisfy_lmc0_symmetry (const conference_column &c, const symmdata &sd, int rowstart = 2) {
        const int verbose = 0;

        if (verbose) {
                printf ("satisfy_lmc0_symmetry: ");
                display_vector (c);
                printf ("\n");
				if (verbose >= 2) {
					printf("satisfy_lmc0_symmetry: sd: ");
					sd.show();
				}
		}

		const std::vector< int > check_indices = sd.checkIdx();

		bool return_value = satisfy_lmc0_symmetry(c, check_indices, rowstart);
        if (verbose >= 2) {
                printf ("satisfy_lmc0_symmetry: return %d\n", return_value);
        }
        return return_value;
}

/// return column of an array in conference_column format
conference_column getColumn (const array_link &al, int c) {
        conference_column cx (al.n_rows);
        std::copy (al.array + c * al.n_rows, al.array + (c + 1) * al.n_rows, cx.begin ());
        return cx;
}

// return true if the extension column satisfies the inner product check
int ipcheck (const conference_column &col, const array_link &al, int cstart = 2, int verbose = 0) {
        for (int c = cstart; c < al.n_columns; c++) {
                if (innerprod (col, al, c) != 0) {
                        if (verbose) {
                                printf ("ipcheck: column %d to %d (inclusive), failed at col %d\n", c, cstart,
                                        al.n_columns + 1);
                        }
                        return false;
                }
        }
        return true;
}

int minz (const array_link &al, int column_idx) {
        const int N = al.n_rows;
        int minzidx = -1;
        const int nr = al.n_rows;
        if (column_idx == -1) {
                for (int k = 0; k < al.n_columns; k++) {
                        for (int r = 0; r < N; r++) {
                                if (al._at (r, k) == 0) {
                                        minzidx = std::min (minzidx, r);
                                }
                        }
                }
                return minzidx;
        } else {
                for (int r = 0; r < N; r++) {
                        if (al._at (r, column_idx) == 0) {
                                return r;
                        }
                }
        }
        return minzidx;
}

int maxz (const array_link &al, int column_index) {
        int maxzidx = -1;
        const int nr = al.n_rows;
        if (column_index == -1) {
                for (int k = 0; k < al.n_columns; k++) {
                        for (int r = nr - 1; r >= maxzidx; r--) {
                                if (al._at (r, k) == 0) {
                                        maxzidx = std::max (maxzidx, r);
                                }
                        }
                }
                return maxzidx;
        } else {
                for (int r = nr - 1; r >= 0; r--) {
                        if (al._at (r, column_index) == 0) {
                                return r;
                        }
                }
        }
        return maxzidx;
}

/** Filter candidate extensions based on symmetry propery
 *
 * The subblock S obtained by deleting the first row or column should be symmetric or anti-symmetric.
 * See "Conference matrices", Peter J. Cameron, CSG, 7 October 2011
 *
 */
std::vector< conference_column > filterCandidatesSymm (const std::vector< conference_column > &extensions, const array_link &als,
                                           int verbose) {
        const int N = als.n_rows;

        int mval = 0;
        if (N % 4 == 2) {
                // S is symmetric
                mval = 1;
        } else {
                // else S is anti-symmetric
                mval = -1;
        }
        if (verbose >= 2)
                printf ("N %d, mval %d\n", N, mval);

        const int k = als.n_columns;
        conference_column comparerow (N);
        for (int i = 0; i < k; i++)
                comparerow[i] = als.at (k, i);

        std::vector< conference_column > e (0);
        for (size_t i = 0; i < extensions.size (); i++) {
                const conference_column &ex = extensions[i];

                int good = 1;
                for (int x = 2; x < k; x++) {
                        if (comparerow[x] != mval * ex[x]) {
                                good = 0;
                                break;
                        }
                }
                if (good)
                        e.push_back (extensions[i]);
        }
        if (verbose >= 1)
                printf ("filterCandidatesSymm: k %d, filter %d/%d\n", k, (int)e.size (), (int)extensions.size ());
        return e;
}

/// create J2 table as intermediate result for J-characteristic calculations for conference matrices
array_link createJ2tableConference(const array_link &confmatrix) {
	const int nr = confmatrix.n_rows;
	const int nc = (confmatrix.n_columns + 1) * confmatrix.n_columns / 2;
	// fill double column table
	array_link dtable(nr, nc, -1);

	// loop over all column pairs
	int idx = 0;
	for (int i = 0; i < confmatrix.n_columns; i++) {
		for (int j = 0; j <= i; j++) {
			// loop over all rows of original array
			const array_t *p1 = confmatrix.array + confmatrix.n_rows * i;
			const array_t *p2 = confmatrix.array + confmatrix.n_rows * j;
			array_t *pout = dtable.array + idx * dtable.n_rows;
			for (int x = 0; x < nr; x++) {
				pout[x] = p1[x] * p2[x];
			}
			idx++;
		}
	}

	return dtable;
}

/// filter candidate extensions on J3 value (only pairs with extension candidate are checked)
std::vector< conference_column > filterJ3 (const std::vector< conference_column > &extensions, const array_link &als, int verbose) {

        const int N = als.n_rows;

        array_link dtable = createJ2tableConference (als);

        if (verbose >= 2) {
                printf ("array + dtable\n");
                als.showarray ();
                printfd ("dtable:\n");
                dtable.showarray ();
        }

        int nc = dtable.n_columns;

        if (verbose >= 1) {
                printf ("filterJ3: array %dx%d, nc %d\n", N, als.n_columns, nc);
        }

        std::vector< conference_column > e2 (0);
        for (size_t i = 0; i < extensions.size (); i++) {
                const conference_column &c = extensions[i];

                int jv = 0;
                for (int idx1 = 0; idx1 < nc; idx1++) {
                        jv = 0;

                        const array_t *o1 = dtable.array + dtable.n_rows * idx1;
                        for (int xr = 0; xr < N; xr++) {
                                jv += (o1[xr]) * (c[xr]);
                        }

                        if (jv != 0)
                                break;
                }

                if (jv == 0)
                        e2.push_back (c);
        }

        if (verbose >= 1) {
                printf ("filterJ3: %ld -> %ld extensions\n", (long)extensions.size (), (long)e2.size ());
        }
        return e2;
}

/** filter conferece matrix extension candidates
 *
 * Filtering is based in symmetry and ip
 */
std::vector< conference_column > filterDconferenceCandidates (const std::vector< conference_column > &extensions, const array_link &als,
                                                  int filtersymm, int filterip, int verbose) {
        //	symmetry_group rs = als.row_symmetry_group();
        symmdata sd (als);
        DconferenceFilter dfilter (als, filtersymm, filterip);
        dfilter.filterfirst = 1;

        if (verbose >= 2)
                sd.show (1);

        std::vector< conference_column > e2 (0);
        for (size_t i = 0; i < extensions.size (); i++) {

                if (dfilter.filter (extensions[i])) {
                        e2.push_back (extensions[i]);
                }
        }
        return e2;
}

bool DconferenceFilter::filterJpartial (const conference_column &c, int maxrow) const {
        const int N = als.n_rows;
        long j = partial_inner_product (c, this->als, als.n_columns - 1, maxrow);
        if (std::abs (j) > (N - maxrow)) {
                return false;
        } else {
                return true;
        }
}

DconferenceFilter::DconferenceFilter(const array_link &_als, int filter_symmetry_, int filterj2_, int filterj3_ )
	: als(_als), filtersymm(filter_symmetry_), filterj2(filterj2_), filterj3(filterj3_), filterfirst(0),
	filterzero(0), ngood(0), sd(als) {

	check_indices = sd.checkIdx();
	dtable = createJ2tableConference(als);

	if (als.n_columns >= 2) {
		inline_dtable = als.selectColumns(0) - als.selectColumns(1);
		inline_dtable = hstack(inline_dtable, als.selectColumns(0) + 1);
		inline_dtable = hstack(inline_dtable, als.selectColumns(0) * als.selectColumns(0) - 1);
		inline_dtable = hstack(inline_dtable, als.selectColumns(1) * als.selectColumns(1) - 1);

		minzvalue = minz(als, als.n_columns - 1);

		inline_row = als.n_rows;
		int br = 0;
		for (int i = als.n_rows - 1; i >= 0; i--) {
			for (int c = 0; c < als.n_columns; c++) {
				if (inline_dtable.at(i, 0) != 0) {
					br = 1;
					break;
				}
			}
			if (br) {
				break;
			}
			inline_row = i;
		}
	}
	else {
		inline_row = -1;
	}
}

void DconferenceFilter::show() const {
	myprintf("DconferenceFilter: filterj1 -, filterj2 %d, filterj3 %d, filter_symmetry %d\n", filterj2,
		filterj3, filtersymm);
}

std::vector< conference_column > DconferenceFilter::filterList(const std::vector< conference_column > &lst, int verbose) const {
	std::vector< conference_column > out;
	for (size_t i = 0; i < lst.size(); i++) {
		if (this->filter(lst[i])) {
			out.push_back(lst[i]);
		}
	}
	if (verbose) {
		printfd("filterList: %d -> %d\n", lst.size(), out.size());
	}
	return out;
}

std::vector< conference_column > DconferenceFilter::filterListJ2last(const std::vector< conference_column > &column_list) const {
	std::vector< conference_column > out;
	for (size_t i = 0; i < column_list.size(); i++) {
		if (this->filterJ2last(column_list[i])) {
			out.push_back(column_list[i]);
		}
	}
	return out;
}

/// filter a list of cperms using the filterZero method
std::vector< conference_column > DconferenceFilter::filterListZero(const std::vector< conference_column > &lst) const {
	std::vector< conference_column > out;
	for (size_t i = 0; i < lst.size(); i++) {
		if (this->filterZero(lst[i])) {
			out.push_back(lst[i]);
		}
	}
	return out;
}

/// return True of the extension satisfies all checks
bool DconferenceFilter::filter (const conference_column &candidate) const {
        if (filterfirst) {
                if (candidate[0] < 0) {
                        return false;
                }
        }
        if (filtersymm) {
                if (!satisfy_lmc0_symmetry (candidate, check_indices, 0)) {
                        return false;
                }
        }
        if (filterj2) {
                // perform inner product check for all columns
                if (!ipcheck (candidate, als, 0)) {
                        return false;
                }
        }
        if (filterj3) {
                // perform inner product check for all columns
                if (!this->filterJ3 (candidate)) {
                        return false;
                }
        }
        ngood++;
        return true;
}

/// return True of the extension satisfies all J-characteristic checks
bool DconferenceFilter::filterJ(const conference_column &column, int j2start) const {
	if (filterj2) {
		// perform inner product check for all columns
		if (!ipcheck(column, als, j2start)) {
			return false;
		}
	}
	if (filterj3) {
		// perform inner product check for all columns
		if (!this->filterJ3(column)) {
			return false;
		}
	}
	ngood++;
	return true;
}

/// return True of the extension satisfies all J-characteristic checks for the last columns
bool DconferenceFilter::filterJlast(const conference_column &c, int j2start) const {
	if (filterj2) {
		// perform inner product check for all columns
		if (!ipcheck(c, als, j2start)) {
			return false;
		}
	}
	int startidx = this->dtable.n_columns - this->als.n_columns;
	if (filterj3) {
		// perform inner product check for all columns
		if (!this->filterJ3s(c, startidx)) {
			return false;
		}
	}
	ngood++;
	return true;
}

/// return True of the candidate satisfies the symmetry check
bool DconferenceFilter::filterSymmetry (const conference_column &c) const { return satisfy_lmc0_symmetry (c, check_indices, 0); }

bool DconferenceFilter::filterJ2 (const conference_column &c) const { return ipcheck (c, als, 0); }
bool DconferenceFilter::filterJ2last (const conference_column &c) const { return ipcheck (c, als, als.n_columns - 1); }

bool DconferenceFilter::filterZero(const conference_column &c) const {
	for (int i = 0; i < minzvalue - 1; i++) {
		if (c[i] == 0) {
			return false;
		}
	}
	return true;
}

bool DconferenceFilter::filterReason (const conference_column &c) const {
        if (filterfirst) {
                if (c[0] < 0) {
                        myprintf ("filterfirst\n");
                        return false;
                }
        }
        if (filtersymm) {
                if (!satisfy_lmc0_symmetry (c, sd, 0)) {
                        myprintf ("symmetry\n");
                        return false;
                }
        }
        if (filterj2) {
                // perform inner product check for all columns
                if (!ipcheck (c, als, 0)) {
                        myprintf ("j2\n");
                        return false;
                }
        }
        if (filterj3) {
                // perform inner product check for all columns
                if (!this->filterJ3 (c)) {
                        myprintf ("j3\n");
                        return false;
                }
        }
        ngood++;
        myprintf ("filter check good\n");

        return true;
}

/// return True if the candidate satisfies the J3 check
bool DconferenceFilter::filterJ3(const conference_column &c) const {
	const int nc = dtable.n_columns;
	const int N = als.n_rows;
	int Jvalue = 0;
	for (int column_idx = 0; column_idx < nc; column_idx++) {
		Jvalue = 0;

		const array_t *o1 = dtable.array + dtable.n_rows * column_idx;
		for (int row = 0; row < N; row++) {
			Jvalue += (o1[row]) * (c[row]);
		}

		if (Jvalue != 0) {
			return false;
		}
	}
	return true;
}

/** filter conference matrix extension candidates
 *
 * Filtering is based in symmetry and ip
 *
 * \param extensions List of extensions
 * \param als Design to be extended
 * \param filtersymm If true, then filter based on row symmetries?
 * \param filterip If true, filter based on J2
 */
std::vector< conference_column > filterCandidates (const std::vector< conference_column > &extensions, const array_link &als, int filtersymm,
                                       int filterip, int verbose) {
        symmdata sd (als);

        if (verbose >= 2)
                sd.show (1);

        std::vector< int > checkidx = sd.checkIdx ();

        std::vector< conference_column > e2 (0);
        for (size_t i = 0; i < extensions.size (); i++) {

                if (filterip) {
                        // perform inner product check for all columns
                        if (!ipcheck (extensions[i], als)) {
                                if (verbose >= 2) {
                                        printf ("   extension ");
                                        display_vector (extensions[i]);
                                        printf ("\n");
                                        printf ("filterCandidates: reject due to innerproduct (extension %d)\n",
                                                (int)i);
                                        ipcheck (extensions[i], als, 2, 1);
                                }
                                continue;
                        }
                }
                if (filtersymm) {
                        if (!satisfy_lmc0_symmetry (extensions[i], checkidx)) {
                                if (verbose >= 2) {
                                        printf ("filterCandidates: reject due to row symm: ");
                                        display_vector (extensions[i]);
                                        printf ("\n");
                                }
                                continue;
                        }
                }
                e2.push_back (extensions[i]);
        }
        return e2;
}

/// return True if the candidate satisfies the J3 check for specified pairs
bool DconferenceFilter::filterJ3s(const conference_column &c, int idxstart) const {
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

/// return True if the candidate satisfies the J3 check
bool DconferenceFilter::filterJ3inline(const conference_column &c) const {
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

std::vector< conference_column > generateConferenceExtensions(const array_link &al, const conference_t &conference_type,
	int zero_index, int verbose, int filtersymm, int filterj2) {
	if (conference_type.ctype == conference_t::DCONFERENCE) {
		return generateDoubleConferenceExtensions(al, conference_type, verbose, filtersymm, filterj2);
	}

	if (conference_type.ctype == conference_t::CONFERENCE_NORMAL || conference_type.ctype == conference_t::CONFERENCE_DIAGONAL) {

		if (conference_type.itype == CONFERENCE_RESTRICTED_ISOMORPHISM) {
			return generateConferenceRestrictedExtensions(al, conference_type, zero_index, verbose, filtersymm,
				filterj2);
		}
		if (conference_type.itype == CONFERENCE_ISOMORPHISM) {
			myassert(filterj2 == 1, "for single conference designs J2 should be one");
			int filterj3 = 0;
			int filtersymminline = 1;
			std::vector< conference_column > ee = generateSingleConferenceExtensions(
				al, conference_type, zero_index, verbose, filtersymm, filterj2, filterj3, filtersymminline);
			return ee;
		}
		else {
			throw_runtime_exception( printfstring("invalid isomorphism type for ctype %d", int(conference_type.ctype)) );
		}
	}
	else {
		throw_runtime_exception(printfstring("invalid type of conference class %d", conference_type.ctype) );
	}
	std::vector< conference_column > empty_list;
	return empty_list;
}

/// return number of +1 values in the first column of an array
int nOnesFirstColumn (const array_link &al) {
        int n = 0;
        for (int i = 0; i < al.n_rows; i++)
                n += al.atfast (i, 0) == 1;
        return n;
}

/// helper function to sort rows of an array
indexsort rowsorter (const array_link &al) {
        std::vector< mvalue_t< int > > rr;
        for (int i = 0; i < al.n_rows; i++) {
                mvalue_t< int > m;
                for (int k = 0; k < al.n_columns; k++)
                        m.values.push_back (al.at (i, k));
                rr.push_back (m);
        }
        indexsort is (rr);
        return rr;
}

/// special structure for branch-and-bound generation of candidates
struct branch_t {
		/// current row
        int row; 
		/// value assigned to current row
        int rval;
        int nvals[3]; /// number of 0, 1, and -1 remaining

        void show () const {
                myprintf ("branch_t: row %d, val %d: %d %d %d\n", row, rval, nvals[0], nvals[1], nvals[2]);
        }
};

const int bvals[3] = {0, 1, -1}; // these are ordered

template < class object_t >
/// stack object with minimal amount of memory allocations
class lightstack_t {
        object_t *stack;

      private:
        int stack_max_size;
        int current_position;

      public:
        lightstack_t (int sz) : stack_max_size (sz), current_position (0) { stack = new object_t[sz]; }
        ~lightstack_t () { delete[] stack; }

        bool empty () const { return current_position == 0; }
        void push (const object_t &o) {
                this->stack[current_position] = o;
                current_position++;
        }
        object_t &top () const {
                assert (current_position > 0);
                return this->stack[current_position - 1];
        }
        void pop () { current_position--; }
};

/// add new branches to a stack of branches
inline void branch (lightstack_t< branch_t > &branches, branch_t &b, const int *const(&bvals), int istart = 0,
                    int iend = 3) {
        for (int i = istart; i < iend; i++) {

                if (b.nvals[i] == 0)
                        continue;

                branch_t bnew (b);
                bnew.row++;
                bnew.rval = bvals[i];
                bnew.nvals[i]--;
                branches.push (bnew);
        }
}

template < class NumType >
size_t countvalues (const std::vector< NumType > &c, size_t start, size_t end, NumType value) {
        size_t n = 0;
        for (size_t i = start; i < end; i++) {
                if (c[i] == value)
                        n++;
        }
        return n;
}

void debug_candidate (const conference_column &candidate, const std::vector< int > &check_indices, const char *str) {
        printfd ("%s:\n", str);
        printf (" : perm ");
        print_column (candidate);
        printf ("\n");
        conference_column xx (check_indices.begin (), check_indices.end ());
        printf (" : chec ");
        print_column (xx);
        printf ("\n");
}

std::vector< conference_column > debug_branch (conference_column candidate, int gstart, int gend, int block, int blocksize,
                                   std::vector< int > check_indices, int showd) {
        std::vector< conference_column > cc;

        printf ("#### debug_branch: range %d %d\n", gstart, gend);
        debug_candidate (candidate, check_indices, "debug_branch");
        std::sort (candidate.begin () + gstart, candidate.begin () + gend);
        conference_column candidatetmp (candidate.begin (), candidate.end ());

        const int N = gend - gstart;
        // special construction for larger blocksizes
        lightstack_t< branch_t > branches (3 * (gend - gstart + 2));

        if (showd)
                printfd ("inflation of large block: block %d, blocksize %d\n", block, blocksize);

        // count items;
        // push initial branches
        branch_t b1 = {-1, -1, {-1, -1, -1}};
        b1.nvals[0] = countvalues< signed char > (candidate, gstart, gend, bvals[0]);
        b1.nvals[1] = countvalues< signed char > (candidate, gstart, gend, bvals[1]);
        b1.nvals[2] = countvalues< signed char > (candidate, gstart, gend, bvals[2]);

        printf ("  counts 1: %d 0: %d -1: %d\n", b1.nvals[0], b1.nvals[1], b1.nvals[2]);
        branches.push (b1);

        for (int x = gstart; x < gend; x++) {
                candidatetmp[x] = -9;
        }

        double t0 = get_time_ms ();

        long n = 0;
        do {

                branch_t b = branches.top ();

                branches.pop ();
                if (b.row >= 0) { // special check for dummy first branch element...
                        candidatetmp[b.row + gstart] = b.rval;
                }

                if (b.row >= 0) { // special check for dummy first branch element...
                        if (check_symm_zero (candidatetmp, check_indices, b.row + gstart)) {
                                // all good
                                debug_candidate (
                                    candidatetmp, check_indices,
                                    printfstring ("good    based on symmetry: block %d, row %d, range %d %d", block,
                                                  b.row + gstart, gstart, gend)
                                        .c_str ());
                        } else {
                                // discard branch

                                debug_candidate (
                                    candidatetmp, check_indices,
                                    printfstring ("discard based on symmetry: block %d, row %d", block, b.row + gstart)
                                        .c_str ());
                                continue;
                        }
                }
                if (b.row == N - 1) {
                        n++;
                        cc.push_back (candidatetmp);
                        continue;
                }

                // perform branching
                branch (branches, b, bvals, 0, 3);

        } while (!branches.empty ());
        if (showd) {
                printfd ("inflation of large block: blocksize %d: %d branch calls\n", blocksize, n);
        }
        return cc;
}

void inflateCandidateExtensionHelper (std::vector< conference_column > &list, const conference_column &basecandidate, conference_column &candidate,
                                      int block, const array_link &al, const symmetry_group &alsg,
                                      const std::vector< int > &check_indices, const conference_t &ct, int verbose,
                                      const DconferenceFilter &filter, long &ntotal) {
        const symmdata &sd = filter.sd;
        int nblocks = alsg.ngroups;

        if (block == nblocks) {
                ntotal++;
                const int j2start = std::max (0, al.n_columns - 2);

                if (filter.filterJlast (candidate, j2start)) {
                        list.push_back (candidate);
                }
                return;
        }

        const int blocksize = alsg.gsize[block];

        if (block > nblocks - 8 && blocksize > 1 ) {
                int r = alsg.gstart[block] - 1;

                bool check = filter.filterJpartial (candidate, r);
                if (verbose >= 2) {
                        int j = partial_inner_product (candidate, filter.als, filter.als.n_columns - 1, r);
                        printfd ("partial check: block %d/%d:  r %d, j %d, check %d\n", block, nblocks, r, j, check);
                }
                if (!check)
                        return;
        }
        if (verbose >= 2)
                printfd ("inflateCandidateExtensionHelper: block %d/%d: blocksize %d\n", block, alsg.gsize.size (),
                         blocksize);

        if (blocksize == 1) {
                // easy case
                inflateCandidateExtensionHelper (list, basecandidate, candidate, block + 1, al, alsg, check_indices,
                                                 ct, verbose, filter, ntotal);
                return;
        }
        int gstart = alsg.gstart[block];
        int gend = alsg.gstart[block + 1];

        if (block <= -1) {
                printfd ("row: %d to %d\n", gstart, gend);
                conference_column tmp (candidate.begin () + gstart, candidate.begin () + gend);
                printf ("   current perm: ");
                print_column (tmp);
                printf ("\n");
        }

        if (verbose >= 2)
                printfd ("  split\n");
        if (blocksize < 3 ) {
                unsigned long iter = 0;
                std::sort (candidate.begin () + gstart, candidate.begin () + gend);
                unsigned long nbc = 0;
                do {
                        iter++;

                        if (satisfy_lmc0_symmetry (candidate, check_indices, gstart, gend)) {
                                nbc++;
                                inflateCandidateExtensionHelper (list, basecandidate, candidate, block + 1, al, alsg,
                                                                 check_indices, ct, verbose, filter, ntotal);
                        } else {
                        }
                } while (std::next_permutation (candidate.begin () + gstart, candidate.begin () + gend));
                if (verbose >= 2)
                        printfd ("nbc block %d: %d/%ld\n", block, nbc, iter);

                if (blocksize > 10 && (verbose >= 4)) {
                        printfd ("block %d: nbc %ld\n", block, (long)nbc);
                        debug_candidate (candidate, check_indices, "...");
                }

        } else {
                const int showd = 0; // show debugging output
                std::sort (candidate.begin () + gstart, candidate.begin () + gend);
                conference_column candidatetmp (candidate.begin (), candidate.end ());

                for (int x = gstart; x < gend; x++) {
                        candidatetmp[x] = -9;
                }
                const int N = gend - gstart;
                // special construction for larger blocksizes
                lightstack_t< branch_t > branches (3 * (gend - gstart + 2));

                if (showd)
                        printfd ("inflation of large block: block %d, blocksize %d\n", block, blocksize);

                // push initial branches
                branch_t b1 = {-1, -1, {-1, -1, -1}};
                b1.nvals[0] = countvalues< signed char > (candidate, gstart, gend, bvals[0]);
                b1.nvals[1] = countvalues< signed char > (candidate, gstart, gend, bvals[1]);
                b1.nvals[2] = countvalues< signed char > (candidate, gstart, gend, bvals[2]);

                branches.push (b1);

                double t0 = get_time_ms ();

                long nbc = 0;
                do {
                        branch_t b = branches.top ();

                        branches.pop ();
                        if (b.row >= 0) { // special check for dummy first branch element...
                                candidatetmp[b.row + gstart] = b.rval;
                        }

                        if (b.row >= 0) { // special check for dummy first branch element...
                                if (check_symm_zero (candidatetmp, check_indices, b.row + gstart)) {
                                        // all good
                                } else {
                                        // discard branch
                                        continue;
                                }
                        }
                        if (b.row == N - 1) {
                                nbc++;
                                // call the inflate function
                                inflateCandidateExtensionHelper (list, basecandidate, candidatetmp, block + 1, al,
                                                                 alsg, check_indices, ct, verbose, filter, ntotal);
                                continue;
                        }

                        // perform branching
                        branch (branches, b, bvals, 0, 3);

                } while (!branches.empty ());
                if (verbose >= 2)
                        printfd ("nbc block %d: %d/%ld\n", block, nbc, nbc);

                if (showd) {
                        printfd ("inflation of large block: blocksize %d: %d branch calls\n", blocksize, nbc);
                }
        }
}


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
                                                const conference_t &ct, int verbose, const DconferenceFilter &filter) {
        long ntotal = 0;

        conference_column candidate = basecandidate;
        int block = 0;
        std::vector< conference_column > cc;
        inflateCandidateExtensionHelper (cc, basecandidate, candidate, block, als, alsg, check_indices, ct, verbose,
                                         filter, ntotal);

        if (verbose >= 2 || 0) {
                printfd ("inflateCandidateExtension: generated %ld/%ld candidates (k %d)\n", (long)cc.size (), ntotal,
                         als.n_columns);
        }
        return cc;
}

/// return true if zero is a specified position
inline bool checkZeroPosition(const conference_column &column, int zero_position) {
	if (zero_position <= 0)
		return false;

	if (column[zero_position] == 0) {
		return true;
	}
	else
		return false;
}

std::vector< conference_column > generateSingleConferenceExtensions (const array_link &al, const conference_t &ct, int zero_index,
                                                         int verbose, int filter_symmetry, int filterj2, int filterj3,
                                                         int filter_symmetry_inline) {
        double t0 = get_time_ms ();
        if (verbose) {
                myprintf ("%s: filters: symmetry %d, symmetry inline %d, j2 %d, j3 %d\n", __FUNCTION__, filter_symmetry,
                          filter_symmetry_inline, filterj2, filterj3);
        }
        myassert (al.n_columns > 1, "need at least 2 columns");
		myassert(ct.j1zero==0, "for single conference designs j1zero needs to be zero");

        if (filter_symmetry_inline) {
                if (!filter_symmetry) {
                        myprintf ("filter_symmetry_inline selected, but no filter_symmetry, is this what you want?");
                }
        }
        const int N = al.n_rows;
        DconferenceFilter dfilter (al, filter_symmetry, filterj2, filterj3);
        if (verbose >= 2) {
                dfilter.show ();
        }

        std::vector< int > sidx = dfilter.sd.checkIdx ();
        if (verbose >= 2) {
                print_perm ("DconferenceFilter sidx: ", sidx);
        }
        std::vector< conference_column > candidate_columns;
        conference_column current_column (al.n_rows);
		current_column[0] = 1;

        lightstack_t< branch_t > branches (3 * N);


        // push initial branches
        branch_t b1 = {0, 1, {1, N / 2 - 1, N / 2 - 1}}; // branches starting with a 1
        branches.push (b1);

        std::vector< long > branch_count (N + 1);

        long n = 0;
        do {

                branch_t b = branches.top ();

                if (verbose >= 3) {
                        b.show ();
                }
                branch_count[b.row]++;
                
                branches.pop ();
                current_column[b.row] = b.rval; // update column vector

                if (b.row == dfilter.inline_row && filterj3) {
                        if (!dfilter.filterJ3inline (current_column))
                                // discard branch
                                continue;
                }
                if (b.row == N - 1) {
                        n++;
                        if (verbose >= 3) {
                                myprintf ("n %d: filter %d: ", (int)n, dfilter.filter (current_column));
                                print_column (current_column);
                                printf ("\n");
                        }

                        if (dfilter.filter (current_column)) {
                                if (zero_index >= 0) {
                                        if (!checkZeroPosition (current_column, zero_index)) {
                                                continue;
                                        }
                                }

                                if (verbose >= 2) {
                                        printfd ("## push candidate   : ");
                                        print_column (current_column);
                                        myprintf ("\n");
                                }
                                candidate_columns.push_back (current_column);
                        } else {
                                // discard candindate
                        }
                        continue;
                }
                int istart = 0;
                const int iend = 3;
                if (sidx[b.row + 1] && filter_symmetry_inline) {
                        if (current_column[b.row] == -1) {
                                istart = 2;
                        }
                        if (current_column[b.row] == 1) {
                                istart = 1;
                        }
                }
                for (int i = istart; i < iend; i++) {

                        if (b.nvals[i] == 0)
                                continue;

                        // checks
                        branch_t bnew (b);
                        bnew.row++;
                        bnew.rval = bvals[i];
                        bnew.nvals[i]--;
                        branches.push (bnew);
                }

        } while (!branches.empty ());

        if (verbose) {
               if (verbose>=2) {
                         printf ("branch count:\n");
                         for (int i = 0; i <= N; i++) {
                              printf ("  %d: %ld\n", i, branch_count[i]);
                         }
               }
                printfd ("%s: %.3f [s]: generated %ld/%ld/%ld perms (len %ld)\n", __FUNCTION__, get_time_ms () - t0,
                         (long)candidate_columns.size (), n, factorial< long > (current_column.size ()), (long)current_column.size ());
                if (verbose>=3) {
                 al.show();
                 al.transposed().showarray(); showCandidates ( candidate_columns );
                }
        }
        return candidate_columns;
}

std::vector< conference_column > generateDoubleConferenceExtensions (const array_link &al, const conference_t &ct, int verbose,
                                                         int filtersymm, int filterj2, int filterj3,
                                                         int filtersymminline) {
        if (verbose)
                myprintf (
                    "generateDoubleConferenceExtensions: filters: symmetry %d, symmetry inline %d, j2 %d, j3 %d\n",
                    filtersymm, filtersymminline, filterj2, filterj3);

        myassert (ct.j1zero == 1, "for double conference designs j1zero should 1");

        const int N = al.n_rows;
        DconferenceFilter dfilter (al, filtersymm, filterj2);

        std::vector< int > sidx = dfilter.sd.checkIdx ();
        if (verbose >= 2) {
                print_perm ("sidx: ", sidx);
        }
        std::vector< conference_column > cc;
        conference_column c (al.n_rows);

        lightstack_t< branch_t > branches (3 * N);

        c[0] = 1;

        // push initial branches
        branch_t b1 = {0, 1, {2, N / 2 - 2, N / 2 - 1}};
        branches.push (b1);
        branch_t b0 = {0, 0, {1, N / 2 - 1, N / 2 - 1}};
        branches.push (b0);
        double t0 = get_time_ms ();

        std::vector< long > branch_count (N + 1);

        long n = 0;
        do {

                branch_t b = branches.top ();
                branch_count[b.row]++;
                branches.pop ();
                c[b.row] = b.rval;

                if (b.row == dfilter.inline_row && filterj3) {
                        if (!dfilter.filterJ3inline (c))
                                // discard branch
                                continue;
                }
                if (b.row == N - 1) {
                        n++;
                        if (dfilter.filter (c)) {
                                if (verbose >= 2) {
                                        printfd ("## push candindate   : ");
                                        print_column (c);
                                        myprintf ("\n");
                                }
                                cc.push_back (c);
                        } else {
                                // printfd ( "## discard candindate: " );
                        }
                        continue;
                }
                int istart = 0;
                const int iend = 3;
                if (sidx[b.row + 1] && filtersymminline) {
                        if (c[b.row] == -1) {
                                istart = 2;
                        }
                        if (c[b.row] == 1 && 1) {
                                istart = 1;
                        }
                }
                for (int i = istart; i < iend; i++) {

                        if (b.nvals[i] == 0)
                                continue;

                        // checks

                        branch_t bnew (b);
                        bnew.row++;
                        bnew.rval = bvals[i];
                        bnew.nvals[i]--;
                        branches.push (bnew);
                }

        } while (!branches.empty ());

        if (verbose) {
               if (verbose>=2) {
                         printf ("branch count:\n");
                         for (int i = 0; i <= N; i++) {
                              printf ("  %d: %ld\n", i, branch_count[i]);
                         }
               }
                printfd ("generateDoubleConferenceExtensions: generated %ld/%ld/%ld perms (len %ld)\n",
                         (long)cc.size (), n, factorial< long > (c.size ()), (long)c.size ());
        }
        return cc;
}

/** generate double conference matrices
 *
 * Old vesion that still can handle the j1zero=0 case
 *
 **/
std::vector< conference_column > generateDoubleConferenceExtensions2 (const array_link &al, const conference_t &ct, int verbose,
                                                          int filtersymm, int filterip) {
        assert (ct.itype == CONFERENCE_RESTRICTED_ISOMORPHISM || ct.itype == CONFERENCE_ISOMORPHISM);

        int j1zero = ct.j1zero;

        std::vector< conference_column > cc;

        const int N = ct.N;
        conference_column c (N);

        DconferenceFilter dfilter (al, filtersymm, filterip);
        dfilter.filterfirst = 1;
        dfilter.filterj3 = ct.j3zero;
        unsigned long n = 0;
        for (int i = 0; i < N - 2; i++) {
                if (j1zero && i != (N - 2) / 2)
                        continue;

                // fill initial permutation
                std::fill (c.begin (), c.end (), -1);
                c[0] = 0;
                c[1] = 0;
                for (int k = 2; k < i + 2; k++)
                        c[k] = 1;

                std::sort (c.begin (), c.end ());

                do {
                        n++;

                        if (dfilter.filter (c)) {
                                cc.push_back (c);
                        }
                } while (std::next_permutation (c.begin (), c.end ()));
        }

        if (verbose) {
                printfd ("generateDoubleConferenceExtensions: generated %ld/%ld/%ld perms (len %ld)\n",
                         (long)cc.size (), n, factorial< long > (c.size ()), (long)c.size ());
        }
        return cc;
}

std::vector< conference_column > generateConferenceRestrictedExtensions (const array_link &al, const conference_t &ct, int kz,
                                                             int verbose, int filtersymm, int filterip) {

        const int extcol = al.n_columns;
        const int N = ct.N;

        // special case
        if (extcol == 1) {
                std::vector< conference_column > ee;

                // we have no values +1 in the first column and
                int no = nOnesFirstColumn (al);
                int n1 = no - 1;
                int n2 = N - n1 - 2;

                // in the second column we start with [1,0]^T and then in block1: k1 values +1, block2: k2 values +1
                // we must have k2 = k1+(n2-n1)/2 = k1 +N -2n1-2

                conference_column cc (N);
                cc[0] = 1;
                cc[1] = 0;
                for (int k1 = 0; k1 <= n1; k1++) {
                        int k2 = k1 + (n2 - n1) / 2;
                        if (k2 > n2)
                                continue;
                        if (k2 < 0)
                                continue;
                        printf ("generateConferenceRestrictedExtensions: column 1: n1 %d, n2 %d, k1 %d, k2 %d\n", n1,
                                n2, k1, k2);

                        std::fill (cc.begin () + 2, cc.end (), -1);

                        for (int i = 2; i < 2 + k1; i++)
                                cc[i] = 1;
                        for (int i = 2 + n1; i < 2 + n1 + k2; i++)
                                cc[i] = 1;

                        ee.push_back (cc);
                }
                return ee;
        }

        std::vector< int > moi; // indices of -1 in first column
        for (int i = 0; i < N; i++) {
                if (al.atfast (i, 0) == -1)
                        moi.push_back (i);
        }
        array_link alx = al.clone ();

        // multiply
        for (size_t i = 0; i < moi.size (); i++) {
                alx.negateRow (moi[i]);
        }

        // sort rows of array
        indexsort is = rowsorter (alx);

        // now get candidate columns for the normal case, afterwards convert then using the rowsorter and row negations

        // loop over all possible first combinations
        std::vector< conference_column > ff = get_first (N, kz, verbose);

        if (verbose >= 2) {
                for (size_t i = 0; i < ff.size (); i++) {
                        printf ("extend1 %d: N %d: ", (int)i, N);
                        display_vector (ff[i]);
                        printf ("\n");
                }
        }

        conference_extend_t ce;
        std::vector< conference_column > extensions (0);

        ce.first = ff;
        ce.second = ff;

        array_link als = al.selectFirstColumns (extcol);

        conference_column c0 = getColumn (al, 0);
        conference_column c1 = getColumn (al, 1);
        for (size_t i = 0; i < ce.first.size (); i++) {
                int ip = innerprod (c0, ce.first[i]);

                int target = -ip;

                std::vector< conference_column > ff2 = get_second (N, kz, target, verbose >= 2);
                ce.second = ff2;


                for (size_t j = 0; j < ff2.size (); j++) {
                        conference_column c = ce.combine (i, j);

                        extensions.push_back (c);
                        continue;
                }
        }

        ce.extensions = extensions;
        if (verbose >= 2)
                printf ("generateConferenceExtensions: after generation: found %d extensions\n",
                        (int)extensions.size ());

        // perform row symmetry check
        std::vector< conference_column > e2 = filterCandidates (extensions, al, filtersymm, filterip, verbose);

        if (verbose >= 1)
                printf ("extend_conference: symmetry check %d + ip filter %d: %d->%d\n", filtersymm, filterip,
                        (int)extensions.size (), (int)e2.size ());

        ce.extensions = e2;

        return e2;
}

/** select maximum possible position for a zero in the column of a design assuming the design is in LMC0 format
 *
 * \param maxzpos Maximum position of zero an specified design
 * \param ctype Type of conference matrix
 * \param al Array containing the design
 * \param extcol Extension column
 *
 */
int selectZmax (int maxzpos, const conference_t::conference_type &ctype, const array_link &al, int extcol) {
        if (maxzpos < 0) {
                switch (ctype) {
                case conference_t::CONFERENCE_NORMAL:
                        maxzpos = al.n_rows - 1;
                        break;
                case conference_t::CONFERENCE_DIAGONAL:
                        maxzpos = extcol;
                        break;
                case conference_t::DCONFERENCE:
                        maxzpos = al.n_rows - 1;
                        break;
                default:
                        printfd ("not implemented...\n");
                        maxzpos = al.n_rows - 1;
                }
        }
        return maxzpos;
}

std::vector< conference_column > generateDoubleConferenceExtensionsInflate (const array_link &al, const conference_t &ct,
                                                                int verbose, int filterj2, int filterj3,
                                                                int kstart = 2);

conference_extend_t extend_double_conference_matrix (const array_link &al, const conference_t &ct,
                                                     CandidateGeneratorDouble &cgenerator, int extcol, int verbose,
                                                     int maxzpos) {
        conference_extend_t ce;
        ce.extensions.resize (0);

        const int N = ct.N;
        const int k = extcol;
        const int maxzval = maxz (al);

        if (verbose)
                printf ("--- extend_double_conference_matrix: extcol %d, maxz %d, itype %d ---\n", extcol, maxzval,
                        ct.itype);

        int filterip = 1;
        int filtersymm = 1;
        std::vector< conference_column > cc;

        if (k >= 3 && filtersymm && filterip && ct.j1zero == 1 && 1) {
                cc = cgenerator.generateCandidates (al);
        } else {
                if (k > 3 && ct.j1zero == 1) {
                        cc = generateDoubleConferenceExtensionsInflate (al, ct, verbose, filterip, 1);
                } else {
                        if (ct.j1zero == 1)
                                cc = generateDoubleConferenceExtensions (al, ct, verbose, filtersymm, filterip);
                        else
                                cc = generateDoubleConferenceExtensions2 (al, ct, verbose, filtersymm, filterip);
                }
        }

        if (ct.j3zero) {
                cc = filterJ3 (cc, al, verbose);
        }

        ce.extensions = cc;
        return ce;
}

/** Extend a single conference design with candidate columns */
conference_extend_t extend_conference_matrix(const array_link &al, const conference_t &ct, int extcol,
	int verbose = 1, int maxzpos = -1);

conference_extend_t extend_conference_matrix (const array_link &al, const conference_t &ct, int extcol, int verbose,
                                              int maxzpos) {
        conference_extend_t ce;
        ce.extensions.resize (0);

        const int N = ct.N;
        const int k = extcol;
        const int maxzval = maxz (al);

        if (verbose)
                printf ("--- extend_conference_matrix: extcol %d, maxz %d, itype %d ---\n", extcol, maxzval, ct.itype);

        const int zstart = maxzval + 1;

        maxzpos = selectZmax (maxzpos, ct.ctype, al, extcol);

        for (int ii = zstart; ii < maxzpos + 1; ii++) {
                if (verbose >= 2)
                        printf ("array: zero_index %d: generate\n", ii);
                std::vector< conference_column > extensionsX = generateConferenceExtensions (al, ct, ii, verbose, 1, 1);

                if (verbose >= 2)
                        printf ("array: zero_index %d: %d extensions\n", ii, (int)extensionsX.size ());
                ce.extensions.insert (ce.extensions.end (), extensionsX.begin (), extensionsX.end ());
        }

        return ce;
}

/// extend a conference matrix using a generator for the candidate extensions
conference_extend_t extend_conference_matrix_generator (const array_link &al, const conference_t &ct, int extcol,
                                                        int verbose, int maxzpos,
                                                        const CandidateGenerator &cgenerator) {
        conference_extend_t ce;

        const int N = ct.N;
        const int k = extcol;
        const int maxzval = maxz (al);

        const int zstart = maxzval + 1;

        maxzpos = selectZmax (maxzpos, ct.ctype, al, extcol);
        if (verbose) {
                printfd ("--- extend_conference_matrix: extcol %d, zstart %d, maxz %d, maxzpos %d ---\n", extcol,
                         zstart, maxzval, maxzpos);

                if (verbose >= 2)
                        al.showarray ();
        }

        // loop over all possible positions of the next zero in the design
        for (int ii = zstart; ii < maxzpos + 1; ii++) {
                if (verbose >= 2)
                        printfd ("array: generate with zero at %d\n", ii);
                std::vector< conference_column > extensionsX;

                {
                        // generate possible extension candidates
                        const std::vector< conference_column > &cl = cgenerator.generateCandidatesZero (al, ii);
                        if (verbose > 2) {
                                printfd ("-- ii %d, %d candidates \n", ii, cl.size ());
                                if (verbose >= 2) {
                                        cgenerator.showCandidates (2);
                                }
                        }
                        if (ct.ctype == conference_t::CONFERENCE_DIAGONAL) {
                                extensionsX = filterCandidatesSymm (
                                    cl, al, verbose); // not needed, but might make calculations faster
                                extensionsX = filterCandidates (extensionsX, al, 1, 1, verbose);
                        } else {
                                extensionsX = filterCandidates (cl, al, 1, 1, verbose);
                        }
                }

                if (verbose >= 2) {
                        printf ("array: zero_index %d: %d extensions\n", ii, (int)extensionsX.size ());
                }

                ce.extensions.insert (ce.extensions.end (), extensionsX.begin (), extensionsX.end ());
        }

        return ce;
}

/// sort rows in an array
array_link sortrows (const array_link &al) {
        const int nc = al.n_columns;

        std::vector< mvalue_t< int > > rr (al.n_rows);
        for (int i = 0; i < al.n_rows; i++) {
                mvalue_t< int > &m = rr[i];
                m.values.resize (nc);

                for (int k = 0; k < nc; k++) {
                        m.values[k] = al.atfast (i, k);
                }
        }
        indexsort sorter (rr);
        sorter.show ();

        array_link out (al.n_rows, al.n_columns, 0);
        for (int r = 0; r < al.n_rows; r++) {
                for (int i = 0; i < al.n_columns; i++) {
                        int newrow = sorter.indices[r];
                        out._setvalue (r, i, al.atfast (newrow, i));
                }
        }
        return out;
}

template < typename T >
/// return size of std::vector in memory (in bytes)
size_t vectorsizeof (const typename std::vector< T > &vec) {
        return sizeof (T) * vec.size ();
}

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

conf_candidates_t generateCandidateExtensions (const conference_t ctype, int verbose = 1, int ncstart = 3,
                                               int ncmax = -1, int root = -1) {

        conf_candidates_t cande;

        cande.ce.resize (ctype.N);

        array_link al2 = ctype.create_root ();
        array_link al3 = ctype.create_root_three_columns ();

        if (ncmax == -1) {
                ncmax = ctype.N;
                if (ctype.ctype == conference_t::CONFERENCE_DIAGONAL) {
                        ncmax = ncstart;
                }
        }

        for (int extcol = ncstart - 1; extcol < ncmax; extcol++) {
                std::vector< conference_column > ee;
                {
                        switch (root) {
                        case -1: {
                                if (extcol == 2)
                                        ee = generateConferenceExtensions (al2, ctype, extcol, 0, 0, 1);
                                else
                                        ee = generateConferenceExtensions (al3, ctype, extcol, 0, 0, 1);
                        } break;
                        case 2:
                                ee = generateConferenceExtensions (al2, ctype, extcol, 0, 0, 1);
                                break;
                                break;
                        case 3:
                                ee = generateConferenceExtensions (al3, ctype, extcol, 0, 0, 1);
                                break;
                        default:
                                printfd ("not implemented!\n");
                                break;
                        }
                }

                if ((long)vectorsizeof (ee) > (long(1) * 1024 * 1024 * 1024) / (long)ctype.N) {
                        printfd ("generateCandidateExtensions: set of generated candidates too large, aborting (root "
                                 "%d, extcol %d, ee.size() %d)",
                                 root, extcol, ee.size ());
                        assert (0);
                        exit (0);
                }
                cande.ce[extcol] = ee;
        }

        cande.info (verbose);
        return cande;
}

arraylist_t extend_double_conference (const arraylist_t &lst, const conference_t ctype, int verbose) {
        arraylist_t outlist;
        if (verbose >= 2) {
                printfd ("extend_double_conference: start with %d arrays\n", (int)lst.size ());
        }
        double t0 = get_time_ms ();

        int vb = std::max (0, verbose - 1);

        CandidateGeneratorDouble cgenerator (array_link (), ctype);

        for (size_t i = 0; i < lst.size (); i++) {
                const array_link &al = lst[i];
                int extcol = al.n_columns;
                conference_extend_t ce = extend_double_conference_matrix (al, ctype, cgenerator, extcol, vb, -1);

                arraylist_t ll = ce.getarrays (al);
                outlist.insert (outlist.end (), ll.begin (), ll.end ());

                if (verbose >= 2 || (verbose >= 1 && (i % 1000 == 0 || i == lst.size () - 1))) {
                        double dt = get_time_ms () - t0;
                        printf ("extend_conference: extended array %d/%d to %d/%d arrays (%.1f [s])\n", (int)i,
                                (int)lst.size (), (int)ll.size (), (int)outlist.size (), dt);
                        fflush (0);
                }
        }
        return outlist;
}

arraylist_t extend_conference_restricted (const arraylist_t &lst, const conference_t ctype, int verbose) {
        arraylist_t outlist;

        if (verbose >= 2) {
                printfd ("extend_conference: start %d\n", (int)lst.size ());
        }

        int vb = std::max (0, verbose - 1);

        for (size_t i = 0; i < lst.size (); i++) {
                const array_link &al = lst[i];
                int extcol = al.n_columns;
                conference_extend_t ce = extend_conference_matrix (al, ctype, extcol, vb, -1);

                arraylist_t ll = ce.getarrays (al);

                outlist.insert (outlist.end (), ll.begin (), ll.end ());

                if (verbose >= 2 || (verbose >= 1 && (i % 200 == 0 || i == lst.size () - 1))) {
					printf ("extend_conference: extended array %d/%d to %d arrays\n", (int)i, (int)lst.size (), (int)ll.size());
                        fflush (0);
                }
        }
        return outlist;
}

/// select the unique arrays from a list of arrays. the indices of the unique arrays are returned
std::vector< int > selectUniqueArrayIndices (const arraylist_t &lstr, int verbose) {
        // perform stable sort
        indexsort sortidx (lstr);

        const std::vector< int > &idx = sortidx.indices;

        const size_t nn = lstr.size ();

        std::vector< int > cidx (nn);
        std::vector< int > ridx;

        array_link prev;

        if (lstr.size () > 0)
                prev = lstr[0];
        prev.setconstant (-10);

        int ci = -1;
        for (size_t i = 0; i < idx.size (); i++) {
                array_link al = lstr[idx[i]];
                if (al != prev) {
                        // new isomorphism class
                        ci++;
                        if (verbose >= 3)
                                printf ("selectConferenceIsomorpismClasses: representative %d: index %d\n", (int)ci,
                                        (int)idx[i]);
                        ridx.push_back (idx[i]);
                        prev = al;
                }
                cidx[i] = ci;
        }
        return ridx;
}

/// reduce a matrix to Nauty normal form
array_link reduceMatrixNauty (const array_link &array, matrix_isomorphism_t itype, int verbose) {
        array_link alx;
        switch (itype) {
        case CONFERENCE_ISOMORPHISM: {
                alx = reduceConference (array, verbose >= 2);
        } break;
        case CONFERENCE_RESTRICTED_ISOMORPHISM: {
                arraydata_t arrayclass (3, array.n_rows, 1, array.n_columns);
                array_transformation_t t = reduceOAnauty (array + 1, verbose >= 2, arrayclass);
                alx = t.apply (array + 1) + (-1);
                break;
        }
        default:
                printfd ("error: isomorphism type not implemented\n");
                break;
        }
        return alx;
}

/** Class to select isomorphism classes
 *
 * The selecting is performed by reducing to Nauty normal form. By performing the isomorphism check incrementally we can save memory.
 */
class ConferenceIsomorphismSelector {
      public:
        matrix_isomorphism_t itype;
        int verbose;
        int select_isomorphism_classes; /// if true then select only a single representative for each isomorphism class

        arraylist_t candidates;
        arraylist_t reductions;

      private:
        int nadd;

      public:
        ConferenceIsomorphismSelector (matrix_isomorphism_t itype, int verbose, int select_isomorphism_classes) {
                this->itype = itype;
                this->verbose = verbose;
                this->select_isomorphism_classes = select_isomorphism_classes;
                this->nadd = 0;
        }

        size_t size () const { return candidates.size (); }

        /** Add a set of arrays to the list of isomorphism classes
         */
        void add (const arraylist_t &lst) {

                long nstart = candidates.size ();
                long nextra = lst.size ();
                candidates.insert (candidates.end (), lst.begin (), lst.end ());

                nadd++;
                if (select_isomorphism_classes) {
                        for (size_t i = 0; i < lst.size (); i++) {
                                array_link alr = reduceMatrixNauty (lst[i], itype, verbose);
                                reductions.push_back (alr);
                        }

                        if (nadd % 1000 == 0) {
                                std::vector< int > ridx = selectUniqueArrayIndices (reductions, verbose);

                                arraylist_t tmp;
                                arraylist_t tmpr;
                                for (size_t i = 0; i < ridx.size (); i++) {
                                        size_t ix = ridx[i];
                                        tmp.push_back (candidates[ix]);
                                        tmpr.push_back (reductions[ix]);
                                }
                                candidates = tmp;
                                reductions = tmpr;
                                if (verbose)
                                        printf ("  reduce ... %ld -> %ld \n", nstart, (long)candidates.size ());
                        }
                }

                long nfinal = candidates.size ();

                if (verbose) {
                        printf ("ConferenceIsomorphismSelector: add %ld -> %ld -> %ld\n", (long)nstart,
                                nstart + nextra, nfinal);
                }
        }
};

arraylist_t extend_conference_plain (const arraylist_t &lst, const conference_t ctype, int verbose,
                                     int select_isomorphism_classes) {
        double t0 = get_time_ms ();
        arraylist_t outlist;

        if (verbose >= 2) {
                printfd ("extend_conference: start %d\n", (int)lst.size ());
        }

        int vb = std::max (0, verbose - 1);

        ConferenceIsomorphismSelector selector (ctype.itype, verbose >= 2, select_isomorphism_classes);

        for (size_t i = 0; i < lst.size (); i++) {
                const array_link &al = lst[i];
                int extcol = al.n_columns;
                conference_extend_t conference_extensions = extend_conference_matrix (al, ctype, extcol, vb, -1);

                arraylist_t ll = conference_extensions.getarrays (al);
                const int nn = ll.size ();

                selector.add (ll);

                if (verbose >= 2 || (verbose >= 1 && (i % 1000 == 0 || i == lst.size () - 1))) {
                        printf ("extend_conference: extended array %d/%d to %d arrays (total %ld, %.1f [s])\n", (int)i,
                                (int)lst.size (), nn, (long)selector.size (), get_time_ms () - t0);
                        fflush (0);
                }
        }

        return selector.candidates;
}

arraylist_t extend_conference (const arraylist_t &lst, const conference_t conference_type, int verbose,
                               int select_isomorphism_classes) {
        double t0 = get_time_ms ();

        if (verbose >= 2) {
                printfd ("extend_conference: start with %d arrays\n", (int)lst.size ());
        }

        int subverbose = std::max (0, verbose - 1);

        CandidateGenerator cgenerator (array_link (), conference_type);

        ConferenceIsomorphismSelector selector (conference_type.itype, verbose >= 2, select_isomorphism_classes);

        for (size_t i = 0; i < lst.size (); i++) {
                if (verbose >= 2)
                        printfd ("extend_conference: extend array %d\n", i);

                const array_link &al = lst[i];
                int extcol = al.n_columns;
                conference_extend_t ce =
                    extend_conference_matrix_generator (al, conference_type, extcol, subverbose, -1, cgenerator);

                arraylist_t ll = ce.getarrays (al);
                selector.add (ll);

                if (verbose >= 2 || (verbose >= 1 && (i % 400 == 0 || i == lst.size () - 1))) {
                        printf ("extend_conference: extended array %d/%d to %d arrays (total %ld, %.1f [s])\n", (int)i,
                                (int)lst.size (), (int)ll.size (), (long)selector.size (), get_time_ms () - t0);
                        fflush (0);
                }
        }

        return selector.candidates;
}

std::pair< arraylist_t, std::vector< int > > selectConferenceIsomorpismHelper (const arraylist_t &lst, int verbose,
                                                                               matrix_isomorphism_t itype) {
        const int nn = lst.size ();

        arraylist_t lstr;
        double t0 = get_time_ms ();

        // safety check
        if (lst.size () > 0) {
                if (!lst[0].is_conference ()) {
                        printfd ("error: arrays should have positive integer values\n");
                        arraylist_t lstgood;
                        std::vector< int > cidx;
                        return std::pair< arraylist_t, std::vector< int > > (lstgood, cidx);
                }
        }
        for (int i = 0; i < (int)lst.size (); i++) {
                if (verbose >= 1 && (i % 20000 == 0 || i == (int)lst.size () - 1))
                        printf ("selectConferenceIsomorpismClasses: reduce %d/%d\n", i, (int)lst.size ());
                array_link alx = reduceMatrixNauty (lst[i], itype, 2 * (verbose >= 3));

                lstr.push_back (alx);
        }

        // perform stable sort
        arraylist_t lstgood;
        indexsort sortidx (lstr);

        const std::vector< int > &idx = sortidx.indices;

        std::vector< int > cidx (nn);

        array_link prev;

        if (lst.size () > 0)
                prev = lst[0];
        prev.setconstant (-10);

        int ci = -1;
        for (size_t i = 0; i < idx.size (); i++) {
                array_link al = lstr[idx[i]];
                if (al != prev) {
                        // new isomorphism class
                        if (verbose >= 3)
                                printf ("selectConferenceIsomorpismClasses: representative %d: index %d\n",
                                        (int)lstgood.size (), (int)idx[i]);

                        lstgood.push_back (lst[idx[i]]);
                        prev = al;
                        ci++;
                }
                cidx[idx[i]] = ci;
        }

        if (verbose) {
                double dt = get_time_ms () - t0;
                myprintf ("selectConferenceIsomorpismClasses: select classes %d->%d (%.1f kArrays/h)\n",
                          (int)lst.size (), (int)lstgood.size (), 3600 * 1e-3 * double(lst.size ()) / dt);
        }
        return std::pair< arraylist_t, std::vector< int > > (lstgood, cidx);
}

std::vector< int > selectConferenceIsomorpismIndices (const arraylist_t &lst, int verbose,
                                                      matrix_isomorphism_t itype) {
        std::pair< arraylist_t, std::vector< int > > pp = selectConferenceIsomorpismHelper (lst, verbose, itype);
        return pp.second;
}

arraylist_t selectConferenceIsomorpismClasses (const arraylist_t &lst, int verbose, matrix_isomorphism_t itype) {
        std::pair< arraylist_t, std::vector< int > > pp = selectConferenceIsomorpismHelper (lst, verbose, itype);
        return pp.first;
}

arraylist_t selectLMC0doubleconference (const arraylist_t &list, int verbose, const conference_t &ctype) {
        double t0 = get_time_ms ();

        arraylist_t out;
        for (size_t i = 0; i < list.size (); i++) {
                lmc_t r = LMC0checkDC (list[i]);

                if (verbose >= 2 || (verbose && i % 10000 == 0))
                        printfd ("selectLMC0: i %d/%d, result %d, total %d (%.1f [s])\n", i, list.size (), r,
                                 out.size (), get_time_ms () - t0);
                if (r == LMC_LESS) {
                        // pass, array is not in LMC0 format
                } else {
                        if (verbose >= 2)
                                list[i].showarray ();
                        out.push_back (list[i]);
                }
        }
        return out;
}

arraylist_t selectLMC0 (const arraylist_t &list, int verbose, const conference_t &ctype) {
        arraylist_t out;
        for (size_t i = 0; i < list.size (); i++) {
                lmc_t r = LMC0check (list[i]);

                if (verbose >= 2 || (verbose && i % 10000 == 0))
                        printfd ("selectLMC0: i %d/%d, r %d, total %d\n", i, list.size (), r, out.size ());
                if (r == LMC_LESS) {
                        // pass, array is not in LMC0 format
                } else {
                        if (verbose >= 2)
                                list[i].showarray ();
                        out.push_back (list[i]);
                }
        }
        return out;
}

/** return true if alL is smaller than alR in LMC-0 ordering
 *
 * In LMC0 ordering a column is smaller than another column if
 *
 * i) The zeros occur in earlier positions
 * ii) For equal zeros, use LMC ordering with order 0, 1, -1
 *
 */
bool compareLMC0 (const array_link &alL, const array_link &alR) {
        assert (alL.n_rows == alR.n_rows);
        assert (alL.n_columns == alR.n_columns);

        for (int c = 0; c < alL.n_columns; c++) {
                const array_t *al = alL.array + c * alL.n_rows;
                const array_t *ar = alR.array + c * alR.n_rows;

                // check position of zero(s) in column c
                for (int r = 0; r < alL.n_rows; r++) {
                        if (al[r] == 0 && ar[r] != 0)
                                return true;
                        if (al[r] != 0 && ar[r] == 0)
                                return false;
                }

                // zero is at same position(s) in column, let LMC ordering decide
                for (int r = 0; r < alL.n_rows; r++) {
                        if (al[r] > ar[r])
                                return true; // note the reversed sign here
                        if (al[r] < ar[r])
                                return false; // note the reversed sign here
                }
        }
        // the arrays are equal
        return false;
}

arraylist_t sortLMC0 (const arraylist_t &lst) {
        arraylist_t outlist = lst;
        std::sort (outlist.begin (), outlist.end (), compareLMC0);
        return outlist;
}

/// inflate a list of extensions
std::vector< conference_column > conferenceReduce (const std::vector< conference_column > &ccX, const array_link &als,
                                       const array_link &alfull, const DconferenceFilter &filter,
                                       const conference_t &ct, int verbose) {
        std::vector< conference_column > cci = filter.filterListJ2last (ccX);

        return cci;
}

/** Inflate a list of extensions to extensions for one column extra
 *
 * \param ccX List of candidate extensions
 * \param als Design to be extended
 * \param alfull ??
 * \param filter Class to filter designs
 * \param ct conferece class
 * \param verbose Verbosity level
 */
std::vector< conference_column > extensionInflate (const std::vector< conference_column > &ccX, const array_link &als,
                                       const array_link &alfull, const DconferenceFilter &filter,
                                       const conference_t &ct, int verbose) {
        std::vector< conference_column > cci;
        std::vector< conference_column > cc;

        symmetry_group alfullsg = alfull.row_symmetry_group ();
        const std::vector< int > check_indices = alfullsg.checkIndices ();
        symmetry_group alsg = als.row_symmetry_group ();

        // loop over all candidinates with k columns and inflate to (k+1)-column candidates
        for (size_t i = 0; i < ccX.size (); i++) {
                const conference_column &basecandidate = ccX[i];

                if (verbose > 2)
                        myprintf ("### inflate candidate %d: (sg ngroups %d, sgfull ngroups %d\n", (int)i,
                                  (int)alsg.ngroups, (int)alfullsg.ngroups);
                cc = inflateCandidateExtension (basecandidate, als, alsg, check_indices, ct, verbose, filter);

                if (verbose >= 2) {
                        myprintf ("### inflate: array %d/%d: generated %ld candidates\n", (int)i, (int)ccX.size (),
                                  (long)cc.size ());
                }
                cci.insert (cci.begin (), cc.begin (), cc.end ());
        }
        return cci;
}

std::vector< conference_column > generateDoubleConferenceExtensionsInflate (const array_link &al, const conference_t &ct,
                                                                int verbose, int filterj2, int filterj3, int kstart) {
        // kstart=1;
        double t00 = get_time_ms ();
        int kfinal = al.n_columns;
        if (kstart < 0)
                kstart = al.n_columns - 1;
        assert (kstart >= 1);

        std::vector< conference_column > cci;

        array_link als = al.selectFirstColumns (kstart);

        double t0 = get_time_ms ();
        std::vector< conference_column > ccX = generateDoubleConferenceExtensions (als, ct, verbose, 1, filterj2, filterj3);
        if (verbose)
                printf ("generateDoubleConferenceExtensionsInflate: extend array with %d columns (kfinal %d, initial "
                        "array %d columns)\n",
                        al.n_columns, kfinal, als.n_columns);
        if (verbose)
                printf ("   initial generation: dt %.1f [ms]\n", 1e3 * (get_time_ms () - t0));

        for (int kx = kstart; kx < kfinal; kx++) {
                als = al.selectFirstColumns (kx);
                array_link alx = al.selectFirstColumns (kx + 1);
                DconferenceFilter filter (alx, 1, filterj2, filterj3);

                if (verbose)
                        printf (
                            "## generateDoubleConferenceExtensionsInflate: at %d columns: start with %d extensions\n",
                            kx + 1, (int)ccX.size ());

                cci = extensionInflate (ccX, als, alx, filter, ct, verbose);

                if (verbose) {
                        printf ("## generateDoubleConferenceExtensionsInflate: at %d columns: total inflated: %ld "
                                "candidates for column %d\n",
                                kx + 1, cci.size (), kx + 1);
                        printf ("   dt %.1f [ms]\n", 1e3 * (get_time_ms () - t00));
                }

                ccX = cci;
        }

        return cci;
}

/// return all candidates for the kth column
conference_column_list CandidateGeneratorBase::candidates (int k) { return this->candidate_list[k]; }

CandidateGeneratorBase::CandidateGeneratorBase (const array_link &al, const conference_t &ct_) {
        this->ct = conference_t (ct_);
        this->al = al;
        this->verbose = 1;
        this->last_valid = 0;
        this->candidate_list.clear ();
        this->candidate_list.resize (ct.N + 1); // set a safe max
}

void CandidateGeneratorBase::showCandidates (int verbose) const {
          myprintf ("CandidateGenerator: N %d\n", this->ct.N);
          for (int i = 2; i <= last_valid; i++) {
               myprintf ("CandidateGenerator: number of candidates for %dth column: %ld\n", i,
                         (long)candidate_list[i].size ());
               if (verbose >= 2) {
                         ::showCandidates (candidate_list[i]);
               }
          }
}
CandidateGeneratorConference::CandidateGeneratorConference (const array_link &al, const conference_t &ct_)
    : CandidateGeneratorBase (al, ct_) {
        if (ct_.j1zero != 0) {
                myprintf ("error: j1zero should be zero for conference designs!\n");
        }
}

CandidateGeneratorDouble::CandidateGeneratorDouble (const array_link &al, const conference_t &ct_)
    : CandidateGeneratorBase (al, ct_) // , filter ( DconferenceFilter ( al, 1, 1, 1 ) )
{}

std::vector< conference_column > CandidateGeneratorConference::generateCandidatesZero (const array_link &al, int kz) const {
        const std::vector< conference_column > &cci = this->generateCandidates (al);

        std::vector< conference_column > cci0 = filterZeroPosition (cci, kz);
        return cci0;
}

const std::vector< conference_column > &CandidateGeneratorConference::generateCandidates (const array_link &al) const {
        // assert we have the right settings
        const char *tag = "generateCandidates (conference, zero fixed, cache)";
        const int filterj2 = 1;
        if (ct.j1zero != 0) {
                myprintf ("error: j1zero should be zero for conference designs!\n");
        }
        const int filterj3 = ct.j3zero;
        const int filtersymminline = 1;
        double t00 = get_time_ms ();

        if (verbose >= 2)
                myprintf ("CandidateGenerator::%s: start\n", tag);

        int startcol = this->startColumn (al, 0);
        int kfinal = al.n_columns;
        int ncfinal = al.n_columns + 1;
        int finalcol = kfinal;

        if (verbose >= 2) {
                myprintf ("\n");
                myprintf ("## %s: startcol %d, ncfinal %d\n", tag, startcol, ncfinal);
        }

        std::vector< conference_column > ccX, cci;
        int kstart = -1;

        /* select initial set of columns */
        if (startcol == -1) {
                array_link als = al.selectFirstColumns (START_COL);
                startcol = START_COL + 1;
                int averbose = 0;

                double t0 = get_time_ms ();
                ccX = generateSingleConferenceExtensions (als, ct, -1, averbose, 1, filterj2, filterj3,
                                                          filtersymminline);
                this->candidate_list[START_COL + 1] = ccX;
                last_valid = START_COL + 1;
                kstart = startcol - 1;
        } else {
                ccX = this->candidate_list[startcol];
                last_valid = startcol;
                kstart = startcol - 1;
        }

        /* inflate the columns untill we have reached the target */
        array_link als;
        for (int kx = kstart; kx < kfinal; kx++) {
                als = al.selectFirstColumns (kx);
                array_link alx = al.selectFirstColumns (kx + 1);
                DconferenceFilter filter (alx, 1, filterj2, filterj3);

                if (verbose >= 2)
                        myprintf ("## %s: at %d columns: start with %d extensions, to generate extensions for column "
                                  "%d (%d column array)\n",
                                  tag, kx + 1, (int)ccX.size (), kx + 1, kx + 2);

                cci = extensionInflate (ccX, als, alx, filter, ct, (verbose >= 2) * (verbose - 1));


                if (verbose >= 2) {
                        myprintf ("## %s: at %d columns: total inflated: %ld\n", tag, kx + 1, cci.size ());
                        myprintf ("   dt %.1f [ms]\n", 1e3 * (get_time_ms () - t00));
                }

                ccX = cci;

                this->candidate_list[kx + 2] = ccX;
                this->last_valid = kx + 2;
        }

        if (verbose >= 2)
                printf ("CandidateGenerator::%s: generated %d candidates with %d columns\n", tag,
                        (int)this->candidate_list[ncfinal].size (), ncfinal);

        this->al = al;
        return this->candidate_list[ncfinal];
}

const std::vector< conference_column > &CandidateGeneratorDouble::generateCandidates (const array_link &al) const {
        // assert we have the right settings

        const char *tag = "generateCandidates (double conf matrices, cache)";
        assert (ct.j1zero == 1); // method only valif for j1 and j2 zero
        const int filterj2 = 1;
        const int filterj3 = ct.j3zero;
        double t00 = get_time_ms ();

        int startcol = this->startColumn (al, 0);
        int kfinal = al.n_columns;
        int ncfinal = al.n_columns + 1;
        int finalcol = kfinal;

        if (verbose >= 2) {
                printf ("\n");
                printf ("## %s: startcol %d, ncfinal %d\n", tag, startcol, ncfinal);
        }

        std::vector< conference_column > ccX, cci;
        int kstart = -1;

        /* select initial set of columns */
        if (startcol == -1) {
                array_link als = al.selectFirstColumns (START_COL);
                startcol = START_COL + 1;
                int averbose = 1;

                double t0 = get_time_ms ();
                ccX = generateDoubleConferenceExtensions (als, ct, averbose, 1, filterj2, filterj3);
                this->candidate_list[START_COL + 1] = ccX;
                last_valid = START_COL + 1;
                kstart = startcol - 1;
        } else {
                // NOTE: check this bound is sharp?
                ccX = this->candidate_list[startcol];
                last_valid = startcol;
                kstart = startcol - 1;
        }

        /* inflate the columns untill we have reached the target */
        array_link als;
        for (int kx = kstart; kx < kfinal; kx++) {
                als = al.selectFirstColumns (kx);
                array_link alx = al.selectFirstColumns (kx + 1);
                DconferenceFilter filter (alx, 1, filterj2, filterj3);

                if (verbose >= 2)
                        printf ("## %s: at %d columns: start with %d extensions, to generate extensions for column %d "
                                "(%d column array)\n",
                                tag, kx + 1, (int)ccX.size (), kx + 1, kx + 2);

                cci = extensionInflate (ccX, als, alx, filter, ct, (verbose >= 2) * (verbose - 1));
                cci = filter.filterListZero (cci);

                if (verbose >= 2) {
                        printf ("## %s: at %d columns: total inflated: %ld\n", tag, kx + 1, cci.size ());
                        printf ("   dt %.1f [ms]\n", 1e3 * (get_time_ms () - t00));
                }
                ccX = cci;

                this->candidate_list[kx + 2] = ccX;
                this->last_valid = kx + 2;
        }

        if (verbose >= 2)
                printf ("CandidateGenerator::%s: generated %d candidates with %d columns\n", tag,
                        (int)this->candidate_list[ncfinal].size (), ncfinal);

        this->al = al;
        return this->candidate_list[ncfinal];
}

/// find row sign permutation such that the specified columns has only values +1
void rowlevel_permutation (const array_link &al, rowsort_t *rowperm, const std::vector< int > &colperm,
                           std::vector< int > &rowsignperm, const rowindex_t n_rows, int column) {

        for (rowindex_t r = 0; r < n_rows; r++) {
                if (al.atfast (r, colperm[column]) < 0)
                        rowsignperm[r] = -1;
        }
}

void init_lmc0_rowsort (const array_link &al, int sutk_col, rowsort_t *rowperm, std::vector< int > &rowsignperm,
                        const rowindex_t n_rows, const colindex_t n_cols) {

        for (int i = 0; i < n_rows; i++) {

                int rx = rowperm[i].r;
                int trans_val = (rowsignperm[rx]) * al.at (rx, sutk_col);
                rowperm[i].val = ((trans_val + 3) % 3);
        }
        std::stable_sort (rowperm, rowperm + al.n_rows);
}

/*** compare zero positions in a block of a design
 *
 *
 */
lmc_t lmc0_compare_zeropos_block (const array_link &al, const int x1, const int x2, rowsort_t *rowperm,
                                  const std::vector< int > &colperm, int column, const std::vector< int > &rowsignperm,
                                  const std::vector< int > &colsignperm) {
        for (int i = x1; i < x2; i++) {
                array_t val0 = al.atfast (i, column);
                array_t val = al.atfast (rowperm[i].r, colperm[column]);

                if (val0 != val) {
                        if (val0 == 0) {
                                return LMC_MORE;
                        }
                        if (val == 0)
                                return LMC_LESS;
                }
        }
        return LMC_EQUAL;
}

/* Compare two columns with the zero elements in the same position */
lmc_t compare_conf_columns (const array_link &al, rowsort_t *rowperm, const std::vector< int > &colperm, int column,
                            const std::vector< int > &rowsignperm, const std::vector< int > &colsignperm,
                            const int nrows) {
        int cp = colperm[column];
        int csp = colsignperm[cp];
        for (int i = 0; i < nrows; i++) {
                int value_cdesign_trans = (csp * rowsignperm[rowperm[i].r]) * al.atfast (rowperm[i].r, cp);
                if (((al.atfast (i, column) + 3) % 3) <
                    ((value_cdesign_trans + 3) % 3)) // Transform the elements from (0, 1, -1) to (0, 1, 2)
                        return LMC_MORE;
                if (((al.atfast (i, column) + 3) % 3) > ((value_cdesign_trans + 3) % 3))
                        return LMC_LESS;
        }
        return LMC_EQUAL;
}

lmc_t init_lmc0_sort_comp (const array_link &al, int column, int sel_col, rowsort_t *rowperm,
                           std::vector< int > &rowsignperm, const std::vector< int > &colperm,
                           const std::vector< int > &colsignperm, const int n_rows) {
        lmc_t r = LMC_NONSENSE;
        for (int i = 0; i < n_rows; i++) {

                int rx = rowperm[i].r; 
                int posit_al = al.at (rx, sel_col);
                int trans_val = (colsignperm[sel_col] * rowsignperm[rx]) * posit_al;
                int m = ((trans_val + 3) % 3);
                rowperm[i].val = m;
        }
        std::stable_sort (rowperm, rowperm + n_rows); 

        // Compare zero position
        r = lmc0_compare_zeropos_block (al, 0, n_rows, rowperm, colperm, column, rowsignperm, colsignperm);
        if (r == LMC_LESS || r == LMC_MORE) {
                return r;
        }
        /* Compare whole Columns */
        r = compare_conf_columns (al, rowperm, colperm, column, rowsignperm, colsignperm, n_rows);
        return r;
}

lmc_t LMC0_sortrows_compare (const array_link &al, int column, rowsort_t *rowperm, const std::vector< int > &colperm,
                             const std::vector< int > &rowsignperm, const std::vector< int > &colsignperm,
                             const rowindex_t n_rows, const symmdata &sd) {

        lmc_t r = LMC_NONSENSE;
        const int cp = colperm[column];

        /* Sort rows of the array in blocks*/
        int scol = column - 1;
        int nb = sd.ft.atfast (sd.ft.n_rows - 1, scol); // number of blocks
        /* we check in blocks determined by the ft */
        for (int j = 0; j < nb; j++) {
                int x1 = sd.ft.atfast (2 * j, scol);
                int x2 = sd.ft.atfast (2 * j + 1, scol);

                if ((x2 - x1) > 1) {

                        for (int i = x1; i < x2; i++) {
                                int rx = rowperm[i].r;
                                int current_val = ((colsignperm[cp] * rowsignperm[rx]) * al.atfast (rx, cp));
                                rowperm[i].val = ((current_val + 3) % 3);
                        }

                        flipSort (rowperm, x1, x2 - 1);
                }

                // Compare blocks wrt zero position
                r = lmc0_compare_zeropos_block (al, x1, x2, rowperm, colperm, column, rowsignperm, colsignperm);
                if (r == LMC_LESS || r == LMC_MORE) {
                        return r;
                }
        }

        /* Compare all elements */
        r = compare_conf_columns (al, rowperm, colperm, column, rowsignperm, colsignperm, n_rows);
        return r;
}

/** Recursive function to perform LMC0 test.
 *
 * The rowsignperm is assumed to be fixed. The colsignperm is determined by setting the value for the first row to +1.
 *
 *
 */
lmc_t LMC0_columns (const array_link &al, rowsort_t *rowperm, std::vector< int > colperm, int column,
                    const std::vector< int > &rowsignperm, std::vector< int > colsignperm, const int ncols,
                    const int nrows, const symmdata &sd) {
        // OPTIMIZE: pass colperm and colsignperm as reference?

        lmc_t r = LMC_NONSENSE;

        for (int c = column; c < ncols; c++) {

                std::swap (colperm[c], colperm[column]);

                /* i. Apply the correct column level permutation to make the element X(1,k) equal to 1*/
                int current_sign_col = colsignperm[colperm[column]];
                int current_val_firstrow =
                    (rowsignperm[rowperm[0].r] * current_sign_col) * (al.atfast (rowperm[0].r, colperm[column]));
                colsignperm[colperm[column]] = colsignperm[colperm[column]] * current_val_firstrow;

                /* ii. Sort rows using the ordering 0, 1, -1 and compare */
                r = LMC0_sortrows_compare (al, column, rowperm, colperm, rowsignperm, colsignperm, nrows, sd);

                if (r == LMC_EQUAL) {
                        r = LMC0_columns (al, rowperm, colperm, column + 1, rowsignperm, colsignperm, ncols, nrows,
                                          sd);
                }
                if (r == LMC_LESS) {
                        break; // array is not minimal form
                }

                colsignperm[colperm[column]] = current_sign_col;

                std::swap (colperm[c], colperm[column]);
        }
        return r;
}

lmc_t LMC0checkDC (const array_link &al, int verbose) {
        /*0. Initialize data */
        lmc_t result = LMC_MORE;

        if (!al.is_conference (2)) {
                printfd ("error: input array is not a conference design");
                return LMC_NONSENSE;
        }
        const int ncols = al.n_columns;
        const int nrows = al.n_rows;

        std::vector< int > colperm (ncols);
        init_perm (colperm);

        std::vector< int > colsignperm (ncols);
        init_signperm (colsignperm);

        // row sign permutations are not allowed for double conference matrices, but we have them to re-use code
        std::vector< int > rowsignperm (nrows);
        init_signperm (rowsignperm);

        rowsorter_t rowsorter(nrows);
        rowsort_t *rowsort = rowsorter.rowsort;

        for (rowindex_t i = 0; i < nrows; i++) {
                rowsort[i].val = i;
        }

        symmdata sd (al);

        for (int sel_col = 0; sel_col < ncols; sel_col++) {
                cprintf (verbose, "\nLMC0checkDC: selected column %d\n", sel_col);

                /*1. Select the first (sel_col) column */
                std::swap (colperm[0], colperm[sel_col]);

                int colsign0 = colsignperm[colperm[0]];

                /* 2. Select a possible sign switch for the selected column */
                for (int k = 0; k < 2; k++) {

                        colsignperm[colperm[0]] = colsign0 * ((k == 0) ? 1 : -1);
                        if (verbose >= 2) {
                                myprintf (" # LMC0checkDC: selected sign switch %d: colsignperm ",
                                          colsignperm[colperm[0]]);
                                display_vector (colsignperm);
                                myprintf ("\n");
                        }

                        /* 3. Select one of the possible zero's in the selected column */
                        for (int zidx = 0; zidx < 2; zidx++) {

                                /* 4. Find permutation to sort the array and compare column */
                                result = init_lmc0_sort_comp (al, 0, sel_col, rowsort, rowsignperm, colperm,
                                                              colsignperm, nrows);

                                std::swap (rowsort[0].r, rowsort[zidx].r);

                                if (result == LMC_LESS) {
                                        cprintf (verbose, "LMC0checkDC: init_lmc0_sort_comp result is %d, abort\n",
                                                 result);
                                        return result;
                                }

                                if (verbose >= 2) {
                                        printf ("--- LMC0checkDC: sel_col %d, zidx %d\n", sel_col, zidx);
                                        print_rowsort (rowsort, al.n_rows);
                                }

                                // select the proper column sign permutations
                                for (int c = 1; c < ncols; c++) {
                                        int current_sign_col = colsignperm[colperm[c]];
                                        int value = (rowsignperm[rowsort[0].r] * current_sign_col) *
                                                    (al.atfast (rowsort[0].r, colperm[c]));

                                        if (current_sign_col * value == 0) {
                                                printf (
                                                    "LMC0checkDC: multiple zeros in a single row not supported...\n");
                                                exit (1);
                                        }

                                        cprintf (verbose >= 2, "   LMC0checkDC:  current_sign_col * value %d\n",
                                                 current_sign_col * value);
                                        colsignperm[colperm[c]] = current_sign_col * value;
                                }
                                if (verbose >= 2) {
                                        cprintf (verbose >= 2, " # LMC0checkDC: colsignperm for zidx %d: ", zidx);
                                        display_vector (colsignperm);
                                        myprintf ("\n");
                                }

                                /* 5. Select the next column */
                                result =
                                    LMC0_columns (al, rowsort, colperm, 1, rowsignperm, colsignperm, ncols, nrows, sd);
                                if (result == LMC_LESS) {
                                        cprintf (verbose, "LMC0checkDC: zidx %d, result is %d, abort\n", zidx, result);

                                        return result;
                                }

                                std::swap (rowsort[0].r, rowsort[zidx].r);
                        }
                }
                std::swap (colperm[0], colperm[sel_col]);
                init_signperm (rowsignperm);
        }

        return result;
}

lmc_t LMC0check (const array_link &al, int verbose) {
        /*0. Initialize data */
        lmc_t result = LMC_MORE;

        if (!al.is_conference ()) {
                printfd ("error: input array is not a conference design");
                return LMC_NONSENSE;
        }
        const int ncols = al.n_columns;
        const int nrows = al.n_rows;

        std::vector< int > colperm (ncols);
        init_perm (colperm);

        std::vector< int > colsignperm (ncols);
        init_signperm (colsignperm);

        std::vector< int > rowsignperm (nrows);
        init_signperm (rowsignperm);

        rowsorter_t rowsorter(nrows);
        rowsort_t *rowsort = rowsorter.rowsort;

        for (rowindex_t i = 0; i < nrows; i++) {
                rowsort[i].val = i;
        }

        symmdata sd (al);

        for (int sel_col = 0; sel_col < ncols; sel_col++) {

                /*1. Select the first (sel_col) column */
                std::swap (colperm[0], colperm[sel_col]);

                /*2. Find row-level permutation such that the first column only contains ones */
                rowlevel_permutation (al, rowsort, colperm, rowsignperm, nrows, 0); //

                /* 3. Find permutation to sort the array and compare column */
                result = init_lmc0_sort_comp (al, 0, sel_col, rowsort, rowsignperm, colperm, colsignperm, nrows);
                if (result == LMC_LESS) {
                        return result;
                }

                /* 4. Select one of two possible sign permutations for the first row */
                for (int r_sign = 0; r_sign < 2; r_sign++) {
                        rowsignperm[rowsort[0].r] = 2 * r_sign - 1;

                        /* 5. Select the next column */
                        result = LMC0_columns (al, rowsort, colperm, 1, rowsignperm, colsignperm, ncols, nrows, sd);
                        if (result == LMC_LESS) {
                                return result;
                        }
                }

                std::swap (colperm[0], colperm[sel_col]);
                init_signperm (rowsignperm);
        }

        return result;
}

