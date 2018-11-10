#include <string>

#include "arraytools.h"
#include "mathtools.h"
#include "tools.h"

#include "graphtools.h"

#ifdef RPACKAGE
#define printf notallowed
#endif

template < class Type >
/// substract minimum from list of values
std::vector< Type > subtract_minimum (std::vector< Type > &v) {
        std::vector< Type > w (v.size ());
        std::copy (v.begin (), v.end (), w.begin ());
        int rmin = w[0];
        for (typename std::vector< Type >::size_type j = 0; j < w.size (); j++) {
                rmin = std::min (rmin, w[j]);
        }

        for (typename std::vector< Type >::size_type i = 0; i < w.size (); i++)
                w[i] -= rmin;
        return w;
}

/* Interface to Nauty code
 *
 */

namespace nauty {
#include "nauty.h"
/* MAXN=0 is defined by nauty.h, which implies dynamic allocation */

template < class Type >
/// return vector with unique elements
std::vector< Type > uniquevec (const std::vector< Type > &v) {
        std::vector< Type > w = v;
        std::sort (w.begin (), w.end ());
        typename std::vector< Type >::iterator last = std::unique (w.begin (), w.end ());
        w.erase (last, w.end ());
        return w;
}

/// set colors in nauty format
void setcolors (std::vector< int > colors, int *lab, int *ptn) {
        const int verbose = 0;
        const int n = colors.size ();
        std::vector< int > ucols = uniquevec (colors);
        if (verbose)
                myprintf ("setcolors: found %d/%d unique colors\n", (int)ucols.size (), n);
        if (verbose) {
                display_vector (colors);
                myprintf ("\n");
        }
        int x = -1;
        for (size_t i = 0; i < ucols.size (); i++) {
                if (verbose)
                        myprintf ("setcolors: color %d: %d\n", (int)i, ucols[i]);
                for (size_t j = 0; j < colors.size (); j++) {
                        if (colors[j] == ucols[i]) {
                                x++;
                                lab[x] = j;
                                ptn[x] = 1;
                        }
                }
                ptn[x] = 0;
        }
        if (verbose) {
                myprintf ("setcolors: lab and ptn:\n");
                print_perm (lab, n);
                print_perm (ptn, n);
        }
}

std::vector< int > reduceNauty (const array_link &G, std::vector< int > colors, int verbose) {
        if (!G.isSymmetric ()) {
                printfd ("reduceNauty: array is not symmetric, operation not well defined\n");
        }

        if (verbose) {
                myprintf ("reduceNauty: %d vertices\n", G.n_rows);
                myprintf ("  colors: ");
                print_perm (colors);
                myprintf ("\n");
        }
        if ((int)colors.size () != G.n_rows || G.n_rows != G.n_columns) {
                myprintf ("reduceNauty: input sizes not valid");
                return std::vector< int > ();
        }

        int nvertices = G.n_rows;

        /* DYNALLSTAT declares a pointer variable (to hold an array when it
           is allocated) and a size variable to remember how big the array is.
           Nothing is allocated yet.  */

        DYNALLSTAT (graph, g, g_sz);
        DYNALLSTAT (graph, canong, canong_sz);
        DYNALLSTAT (int, lab, lab_sz);
        DYNALLSTAT (int, ptn, ptn_sz);
        DYNALLSTAT (int, orbits, orbits_sz);
        static DEFAULTOPTIONS_GRAPH (options);
        statsblk stats;

        int m;
        // set *gv;

        /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
           Here we change those options that we want to be different from the
           defaults.  writeautoms=TRUE causes automorphisms to be written. */

        options.writeautoms = TRUE;
        options.writeautoms = FALSE;

        options.getcanon = true;

        int n = nvertices;

        /* The nauty parameter m is a value such that an array of
           m setwords is sufficient to hold n bits.  The type setword
           is defined in nauty.h.  The number of bits in a setword is
           WORDSIZE, which is 16, 32 or 64.  Here we calculate
           m = ceiling(n/WORDSIZE). */

        m = SETWORDSNEEDED (n);

        /* The following optional call verifies that we are linking
           to compatible versions of the nauty routines. */

        nauty_check (WORDSIZE, m, n, NAUTYVERSIONID);

        /* Now that we know how big the graph will be, we allocate
         * space for the graph and the other arrays we need. */

        DYNALLOC2 (graph, g, g_sz, m, n, "malloc");
        DYNALLOC2 (graph, canong, canong_sz, m, n, "malloc");

        DYNALLOC1 (int, lab, lab_sz, n, "malloc");
        DYNALLOC1 (int, ptn, ptn_sz, n, "malloc");
        DYNALLOC1 (int, orbits, orbits_sz, n, "malloc");

        EMPTYGRAPH (g, m, n);

        for (int ix = 0; ix < nvertices; ix++) {
                for (int iy = 0; iy < nvertices; iy++) {
                        if (G.atfast (ix, iy) > 0) {
                                if (verbose >= 3) {
                                        myprintf ("adding edge: %d->%d: %d\n", ix, iy, m);
                                }
                                ADDONEEDGE (g, ix, iy, m);
                        }
                }
        }

        setcolors (colors, lab, ptn);
        options.defaultptn = false; // use the coloring!

        if (verbose >= 2) {
                myprintf ("options.defaultptn: %d\n", options.defaultptn);
                myprintf (" lab: \n ");
                print_perm (lab, n);
                myprintf (" ptn: \n ");
                print_perm (ptn, n);
        }

        if (verbose)
                myprintf ("reduceNauty: calling densenauty\n");
        if (verbose >= 2) {
                myprintf ("reduceNauty: options.defaultptn %d, options.getcanon %d, options.digraph %d\n",
                          options.defaultptn, options.getcanon, options.digraph);
        }
        densenauty (g, lab, ptn, orbits, &options, &stats, m, n, canong);
        if (verbose >= 3) {
                myprintf ("Generators for Aut(C[%d]):\n", n);

                myprintf ("order = ");
                writegroupsize (stdout, stats.grpsize1, stats.grpsize2);
                myprintf ("\n");
        }

        if (verbose >= 2) {
                myprintf ("options.defaultptn: %d\n", options.defaultptn);
                myprintf (" lab: \n ");
                print_perm (lab, n);
                myprintf (" ptn: \n ");
                print_perm (ptn, n);
        }

        std::vector< int > tr (nvertices);
        std::copy (lab, lab + nvertices, tr.begin ());

        //  delete allocated variables
        DYNFREE (g, g_sz);
        DYNFREE (canong, canong_sz);
        DYNFREE (lab, lab_sz);
        DYNFREE (ptn, ptn_sz);
        DYNFREE (orbits, orbits_sz);

        return tr;
}

} // end of nauty namespace

array_transformation_t reduceOAnauty (const array_link &al, int verbose) {
        arraydata_t ad = arraylink2arraydata (al);
        return reduceOAnauty (al, verbose, ad);
}

array_transformation_t reduceOAnauty (const array_link &al, int verbose, const arraydata_t &arrayclass) {
        if (verbose >= 2) {
			myprintf ("reduceOAnauty: running on class:\n"); arrayclass.show();
		}
        std::pair< array_link, std::vector< int > > Gc = array2graph (al, verbose, arrayclass);

        array_link &G = Gc.first;
        std::vector< int > &colors = Gc.second;

        if (verbose >= 2) {
                myprintf ("graph:\n");
                G.showarray ();
        }

        std::vector< int > tr = nauty::reduceNauty (G, colors);
        tr = invert_permutation (tr);

        if (verbose >= 2) {
                myprintf ("reduceOAnauty: calculate array_transformation_t from nauty relabelling\n");
        }
        array_transformation_t ttm = oagraph2transformation (tr, arrayclass, verbose >= 2);
        if (verbose >= 2) {
			ttm.show();
			myprintf ("reduceOAnauty: returning array_transformation_t\n");
        }

        return ttm;
}

template < class IntType >
/// helper function
std::vector< int > indexvector (const std::vector< IntType > s) {
        int n = std::accumulate (s.begin (), s.end (), 0);

        std::vector< int > v (n);

        int x = 0;
        for (int i = 0; i < s.size (); i++) {
                for (int j = 0; j < s[i]; j++) {
                        v[x] = i;
                        x++;
                }
        }
        return v;
}

std::pair< array_link, std::vector< int > > array2graph (const array_link &al, int verbose) {
        arraydata_t arrayclass = arraylink2arraydata (al);
        return array2graph (al, verbose, arrayclass);
}

/**  Convert orthogonal array to graph representation
 *
 *   The conversion method is as in Ryan and Bulutoglu.
 *   The resulting graph is bi-partite.
 *   The graph representation can be used for isomorphism testing.
 *
 * For each array the graph consists of:
 *
 * i) For each row a vertex
 * ii) For each column a vertex
 * iii) For each column j there are s_j vertices corresponding to the levels
 *
 * The vertices are labelled according to these 3 groups. Next edges are created:
 *
 * i) For each entry of the array the corresponding row is connected to the column-level vertex
 * ii) Each column vertex is connected to all the column-level vertices with the same column index
 */
std::pair< array_link, std::vector< int > > array2graph (const array_link &al, int verbose,
                                                         const arraydata_t &arrayclass) {
        int nrows = al.n_rows;
        int ncols = al.n_columns;
        const std::vector< int > s = arrayclass.factor_levels ();
		
        int nRowVertices = nrows;
        int nColVertices = ncols;
        int nColumnLevelVertices = std::accumulate (s.begin (), s.end (), 0);

        int nVertices = nrows + ncols + nColumnLevelVertices; // number of vertices in the graph
        int colOffset = nrows;

        std::vector< int > vertexOffsets (s.size () + 1);

        std::vector< int > cs = cumsum0 (s);

        for (size_t j = 0; j < s.size (); j++)
                vertexOffsets[j] = colOffset + ncols + cs[j];

        std::vector< int > colors (nVertices);

        symmetry_group sg (arrayclass.factor_levels (), 0);

        // row vertices: color 0
        for (int i = 0; i < nrows; i++)
                colors[i] = 0;
        // column vertices: color 1 + column group index
        for (int i = 0; i < ncols; i++)
                colors[i + nrows] = 1 + sg.gidx[i];

        // other vertices: color 2
        for (int i = 0; i < nColumnLevelVertices; i++)
                colors[i + nrows + ncols] = (arrayclass.ncolgroups - 1) + 2;

        if (verbose)
                myprintf ("array2graph: generating graph of size %d=%d+%d+%d\n", nVertices, nrows, ncols,
                          nColumnLevelVertices);
        array_link G (nVertices, nVertices, 0); // graph
        G.setconstant (0);

        for (int r = 0; r < nrows; r++) {
                for (int c = 0; c < ncols; c++) {
                        int idx = al.at (r, c) + vertexOffsets[c];
                        G.at (r, idx) = 1;
                        G.at (idx, r) = 1;
                }
        }

        if (nColVertices > 0) {
                int colidx = 2;
                for (int col = 0; col < ncols; col++) {
                        for (int sx = 0; sx < s[col]; sx++) {
                                int sel = vertexOffsets[col] + sx;
                                G.at (colOffset + col, sel) = colidx;
                                G.at (sel, colOffset + col) = colidx;
                        }
                }
        }

        return std::pair< array_link, std::vector< int > > (G, colors);
}

/// apply a vertex permutation to a graph
array_link transformGraph (const array_link &G, const std::vector< int > tr, int verbose) {
        array_link H = G;

        for (int i = 0; i < H.n_rows; i++) {
                for (int j = 0; j < H.n_columns; j++) {
                        int ix = tr[i];
                        int jx = tr[j];
                        H.at (ix, jx) = G._at (i, j);
                        ;
                }
        }
        return H;
}

/// From a relabelling of the graph return the corresponding array transformation
array_transformation_t oagraph2transformation (const std::vector< int > &pp, const arraydata_t &arrayclass,
                                               int verbose) {
        if (arrayclass.ismixed ()) {
                printfd ("note: oagraph2transformation not tested for mixed-level designs\n");
                arrayclass.show ();
        }
        /// invert the labelling transformation
        std::vector< int > ppi = invert_permutation (pp);

        // extract colperms and rowperm and levelperms from pp
        array_transformation_t ttr (arrayclass);

        if (verbose) {
                myprintf ("labelling2transformation: class %s\n", arrayclass.idstr ().c_str ());
                myprintf ("labelling2transformation: pp  ");
                display_vector (pp);
                myprintf ("\n");
                myprintf ("labelling2transformation: ppi ");
                display_vector (ppi);
                myprintf ("\n");
        }
        const int N = arrayclass.N;
        std::copy (pp.begin (), pp.begin () + N, ttr.rperm);
        int rmin = pp.size ();
        for (int j = 0; j < N; j++)
                rmin = std::min (rmin, (int)ttr.rperm[j]);
        for (int i = 0; i < N; i++)
                ttr.rperm[i] -= rmin;
        ttr = ttr.inverse ();

        if (verbose) {
                myprintf ("labelling2transformation: rowperm ");
                print_perm (ttr.rperm, N);
        }

        int ncols = arrayclass.ncols;
        array_transformation_t ttc (arrayclass);
        std::vector< int > colperm (arrayclass.ncols);
        std::copy (pp.begin () + N, pp.begin () + N + ncols, colperm.begin ());
        if (verbose >= 2) {
                printfd ("colperm: ");
                display_vector (colperm);
                myprintf ("\n");
        }
        colperm = subtract_minimum (colperm);
        colperm = invert_permutation (colperm);
        ttc.setcolperm (colperm);

        if (verbose) {
                printfd ("labelling2transformation: colperm ");
                display_vector (colperm);
                myprintf ("\n");
        }

        std::vector< int > s = arrayclass.factor_levels ();

        int ns = std::accumulate (s.begin (), s.end (), 0);
        array_transformation_t ttl (arrayclass);

        std::vector< int > lvlperm (ns);
        std::copy (pp.begin () + N + ncols, pp.begin () + N + ncols + ns, lvlperm.begin ());
        lvlperm = subtract_minimum (lvlperm);

        std::vector< int > cs = cumsum0 (s);

        for (int ii = 0; ii < ncols; ii++) {
                std::vector< int > ww (lvlperm.begin () + cs[ii],
                                       lvlperm.begin () + cs[ii + 1]); //  = lvlperm[cs[ii]:cs[ii + 1]]

                indexsort is (ww);

                if (verbose >= 2) {
                        printfd ("index sorted: ");
                        display_vector (is.indices);
                        myprintf ("\n");
                }
                ww = is.indices;

                if (verbose >= 1) {
                        printfd ("oagraph2transformation: lvlperm %d: ", ii);
                        display_vector (ww);
                        myprintf ("\n");
                        fflush (0);
                }
                ttl.setlevelperm (ii, ww);
        }

        ttl = ttl.inverse ();

        array_transformation_t tt = ttr * ttc * ttl;
        return tt;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
