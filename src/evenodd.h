/*! \file evenodd.h
 *  \brief Contains functions to generate even-odd designs
 *
 * The generation is done by defining a special ordering in the set of designs.
 * The primary ordering is based in the J5 value of 5-column designs, the secondary ordering is the regular LMC
 * ordering.
 */

#pragma once

#include <algorithm>
#include <map>

#include "arraytools.h"
#include "extend.h"
#include "lmc.h"

/// structure containing current position in search tree
struct depth_path_t {
        /// vector with current position
        std::vector< int > ncurr;  
        /// vector with target
        std::vector< int > nmax;    
        /// number of extension columns
        std::vector< int > necols;  
        /// number of good extension columns
        std::vector< int > ngecols; 
        int depthstart;

        depth_path_t () {}

        void updatePositionGEC (int k, int goodextensioncols) { ngecols[k] = goodextensioncols; }
        void updatePosition (int k, int c, int m, int extensioncols, int goodextensioncols) {
                ncurr[k] = c;
                nmax[k] = m;
                necols[k] = extensioncols;
                ngecols[k] = goodextensioncols;
        }
        void show (int depth, int maxentries = 8) const {
                for (int i = depthstart; i <= depth; i++) {

                        if ((i - depthstart) == maxentries) {
                                myprintf ("... ");
                                break;
                        }

                        myprintf ("%d: %d/%d (%d->%d) ", i, ncurr[i], nmax[i], necols[i], ngecols[i]);
                }
                myprintf ("\n");
        }
        void init (int ncols, int _depthstart = 9) {
                ncurr.resize (ncols + 1);
                nmax.resize (ncols + 1);
                necols.resize (ncols + 1);
                ngecols.resize (ncols + 1);
                depthstart = _depthstart;
        }
};

/// structure to count and show number of arrays generated, the structure is thread safe
struct counter_t {

	// vector with number of arrays found for each column
        std::vector< int > nfound; 

        counter_t (int n);

        void addNfound (int col, int num);

        long nArrays () const ;
        void addNumberFound (int n, int k) ;

        void clearNumberFound ();

        void addNumberFound (const counter_t &de);

        /// show information about the number of arrays found
         void showcountscompact () const;

        /// show information about the number of arrays found
        void showcounts (const arraydata_t &ad) const;

        /// show information about the number of arrays found
        void showcounts (const char *str, int first, int last) const;
};

/** Helper structure for dynamic extension
 *
 * In this structure we keep track of pointers to valid column extensions
 *
 */
struct depth_extend_sub_t {
      public:
        std::vector< int > lmctype;
        /// last column changed in lmc check
        std::vector< int > lastcol;

        std::vector< double > strengthcheck;
        std::vector< int > valididx;
        // param
        int verbose;

        depth_extend_sub_t (int nn = 0) : verbose (0) { resize (nn); };

        void resize (int nn) {
                this->lmctype.resize (nn);
                this->lastcol.resize (nn);
                this->strengthcheck.resize (nn);

                std::fill (this->lmctype.begin (), this->lmctype.begin () + nn, LMC_MORE);
                std::fill (this->lastcol.begin (), this->lastcol.begin () + nn, -1);
        }

        inline size_t n () const { return lmctype.size (); }

        std::vector< int > updateExtensionPointers (int extcol) {
                if (verbose >= 3)
                        myprintf (
                            "updateExtensionPointers: determine extensions that can be used at the next stage\n");

                std::vector< int > pointers;
                //
                for (size_t i = 0; i < lmctype.size (); i++) {
                        // we need proper strength
                        if (strengthcheck[i]) {
                                if (lastcol[i] >= extcol || lastcol[i] == -1 || extcol < 5) {
                                        pointers.push_back (valididx[i]);
                                }
                        }
                }
                if (verbose >= 2)
                        myprintf ("updateExtensionPointers: extcol %d, kept %ld/%ld pointers\n", extcol,
                                  pointers.size (), lmctype.size ());
                return pointers;
        }

        /// initialize the new list of extension columns
        arraylist_t initialize (const arraylist_t &alist, const arraydata_t &adf, const OAextend &oaextend);

        /// select the arrays with are LMC and hence need to be written to disk
        inline arraylist_t selectArraysZ (const arraylist_t &alist) const {
                if (verbose >= 2)
                        myprintf ("depth_extend_sub_t: selectArrays: alist.size() %ld, lmctype %ld\n", alist.size (),
                                  lmctype.size ());
                arraylist_t ga;
                for (size_t i = 0; i < lmctype.size (); i++) {
                        if (verbose >= 3)
                                myprintf ("  depth_extend_sub_t.selectArraysZ: array %ld: lmctype %d\n", i,
                                          lmctype[i]);
                        if (lmctype[i] >= LMC_EQUAL) {
                                array_link ee = alist[i];

                                ga.push_back (ee);
                        }
                }
                if (verbose)
                        myprintf ("dextend_sub_t: selected %d/%d arrays\n", (int)ga.size (), (int)alist.size ());
                return ga;
        }

	inline arraylist_t selectArraysXX (const array_link &al, const arraylist_t &elist) const {
                if (verbose >= 2)
                        myprintf ("depth_extend_sub_t: selectArraysXX: alist.size() %ld, lmctype %ld\n", elist.size (),
                                  lmctype.size ());
                arraylist_t ga;
                for (size_t i = 0; i < n (); i++) {
                        if (verbose >= 3) {
                                myprintf ("  selectArraysXX lmctype %d\n", lmctype[i]);
                        }
                        if (lmctype[i] >= LMC_EQUAL) {
                                array_link ee = hstacklastcol (al, elist[valididx[i]]);

                                ga.push_back (ee);
                        }
                }
                if (verbose >= 1)
                        myprintf ("dextend_sub_t: selected %d/%d arrays\n", (int)ga.size (), (int)elist.size ());
                return ga;
        }


        void info () const {
                size_t number_lmc = 0;
                for (size_t t = 0; t < lmctype.size (); t++) {
                        number_lmc += (lmctype[t] > +LMC_EQUAL);
                }
                if (verbose) {
                        myprintf ("lmc %ld/%d\n", number_lmc, (int)lmctype.size ());
                        myprintf ("valididx size %ld\n", valididx.size ());
                }
        }
};

/** @brief Helper structure for dynamic extension
 *
 * This structure allows for writing the generated arrays to disk.
 * It also contains functions to print progress of the extension.
 *
 * Multiple copies of this class are made, but they all share the same counter_t and arraywriter_t object. Also t0 and
 * tp are shared
 *
 */
struct depth_extend_t {
      public:
        int verbose;
        OAextend oaextend;
        const arraydata_t *ad;

        int loglevelcol;
        /// print progress every x seconds
        double logtime;

        arraylist_t extension_column_list; // list of possible extensions

        /// if set to true write arrays to disk
        int writearrays;

        int discardJ5; /// if true, then we discard the designs which have J5 maximal
        long discardJ5number;

        // shared by mutiple instances of dextend_t (could be made static)
        arraywriter_t *arraywriter;
        counter_t *counter;

        static double t0; // time since start of calculation
        static double tp; // time since last progress report

      private:
        long narraysmax;
        depth_path_t searchpath;

      public:
        // constructure function
        depth_extend_t (const arraydata_t *ad_, double _logtime = 10000000, int _discardJ5 = -1)
            : verbose (1), ad (ad_), discardJ5 (_discardJ5) {
                loglevelcol = -1;
                t0 = get_time_ms ();
                tp = get_time_ms ();

                logtime = _logtime;
                if (ad == 0) {
                        myprintf ("depth_extend_t: pointer to arraydata_t is zero!");
                }

                writearrays = 1;
                narraysmax = LONG_MAX - 1;

                arraywriter = 0;
                counter = 0;

                discardJ5number = 0;

                searchpath.init (ad->ncols);
        };

        depth_extend_t (const depth_extend_t &de) {
                verbose = de.verbose;
                oaextend = de.oaextend;
                ad = de.ad;
                loglevelcol = de.loglevelcol;
                logtime = de.logtime;
                extension_column_list = de.extension_column_list;
                arraywriter = de.arraywriter;
                writearrays = de.writearrays;
                narraysmax = de.narraysmax;
                searchpath = de.searchpath;
                discardJ5 = de.discardJ5;
                discardJ5number = de.discardJ5number;

                counter = de.counter;
        }
        ~depth_extend_t () {
        }

      public:
        void show () { myprintf ("depth_extend_t: logtime %.1f [s]\n", logtime); }

        void setNarraysMax (long n) { this->narraysmax = n; }

        // helper function, thread safe
        void maxArrayCheck () {
                if (arraywriter == 0) {
                        return;
                }
#ifdef DOOPENMP
#pragma omp critical
#endif
                {
                        if (arraywriter->nwritten > this->narraysmax) {
                                myprintf ("dextend_t: number of arrays written: %d, quitting\n",
                                          arraywriter->nwritten);
                                this->counter->showcounts (*this->ad);
                                this->arraywriter->closeafiles ();
                                exit (0);
                        }
                }
        }

      public:
        void showsearchpath (int depth) const { searchpath.show (depth); }

        /// show information about the progress of the loop
        bool showprogress (int showtime = 1, int depth = 0, int forcelog = 0) {
                {
                        double currenttime = get_time_ms ();
                        double dt = currenttime - tp;
                        if ((dt > logtime) || forcelog) {
#ifdef DOOPENMP
#pragma omp critical
#endif
                                tp = get_time_ms ();
                                this->arraywriter->flush ();
                                double dt0 = currenttime - t0;
                                if (showtime) {
                                        int na = this->counter->nArrays ();

#ifdef DOOPENMP
                                        myprintf ("-- depth_extend: progress: %.1f [s], narrays %d (%.1f arrays/s), "
                                                  "thread %d/%d\n",
                                                  dt0, na, na / dt0, omp_get_thread_num (), omp_get_num_threads ());
#else
                                        myprintf ("-- depth_extend: progress: %.1f [s], narrays %d (%.1f arrays/s)\n",
                                                  dt0, na, na / dt0);
#endif

                                        if (depth > 0) {
                                                myprintf ("-- depth %2d: %.1f [s]: ", depth, dt0);
                                                searchpath.show (depth);
                                        }
                                }
                                return true;
                        } else {
                                return false;
                        }
                }
        }
        inline void info () const {
                myprintf ("depth_extend: ");
                ad->show ();
        }

        /// set the position in the dextend structure
        void setposition (int k, int c, int m, int extensioncols = -1, int goodextensioncols = -1) {
#ifdef DOOPENMP
#pragma omp critical
#endif
                { searchpath.updatePosition (k, c, m, extensioncols, goodextensioncols); }
        }
        /// set the position in the dextend structure
        void setpositionGEC (int k, int goodextensioncols) {
#ifdef DOOPENMP
#pragma omp critical
#endif
                { searchpath.updatePositionGEC (k, goodextensioncols); }
        }
};

enum depth_alg_t { DEPTH_DIRECT, DEPTH_EXTENSIONS };

/// Helper structure for the even-odd depth extension
struct depth_extensions_storage_t {

        void resize (size_t s) {
                columnextensionsList.resize (s);
                goodarrayslist.resize (s);
                depthalglist.resize (s);
                dextendsubList.resize (s);
        }

        void set (int ai, const arraylist_t &goodarrays, const arraylist_t &extension_column_list,
                  depth_alg_t depthalg, const depth_extend_sub_t &dextendsub) {
                this->goodarrayslist[ai] = (goodarrays);
                this->columnextensionsList[ai] = (extension_column_list);
                this->depthalglist[ai] = (depthalg);
                this->dextendsubList[ai] = (dextendsub);
        }
        std::vector< arraylist_t > columnextensionsList;
        std::vector< arraylist_t > goodarrayslist;
        std::vector< depth_alg_t > depthalglist;
        std::vector< depth_extend_sub_t > dextendsubList;
};

/** Extend arrays using a depth-first or breadth-first approach
 *
 * @param goodarrays List of arrays to extend
 * @param depthalg Extend using depth-first or breadth-first
 * @param dextend Option structure for the extension
 * @param dextendsublight Data structure for the extensions
 * @param extensioncol Column to extend
 * @param verbose Verbosity level
 */
void processDepth (const arraylist_t &goodarrays, depth_alg_t depthalg, depth_extend_t &dextend,
                   depth_extend_sub_t &dextendsublight, int extensioncol, int verbose = 0);

/// depth-first extension of arrays. depending on the symmetry group of the array to be extended a direct method is
/// used or a method with caching of candidate columns
void depth_extend_hybrid (const arraylist_t &alist, depth_extend_t &dextend, int extcol, const OAextend &oaextendx,
                          int verbose);

/// variation of depth_extend for arrays with large symmetry groups
void depth_extend_direct (const arraylist_t &alist, depth_extend_t &dextend, int extcol, const OAextend &oaextendx,
                          int verbose);

/// depth extend a single array
void depth_extend_array (const array_link &al, depth_extend_t &dextend, const arraydata_t &adfull, int verbose,
                         depth_extensions_storage_t *ds = 0, int = 0);

/// callback function for Pareto calculations
typedef Pareto< mvalue_t< long >, array_link >::pValue (*pareto_cb) (const array_link &, int);
/// callback function for Pareto calculations with cache
typedef Pareto< mvalue_t< long >, array_link >::pValue (*pareto_cb_cache) (const array_link &, int, rankStructure &rs);

template < class IndexType >
inline typename Pareto< mvalue_t< long >, IndexType >::pValue
calculateArrayParetoJ5Cache (const array_link &al, int verbose, rankStructure &rs) {

        const int N = al.n_rows;
        int rank_modelmatrix = rs.rankxf (al);
        mvalue_t< long > a3a4 = A3A4 (al);
        mvalue_t< long > f4 = F4 (al);

        typename Pareto< mvalue_t< long >, IndexType >::pValue p;
        p.push_back (rank_modelmatrix);  
        p.push_back (a3a4); 
        p.push_back (f4); 
        addJmax< IndexType > (al, p, verbose);

        return p;
}

/// add arrays to set of Pareto results
void addArraysToPareto (Pareto< mvalue_t< long >, array_link > &pset, pareto_cb paretofunction,
                        const arraylist_t &arraylist, int jj, int verbose);
/// add arrays to set of Pareto results
void addArraysToPareto (Pareto< mvalue_t< long >, array_link > &pset, pareto_cb_cache paretofunction,
                        const arraylist_t &arraylist, int jj, int verbose);

/** helper class for indexing statistics of designs
 *
 * The index consists of the number of columns and the value for the J-characteristic
 */
struct jindex_t {
		/// number of columns
        int k; 
		/// J-value
        int j; 

        jindex_t (int colindex, int jvalue) : k (colindex), j (jvalue) {}

      public:
        bool operator< (const jindex_t &rhs) const {
                if (this->k < rhs.k) {
                        return true;
                }
                if (this->k > rhs.k) {
                        return false;
                }

                return (this->j < rhs.j);
        }
        std::string toString () const { return printfstring ("k%d-j%d", this->k, this->j); }
};

/// object to hold counts of maximum J_k-values
class Jcounter {
      public:
        /// number of rows
        int N; 
        int jj;
        std::vector< int > fvals;
        std::map< jindex_t, long > maxJcounts;
        /// time needed for calculation
        double dt; 

        Jcounter () : N (-1), jj (-1), dt (0) { }

        Jcounter (int N, int jj = 5, int k = -1) { this->init (N, jj, k); }

		bool validData();

        /// return true if specified column is in the data
        bool hasColumn (int col) const;

        bool isOpen () const { return N > 0; }
        void showPerformance () const {
                myprintf ("Jcounter: %.1f Marrays/h\n", (1e-6 * 3600.) * double(narrays ()) / this->dt);
        }
        long narrays () const;

        /// show statistics of the object
        void show () const;

        int maxCols () const;

        long getCount (int k, int j) const;

        std::vector< long > getTotalsJvalue (int jval) const;
        std::vector< long > getTotals () const;
        
        /// show statistics of the object
        void showcompact () const;

        Jcounter &operator+= (Jcounter &jc);

        /// add list of arrays to object
        void addArrays (const arraylist_t &arraylist, int verbose = 0);

        /// add single array to statistics object
		void addArray(const array_link &al, int verbose = 0);

      private:
		  void init(int N, int jj, int k = -1);
};

/// read statistics object from disk
Jcounter readStatisticsFile (const char *numbersfile, int verbose);

/// write statistics object to disk
void writeStatisticsFile (const char *numbersfile, const Jcounter &jc, int verbose);

/// calculate J-value statistics
Jcounter calculateJstatistics (const char *afile, int jj = 5, int verbose = 1);

/** Return -1 if the first array is smaller in J54 ordering than the second array, 0 if equal and 1 otherwise **/
int compareJ54(const array_link &lhs, const array_link &rhs);
