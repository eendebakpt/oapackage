/*! \file extend.h
 *  \brief Contains functions to generate and extend orthogonal arrays.
 *
 */

#ifndef EXTEND_H
#define EXTEND_H

#ifdef OAEXTEND_MULTICORE
#include <mpi.h>
#endif
#include "arraytools.h"
#include "lmc.h"
#include "oaoptions.h"
#include "tools.h"

#ifdef SWIG
%ignore extendpos;
#endif

/** @brief Options for the extend code
 *
 * class containing parameters of the extension and LMC algorithm
 *
 */
class OAextend {
      public:
        //! time before printing progress of single extension, [seconds]
        double singleExtendTime;
        //! number of arrays LMC tested before printing progress of single extension
        int nLMC;

        //! perform LMC test after generation of array
        int checkarrays;

        //! if true then return at once if a single extension has been found
        int check_maximal;

        //! adds a symmetry check to the extension algorithm based in symmetry of row permutations
        int use_row_symmetry;

        /// init column with previous column in extension (if in the same column group)
        int init_column_previous;

        /// Specification of how to use the generated extensions
        enum extendarray_mode_t {
			/// append extension column to extension list
			APPENDEXTENSION,
			/// append full array to extension list
			APPENDFULL,
			/// store extension to disk
			STOREARRAY,
			/// do not store generated extensions
			NONE };

        /// determines how the extension arrays are stored
		extendarray_mode_t extendarraymode;

        arrayfile_t storefile; 

        // special cases
        j5structure_t j5structure;

      private:
        /// Algorithm mode
        algorithm_t algmode; 

      public:
	/** Options for the extension algorithm
	 * 
	 * 
	 */
        OAextend ()
            : singleExtendTime (10.0), nLMC (40000), checkarrays (1), check_maximal (0), use_row_symmetry (1),
              init_column_previous (1), extendarraymode (APPENDFULL), j5structure (J5_45), algmode (MODE_AUTOSELECT){};
	/// @copydoc OAextend()
        OAextend (const OAextend &o) : singleExtendTime (o.singleExtendTime) {
                this->nLMC = o.nLMC;
                this->checkarrays = o.checkarrays;
                this->use_row_symmetry = o.use_row_symmetry;
                this->init_column_previous = o.init_column_previous;
                this->extendarraymode = o.extendarraymode;
                this->j5structure = o.j5structure;
                this->check_maximal = o.check_maximal;

                this->algmode = o.algmode;
                // we do not copy the storefile: this->storefile = o.storefile;
        };
	/** @copydoc OAextend()
	 * 
	 * The algorithm is automatically determined from the specified arrayclass.
	 * 
	 */
        OAextend (arraydata_t &arrayclass)
            : singleExtendTime (10.0), nLMC (40000), checkarrays (1), check_maximal (0), use_row_symmetry (1),
              init_column_previous (1), extendarraymode (APPENDFULL), j5structure (J5_45), algmode (MODE_AUTOSELECT) {
                setAlgorithmAuto (&arrayclass);
        };
        /// Set the algorithm to use for LMC checks
        void setAlgorithm (algorithm_t algorithm, arraydata_t *ad = 0);
        /// Set the algorithm automatically
        void setAlgorithmAuto (arraydata_t *ad = 0);

        /// Return algorithm used
        algorithm_t getAlgorithm () const { return this->algmode; };
        /// Return algorithm used (as string)
        std::string getAlgorithmName () const { return algnames (this->algmode); };

		/// update the options structuer with the specified class of designs
        void updateArraydata (arraydata_t *arrayclass = 0) const;

		/** Return preferred extension algorithm
		 *
		 * \param arrayclass Class of designs to extend
		 * \param verbose Verbosity level
		 * \return Algorithm selected to be used for this class
		 *
		 **/
        static algorithm_t getPreferredAlgorithm (const arraydata_t &arrayclass, int verbose = 0) {
                if (verbose)
                        myprintf ("getPreferredAlgorithm: ad.ncolgroups %d, ad.s[0] %d\n", arrayclass.ncolgroups, arrayclass.s[0]);
                if (arrayclass.ncolgroups == 1 && arrayclass.s[0] == 2 && (arrayclass.strength == 3)) {
                        return MODE_J4;
                } else
                        return MODE_ORIGINAL;
        }

        /// print configuration to stdout
		void info(int verbose = 1) const;

        std::string __repr__ () const;
};

struct extend_data_t;

/*!
 *	Struct holds information abount the current position in the array and some characteristics of the array. This
 reduces
 *	number of parameters for a function.
 \brief Position Information
 */
struct extendpos {
        /// Current row of cell being filled
        rowindex_t row;
        /// Current column being filled
        colindex_t col;
        /// Used to test the possibility of a value, leaves the OA in tact
        array_t value;
        /// Array data
        arraydata_t *ad;

        extendpos (colindex_t extensioncol, arraydata_t *adp) : row (0), col (extensioncol), ad (adp){};

        void show () {
                myprintf ("extendpos struct: N %d, col %d, ncols %d\n", this->ad->N, this->col, this->ad->ncols);
        }
};

/** Extend a list of orthogonal arrays
*
* \param array_list The list of arrays to be extended
* \param array_class Class of arrays to generate
* \param oaextend_options Parameters for the extension algorithm
* \return List of all generated arrays
* 
* @see extend_array(const array_link &, arraydata_t &, OAextend const &)
*/
arraylist_t extend_arraylist (const arraylist_t &array_list, arraydata_t &array_class, OAextend const &oaextend_options);

/** Extend a list of arrays with default options
*
* @see extend_array(const array_link &, arraydata_t &, OAextend const &)
*/
arraylist_t extend_arraylist (const arraylist_t &array_list, const arraydata_t &arrayclass);

/** @copydoc extend_arraylist(const arraylist_t &, arraydata_t &, OAextend const &)
 * 
 * \param extensioncol Index of column to be added to the designs
 * \param extensions List to append generated designs to
 * \return Number of candidate arrays generated
 */
int extend_arraylist (const arraylist_t &array_list, arraydata_t &array_class, OAextend const &oaextend_options, colindex_t extensioncol,
                      arraylist_t &extensions);

/** Extend a single orthogonal array
 *
 * \param array The array to be extended
 * \param array_class Class of arrays to generate
 * \param oaextend Parameters for the extension algorithm
 *
 */
arraylist_t extend_array (const array_link &array, arraydata_t &array_class, OAextend const &oaextend);

/** Extend a single orthogonal array with the default LMC algorithm
 *
 * @see extend_array(const array_link &, arraydata_t &, OAextend const &)
 */
arraylist_t extend_array (const array_link &array, arraydata_t &arrayclass);

/** Extend an orthogonal array with a single column
 *
 * @see extend_array(const array_link &, arraydata_t &, OAextend const &)
 *
 * @param array Array to extend
 * @param arrayclass Array data for the full array
 * @param extension_column Column to extend
 * @param extensions List to which generated valid extensions are added
 * @param oaextend Structure with options
 * @return Number of candidate extensions generated
 */
int extend_array (const array_link &array, const arraydata_t *arrayclass, const colindex_t extension_column, arraylist_t &extensions,
                  OAextend const &oaextend);

/** Run the LMC extension algorithm starting with the root array 
 *
 * @see extend_array(const array_link &, arraydata_t &, OAextend const &)
 */
arraylist_t runExtendRoot (arraydata_t arrayclass, int max_number_columns, int verbose = 0);


enum dfilter_t {
  /// no filtering on D-efficiency
  DFILTER_NONE,
  /// filtering on D-efficiency
  DFILTER_BASIC,
  /// filtering on D-efficiency with multi column prediction
  DFILTER_MULTI };

enum dcalc_mode { 
  /// always calculate efficiency
  DCALC_ALWAYS,
  /// only calculate efficiency for LMC_LESS
  DCALC_COND };

/** @brief Structure for dynamic extension of arrays based on D-efficiencies
 *
 */
struct dextend_t {

        static const int NO_VALUE = 0;

		/// results of minimal form calculations
        std::vector< lmc_t > lmctype;
        /// last column changed in lmc check
        std::vector< int > lastcol;

        /// calculated efficiency values
        std::vector< double > Deff;

        /// indices of filtered arrays
        std::vector< int > filter;

		dfilter_t filtermode;

		dcalc_mode Dcheck;

        /// perform immediate LMC check in extension
        int directcheck;

        dextend_t () : filtermode (DFILTER_MULTI), Dcheck (DCALC_COND), directcheck (1){};

        void resize (int nn) {
                this->lmctype.resize (nn);
                this->lastcol.resize (nn);
                this->filter.resize (nn);
                this->Deff.resize (nn);

                std::fill (this->lmctype.begin (), this->lmctype.begin () + nn, LMC_MORE);
                std::fill (this->lastcol.begin (), this->lastcol.begin () + nn, -1);
        }

        /// perform filtering using D-efficiency
        void DefficiencyFilter (double Dfinal, int k, int kfinal, double Lmax, int verbose = 1);

        /// filter the arrays based on values in filter
        std::vector< int > filterArrays (const array_link &al, const arraylist_t &earrays, arraylist_t &earraysout,
                                         std::vector< std::vector< double > > &edata, int verbose = 1);

        /// total number of arrays found
        long ntotal;
        /// total number of arrays found in LMC form
        long nlmc;
        /// total number of arrays found passing all tests
        long n;

        double DmaxDiscard;
        long nmaxrnktotal; // total number of arrays found with max rank
};


#endif
