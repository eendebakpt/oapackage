/*! \file extend.h
 *  \brief Contains functions to generate and extend orthogonal arrays.
 *
 */

#ifndef EXTEND_H
#define EXTEND_H

#ifdef OAEXTEND_MULTICORE
#include <mpi.h>
#endif
#include "oaoptions.h"
#include "arraytools.h"
#include "tools.h"
#include "lmc.h"

#define stringify( name ) # name


#ifdef SWIG
%ignore extendpos;
%ignore progress_column;
%ignore check_block_exchange;
%ignore init_column;
#endif



/** @brief Options for the extend code
 *
 * class containing parameters of the extension and LMC algorithm
 *
 */
class OAextend
{
public:
    //! time before printing progress of single extension, [seconds]
    double singleExtendTime;
    //! number of arrays LMC tested before printing progress of single extension
    int nLMC ;
    //! perform LMC test after generation of array
    int checkarrays;
    //! adds a symmetry check to the extension algorithm based in symmetry of row permutations
    int use_row_symmetry;

    /// init column with previous column in extension (if in the same column group)
    int init_column_previous;

    /// append full array, append only extension column, store array to disk, or do nothing
    enum {APPENDEXTENSION, APPENDFULL, STOREARRAY, NONE};
    /// determined how the arrays are stored
    int extendarraymode;
    arrayfile_t storefile;	 //NOTE: we should make a copy constructor and assignment operator

    // special cases
    j5structure_t j5structure;

private:
    /// Algorithm mode
    algorithm_t algmode;	// MODE_ORIGINAL: original, MODE_J4: j4 check, ...

public:
    OAextend() : singleExtendTime(10.0), nLMC(40000), checkarrays(1), use_row_symmetry(1), init_column_previous(1), extendarraymode(APPENDFULL), j5structure(J5_45), algmode(MODE_ORIGINAL)
    {
#ifdef OADEV
        algmode = MODE_AUTOSELECT;
#endif
    };
    OAextend( const OAextend &o) : singleExtendTime(o.singleExtendTime)
    {
        this->nLMC = o.nLMC;
        this->checkarrays = o.checkarrays;
        this->use_row_symmetry = o.use_row_symmetry;
        this->init_column_previous = o.init_column_previous;
        this->extendarraymode = o.extendarraymode;
        this->j5structure = o.j5structure;

        this->algmode=o.algmode;
        // we do NOT copy the storefile: this->storefile = o.storefile;

    };
    OAextend( arraydata_t &ad) : singleExtendTime(10.0), nLMC(40000), checkarrays(1), use_row_symmetry(1), init_column_previous(1), extendarraymode(APPENDFULL), j5structure(J5_45), algmode(MODE_ORIGINAL)
    {

#ifdef OADEV
        algmode = MODE_AUTOSELECT;
#endif
        setAlgorithmAuto(&ad);
    };
    /// Set the algorithm to use for LMC checks
    void setAlgorithm(algorithm_t algorithm,  arraydata_t *ad = 0);
    /// Set the algorithm automatically
    void setAlgorithmAuto( arraydata_t *ad = 0);

    /// Return algorithm used
    algorithm_t getAlgorithm() const
    {
        return this->algmode;
    };
    /// Return algorithm used (as string)
    std::string getAlgorithmName() const
    {
        return algnames(this->algmode);
    };

    void updateArraydata(arraydata_t *ad = 0) const;

    /// return preferred extension algorithm
    static inline algorithm_t getPreferredAlgorithm(const arraydata_t &ad, int verbose=0)
    {
        //verbose=2;
        if (verbose)
            myprintf("getPreferredAlgorithm: ad.ncolgroups %d, ad.s[0] %d\n", ad.ncolgroups, ad.s[0]);
        if (ad.ncolgroups==1 && ad.s[0]==2 && (ad.strength==3) ) {
            //printf(" using MODE_J4\n");
            return MODE_J4;
        } else return MODE_ORIGINAL;
    }

    /// print configuration to stdout
    void info(int vb=1) const
    {
        std::cout << __repr__();

        if (vb>1) {
            myprintf("OAextend: use_row_symmetry: %d\n", this->use_row_symmetry );
        }
        std::cout << std::endl;
    }
    std::string __repr__() const;
};

struct  extend_data_t;


/*!
 *	Struct holds information abount the current position in the array and some characteristics of the array. This reduces
 *	number of parameters for a function.
 \brief Position Information
 */
struct extendpos {
    /// Current row of cell being filled
    rowindex_t	row;
    /// Current column being filled
    colindex_t	col;
    /// Used to test the possibility of a value, leaves the OA in tact
    array_t	value;
    /// Array data
    arraydata_t *ad;

    extendpos(colindex_t extensioncol, arraydata_t *adp): row(0), col(extensioncol), ad(adp) {};

#ifdef OADEBUG
    void show()
    {
        myprintf("extendpos struct: N %d, col %d, ncols %d\n", this->ad->N, this->col, this->ad->ncols);
    }
#endif
};

/* functions */

double progress_column(array_t *column, extendpos *p);


/* Public part of interface */

/// extend a list of arrays
int extend_arraylist(const arraylist_t & alist, arraydata_t &fullad,   OAextend const &oaextend, colindex_t extensioncol, arraylist_t &extensions);

/// extend a list of arrays
arraylist_t extend_arraylist(const arraylist_t & alist, arraydata_t &fullad,   OAextend const &oaextend);

/// extend a list of arrays with default options
arraylist_t extend_arraylist(const arraylist_t & alist, const arraydata_t &fullad);

/// extend a single array
arraylist_t extend_array(const array_link &al, arraydata_t &fullad,   OAextend const &oaextend);

/// extend a single array with the default LMC algorithm
arraylist_t extend_array(const array_link &al, arraydata_t &arrayclass);

/// extend an array with a single column
int extend_array(carray_t *array, const arraydata_t *, const colindex_t extensioncol, arraylist_t &solutions, OAextend const &oaextend );

/// simple wrapper function
arraylist_t runExtendRoot(arraydata_t adata, int nmax, int verbose=0);


/* Helper functions and development code */

enum {DFILTER_NONE, DFILTER_BASIC, DFILTER_MULTI};

enum {DCALC_ALWAYS, DCALC_COND};


/** @brief Helper structure for dynamic extension
 *
 *
 *
 */
struct dextend_t {

    static const int NO_VALUE = 0;

    std::vector<lmc_t> lmctype;
    /// last column changed in lmc check
    std::vector<int> lastcol;

    /// A values
    std::vector<double> Deff;

    /// indices of filtered arrays
    std::vector<int> filter;


    /// check mode: 0: no filter, 1: classic, 2: predict
    int filtermode;
    /// DCALC_ALWAYS: always calculate A, DCALC_COND: only for LMC_LESS
    int Dcheck;

    /// perform immediate LMC check in extension
    int directcheck;

    dextend_t() : filtermode(DFILTER_MULTI), Dcheck(DCALC_COND), directcheck(1) {};

    void resize(int nn)
    {
        this->lmctype.resize(nn);
        this->lastcol.resize(nn);
        this->filter.resize(nn);
        this->Deff.resize(nn);

        std::fill(this->lmctype.begin(), this->lmctype.begin() + nn, LMC_MORE);
        std::fill(this->lastcol.begin(), this->lastcol.begin() + nn, -1);
    }

    /// perform filtering using D-efficiency
    void DefficiencyFilter(double Dfinal, int k,int kfinal, double Lmax, int verbose=1);

    /// filter the arrays based on values in filter
    std::vector<int> filterArrays(const array_link &al, const arraylist_t &earrays, arraylist_t &earraysout, std::vector<std::vector<double> > &edata, int verbose=1);

    // output
    long ntotal; /// total number of arrays found
    long nlmc; /// total number of arrays found in LMC form
    long n; /// total number of arrays found passing all tests

    double DmaxDiscard;
    long nmaxrnktotal;	// total number of arrays found with max rank
};



/*!  \brief Branche Information for extension algorithm
 *
 * This struct is used to keep track of the branches of different possibilities. If more than one option is available,
 * this location should be stored at st[count] and count should be increased. For short, this could be done by
 * st[count++]. If a column if filled, the highest element in use, gives the last branching point. Struct is used to
 * represent a stack, since no random access is needed.
 */
struct split {
    ///	Pointer to the stack
    rowindex_t	*st;		//alternative for a stack object
    ///	Counter for number of elements on the stack and the number of options at the last spot, is obsolete?
    rowindex_t	count;
    /// size of stack
    rowindex_t stacksize;

    /// number of valid positions
    array_t *nvalid;
    array_t *cvalidpos;
    ///  valid positions
    array_t **valid;


    /**
     * @brief Constructor for stack structure
     * @param stacksize
     * @param maxs Maximum value in array
     */
    split(rowindex_t stacksize, int maxs = 26)
    {
        log_print(DEBUG+1, "split constructor: size %d\n", stacksize);
        this->stacksize = stacksize;
        this->st = new rowindex_t [stacksize];
        this->nvalid = new array_t [stacksize];
        this->cvalidpos = new array_t [stacksize];
        this->valid = malloc2d<array_t>(stacksize, maxs+1);
        this->count = 0;		//initialize stack
    }

    ~split()
    {
        delete [] this->st;
        delete [] this->nvalid;
        delete [] this->cvalidpos;
        free2d(this->valid); //, this->stacksize);
    }

    inline rowindex_t position()
    {
        return this->count - 1;
    };

    /// show information about stack
    void show()
    {
        myprintf("stack: %d elements\n", this->count);
        for(int i=0; i<this->count; i++) {
            myprintf("row %d (pos %d): ", this->st[i], this->cvalidpos[i]);
            print_perm(this->valid[i], this->nvalid[i]);
        }
    }
};



#endif

