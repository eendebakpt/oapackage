/** \file arraytools.h

 \brief Contains the array_link class and related classes.

 C++ Interface: arraytools

 This file contains definitions are functions to work with (orthogonal) arrays.
 The code is generic (with templates) and partly inlined for speed.

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#ifndef ARRAYTOOLS_H
#define ARRAYTOOLS_H


#ifdef WIN32
#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable: 4996)
#pragma warning(disable: 4018)
#pragma warning(disable: 4244)
#endif


#ifdef WIN32
#ifdef FULLPACKAGE
#include "msstdint.h"
#endif
#else
#ifdef _WIN32			// || __CYGWIN__
// No visual studio!

#ifdef FULLPACKAGE
#ifndef int32_t
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#endif
#endif

#else
// assume zlib is present on unix
#ifdef NOZLIB
#else
#ifdef FULLPACKAGE
#ifndef USEZLIB
#define USEZLIB 1
#endif
#endif
#endif
#endif
#endif

#ifdef FULLPACKAGE
#include <iostream>
#endif

#include "printfheader.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <deque>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdarg.h>

#include <stdexcept>

#include <Eigen/Core>

namespace Eigen
{
typedef Eigen::Matrix < long double, Eigen::Dynamic,
        Eigen::Dynamic > MatrixXld;
}

/// default float matrix type used
//typedef Eigen::MatrixXf MatrixFloat; typedef Eigen::ArrayXf ArrayFloat; typedef float eigenFloat;


// use double
typedef
Eigen::MatrixXd
MatrixFloat;
typedef
Eigen::ArrayXd
ArrayFloat;
typedef
Eigen::VectorXd
VectorFloat;
typedef double
eigenFloat;

// use long double
//typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixFloat; typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayFloat; typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> VectorFloat; typedef long double eigenFloat;




/// show information about Eigen matrix
//void eigenInfo ( const Eigen::MatrixXd m, const char *str="eigen", int verbose=1 );
void
eigenInfo ( const MatrixFloat m, const char *str = "eigen", int verbose = 1 );


// helper function for Python interface
void
eigen2numpyHelper ( double *pymat1, int n, const MatrixFloat & m );


Eigen::VectorXd dummy ();
Eigen::MatrixXd dummy2 ();

#ifdef USEZLIB
#include <zlib.h>
#endif

#include "oaoptions.h"
#include "mathtools.h"
#ifdef FULLPACKAGE
#include "md5.h"
#include "bitarray/bit_array.h"
#endif



#ifdef SWIG
// only export high level IO functions
%ignore::array_diff;
%ignore::write_array;
%ignore::write_array_latex;
%ignore::finish_arrayfile;
//%ignore append_arrays;
%ignore arrayfile_t::arrayNbits;
%ignore foldtest;
%ignore arraydata_t::complete_arraydata_splitn;
%ignore::writebinheader;
#endif

extern
"C" {
}

//typedef int int32_t;

#ifdef OADEBUG

typedef int
array_t;		/** type of array elements,  should be signed! */

#ifdef SWIGR
typedef int
carray_t;			/* array_t should be signed! */
#else
typedef const int
carray_t;			/* array_t should be signed! */
//typedef int carray_t; /* array_t should be signed! */
#endif

/* change definition below together with array_t !!!! */
#define MPI_ARRAY_T MPI_INT
/*other options for MPI_ARRAY_T are: char: MPI_CHAR, short: MPI_SHORT, int: MPI_INT, long: MPI_LONG */

typedef int
rowindex_t;			/** type used for row indexing */
typedef int
colindex_t;		/** type used for column indexing */
typedef const int
const_colindex_t;		    /** constant version of type used for column indexing */

//typedef long colindex_t; /** type used for column indexing */
//typedef const long const_colindex_t; /** constant version of type used for column indexing */

#else

typedef short int
array_t;			/** type of elements in an orthogonal array *//* array_t should be signed! */
typedef const short int
carray_t;				/** constant version of array_t */

/* change definition below together with array_t !!!! */
#define MPI_ARRAY_T MPI_SHORT
/*other options for MPI_ARRAY_T are: char: MPI_CHAR, short: MPI_SHORT, int: MPI_INT, long: MPI_LONG */

typedef short int
rowindex_t;				/** type used for row indexing */
typedef int
colindex_t;		/** type used for column indexing */
typedef const int
const_colindex_t;		    /** constant version of type used for column indexing */

#endif /* OADEBUG */

typedef array_t *
array_p;			/** pointer to array */
typedef carray_t *
carray_p;			/** point to constant array */

//#define XX
#ifdef XX
  typedef std::vector < int > rowperm_t;				/** type of row permutation */
#else
  typedef rowindex_t * rowperm_t;			/** type of row permutation */
#endif
typedef colindex_t * colperm_t;		       /** type of column permutation */
typedef array_t * levelperm_t;		       /** type of level permutation */


// used to calculate the value (index) of values in a column combination
// this index is used in the strength calculations
// maximum value if of order max(s)*t
typedef int
vindex_t;			/* value index type */

// forward declarations
struct array_link;
struct arraydata_t;


/// possible values for J-values of 2-level design
inline
std::vector < int >
Fval ( int N, int strength )
{
    int x = pow ( ( double ) 2, strength + 1 );	// TODO: replace by integer power
    int nn = floor ( ( double ) N / x ) + 1;
    std::vector < int >Fv ( nn );
    for ( int i = 0; i < nn; i++ ) {
        Fv[i] = N - x * i;
    }
    return Fv;
}


/// return true if the specified file exists
bool file_exists ( const std::string filename );
/// return true if the specified file exists
bool file_exists ( const char *filename );

/// return true if the specified oa file exists
bool oa_file_exists ( const char *filename );

/// return true if the specified oa file exists
bool oa_file_exists ( const std::string filename );

/// create J2 table as intermediate result for J-characteristic calculations for conference matrices
array_link createJ2tableConference ( const array_link & confmatrix );

/// create J2 table as intermediate result for J-characteristic calculations
array_link createJdtable ( const array_link & al );

enum ordering_t {
    ORDER_LEX,
    ORDER_J5
};

/** @brief Contains properties of the design (number of rows, columns, levels)
 *
 * Constructor: arrayclass = arraydata_t(s, N, strength,ncolumns)
 */
struct arraydata_t {
    rowindex_t N;		/** number of runs */
    colindex_t ncols;		/** total number of columns (factors) in the design */
    colindex_t strength;		/** strength of the design */
    array_t *s;	    /** pointer to levels of the array */

    ordering_t order;		/** Ordering used for arrays */

    /* derived data */
    colindex_t ncolgroups;	/// number of groups of columns with the same number of levels
    colindex_t *colgroupindex;	/// specifies for each column the index of the column group
    colindex_t *colgroupsize;
    int oaindex;			/* index of the array */

public:
    /// create new arraydata_t object
    arraydata_t ( array_t s, rowindex_t N, colindex_t strength,
                  colindex_t ncols );
    arraydata_t ( const std::vector < int >s, rowindex_t N,
                  colindex_t strength, colindex_t ncols );
    arraydata_t ( const array_t * s_, rowindex_t N, colindex_t strength,
                  colindex_t ncols );
    arraydata_t ( const arraydata_t & adp );	/// copy constructor

    arraydata_t ( const arraydata_t * adp, colindex_t newncols );	/// copy constructor
    arraydata_t ();		/// dummy constructor

    ~arraydata_t ();		/// destructor

    /// return true if the array is of mixed type
    bool ismixed () const;

    /// return true if the array is a 2-level array
    bool is2level () const;

    /// return random array from the class. this operation is only valid for strength 0 or 1
    array_link randomarray ( int strength = 0, int ncols = -1 ) const;

    /**
     * @brief Write file with design of OA
     * @param file
     * @param ad Arraydata structure to write
     * @return
     */
    void writeConfigFile ( const char *filename ) const;

    /// @brief assignment operator
    inline arraydata_t & operator= ( const arraydata_t & ad2 ) {
        this->N = ad2.N;
        this->strength = ad2.strength;
        this->ncols = ad2.ncols;
        this->order = ad2.order;
        if ( s != 0 ) {
            delete[] s;
        }
        this->s = new array_t[this->ncols];

        if ( ad2.s == 0 ) {
            printf ( "error: invalid arraydata_t structure\n" );
        }
        std::copy ( ad2.s, ad2.s + this->ncols, s );
        return *this;
    }

    /// @brief Comparison operator
    inline int operator== ( const arraydata_t & ad2 ) {
        if ( this->N != ad2.N ) {
            return 0;
        }

        if ( this->ncols != ad2.ncols ) {
            return 0;
        }
        if ( !std::equal ( this->s, this->s + this->ncols, ad2.s ) ) {
            return 0;
        }
        if ( this->strength != ad2.strength ) {
            return 0;
        }
        if ( this->order != ad2.order ) {
            return 0;
        }

        return 1;
    };


    std::string idstr () const;
    std::string idstrseriesfull () const;
    std::string fullidstr ( int series = 0 ) const;
    std::string latexstr ( int cmd = 0, int series = 0 ) const;

public:
    arraydata_t reduceColumns ( int k ) {
        arraydata_t adata ( this, k );
        return adata;
    }
    std::string showstr () const;
    void show ( int verbose = 1 ) const;
    void complete_arraydata ();
    
    /// check whether the LMC calculation will overflow
    void lmc_overflow_check () const;
    
    // complete arraydata but split the column groups at the last column
    void complete_arraydata_fixlast ();

    // complete arraydata but split the column groups at ns
    void complete_arraydata_splitn ( int ns );

    // set column groups at positions given by argument vector
    void set_colgroups ( const std::vector < int >splits );

    // set column group equal to that of a symmetry group + one remaining component
    void set_colgroups_jj ( const symmetry_group & sg, int jj );
    /// set column group equal to that of a symmetry group
    void set_colgroups ( const symmetry_group & sg );
    void show_colgroups () const {
        myprintf ( "arraydata_t: colgroups: " );
        print_perm ( this->colgroupindex, this->ncolgroups );
        myprintf ( "                  size: " );
        print_perm ( this->colgroupsize, this->ncolgroups );
    }

    void calcoaindex ( colindex_t strength ) {
        int combs = 1;
        for ( int i = 0; i < this->strength; i++ ) {
            combs *= this->s[i];
        }

        if ( combs == 0 ) {
            this->oaindex = 0;
        } else {
            this->oaindex = this->N / combs;
        }
    }

    /// return the root array for the class
    array_link create_root () const;

    int getfactorlevel ( int idx ) const {
        if ( idx < 0 ) {
            return -1;
        }
        if ( idx >= this->ncols ) {
            return -1;
        }
        return this->s[idx];
    }

    std::vector < int >getS () const {
        std::vector < int >s ( this->ncols );
        for ( int i = 0; i < this->ncols; i++ ) {
            s[i] = this->s[i];
        }
        return s;
    }

    /**
     * @brief Reset strength of arraydata
     * @param t
     */
    void reset_strength ( colindex_t t ) {
        strength = t;
        delete[]colgroupindex;
        delete[]colgroupsize;
        complete_arraydata ();
    }

    /// Return index of the column group for a column
    colindex_t get_col_group ( const colindex_t col ) const {
        colindex_t j = 0;
        for ( colindex_t i = 0; i < ncolgroups; i++ ) {
            if ( colgroupindex[i] <= col ) {
                j = i;
            } else {
                break;
            }
        }
        return j;
    }

};

/// Read array configuration from file
arraydata_t *readConfigFile ( const char *file );


/**
 * @brief Function similar to printf returning C++ style string
 * @param message
 * @return
 */
inline
std::string
printfstring ( const char *message, ... )
{
    char buf[8 * 1024];

    va_list va;
    va_start ( va, message );
#ifdef RPACKAGE
    myprintf ( "printfstring: not implemented in R\n" );
#else
    vsprintf ( buf, message, va );
#endif
    va_end ( va );

    std::string str ( buf );
    return str;
}


/**
 * @brief Make a copy of an array
 */
inline void
copy_array ( const array_t * src, array_t * const dst, const int nrows,
             const int ncols )
{
    memcpy ( dst, src, sizeof ( array_t ) * nrows * ncols );
}



/**
 * @brief Delete an array
 * @param array
 * @return
 */
inline int
destroy_array ( array_t * array )
{
    free ( array );
    return 0;
}

/**
 * @brief Create an array
 * @param nrows Number of rows
 * @param ncols Number of columns
 * @return
 */
static inline array_t *
create_array ( const int nrows, const int ncols )
{
    //myprintf("  create_array: size %d (%d %d)\n", nrows*ncols, nrows, ncols);
    array_t *array = ( array_t * ) malloc ( nrows * ncols * sizeof ( array_t ) );

#ifdef OADEBUG
    if ( array == NULL ) {
        myprintf ( "problem with malloc %d %d, exiting!!!\n", nrows, ncols );
        throw;
    }
#endif
    return array;
}

/**
 * @brief Create an array from an arraydata_t structure
 */
inline array_t *
create_array ( const arraydata_t * ad )
{
    return create_array ( ad->N, ad->ncols );
}


/**
 * @brief Compare 2 columns of an array
 * @param A
 * @param col
 * @param col2
 * @param nrows
 * @param rstart
 * @param rend
 * @return
 */
inline int
equal_array_cols ( carray_t * A, colindex_t col, colindex_t col2,
                   rowindex_t nrows, rowindex_t rstart, rowindex_t rend )
{
    return std::equal ( A + col * nrows + rstart, A + col * nrows + rend,
                        ( A + col2 * nrows + rstart ) );
}

/**
 * @brief Clone an array
 */
inline array_t *
clone_array ( const array_t * const array, const rowindex_t nrows,
              const colindex_t ncols )
{
    array_t *clone = create_array ( nrows, ncols );
    copy_array ( array, clone, nrows, ncols );

    return clone;
}


/**
 * @brief Perform inverse column permutation on an array
 * @param source
 * @param target
 * @param perm
 * @param nrows
 * @param ncols
 */
inline void
perform_inv_column_permutation ( const array_t * source, array_t * target,
                                 colperm_t perm, int nrows, int ncols )
{
    for ( int i = 0; i < ncols; i++ ) {
        memcpy ( &target[i * nrows], &source[perm[i] * nrows],
                 nrows * sizeof ( array_t ) );
    }
}

inline void
perform_column_permutation ( carray_t * source, array_t * target,
                             colperm_t perm, int nrows, int ncols )
{
    for ( int i = 0; i < ncols; i++ ) {
        memcpy ( &target[perm[i] * nrows], &source[i * nrows],
                 nrows * sizeof ( array_t ) );
    }
}


/**
 * @brief Perform a row permutation
 * @param source Source array
 * @param target Target array
 * @param perm Permutation to perform
 * @param nrows Number of rows
 * @param ncols Numer of columns
 */
inline void
perform_row_permutation ( const array_t * source, array_t * target,
                          rowperm_t perm, int nrows, int ncols )
{
    for ( int i = 0; i < ncols; i++ )
        for ( int j = 0; j < nrows; j++ ) {
            target[nrows * i + perm[j]] = source[nrows * i + j];
        }
}

/// apply inverse row permutation
inline void
perform_inv_row_permutation ( const array_t * source, array_t * target,
                              rowperm_t perm, int nrows, int ncols )
{
    for ( int i = 0; i < ncols; i++ )
        for ( int j = 0; j < nrows; j++ ) {
            target[nrows * i + j] = source[nrows * i + perm[j]];
        }
}

/** Return example array */
array_link exampleArray ( int idx = 0, int verbose = 0 );

/// calculate J-characteristics for a conference design
std::vector<int> Jcharacteristics_conference ( const array_link &al, int jj, int verbose = 0 );

/*! \brief Wrapper class for an array

 The array_link struct is a struct that represents an arrays. Copying of array links
 is done with shallow copy or deep copy depending on compile time options!
  */
struct array_link {
    //! Number of rows in array
    rowindex_t n_rows;
    //! Number of columns in array
    colindex_t n_columns;
    //! Index number
    int index;
    //! Pointer to an array data
    array_t *array;

    static const int INDEX_NONE = 0;
    static const int INDEX_ERROR = -1;
    static const int INDEX_DEFAULT = 0;

    /// Constructor functions
    array_link ();
    array_link ( rowindex_t nrows, colindex_t ncols, int index );
    array_link ( rowindex_t nrows, colindex_t ncols, int index,
                 carray_t * data );
    ~array_link ();
    array_link ( const array_link & );
    array_link ( Eigen::MatrixXd & m );

    array_link clone () const;

public:
    /// print an array to output stream
    friend std::ostream & operator<< ( std::ostream &, const array_link & A );

    /// print array to stdout
    void showarray () const;

    /// print array to stdout
    void showarraycompact () const;

    /// print array properties to stdout
    void showproperties () const;

    /// return true if the arra is a 2-level array (e.g. only contains 0 and 1)
    bool is2level () const;

    /// return true if the array is a +1,0, -1 valued array
    bool is_conference () const;

    /// return true if the array is symmetric
    bool isSymmetric() const;
    
    // manipulation of arrays

    /// make the array symmetric by copying the upper-right to the lower-left
    void makeSymmetric();
    
    /// return array with selected column removed
    array_link deleteColumn ( int index ) const;

    /// return array with first n rows
    array_link selectFirstRows ( int n ) const;

    /// return array with first n columns selected
    array_link selectFirstColumns ( int n ) const;

    /// return array with last n columns selected
    array_link selectLastColumns ( int n ) const;

    /// select columns from an array
    array_link selectColumns ( const std::vector < int >c ) const;

    /// select single column from an array
    array_link selectColumns ( int c ) const;

    /// set a column of the array to the given vector
    void setColumn ( int c, const std::vector < int >v ) {
        std::copy ( v.begin (), v.end (), this->array + c * this->n_rows );
    }
    void setColumn ( int c, const std::vector < signed char >v ) {
        std::copy ( v.begin (), v.end (), this->array + c * this->n_rows );
    }

    /// return transposed array
    array_link transposed () const;

    // statistical properties of the array

    /// calculate D-efficiency
    double Defficiency () const;

    /// calculate main effect robustness (or Ds-optimality)
    double DsEfficiency ( int verbose = 0 ) const;

    /// calculate D-efficiency, calculate main effect robustness (or Ds-optimality) and D1-efficiency
    std::vector < double >Defficiencies ( int verbose = 0, int addDs0 = 0 ) const;

    /*** calculate average variation inflation factor
     *
     * If the VIF is infinite, the value 0 is return. The VIF takes values between 1 and infinity.
     */
    double VIFefficiency () const;

    /// calculate A-efficiency
    double Aefficiency () const;

    /// calculate E-efficiency
    double Eefficiency () const;

    /// Calculate F-values of a 2-level matrix
    std::vector < int >Fvalues ( int jj ) const;

    /// Calculate F-values of a conference design
    std::vector < int >FvaluesConference ( int jj ) const;

    /// Calculate J-characteristics of matrix (the values are signed)
    std::vector < int >Jcharacteristics ( int jj = 4 ) const;

    /// Calculate the projective estimation capacity sequence
    std::vector < double >PECsequence () const;

    /// calculate rank of array
    int rank () const;

    /// calculate generalized wordlength pattern
    std::vector < double >GWLP ( int truncate = 1, int verbose = 0 ) const;

    /// calculate strength of an array
    int strength () const;

    /// return true if the array is a foldover array
    bool foldover () const;

    // return value of minimum element in array
    array_t min () const;
    // return value of maximum element in array
    array_t max () const;


    /** calculate centered L2 discrepancy
     *
     * The method is from "A connection between uniformity and aberration in regular fractions of two-level factorials", Fang and Mukerjee, 2000
     */
    double CL2discrepancy () const;

    /// apply a random permutation of rows, columns and levels
    array_link randomperm () const;
    /// apply a random permutation of columns
    array_link randomcolperm () const;
    /// apply a random permutation of row
    array_link randomrowperm () const;

    /// This function calculates Helmert contrasts for the factors of an input design.
    /// implementation from code written by Eric Schoen, Dept. of Applied Economics, University of Antwerp, Belgium
    MatrixFloat getModelMatrix ( int order, int intercept = 1 ) const;

    /* Interal function (public, but not in documentation */

    array_link & operator= ( const array_link & rhs );	// assignment
    array_link & deepcopy ( const array_link & rhs );	// assignment
    array_link & shallowcopy ( const array_link & rhs );	// assignment
    int operator== ( const array_link & rhs ) const;
    int operator!= ( const array_link & rhs ) const;
    int operator< ( const array_link & rhs ) const;
    int operator> ( const array_link & rhs ) const;

    int equalsize ( const array_link & rhs ) const {
        return ( this->n_rows == rhs.n_rows && this->n_columns == rhs.n_columns );
    }

    array_link operator + ( const array_link & ) const;
    array_link operator + ( array_t v ) const;
    array_link operator - ( const array_link & ) const;
    array_link operator - ( array_t v ) const;

    /// elementwise multiplication
    array_link operator * ( const array_link & rhs ) const;

    array_link operator * ( array_t val ) const {
        array_link al ( *this );
        int NN = this->n_rows * this->n_columns;
        for ( int i = 0; i < NN; i++ ) {
            al.array[i] *= val;
        }
        return al;
    }

    array_link operator *= ( array_t val ) {
        int NN = this->n_rows * this->n_columns;
        for ( int i = 0; i < NN; i++ ) {
            this->array[i] *= val;
        }
        return *this;
    }

    array_link operator += ( array_t val ) {
        int NN = this->n_rows * this->n_columns;
        for ( int i = 0; i < NN; i++ ) {
            this->array[i] += val;
        }
        return *this;


    }
    array_link operator -= ( array_t val ) {
        int NN = this->n_rows * this->n_columns;
        for ( int i = 0; i < NN; i++ ) {
            this->array[i] -= val;
        }
        return *this;
    }

    /// get element from array, no error checking, inline version
    inline const array_t & atfast ( const rowindex_t r, const colindex_t c ) const {
        return this->array[r + this->n_rows * c];
    }
    /// get element from array, no error checking, inline version
    inline array_t & atfast ( const rowindex_t r, const colindex_t c ) {
        return this->array[r + this->n_rows * c];
    }
    array_t _at ( const rowindex_t, const colindex_t ) const;	/// get element at specified position, no error checking
    array_t _at ( const int index ) const;	/// get element at specified position, no error checking

    array_t at ( const rowindex_t, const colindex_t ) const;	/// get element at specified position
    array_t at ( const int index ) const;	/// get element at specified position
    array_t & at ( const rowindex_t, const colindex_t );	/// get element at specified position

    /// set all elements in the array to a value
    void setconstant ( array_t val );

    /// set value of an array
    void setvalue ( int row, int col, int val );
    void setvalue ( int row, int col, double val );	/// set value of an array

    void _setvalue ( int row, int col, int val );	/// set value of an array, no error checking!

    /// multiply a row by -1
    void negateRow ( rowindex_t r ) {
        for ( int c = 0; c < this->n_columns; c++ ) {
            this->atfast ( r, c ) *= -1;
        }
    }
    /// print information about array
    void show () const {
        myprintf ( "index: %d, (%d, %d), array %p\n", index, n_rows, n_columns,
                   ( void * ) array );
    }
    std::string showstr () const {
        std::stringstream s;
        s << "array_link: " << n_rows << ", " << n_columns << "";
        std::string rs = s.str ();
        return rs;			//printfstring("index: %d, (%d, %d), array %p\n", index, n_rows, n_columns, array);
    }

    /// return md5 sum of array representation (as represented with 32bit int datatype in memory)
    std::string md5 () const;

    /// return true if two columns are equal
    bool columnEqual(int rl, const array_link &rhs, int rr) const;
    
    /// return index of first different column
    int firstColumnDifference ( const array_link & A ) const;

    /** calculate row and column index of first difference between two arrays
     *
     * The difference is according to the column-major ordering.
     */
    bool firstDiff ( const array_link & A, int &r, int &c, int verbose = 1 ) const {
        r = 0;
        c = 0;
        for ( c = 0; c < this->n_columns; c++ ) {
            for ( r = 0; r < this->n_rows; r++ ) {
                if ( this->at ( r, c ) != A.at ( r, c ) ) {
                    if ( verbose ) {
                        myprintf ( "first difference of array at %d, %d\n", r, c );
                    }
                    return true;
                }
            }
        }
        return false;
    }

    /// create root in arraylink
    void create_root ( const arraydata_t & ad );

    double nonzero_fraction () const;

    /// fill array with zeros
    void clear () {
        std::fill ( array, array + n_rows * n_columns, 0 );
    };

    // getarraydata (Python interface). this needs to be of type int32 (default python int type)
    void getarraydata ( int *pymat1, int n ) {
        std::copy ( this->array, this->array + n, pymat1 );
    }

    /// internal function
    template < class numtype > void setarraydata ( const numtype * tmp, int n ) {
        //for(size_t i=0; i<n; i++) {
        //    this->array[i]=tmp[i];
        //}
        if ( n != this->n_rows * this->n_columns )
            myprintf
            ( "array_link:setarraydata: warning: number of elements incorrect: n %d, %d %d\n",
              n, this->n_rows, this->n_columns );
        std::copy ( tmp, tmp + n, this->array );
    }
    /// special method for SWIG interface
    void setarraydata ( std::vector < int >tmp, int n ) {
        std::copy ( tmp.begin (), tmp.begin () + n, this->array );
    }
    /// internal function
    template < class numtype > void
    setarraydata ( std::vector < numtype > tmp, int n ) {
        std::copy ( tmp.begin (), tmp.begin () + n, this->array );
    }

    array_t maxelement () const {
        array_t maxelem = 0;
        for ( int i = 0; i < this->n_rows * this->n_columns; i++ ) {
            if ( this->array[i] > maxelem ) {
                maxelem = array[i];
            }
        }
        return maxelem;
    }

    /// set column to values
    void setcolumn ( int c, const array_link & al, int sc = 0 ) {
        assert ( c >= 0 );
        assert ( c <= this->n_columns );
        assert ( this->n_rows == al.n_rows );
        std::copy ( al.array + sc * this->n_rows,
                    al.array + ( sc + 1 ) * this->n_rows,
                    this->array + this->n_rows * c );
    }


public:
    array_link ( const array_link &, const std::vector < int >&colperm );
    array_link ( const array_t * array, rowindex_t nrows, colindex_t ncols,
                 int index = 0 );
    array_link ( const array_t * array, rowindex_t nrows, colindex_t ncolsorig,
                 colindex_t ncols, int index );
    array_link ( const std::vector < int >&v, rowindex_t nrows, colindex_t ncols,
                 int index = 0 );

public:
    void init ( rowindex_t r, colindex_t c );	// made public for python interface

    /// return the row_symmetry group of an array
    symmetry_group row_symmetry_group () const;

#ifdef FULLPACKAGE
    /// return the LMC form of the array
    array_link reduceLMC () const;
    /// return the delete-one-factor-projection form of the array
    array_link reduceDOP () const;
#endif

    /// return the array as an Eigen matrix
    inline MatrixFloat getEigenMatrix () const {
        int k = this->n_columns;
        int n = this->n_rows;
        MatrixFloat mymatrix = MatrixFloat::Zero ( n, k );

        // init array
        for ( int c = 0; c < k; ++c ) {
            int ci = c * n;
            for ( int r = 0; r < n; ++r ) {
                mymatrix ( r, c ) = this->array[r + ci];
            }
        }
        return mymatrix;
    }

    /// return true of specified column is smaller than column in another array
    inline int columnGreater ( int c1, const array_link & rhs, int c2 ) const {

#ifdef OADEBUG
        if ( ( this->n_rows != rhs.n_rows ) || c1 < 0 || c2 < 0
                || ( c1 > this->n_columns - 1 ) ) {
            myprintf
            ( "array_link::columnGreater: warning: comparing arrays with different sizes\n" );
            return 0;
        }
#endif

        int n_rows = this->n_rows;
        return std::lexicographical_compare ( rhs.array + c2 * n_rows,
                                              rhs.array + c2 * n_rows + n_rows,
                                              array + c1 * n_rows,
                                              array + c1 * n_rows + n_rows );

    }

//private:
    std::string showarrayS () const;

#ifdef SWIGCODE
    long data ();			/// return pointer to data, needed for swig interface
#endif

};

// simple permutation type
typedef std::vector < signed char > cperm;
typedef std::vector <cperm > cperm_list;


// concatenate 2 arrays in vertical direction
array_link
hstack ( const array_link & al, const array_link & b );

// concatenate 2 arrays in vertical direction
array_link
hstack ( const array_link & al, const cperm & b );

// concatenate 2 arrays in horizontal direction
array_link
hstack ( const array_link & al, const array_link & b );
// concatenate the last column of array B to array A
array_link
hstacklastcol ( const array_link & A, const array_link & B );


inline cperm vstack ( const cperm & A, const cperm & B )
{
    cperm
    c ( A.size () + B.size () );

    std::copy ( A.begin (), A.end (), c.begin () );
    std::copy ( B.begin (), B.end (), c.begin () + A.size () );
    return c;
}

/// perform column permutation for an array
void perform_column_permutation ( const array_link source, array_link & target,
                                  const std::vector < int >perm );

/// perform row permutation for an array
void perform_row_permutation ( const array_link source, array_link & target,
                               const std::vector < int >perm );

/// create arraydata_t structure from array
arraydata_t arraylink2arraydata ( const array_link & al, int extracols = 0, int strength =  2 );


/// container with arrays
typedef std::deque <array_link > arraylist_t;
/* // //typedef std::vector<array_link> arraylist_t; */

/// add a constant value to all arrays in a list
inline arraylist_t
addConstant ( const arraylist_t & lst, int v )
{
    arraylist_t out ( lst.size () );

    for ( size_t i = 0; i < lst.size (); i++ ) {
        out[i] = lst[i] + v;
    }

    return out;
}

/** Return number of arrays with j_{2n+1}=0 for n<m */
std::vector < int > getJcounts ( arraylist_t * arraylist, int N, int k, int verbose = 1 );


/** @brief Predict j4(1,2,3,k) using the theorem from Deng
 * This works only for 2-level arrays. The 0 corresponds to a +
 *
 */
inline int predictJ ( const array_t * array, const int N, const int k )
{
    int t = N / 4;
    int tt = t / 2;

    // x1 is the number of + entries in (s)
    int x1 = 0;
    for ( int i = 0; i < tt; i++ ) {
        if ( array[k * N + i] == 0 ) {
            x1++;
        }
    }
    for ( int i = tt; i < t; i++ ) {
        // TODO: value for second loop can be calculated from result of first loop
        if ( array[k * N + i] == 1 ) {
            x1++;
        }
    }

    return 8 * x1 - N;
}

#include <map>

/**
 * @brief struct to hold data of an array, e.g. J-characteristic. Abstract base class
 *
 */
class jstructbase_t
{
public:
    std::vector < int > values; // calculated J-characteristics
    std::vector < int > jvalues; // possible values for J-characteristics
    std::map < int, int > jvalue2index; // map from j-value to index
    int jj;

public:
private:

public:
    /// calculate maximum J value
    int maxJ () const;

    /// calculate possible values in F vector
    std::vector < int > Jvalues ( ) const {
        return this->jvalues;
    }

    /// calculate histogram of J values
    std::vector < int >calculateF () const;

    virtual void calc ( const array_link & al ) = 0;

    /// Show contents of structure
    void show ();
    void showdata ( int verbose=1 );
    std::string showstr ();

    /// return 1 if all vals are zero
    int allzero () {
        for ( size_t i = 0; i < this->jvalues.size(); ++i ) {
            if ( this->jvalues[i] != 0 ) {
                return 0;
            }
        }
        return 1;

    }
};


/// structure containing data related to symmetries of arrays
struct symmdata {
public:
    array_link rowvalue;
    array_link orig;

    array_link ft;

    symmdata ( const array_link  &al, int minlen=1 );
    void show ( int verbose=1 ) const {
        printf ( "symmdata: rowvalues\n" );
        this->rowvalue.showarray();
        if ( verbose>=2 ) {
            printf ( "symmdata: ft:" );
            this->ft.show();
            this->ft.showarray();
        }
    }

    /// list with indices set to check for symmetry reductions
    std::vector<int> checkIdx ( int col=-1 ) const {
        const int N = this->orig.n_rows;
        if ( col<0 ) {
            col = orig.n_columns-1;
        }

        std::vector<int> idx ( N );

        // never check first index
        for ( int row=1; row<N; row++ ) {
            if ( this->rowvalue._at ( row, col ) ==this->rowvalue._at ( row-1, col ) ) {
                idx[row]=1;
            }
        }
        return idx;
    }
};
/**
 * @brief struct to hold data of an array, e.g. J-characteristic, rank
 *
 * See papers: Minimum G2-aberration properties of two-level foldover designs, Butler, 2004
 *  Design Selection and Classification for Hadamard Matrices Using Generalized Minimum Aberration Criteria, Deng and Tang
 *
 */
class jstruct_t
{
public:
    /// number of rows
    int N;
    int k;
    int jj;
    int nc;
    std::vector < int >values;
    double A;                   // abberation

public:
    jstruct_t ();
    jstruct_t ( const int N, const int K, const int jj = 4 );
    jstruct_t ( const jstruct_t & js );
    jstruct_t ( const array_link & al, int jj = 4 );
    ~jstruct_t ();

private:
    /// init data structures
    void init ( int N, int k, int jj );
    /// calculate J-characteristics of a 2-level array
    void calc ( const array_link & al );
    /// calculate J-characteristics of a 2-level array, special function for jj=4
    void calcj4 ( const array_link & al );
    /// calculate J-characteristics of a 2-level array, special function for jj=5
    void calcj5 ( const array_link & al );

public:
    jstruct_t & operator= ( const jstruct_t & rhs );    // assignment

    /// calculate maximum J value
    int maxJ () const;

    /// calculate possible values in F vector
    std::vector < int >Fval ( int strength = 3 ) const;

    /// calculate histogram of J values for a 2-level array
    std::vector < int >calculateF ( int strength = 3 ) const;

    // calculate aberration value
    void calculateAberration () {
        // TODO: find reference
        jstruct_t *js = this;
        js->A = 0;
        for ( int i = 0; i < js->nc; i++ ) {
            js->A += js->values[i] * js->values[i];
        }
        js->A /= N * N;
    }
    /// Show contents of structure
    void show ();
    void showdata ();
    std::string showstr ();

    /// return 1 if all vals are zero
    int allzero () {
        for ( int i = 0; i < this->nc; ++i ) {
            if ( this->values[i] != 0 ) {
                return 0;
            }

        }
        return 1;

    }
};

class jstructconference_t : public jstructbase_t
{
public:
    jstructconference_t ( int N, int jj = 4 ) {
        this->jj = jj;
        calcJvalues ( N, 4 );
    }
    jstructconference_t ( const array_link & al, int jj = 4 ) {
        this->jj = jj;
        const int N = al.n_rows;
        calcJvalues ( N, 4 );
        calc ( al );
    }
    //~jstruct_t ();
private:
    void calcJvalues ( int N, int jj ) {
        assert ( jj==4 );
        int nn = floor ( double(  int( ( N-jj+1 ) /4 ) ) ) +1;
        this->jvalues = std::vector<int> ( nn );
        this->jvalue2index.clear();
        for ( size_t i=0; i<jvalues.size(); i++ ) {
            int jval= ( N-jj ) - i*4;
            jvalues[i] = jval;
            jvalue2index[jval] = i;
        }
    }

    void calc ( const array_link &al ) {
        values = Jcharacteristics_conference ( al, this->jj );
    }
};


/// set first columns of an array to root form
void create_root ( array_t * array, const arraydata_t * ad );
/// Creates the root of an OA. The root is appended to the current list of arrays
void create_root ( const arraydata_t * ad, arraylist_t & solutions );

/**
 * @brief Assignment operator
 */
inline array_link &
array_link::operator= ( const array_link & rhs )
{

#ifdef CLEAN_ARRAY_LINK
    return deepcopy ( rhs );
#else
    not implemented
#endif
}

inline array_link &
array_link::shallowcopy ( const array_link & rhs )
{
    //myprintf("array_link::operator= (index %d)\n", rhs.index);
    this->n_rows = rhs.n_rows;
    this->n_columns = rhs.n_columns;
    this->index = rhs.index;
    this->array = rhs.array;
    return *this;
}


inline array_link & array_link::deepcopy ( const array_link & rhs )
{
    this->n_rows = rhs.n_rows;
    this->n_columns = rhs.n_columns;
    this->index = rhs.index;
    // perform deep copy
    if ( this->array ) {
        //  myprintf("  destroy array %d\n", this->array);
        destroy_array ( this->array );
    }
    if ( rhs.array == 0 ) {
        this->array = create_array ( this->n_rows, this->n_columns );
    } else {
        this->array = clone_array ( rhs.array, this->n_rows, this->n_columns );
    }
    return *this;
}

/**
 * @brief Comparision operator for the array link
 */
inline int
array_link::operator< ( const array_link & rhs ) const
{
#ifdef OADEBUG
    if ( ( this->n_rows != rhs.n_rows ) || ( this->n_columns != rhs.n_columns ) ) {
        myprintf
        ( "array_link::operator< comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
          this->index, rhs.index, this->n_rows, this->n_columns,
          rhs.n_rows, rhs.n_columns );
        return 0;
    }
#endif

    return std::lexicographical_compare ( array, array + n_rows * n_columns,
                                          rhs.array,
                                          rhs.array + n_rows * n_columns );
}

/**
 * @brief Comparision operator for the array link
 */
inline int
array_link::operator> ( const array_link & rhs ) const
{
#ifdef OADEBUG
    if ( ( this->n_rows != rhs.n_rows ) || ( this->n_columns != rhs.n_columns ) ) {
        myprintf
        ( "array_link::operator< comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
          this->index, rhs.index, this->n_rows, this->n_columns,
          rhs.n_rows, rhs.n_columns );
        return 0;
    }
#endif

    return std::lexicographical_compare ( rhs.array,
                                          rhs.array + n_rows * n_columns,
                                          array, array + n_rows * n_columns );
}

/**
 * @brief Comparision operator for the array link
 */
inline int
array_link::operator== ( const array_link & b ) const
{
    if ( ( this->n_rows != b.n_rows ) || ( this->n_columns != b.n_columns ) ) {
#ifdef OADEBUG
        myprintf
        ( "array_link::operator== comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
          this->index, b.index, this->n_rows, this->n_columns, b.n_rows,
          b.n_columns );
#endif
        return 0;
    }
    return std::equal ( array, array + n_rows * n_columns, b.array );
}

/**
 * @brief Comparision operator for the array link
 */
inline int array_link::operator!= ( const array_link & b ) const
{
    if ( ( this->n_rows != b.n_rows ) || ( this->n_columns != b.n_columns ) ) {
#ifdef OADEBUG
        myprintf
        ( "array_link::operator== comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
          this->index, b.index, this->n_rows, this->n_columns, b.n_rows,
          b.n_columns );
#endif
        return 0;
    }

    return ( !std::equal ( array, array + n_rows * n_columns, b.array ) );
}

/// Compare 2 arrays and return position of first difference
int array_diff ( carray_p A, carray_p B, const rowindex_t r,
                 const colindex_t c, rowindex_t & rpos, colindex_t & cpos );

/// Compare 2 arrays and return 0 if equal
int array_cmp ( carray_p A, carray_p B, const rowindex_t r,
                const colindex_t c );


/* analyse arrays */

/// helper function to calculate J-values
inline int fastJupdateValue ( rowindex_t N, carray_t * tmpval )
{
    int jval = 0;
    for ( rowindex_t r = 0; r < N; r++ ) {
        jval += tmpval[r] % 2;
    }
    jval = 2 * jval - N;
    return ( jval );
}

/// helper function to calculate J-values
inline void fastJupdate ( const array_t * array, rowindex_t N, const int J, const colindex_t * pp, array_t * tmp )
{
    //int jval=0;

    for ( int i = 0; i < J; i++ ) {
        carray_t *cp = array + N * pp[i];
        for ( rowindex_t r = 0; r < N; r++ ) {
            tmp[r] += cp[r];
        }
    }
    return;
}


/** Calculate J-value for a 2-level array
 */
int jvalue ( const array_link & ar, const int J, const int *pp );

/// calculate J-value for a 2-level array
int jvaluefast ( const array_t * array, rowindex_t N, const int J, const colindex_t * pp );

/// Analyse a list of arrays
std::vector < jstruct_t > analyseArrays ( const arraylist_t & arraylist, const int verbose, const int jj = 4 );

/** \brief Contains a transformation of an array
 *
 * Contains an array transformation. The transformation consists of column, row and
 * level permutations. The level and column permutations are not commutative (since the level permutations
 * are tied to a particular column). We apply the column permutations first.
 *
 */
class array_transformation_t
{
public:
    rowperm_t rperm;		/// row permutation
    colperm_t cperm;		/// column permutation
    levelperm_t *lperms;		/// level permutations
    const arraydata_t *ad;	/// type of array

public:
    array_transformation_t ( const arraydata_t * ad );
    array_transformation_t ( const arraydata_t & ad );
    array_transformation_t ();	/// default constructor
    array_transformation_t ( const array_transformation_t & at );	/// copy constructor
    array_transformation_t & operator= ( const array_transformation_t & at );	/// assignment operator
    ~array_transformation_t ();	/// destructor

    /// show the array transformation
    void show () const;

    /// return true if the transformation is equal to the identity
    bool isIdentity () const;

    /// return the inverse transformation
    array_transformation_t inverse () const;

    /// return the transformation to the identity transformation
    void reset ();

    /// initialize to a random transformation
    void randomize ();

    /// initialize with a random column permutation
    void randomizecolperm ();
    /// initialize with a random row permutation
    void randomizerowperm ();

    /// apply transformation to an array_link object
    array_link apply ( const array_link & al ) const {
        array_link trx ( al );
        this->apply ( al.array, trx.array );
        return trx;
    }
    /// apply transformation to an array_link object
    array_link applygeneric ( const array_link & al ) const {
        // transform
        array_link tmp = al + 1;

        array_link trx ( tmp );
        this->apply ( al.array, trx.array );

        // reverse transformation
        return trx + ( -1 );
    }

    /// composition operator. the transformations are applied from the left
    array_transformation_t operator* ( const array_transformation_t b ) {

        array_transformation_t c ( this->ad );

        const array_transformation_t & a = *this;

        const int nc = this->ad->ncols;

        // perform the rows permutations
        perform_inv_perm ( b.rperm, c.rperm, this->ad->N, a.rperm );

        // perform the column permutations
        perform_inv_perm ( b.cperm, c.cperm, nc, a.cperm );

        /* level permutations */
        for ( colindex_t ci = 0; ci < ad->ncols; ci++ ) {
            levelperm_t l1 = b.lperms[a.cperm[ci]];
            levelperm_t l2 = a.lperms[ci];

            composition_perm ( l1, l2, this->ad->s[ci], c.lperms[ci] );
        }

        return c;
    }

    /// apply transformation to an array
    void apply ( array_t * sourcetarget );

    /// apply transformation to an array
    void apply ( const array_t * source, array_t * target ) const {

        array_t *tmp = create_array ( ad );

        /* column permutations */
        perform_inv_column_permutation ( source, tmp, cperm, ad->N, ad->ncols );

        /* level permutations */
        for ( colindex_t c = 0; c < ad->ncols; c++ ) {
#ifdef SAFELPERM
            safe_perform_level_perm ( tmp + c * ad->N, ad->N, lperms[c], ad->s[c] );
#else
            perform_level_perm ( tmp + c * ad->N, ad->N, lperms[c] );
#endif
        }

        /* row permutations */
        perform_inv_row_permutation ( tmp, target, rperm, ad->N, ad->ncols );

        destroy_array ( tmp );
    }

    /// apply transformation and show resulting array
    void print_transformed ( carray_t * source ) const;

    void show ( std::ostream & out ) const;

    std::vector < int >rowperm () const;        /// return the row permutation of the transformation
    std::vector < int >colperm () const;	/// return the column permutation of the transformation
    std::vector < int >lvlperm ( int c ) const; /// return the level permutations of the transformation


    void setrowperm ( std::vector < int > rp );
    void setcolperm ( std::vector < int > colperm );
    void setlevelperm ( int colindex, std::vector < int > lvlperm );

private:
    void init ();			/// initialize permutation structures
    void free ();			/// free permutation structures and arraydata_t structure
};


/** \brief Contains a transformation of a conference matrix
 *
 * Contains an array transformation. The transformation consists of column permutations, row permutations and sign switches for both the rows and columns.
 *
 * The sign switches and the permutations are not commutative. We apply the permutations first and then the sign flips.
 *
 */
class conference_transformation_t
{
public:
    std::vector < int > rperm;	/// row permutation
    std::vector < int > cperm;	/// column permutation

    std::vector < int > cswitch;	/// sign flips for the columns
    std::vector < int > rswitch;	/// sign flips for the columns

    int nrows;
    int ncols;

public:
    conference_transformation_t ();	/// default constructor
    conference_transformation_t ( int nrows, int ncols );
    conference_transformation_t ( const array_link & al );
    conference_transformation_t ( const conference_transformation_t & T );

    /// show the array transformation
    void show ( int verbose = 1 ) const;

    /// return true if the transformation is equal to the identity
    bool isIdentity () const;

    /// return the inverse transformation
    conference_transformation_t inverse () const;

    /// return the transformation to the identity transformation
    void reset ();

    /// initialize to a random transformation
    void randomize ();

    /// initialize with a random column permutation
    void randomizecolperm ();
    /// initialize with a random row permutation
    void randomizerowperm ();
    /// initialize with random col switches
    void randomizecolflips ();
    /// initialize with random row switches
    void randomizerowflips ();

    /// apply transformation to an array_link object
    array_link apply ( const array_link & al ) const;

    int operator== ( const conference_transformation_t & rhs ) const;

    /** composition operator. the transformations are applied from the left
     * 
     * E.g. (T1*T2)(x) = T1(T2(x))
     * 
     */
    conference_transformation_t operator* ( const conference_transformation_t &rhs ) const {
        const int N = this->nrows;
        const int ncols = this->ncols;

        conference_transformation_t c ( N, ncols );

        const conference_transformation_t & lhs = *this;


        // perform the rows permutations       
        composition_perm ( rhs.rperm, lhs.rperm, c.rperm );
        
        // perform the column permutations
        composition_perm ( rhs.cperm, lhs.cperm, c.cperm );

        /* rowsign switches */
        for ( rowindex_t ri = 0; ri < N; ri++ ) {
            int riz = rhs.rperm[ri];
            int rix = c.rperm[ri];
            c.rswitch[rix] = lhs.rswitch[rix] * rhs.rswitch[riz];
            //myprintf(" ri %d, riz %d, rix %d: lhs.rswitch[rix] %d rhs.rswitch[riz] %d ...\n", ri, riz, rix, lhs.rswitch[rix], rhs.rswitch[riz]);
        }

        /* column sign switches */
        for ( colindex_t ci = 0; ci < ncols; ci++ ) {
            int ciz = rhs.cperm[ci];
            int cix = c.cperm[ci];
            c.cswitch[cix] =  lhs.cswitch[cix] * rhs.cswitch[ciz];
        }

        return c;
    }

    //void show ( std::ostream &out ) const;

    void setrowperm ( std::vector < int >rp ) {
        rperm = rp;
    };
    void setcolperm ( std::vector < int >cp ) {
        cperm = cp;
    };

private:
    void init ( int nr, int nc );	/// initialize permutation structures
};



/** functions for working with array files*/

#ifdef FULLPACKAGE

void showArrayList ( const arraylist_t & lst );

namespace arrayfile
{

/// format mode
enum arrayfilemode_t
{ ATEXT, ALATEX, ABINARY, ABINARY_DIFF, ABINARY_DIFFZERO, AERROR, A_AUTOMATIC };
enum afilerw_t
{ READ, WRITE, READWRITE };

/** @brief Structure for reading or writing a file with arrays
 *
 * The format of array files is described in the file FORMAT.txt
 *
 */
struct arrayfile_t {

public:
    std::string filename;
    int iscompressed;
    int nrows;
    int ncols;

    /// number of bits used when storing an array
    int nbits;

    /// file mode, can be ATEXT or ABINARY, ABINARY_DIFF, ABINARY_DIFFZERO
    arrayfilemode_t mode;
    /// file opened for reading or writing
    afilerw_t rwmode;

    // we cannot define SWIG variables as int32_t, we get errors in the Python module for some reason
    int narrays;

    int narraycounter;

    static const int NARRAYS_MAX = 2 * 1000 * 1000 * 1000;	/* maximum number of arrays in structure */

public:

    /// default constructor
    arrayfile_t ();

    /// open existing array file
    arrayfile_t ( const std::string fname, int verbose = 1 );
    /// open new array file for writing
    arrayfile_t ( const std::string fname, int nrows, int ncols,
                  int narrays = -1, arrayfilemode_t m = ATEXT, int nb = 8 );
    /// destructor function, closes all filehandles
    ~arrayfile_t ();

    /// close current file and open a new file for writing
    void createfile ( const std::string fname, int nrows, int ncols,
                      int narrays = -1, arrayfilemode_t m = ATEXT, int nb = 8 );

    /// close the array file
    void closefile ();
    /// return true if file is open
    int isopen () const;
    /// seek to specified array position
    int seek ( int pos );
    /// read array and return index
    int read_array ( array_link & a );

    /// read next array from the file
    array_link readnext ();

    /// read set of array from the file
    arraylist_t readarrays ( int nmax = NARRAYS_MAX, int verbose = 1 );

    /// flush any open file pointer
    void flush ();

    /// return true if the file has binary format
    bool isbinary () const;

    /// append list of arrays to the file
    int append_arrays ( const arraylist_t & arrays, int startidx = -1 );

    /// append a single array to the file
    void append_array ( const array_link & a, int specialindex = -1 );

    int swigcheck () const {
#ifdef SWIGCODE
        if ( sizeof ( int ) != 4 ) {
            fprintf ( stderr, "arrayfile_t: error int is not 32-bit?!\n" );
        }
        return 1;
#else
        return 0;
#endif
    }

    std::string showstr () const {
        if ( this->isopen () ) {
            std::string modestr;
            switch ( mode ) {
            case ALATEX:
                modestr = "latex";
                break;
            case ATEXT:
                modestr = "text";
                break;
            case ABINARY:
                modestr = "binary";
                break;
            case ABINARY_DIFF:
                modestr = "binary_diff";
                break;
            case ABINARY_DIFFZERO:
                modestr = "binary_diffzero";
                break;
            case AERROR:
                modestr = "invalid";
                break;
            default:
                modestr = "error";
                myprintf ( "arrayfile_t: showstr(): no such mode\n" );
                break;
            }

            int na = narrays;
            if ( this->rwmode == WRITE ) {
                na = narraycounter;
            }

            std::string s =
                printfstring ( "file %s: %d rows, %d columns, %d arrays",
                               filename.c_str (), nrows, ncols, na );
            s +=
                printfstring ( ", mode %s, nbits %d", modestr.c_str (), nbits );
            return s;
        } else {
            std::string s = "file " + filename + ": invalid file";
            return s;
        }
    }

    size_t pos () const {
        return narraycounter;
    }
    /// return true of the file format has random access mode
    bool hasrandomaccess () const {
        return ( this->mode == ABINARY );
    }


private:

public:			// hack
    FILE * nfid;
#ifdef USEZLIB
    gzFile gzfid;
#else
    int gzfid;			// dummy
#endif

    int verbose;

private:
    array_link diffarray;

    /// return header size for binary format array
    int headersize () const {
        return 8 * sizeof ( int32_t );
    };
    /// return size of bit array
    int barraysize () const {
        int num = sizeof ( int32_t );

        switch ( this->nbits ) {
        case 8:
            num += nrows * ncols;
            break;
        case 32:
            num += nrows * ncols * 4;
            break;
        case 1: {
            word_addr_t num_of_words = nwords ( nrows * ncols );	//myprintf("num_of_words: %d\n", (int)num_of_words);
            num += sizeof ( word_t ) * num_of_words;
        }
        break;
        default:
            myprintf ( "error: number of bits undefined\n" );
            break;
        }
        return num;

    };

    /// wrapper function for fwrite or gzwrite
    size_t afwrite ( void *ptr, size_t t, size_t n ) {
        if ( this->nfid == 0 ) {
            myprintf
            ( "afwrite: not implemented, we cannot write compressed files\n" );
            return 0;
        }
        return fwrite ( ptr, t, n, nfid );
    }

    /// wrapper function for fread or gzread
    size_t afread ( void *ptr, size_t sz, size_t cnt );


public:
    // update numbers count for a file structure
    void updatenumbers ();


    /// read array and return index
    int read_array ( array_t * array, const int nrows, const int ncols );

    void finisharrayfile () {
        if ( this->mode == ATEXT ) {
            fprintf ( this->nfid, "-1\n" );
        }
        this->closefile ();	//afile);
    }

    void setVerbose ( int v ) {
        this->verbose = v;
    }

private:
    int read_array_binary_zero ( array_link & a );
    void write_array_binary ( carray_t * array, const int nrows,
                              const int ncols );
    void write_array_binary ( const array_link & A );           // write an array in binary mode to a file
    void write_array_binary_diff ( const array_link & A );  // write an array in binary mode to a file
    void write_array_binary_diffzero ( const array_link & A ); // write an array in binary mode to a file

public:
    int getnbits () {
        return nbits;
    }

    /// parse string to determine the file mode
    static arrayfile::arrayfilemode_t parseModeString ( const std::
            string format ) {
        arrayfile::arrayfilemode_t mode = arrayfile::ATEXT;
        if(format=="AUTO" || format == "A") {
            mode = arrayfile::A_AUTOMATIC;
            
        } else {
        if ( format == "BINARY" || format == "B" ) {
            mode = arrayfile::ABINARY;
        } else {
            if ( format == "D" || format == "DIFF" ) {
                mode = arrayfile::ABINARY_DIFF;
            } else {
                if ( format == "Z" || format == "DIFFZERO" ) {
                    mode = arrayfile::ABINARY_DIFFZERO;
                } else {
                    mode = arrayfile::ATEXT;
                }
            }
        }
        }
        return mode;
    }

    /// return number of bits necessary to store an array
    static int arrayNbits ( const arraydata_t & ad ) {
        int m = 0;
        for ( int i = 0; i < ad.ncols; ++i ) {
            //myprintf("s[i]: %d\n", ad.s[i]);
            if ( ad.s[i] > m ) {
                m = ad.s[i];
            }
        }
        //myprintf("m %d\n", m);
        if ( m == 2 ) {
            return 1;    // bit
        } else if ( m < 120 ) {
            return 8;    // char
        } else {
            return 32;    // int32_t
        }
    }

    /// return number of bits necessary to store an array
    static int arrayNbits ( const array_link & A ) {
        int m = A.max ();
        int amin = A.min ();
        m = std::max ( m, -amin + 1 );
        if ( m == 1 ) {
            return 1;    // bit
        } else if ( m < 124 ) {
            return 8;    // char
        } else {
            return 32;    // int32_t
        }

    };

protected:
    void writeheader ();
    void read_array_binary ( array_t * array, const int nrows,
                             const int ncols );

};

}

using namespace arrayfile;

/// return number of arrays in an array file
long nArrays ( const char *fname );

/// return number of arrays in an array file
inline void
arrayfileinfo ( const char *fname, int &n, int &nr, int &nc )
{
    arrayfile_t af ( fname, 0 );
    n = af.narrays;
    nr = af.nrows;
    nc = af.ncols;
    af.closefile ();
}

/// read list of arrays from file and append to list
int readarrayfile ( const char *fname, arraylist_t * arraylist, int verbose =
                        1, int *setcols = 0, rowindex_t * setrows =
                        0, int *setbits = 0 );
/// read list of arrays from file
arraylist_t readarrayfile ( const char *fname, int verbose = 1, int *setcols =
                                0 );

const int NRAUTO = 0;
/// write a list of arrays to file on disk
int writearrayfile ( const char *fname, const arraylist_t * arraylist,
                     arrayfile::arrayfilemode_t mode =
                         arrayfile::ATEXT, int nrows = NRAUTO, int ncols = NRAUTO );

/// write a list of arrays to file on disk
int writearrayfile ( const char *fname, const arraylist_t arraylist,
                     arrayfile::arrayfilemode_t mode =
                         arrayfile::ATEXT, int nrows = NRAUTO, int ncols = NRAUTO );

/// write a single array to file
int writearrayfile ( const char *fname, const array_link & al,
                     arrayfile::arrayfilemode_t mode = arrayfile::ATEXT );

/// append a single array to an array file. creates a new file if no file exists
int appendarrayfile ( const char *fname, const array_link al );

/// Make a selection of arrays from binary array file, append to list
void selectArrays ( const std::string filename, std::vector < int >&idx,
                    arraylist_t & fl, int verbose = 0 );

/// Select a single array from a file
array_link selectArrays ( std::string filename, int ii );

arrayfile_t *create_arrayfile ( const char *fname, int rows, int cols,
                                int narrays, arrayfile::arrayfilemode_t mode =
                                    arrayfile::ATEXT, int nbits = 8 );

int save_arrays ( arraylist_t & solutions, const arraydata_t * ad,
                  const int n_arrays, const int n_procs,
                  const char *resultprefix, arrayfile::arrayfilemode_t mode =
                      ATEXT );

#endif // FULLPACKAGE

template < class atype >
/// write array to output stream
void
write_array_format ( std::ostream & ss, const atype * array, const int nrows,
                     const int ncols, int width = 3 )
{
    assert ( array != 0 || ncols == 0 );

    int count;
    for ( int j = 0; j < nrows; j++ ) {
        count = j;
        for ( int k = 0; k < ncols; k++ ) {
            const char *s = ( k < ncols - 1 ) ? " " : "\n";
            ss << std::setw ( width ) << array[count] << s;
            count += nrows;
        }
    }
}

/// Make a selection of arrays
arraylist_t selectArrays ( const arraylist_t & al, std::vector < int >&idx );
/// Make a selection of arrays
arraylist_t selectArrays ( const arraylist_t & al, std::vector < long >&idx );

/// Make a selection of arrays, append to list
void selectArrays ( const arraylist_t & al, std::vector < int >&idx,
                    arraylist_t & fl );
void selectArrays ( const arraylist_t & al, std::vector < long >&idx,
                    arraylist_t & fl );

/// Make a selection of arrays, keep
template < class Container, class IntType > void
keepElements ( Container & al, std::vector < IntType > &idx )
{
    for ( int jj = idx.size () - 1; jj >= 0; jj-- ) {
        //myprintf("keepElements: %d\n", jj);
        if ( !idx[jj] ) {
            al.erase ( al.begin () + jj );
        }
    }
}

/// Make a selection of arrays, remove
template < class Container, class IntType > void
removeElements ( Container & al, std::vector < IntType > &idx )
{
    for ( int jj = idx.size () - 1; jj >= 0; jj-- ) {
        if ( idx[jj] ) {
            al.erase ( al.begin () + jj );
        }
    }
}

template < class MType >
/// Make a selection of arrays from a list, append to list
void
selectArraysMask ( const arraylist_t & al, std::vector < MType > &mask,
                   arraylist_t & rl )
{
    assert ( al.size () == mask.size () );
    //myprintf("selectArraysMask: al.size() %d mask.size() %d\n", al.size(), mask.size());
    for ( int idx = 0; idx < al.size (); idx++ ) {
        if ( mask[idx] ) {
            //myprintf("selectArraysMask: idx %d %d\n", idx, mask[idx]);
            rl.push_back ( al.at ( idx ) );
        }
    }
}


template < class IndexType >
/// Append selection of arrays to existing list
void
appendArrays ( const arraylist_t & al,
               const typename std::vector < IndexType > &idx,
               arraylist_t & lst )
{
    for ( typename std::vector < IndexType >::const_iterator it = idx.begin ();
            it < idx.end (); ++it ) {
        lst.push_back ( al.at ( *it ) );
    }
}

/// Append selection of arrays to existing list
inline void
appendArrays ( const arraylist_t & al, arraylist_t & dst )
{
    for ( arraylist_t::const_iterator it = al.begin (); it < al.end (); ++it ) {
        dst.push_back ( *it );
    }
}


/** Write a formatted array
 */
template < class atype > void
write_array_format ( const atype * array, const int nrows, const int ncols,
                     int width = 3 )
{
    int count;
    for ( int j = 0; j < nrows; j++ ) {
        count = j;
        for ( int k = 0; k < ncols; k++ ) {
            const char *s = ( k < ncols - 1 ) ? " " : "\n";
            myprintf ( "%3i%s", static_cast < int > ( array[count] ), s );
            count += nrows;
        }
    }
#ifdef RPACKAGE
#else
#ifdef FULLPACKAGE
    fflush ( stdout );
    setbuf ( stdout, NULL );
#endif
#endif
}

/** @brief Write a formatted array
 */
template < class atype > void
write_array_format ( FILE * fid, const atype * array, const int nrows,
                     const int ncols )
{
    int count;
    for ( int j = 0; j < nrows; j++ ) {
        count = j;
        for ( int k = 0; k < ncols; k++ ) {
            const char *s = ( k < ncols - 1 ) ? " " : "\n";
            fprintf ( fid, "%3i%s", static_cast < int > ( array[count] ), s );
            count += nrows;
        }
    }
}


template < class atype >
/// write an array in latex style
void
write_array_latex ( std::ostream & ss, const atype * array, const int nrows,
                    const int ncols )
{
    int count;

    ss << "\\begin{tabular}{";
    for ( int x = 0; x < ncols; x++ ) {
        ss << 'c';
    }
    ss << "}" << std::endl;

    for ( int j = 0; j < nrows; j++ ) {
        count = j;
        for ( int k = 0; k < ncols; k++ ) {
            const char *s = ( k < ncols - 1 ) ? " & " : " \\\\ \n";
            ss << array[count] << s;
            count += nrows;
        }
    }
    ss << "\\end{tabular}" << std::endl;
}

/// structure to write arrays to disk, thread safe
struct arraywriter_t {
public:

    // since depth_extend is a depth first approach we need to store arrays with a different number of columns
    std::vector < arrayfile_t * >afiles;

    bool writearrays; /// only write arrays if this variable is true
    int nwritten;     /// number of arrays written to disk
    int verbose; /// verbosity level

public:
    arraywriter_t () {
        writearrays = true;
        verbose = 1;
    };

    ~arraywriter_t () {
        flush ();
        closeafiles ();
    }

    void flush () {
        for ( size_t i = 0; i < afiles.size (); i++ ) {
            arrayfile_t *af = afiles[i];
            if ( af != 0 ) {
                #pragma omp critical
                af->updatenumbers ();
                af->flush ();
            }
        }
    }

    /// write a single array to disk
    void writeArray ( const array_link & A ) {
        // writing arrays with multiple threads at the same time is not supported
#ifdef DOOPENMP
        #pragma omp critical
#endif
        {
            size_t i = A.n_columns;
            if ( writearrays ) {
                if ( i < afiles.size () && i>=0 ) {
                    afiles[i]->append_array ( A );
                } else {
                    fprintf ( stderr,
                              "depth_extend_t: writeArray: problem: array file for %d columns was not opened\n",
                              ( int ) i );
                }
                nwritten++;
            }
        }
    }

    // write a list of arrays to disk
    void writeArray ( const arraylist_t & lst ) {
        for ( size_t j = 0; j < lst.size (); j++ ) {
            const array_link & A = lst[j];
            writeArray ( A );
        }
    }

    // initialize the result files
    void initArrayFiles ( const arraydata_t & ad, int kstart,
                          const std::string prefix, arrayfilemode_t mode =
                              ABINARY_DIFF ) {
        afiles.clear ();
        afiles.resize ( ad.ncols + 1 );
        nwritten = 0;

        for ( size_t i = kstart; i <= ( size_t ) ad.ncols; i++ ) {
            arraydata_t ad0 ( &ad, i );
            std::string afile = prefix + "-" + ad0.idstr () + ".oa";
            if ( verbose >= 3 )
                printf ( "depth_extend_t: creating output file %s\n",
                         afile.c_str () );

            int nb = arrayfile_t::arrayNbits ( ad );
            afiles[i] = new arrayfile_t ( afile, ad.N, i, -1, mode, nb );
        }
    }

    /// return the total number arrays written to disk
    int nArraysWritten () const {
        return nwritten;
    }

public:
    void closeafiles () {
        for ( size_t i = 0; i < afiles.size (); i++ ) {
            delete afiles[i];
        }
        afiles.clear ();

    }

};


/// Read header for binary data file. Return true if valid header file
inline bool
readbinheader ( FILE * fid, int &nr, int &nc )
{
    if ( fid == 0 ) {
        return false;
    }

    double h[4];
    int nn = fread ( h, sizeof ( double ), 4, fid );
    nr = ( int ) h[2];
    nc = ( int ) h[3];

    //myprintf("readbinheader: nn %d magic %f %f %f %f check %d %d\n", nn, h[0], h[1], h[2], h[3],  h[0]==30397995, h[1]==12224883);
    bool valid = false;

    // check 2 numbers of the magic header
    if ( nn == 4 && h[0] == 30397995 && h[1] == 12224883 ) {
        return true;
    }

    return valid;
}

/// Write header for binary data file
inline void
writebinheader ( FILE * fid, int nr, int nc )
{
    double h[4];
    // write 2 numbers of the magic header
    h[0] = 30397995;
    h[1] = 12224883;
    h[2] = nr;
    h[3] = nc;
    fwrite ( h, sizeof ( double ), 4, fid );
}

template < class Type >
/// Write a vector of integer elements to file
void
doublevector2binfile ( const std::string fname, std::vector < Type > vals,
                       int writeheader = 1 )
{
    FILE *fid = fopen ( fname.c_str (), "wb" );
    if ( fid == 0 ) {
        fprintf ( stderr, "doublevector2binfile: error with file %s\n",
                  fname.c_str () );
        throw;
    }
    if ( writeheader ) {
        writebinheader ( fid, vals.size (), 1 );
    }
    for ( unsigned int i = 0; i < vals.size (); i++ ) {
        double x = vals[i];
        fwrite ( &x, sizeof ( double ), 1, fid );
    }
    fclose ( fid );

}

/// Write a vector of vector  elements to binary file
inline void
vectorvector2binfile ( const std::string fname,
                       const std::vector < std::vector < double > >vals,
                       int writeheader, int na )
{
    FILE *fid = fopen ( fname.c_str (), "wb" );

    if ( fid == 0 ) {
        fprintf ( stderr, "vectorvector2binfile: error with file %s\n",
                  fname.c_str () );

        throw;
    }

    if ( na == -1 ) {
        if ( vals.size () > 0 ) {
            na = vals[0].size ();
        }
    }
    if ( writeheader ) {
        writebinheader ( fid, vals.size (), na );
    } else {
        myprintf ( "warning: legacy file format\n" );
    }
    for ( unsigned int i = 0; i < vals.size (); i++ ) {
        const std::vector < double >x = vals[i];
        if ( ( int ) x.size () != na ) {
            myprintf
            ( "error: writing incorrect number of elements to binary file\n" );
        }
        for ( unsigned int j = 0; j < x.size (); j++ ) {
            fwrite ( & ( x[j] ), sizeof ( double ), 1, fid );
        }
    }
    fclose ( fid );
}

/* Convertion to Eigen matrices */


/// convert 2-level array to second order interaction matrix in Eigen format
MatrixFloat array2eigenX2 ( const array_link & al );
MatrixFloat array2eigenX1 ( const array_link & al, int intercept = 1 );

/// convert 2-level array to second order model matrix (intercept, X1, X2)
MatrixFloat array2eigenModelMatrix ( const array_link & al );


/// convert array to model matrix in Eigen format
MatrixFloat array2eigenME ( const array_link & al, int verbose = 1 );

/// create first and second order model matrix for mixed-level array
std::pair < MatrixFloat,
    MatrixFloat > array2eigenModelMatrixMixed ( const array_link & al,
            int verbose = 1 );


/// return index of specified array in a file. returns -1 if array is not found
int arrayInFile ( const array_link &al, const char *afile, int verbose=1 );

/// return index of specified array in a list. returns -1 if array is not found
int arrayInList ( const array_link &al, const arraylist_t &ll, int verbose=1 );

#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;


