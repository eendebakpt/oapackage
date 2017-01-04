/*! \file tools.h
 *  \brief Contains various functions
 *
 * Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 * Copyright: See LICENSE.txt file that comes with this distribution
 */
#ifndef TOOLS_H
#define TOOLS_H

#include <stdarg.h>
#include <list>
#include <ostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <time.h>
#include <algorithm>
#include <iterator>

#ifdef FULLPACKAGE
#endif
#include "printfheader.h"

inline std::string base_name ( std::string const & path )
{
    return path.substr ( path.find_last_of ( "/\\" ) + 1 );
}

inline void printfd_handler ( const char *file, const char* func, int line, const char* message, ... )
{
    std::string s = file;
    s=base_name ( s );

    const char *fileshort = s.c_str();
    myprintf ( "file %s: function %s: line %d: ", fileshort, func, line );
#ifdef FULLPACKAGE
	char buf[64 * 1024];

    va_list va;
    va_start ( va, message );
    //vprintf ( message, va );
    vsprintf ( buf, message, va );
    va_end ( va );
	myprintf("%s", buf);
#else
    myprintf("printfd_handler not implemented");
#endif
}

//#define printfd(MESSAGE) printfd_handler(__FILE__, __LINE__, MESSAGE)
#define printfd(...) printfd_handler(__FILE__,__FUNCTION__, __LINE__, __VA_ARGS__)

#include "arraytools.h"
#include "mathtools.h"

#if (_MSC_VER >= 100)
#pragma warning(disable: 4018)
#pragma warning(disable: 4996)
#endif

#ifdef SWIG
%ignore nullStream;
%ignore logstream;
%ignore next_comb;
%ignore next_comb2;
%ignore init_restart;
%ignore set_colcombs_fixed;
%ignore addelement;
%ignore row_rank_partial;
%ignore row_rank;
#endif

#ifdef _WIN32
//#define prefetch(x)
#else
#ifndef ARCH_HAS_PREFETCH
#ifndef prefetch
#ifdef DOPREFETCH
// has some conflict with Eigen
#define prefetch(x) __builtin_prefetch(x)
#endif
#endif
#endif
#endif


/* do these work with recent GCC? */
//#pragma GCC diagnostic warning "-Wformat"
//#pragma warning (disable : 4018) // warning C4018: '==' : signed/unsigned mismatch
//#if defined __GNUC__
//#pragma DISABLE_WARNINGS
//#endif


//! Loglevel definitions. The loglevel determines to amount of output to stdout
enum loglevel_t {LOGERROR, SYSTEM, QUIET, NORMAL, DEBUG, EXTRADEBUG};


int log_print ( const int level, const char *message, ... );

//struct split;
int getloglevel();
void setloglevel ( int n );
bool checkloglevel ( int l );


/** \brief Null stream */
class nullStream : public std::ostream
{
public:
    nullStream () : std::ostream ( NULL ) {}
};

#ifdef FULLPACKAGE
std::ostream& logstream ( int level );
#endif
//ostream &streamloglevel(ostream &stream);

std::string system_uname();

inline void mycheck_handler ( const char *file, const char* func, int line, int condition, const char* message, ... )
{
    if ( condition==0 ) {
        va_list		va;
        va_start ( va, message );
        myprintf ( "mycheck: %s: %s (line %d): ", file,func, line );
#ifdef RPACKAGE
        myprintf("(not implemented) %s", message);
#else
        vprintf ( message, va );
#endif
        va_end ( va );
//  myprintf ( "mycheck %d: %s", condition, str);
#ifdef RPACKAGE
        throw;
#else
        exit ( 1 );
#endif
    }


}


#define mycheck(...) mycheck_handler(__FILE__,__FUNCTION__, __LINE__, __VA_ARGS__)

inline void myassert ( int condition, const char *str = 0 )
{
    if ( condition==0 ) {
        if (str==0)
            myprintf ( "myassert: error\n" );
        else
            myprintf ( "myassert: %s", str );
#ifdef RPACKAGE
        throw;
#else
        exit ( 1 );
#endif
    }
}

#ifdef OADEBUG
inline void myassertdebug ( int condition, const char *str )
{
    if ( condition==0 ) {
        myprintf ( "myassert: %s", str );
        myprintf ( "... aborting\n" );
        exit ( 1 );
    }
}
#else
#define myassertdebug(a,b)
inline void myassertdebug2 ( int condition, const char *str ) {}
#endif

inline int cprintf ( int check, const char *message, ... )
{
    int n=0;
    if ( check ) {
        va_list va;
        va_start ( va, message );
#ifdef RPACKAGE
        n = -1;
        myprintf("cprintf: not implemented\n");
#else
        n = vprintf ( message, va );
#endif
        va_end ( va );
    }
    return n;
}

/// flush to stdout
inline void ff()
{
#ifdef RPACKAGE
#else
    fflush ( stdout );
#endif
}

/******************/

template<class A>
/**
 * Delete a pointer and set to zero.
 */
inline void safedelete ( A *p )
{
    if ( p!=0 )
        delete p;
    p=0;
}

template<class A>
/**
 * Delete array and set pointer to zero
 * @param p
 */
inline void safedeletearray ( A *p )
{
    if ( p!=0 )
        delete [] p;
    p=0;
}

template <class Type>
/*!
  Gives next combination for k elements out of n based on an algorithm from wikipedia.
  The generate is sorted.
  \brief Gives next combination
  \param comb Pointer to combination
  \param k Number of the current combination
  \param n Number of elements in combination
  */
int next_comb ( std::vector<Type> &comb, int k, int n )
{
    int             i;// = k - 1;
    const int       offset = n - k + 1;
    //comb[k - 1] = n - 1;
    i = k - 1;
    comb[i]++;
    while ( ( comb[i] >= offset + i ) && ( i > 0 ) ) {
        i--;
        comb[i]++;
    }

    if ( comb[0] > n - k )
        return 0; /* No more combinations can be generated */

    /* comb now looks like (…, x, n, n, n, …, n).
       Turn it into (…, x, x + 1, x + 2, …) */
    for ( i++; i < k; i++ )
        comb[i] = comb[i - 1] + 1;

    return 1;
}


/** Go to next combination in sequence */
int next_comb ( int *comb, int k, int n );
/** Go to next combination in sequence */
int next_comb_s ( int *comb, int k, int n );

/*!
  Calculates the rank of a row, in order to determine the sorting order.
  \brief Calculates row rank
  \todo Combine n_columns and n_rows into *nrs, to reduce overhead
  \param array
  \param n_columns
  \param n_rows
  \param index
  */
static inline int row_rank ( array_t *array, const int n_columns, const int n_rows, const int *index )
{
    register int	i, sum = 0, j = 0;
    for ( i = 0; i < n_columns; i++ ) {
        sum += index[i] * array[j];
        j += n_rows;
    }
    return sum;
}

/**
 * @brief Returns the value of (part of) a row
 * @param array
 * @param start_idx
 * @param end_idx
 * @param n_rows
 * @param index Value index for each of the columns of the array
 * @return
 */
static inline int row_rank_partial ( const carray_t *array, rowindex_t n_rows, const vindex_t *index, colindex_t start_idx, colindex_t end_idx, const colperm_t &colperm, rowindex_t row )
{
    int	sum = 0;
    const array_t *ar = array+row;
    for ( colindex_t i = start_idx; i < end_idx; i++ ) {
        sum += index[i] * ar[n_rows*colperm[i]];
    }
    return sum;
}

/**
 * @brief Returns the value of (part of) a row
 * @param array
 * @param start_idx
 * @param end_idx
 * @param n_rows
 * @param index Value index for each of the columns of the array
 * @return
 */
static inline array_t row_rank_partial ( carray_t *array, colindex_t start_idx, colindex_t end_idx, rowindex_t n_rows, const vindex_t *index )
{
    array_t	sum = 0;
    int j = 0;
    j += n_rows*start_idx;
    for ( colindex_t i = start_idx; i < end_idx; i++ ) {
        sum += index[i] * array[j];
        j += n_rows;
    }
    return sum;
}


/**
 * @brief Returns the value of (part of) a row
 * @param array
 * @param start_idx
 * @param end_idx
 * @param n_rows
 * @param index Value index for each of the columns of the array
 * @return
 */
static inline array_t row_rank_partial ( carray_t *array, const colindex_t start_idx, const colindex_t end_idx, const rowindex_t row, const rowindex_t n_rows, const int *index )
{
    register int	i, sum = 0, j = row;
    j += n_rows*start_idx;
    for ( i = start_idx; i < end_idx; i++ ) {
        sum += index[i] * array[j];
        j += n_rows;
    }
    return sum;
}


template <class Object>
/**
 * @brief Template to swap two objects of arbitrary datatype
 * Please use std::swap instead
 * @param a
 * @param b
 */
void swap_object ( Object &a, Object &b )
{
    Object tmp;
    tmp = a;
    a = b;
    b = tmp;
}

template <class rtype>
/** Calculate the number of elements in a 2D table with rows with different sizes */
inline int malloc2d_nelements ( const int nrows, const rtype *rowsizes )
{
    int nelements = 0;
    for ( int i=0; i<nrows; i++ )
        nelements += rowsizes[i];
    return nelements;
}


template <class DataType, class rtype>
/**
 * @brief Allocate a 2-dimensional array with non-uniform rows
 * @param nrows Number of rows in the table
 * @param rowsizes Size of each row
 * @return
 */
DataType **malloc2d_irr ( const int nrows, const rtype *rowsizes )
{
    //Create a 2D array, but with unequal rows (hence irregular -> irr)
    register int i;
    DataType **data;

    int nelements = malloc2d_nelements ( nrows, rowsizes );

    data = new DataType* [nrows];
    data[0] = new DataType [nelements];
    memset ( data[0] , 0, sizeof ( DataType ) * nelements );

    int offset = 0;
    for ( i = 0; i < nrows; i++ ) {
        data[i] = data[0]+offset;
        offset += rowsizes[i];
    }

    return data;
}

template <class DataType, class rtype>
/**
 * @brief Allocate a 2-dimensional array with non-uniform rows, return size of allocated space
 * @param nrows Number of rows in the table
 * @param rowsizes Size of each row
 * @param nelements This parameter is initialized with the size of the array allocated
 * @return
 */
DataType **malloc2d_irr ( const int nrows, const rtype *rowsizes, int &nelements )
{
    nelements = malloc2d_nelements ( nrows, rowsizes );
    return malloc2d_irr<DataType> ( nrows, rowsizes );
}

template <class DataType, class numtype>
/**
 * @brief Allocate a 2-dimensional array of specified size
 * @param nrows
 * @param rowsize
 * @return
 */
DataType **malloc2d ( const numtype nrows, const int rowsize )
{
    DataType **data;

    data = new DataType* [nrows];
    // if debugging, check for memory allocation
    if ( data==0 ) {
        myprintf ( "malloc2d: error with memory allocation\n" );
        throw;
        //exit(0);
    }

    //myprintf("nrows*rowsize: %d * %d = %d\n", nrows, rowsize, nrows*rowsize);
    data[0] = new DataType [nrows*rowsize];

    int offset = 0;
    for ( int i = 0; i < nrows; i++ ) {
        data[i] = data[0]+offset;
        offset += rowsize;
    }

    return data;
}

template <class DataType>
/**
 * @brief Release a 2-dimensional array
 * @param data
 * @param nrows
 */
void free2d ( DataType **data, const int nrows )
{
    delete [] data[0];
    delete [] data;
    data = 0;
}

template <class DataType>
/**
 * @brief Release a 2-dimensional array
 * @param data
 */
void free2d ( DataType **data )
{
    delete [] data[0];
    delete [] data;
    data = 0;
}

template <class DataType>
/**
 * @brief Release a 2-dimensional non-uniform array
 * @param data
 */
void free2d_irr ( DataType **data )
{
    free2d ( data );
}
template <class DataType>
/**
 * @brief Release a 2-dimensional non-uniform array
 * @param data
 * @param nrows
 */
void free2d_irr ( DataType **data, const int nrows )
{
    free2d ( data );
}


//void show_array(carray_t *array, const int x, const int y);
void print_array ( const char *str, const array_t *array, const rowindex_t r, const colindex_t c );
void print_array ( const array_t *array, const rowindex_t r, const colindex_t c );
/// Print array to stdout
void print_array ( const array_link &A );

#ifdef FULLPACKAGE
template <class atype>
/// print vector
void display_vector ( const std::vector<atype> &v )
{
    const char *sep = " ";
    std::copy ( v.begin(), v.end(), std::ostream_iterator<atype> ( std::cout, sep ) );
}
#else
template <class atype>
void display_vector ( const std::vector<atype> &v )
{
// dummy
}
#endif

template <class atype>
/// print vector
void printf_vector ( const std::vector<atype> &v, const char *format )
{
    for ( unsigned int i=0; i<v.size(); i++ )
        myprintf ( format, v[i] );
}

#ifdef FULLPACKAGE
template <class atype>
void show_array_dyn ( const atype *array, const int x, const int y )
{
    register int	i,j,k;

    for ( i = 0; i < y; i++ ) {
        k = i;
        for ( j = 0; j < x; j++ ) {
            std::cout << std::setw ( 3 ) <<  array[k];
            //log_print(NORMAL, "%3i",(int) array[k]);
            k += y;
        }
        std::cout << "\n";
    }
}
#endif

/// Counts the number of occurences of each value in an array
void countelements ( carray_t* array, const int nelements, const int maxval, int* elements );

/**
 * @brief Add element to element counter
 * @param elem
 * @param elements
 */
inline void addelement ( const array_t elem, int* elements )
{
    //myprintf("adding element as position %d\n", elem);
    elements[elem]++;
}


#ifdef OAANALYZE_DISCR
enum {ANALYSIS_NONE, ANALYSIS_DISCRIMINATOR};
void analysis_print ( const int level, const char *message, ... );

void analyse_discriminant ( int row, int col, lmc_t, int nr, int nc );
void print_discriminant ( int nr, int nc );
void clear_discriminant ( int nr, int nc );
void init_discriminant ( int nr, int nc );

void analysis_init_values();
void analysis_increase_counter ( const std::string &p );
void analysis_show_counter ( const std::string &p );

#endif






/// return time with milisecond precision
double get_time_ms();

/// return time difference with milisecond precision
double get_time_ms ( double t0 );

/// trim a string by removing the specified characters from the left and right
void trim ( std::string& str, const std::string& trimChars = "" );

/// return the current time as a string
inline std::string currenttime()
{
    time_t seconds;
    struct tm *tminfo;
    time ( &seconds );
    tminfo = localtime ( &seconds );
    std::string ts = asctime ( tminfo );
    trim ( ts );
    return ts;
}



/// return string describing array
std::string oafilestring ( const arraydata_t *ad );
/// return string describing array
std::string oafilestring ( rowindex_t rows, colindex_t cols, array_t *s );



template <class numtype>
/** @brief Convert integer to C++ string
 *
 * @param i Integer
 * @return String representation of the integer
 */
inline std::string itos ( numtype i )
{
    std::stringstream s;
    s << i;
    return s.str();
}

/// printf-style function that returns std::string
std::string printfstring ( const char *message, ... );


inline std::string printtime()
{
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    return printfstring ( "%s", asctime ( timeinfo ) );
}

#ifdef OAANALYZE_DISCR
void init_discriminant ( int nr, int nc );
//void print_discriminant(int nr, int nc);
//void clear_discriminant(int nr, int nc);
//void analyse_discriminant(int row, int col, lmc_t lmc, int nr, int nc);

#endif

/* sorting templates */

template <class Object>
/**
 * @brief Tempalate for insertionSort
 * @param x[]
 * @param length
 */
inline void insertionSort ( Object x[],int length )
{
    Object key;
    int i;
    for ( int j=1; j<length; j++ ) {
        key=x[j];
        i=j-1;
        while ( x[i]>key && i>=0 ) {
            x[i+1]=x[i];
            i--;
        }
        x[i+1]=key;
    }
}

template <class itemType, class indexType>
/// sort arrays using bubbleSort
inline void bubbleSort ( itemType a[], indexType l, indexType r )
{
    indexType i, j;

    for ( i=r; i>l; --i )
        for ( j=l; j<i; ++j )
            if ( a[j] > a[j+1] )
                std::swap ( a[j], a[j+1] );
}

template <class itemType, class indexType>
/** sorting similar to bubblesort but fast for sorted arrays
 * 
 * The indices l and r are inclusive.
 */
inline void flipSort ( itemType a[], indexType l, indexType r )
{
    indexType i, j;

    i = r;
    while ( i>l ) {
        indexType ii = i;
        i=l;
        for ( j=l; j<ii; ++j ) {
            if ( a[j] > a[j+1] ) {
                std::swap ( a[j], a[j+1] );
                i=j;
            }
        }
    }
}

template <class Object, class indexType>
/**
 * @brief Template for bubble sort
 * @param obj[]
 * @param array_size
 */
inline void bubbleSort2 ( Object obj[], indexType array_size )
{
    Object temp;

    for ( indexType i = ( array_size - 1 ); i >= 0; i-- ) {
        for ( indexType j = 1; j <= i; j++ ) {
            if ( obj[j] < obj[j-1] ) {
                temp = obj[j-1];
                obj[j-1] = obj[j];
                obj[j] = temp;
            }
        }
    }
}

template<class T>
/// sort list using quickSort
void quickSort ( T a[], const int& leftarg, const int& rightarg )
{
    if ( leftarg < rightarg ) {

        T pivotvalue = a[leftarg];
        int left = leftarg - 1;
        int right = rightarg + 1;

        for ( ;; ) {

            while ( a[--right] > pivotvalue ) {
            };
            while ( a[++left] < pivotvalue ) {
            };

            if ( left >= right )
                break;

            T temp = a[right];
            a[right] = a[left];
            a[left] = temp;
        }

        int pivot = right;
        quickSort ( a, leftarg, pivot );
        quickSort ( a, pivot + 1, rightarg );
    }
}

template <class itemType, class indexType>
/*** sort list using shellSort
 * The indices l and r are inclusive.
 */
void shellSort ( itemType a[], indexType l, indexType r )
{
    static indexType i, j, h;
    static itemType v;

    for ( h=1; h<= ( r-l ) /9; h=3*h+1 ) {
    };
    for ( ; h>0; h/=3 ) {
        for ( i=l+h; i<=r; ++i ) {
            for ( j=i-h, v=a[i]; j>=l && a[j]>v; a[j+h]=a[j], j-=h ) {
            };
            a[j+h] = v;
        }
    }
}

/// replace all occurces of a substring in a string
inline std::string replaceString ( std::string subject, const std::string& search,
                                   const std::string& replace )
{
    size_t pos = 0;
    while ( ( pos = subject.find ( search, pos ) ) != std::string::npos ) {
        subject.replace ( pos, search.length(), replace );
        pos += replace.length();
    }
    return subject;
}

/// print a double value as bits
inline void printdoubleasbits ( double decker )
{
    unsigned char * desmond = ( unsigned char * ) & decker;
    for ( size_t i = 0; i < sizeof ( double ); i++ ) {
        myprintf ( "%02X ", desmond[i] );
    }
    myprintf ( "\n" );
}

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
