/*! \file tools.h
 *  \brief Contains various functions
 *
 * Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 * Copyright: See LICENSE.txt file that comes with this distribution
 */
#pragma once

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <list>
#include <ostream>
#include <sstream>
#include <stdarg.h>
#include <string>
#include <time.h>


#include "printfheader.h"

/// function to print debugging messages
void printfd_handler (const char *file, const char *func, int line, const char *message, ...);

#define printfd(...) printfd_handler (__FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)

#include "arraytools.h"
#include "mathtools.h"

#if (_MSC_VER >= 100)
#pragma warning(disable : 4018)
#pragma warning(disable : 4996)
#endif

#ifdef SWIG
%ignore nullStream;
%ignore logstream;
%ignore next_comb;
%ignore addelement;
#endif

//! Loglevel definitions. The loglevel determines to amount of output to stdout
enum loglevel_t { LOGERROR, SYSTEM, QUIET, NORMAL, DEBUG, EXTRADEBUG };

int log_print (const int level, const char *message, ...);

/// return current level of logging
int getloglevel ();

/// reset the level of logging
void setloglevel (int n);

/// return True if the current logging level is smaller or equal than the specified level
bool checkloglevel (int level);

/** \brief Null stream */
class nullStream : public std::ostream {
      public:
        nullStream () : std::ostream (NULL) {}
};

#ifdef FULLPACKAGE
/// log a stream to stdout if level specied is higher than the current logging level
std::ostream &logstream (int level);
#endif

/// Return string describing the system
std::string system_uname ();

/// return path separator symbol for the current platform
char path_separator ();

/// handler for error messages. throws an std::runtime_error exception
void mycheck_handler (const char *file, const char *func, int line, int condition, const char *error_message, ...);

#define mycheck(...) mycheck_handler (__FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)

/// Check whether the condition is true and throw an expception otherwise
void myassert (int condition, const char *error_message );

/** Throw a runtime_error exception with specified error message
 * 
 * This exception is caught in the SWIG interface.
 */
void throw_runtime_exception (const std::string exception_message);

/// conditional printf
inline int cprintf (int check, const char *message, ...) {
        int n = 0;
        if (check) {
                va_list va;
                va_start (va, message);
                n = vprintf (message, va);
                va_end (va);
        }
        return n;
}

/// flush to stdout
inline void flush_stdout () {
        fflush (stdout);
}

template < class A >
/**
 * Delete object given by a pointer and set to zero.
 */
inline void safedelete (A *pointer) {
        if (pointer != 0)
                delete pointer;
        pointer = 0;
}

template < class A >
/**
 * Delete array and set pointer to zero
 * @param pointer Pointer to allocated array 
 */
inline void safedeletearray (A *pointer) {
        if (pointer != 0)
                delete[] pointer;
        pointer = 0;
}

template < class Type >
/*!
  Gives next combination for k elements out of n based on an algorithm from wikipedia.
  The generate is sorted.
  \brief Gives next combination
  \param comb Pointer to combination
  \param k Number of the current combination
  \param n Number of elements in combination
  */
int next_comb (std::vector< Type > &comb, int k, int n) {
        int i; // = k - 1;
        const int offset = n - k + 1;
        i = k - 1;
        comb[i]++;
        while ((comb[i] >= offset + i) && (i > 0)) {
                i--;
                comb[i]++;
        }

        if (comb[0] > n - k)
                return 0; /* No more combinations can be generated */

        /* comb now looks like (…, x, n, n, n, …, n).
           Turn it into (…, x, x + 1, x + 2, …) */
        for (i++; i < k; i++)
                comb[i] = comb[i - 1] + 1;

        return 1;
}

/** Go to next combination in sequence */
int next_comb (int *comb, int k, int n);
/** Go to next combination in sequence */
int next_comb_s (int *comb, int k, int n);

template < class Object >
/**
 * @brief Template to swap two objects of arbitrary datatype
 * Please use std::swap instead
 * @param a
 * @param b
 */
void swap_object (Object &a, Object &b) {
        Object tmp;
        tmp = a;
        a = b;
        b = tmp;
}

template < class rtype >
/** Calculate the number of elements in a 2D table with rows with different sizes */
inline int malloc2d_nelements (const int nrows, const rtype *rowsizes) {
        int nelements = 0;
        for (int i = 0; i < nrows; i++)
                nelements += rowsizes[i];
        return nelements;
}

template < class DataType, class rtype >
/**
 * @brief Allocate a 2-dimensional array with non-uniform rows
 * @param nrows Number of rows in the table
 * @param rowsizes Size of each row
 * @return
 */
DataType **malloc2d_irr (const int nrows, const rtype *rowsizes) {
        // Create a 2D array, but with unequal rows (hence irregular -> irr)
        register int i;
        DataType **data;

        int nelements = malloc2d_nelements (nrows, rowsizes);

        data = new DataType *[nrows];
        data[0] = new DataType[nelements];
        memset (data[0], 0, sizeof (DataType) * nelements);

        int offset = 0;
        for (i = 0; i < nrows; i++) {
                data[i] = data[0] + offset;
                offset += rowsizes[i];
        }

        return data;
}

template < class DataType, class rtype >
/**
 * @brief Allocate a 2-dimensional array with non-uniform rows, return size of allocated space
 * @param nrows Number of rows in the table
 * @param rowsizes Size of each row
 * @param nelements This parameter is initialized with the size of the array allocated
 * @return
 */
DataType **malloc2d_irr (const int nrows, const rtype *rowsizes, int &nelements) {
        nelements = malloc2d_nelements (nrows, rowsizes);
        return malloc2d_irr< DataType > (nrows, rowsizes);
}

template < class DataType, class numtype >
/**
 * @brief Allocate a 2-dimensional array of specified size
 * @param nrows Number of rows
 * @param rowsize Size of each row
 * @return
 */
DataType **malloc2d (const numtype nrows, const int rowsize) {
        DataType **data;

        data = new DataType *[nrows];
        if (data == 0) {
                throw_runtime_exception ("malloc2d: error with memory allocation");
        }

        data[0] = new DataType[nrows * rowsize];

        int offset = 0;
        for (int i = 0; i < nrows; i++) {
                data[i] = data[0] + offset;
                offset += rowsize;
        }

        return data;
}

template < class DataType >
/**
 * @brief Release a 2-dimensional array
 * @param data
 * @param nrows
 */
void free2d (DataType **data, const int nrows) {
        delete[] data[0];
        delete[] data;
        data = 0;
}

template < class DataType >
/**
 * @brief Release a 2-dimensional array
 * @param data
 */
void free2d (DataType **data) {
        delete[] data[0];
        delete[] data;
        data = 0;
}

template < class DataType >
/**
 * @brief Release a 2-dimensional non-uniform array
 * @param data Pointer to allocated array
 * @param nrows Not used at the moment
 */
void free2d_irr (DataType **data, const int nrows = -1) {
        free2d (data);
}

void print_array (const char *str, const array_t *array, const int nrows, const int ncols);
void print_array (const array_t *array, const rowindex_t nrows, const colindex_t ncols);

#ifdef FULLPACKAGE
template < class atype >
/// print vector using generic std::cout print functionality
void display_vector (const std::vector< atype > &v) {
        const char *sep = " ";
        std::copy (v.begin (), v.end (), std::ostream_iterator< atype > (std::cout, sep));
}
#else
template < class atype > void display_vector (const std::vector< atype > &v) {
        // dummy implementation
        myprintf ("vector(...)");
}
template < int atype > void display_vector (const std::vector< atype > &v) {
        for (int i = 0; i < v.size (); i++) {
                myprintf ("%d", v[i]);
                if (i < v.size () - 1)
                        myprintf ("%s", sep);
        }
}
#endif

template < class atype >
/** print vector using printf function
 *
 * \param vector Vector to be displayed
 * \param format Format to use in printf
 * \param separator Separator symbol to use
 */
void printf_vector (const std::vector< atype > &vector, const char *format, const char *separator = "") {
        for (unsigned int i = 0; i < vector.size (); i++) {
                myprintf (format, vector[i]);
                if (i < vector.size () - 1)
                        myprintf ("%s", separator);
        }
}

/// return time with milisecond precision
double get_time_ms ();

/// return time difference with milisecond precision
double get_time_ms (double t0);

/// trim a string by removing the specified characters from the left and right
void trim (std::string &str, const std::string &trimChars = "");

/// return the current time as a string
std::string currenttime ();

/// return string describing array
std::string oafilestring (const arraydata_t *arrayclass);

template < class numtype >
/** @brief Convert integer to C++ string
 *
 * @param integer_value Integer
 * @return String representation of the integer
 */
inline std::string itos (numtype integer_value) {
        std::stringstream s;
        s << integer_value;
        return s.str ();
}

/// printf-style function that returns std::string
std::string printfstring (const char *message, ...);

template < class Object >
/**
 * @brief Template for insertionSort
 * @param array Data to be sorted
 * @param length Length of array
 */
inline void insertionSort (Object array[], int length) {
        Object key;
        int i;
        for (int j = 1; j < length; j++) {
                key = array[j];
                i = j - 1;
                while (array[i] > key && i >= 0) {
                        array[i + 1] = array[i];
                        i--;
                }
                array[i + 1] = key;
        }
}

template < class itemType, class indexType >
/// sort arrays using bubbleSort
inline void bubbleSort (itemType array[], indexType left, indexType right) {
        indexType i, j;

        for (i = right; i > left; --i)
                for (j = left; j < i; ++j)
                        if (array[j] > array[j + 1])
                                std::swap (array[j], array[j + 1]);
}

template < class itemType, class indexType >
/** Sorting similar to bubblesort but fast for sorted arrays
 *
 * The indices left and right are inclusive.
 */
inline void flipSort (itemType array[], indexType left, indexType right) {
        indexType i, j;

        i = right;
        while (i > left) {
                indexType ii = i;
                i = left;
                for (j = left; j < ii; ++j) {
                        if (array[j] > array[j + 1]) {
                                std::swap (array[j], array[j + 1]);
                                i = j;
                        }
                }
        }
}

template < class Object, class indexType >
/**
 * @brief Template for bubble sort
 * @param array Array to be sorted
 * @param array_size Size of the array
 */
inline void bubbleSort2 (Object array[], indexType array_size) {
        Object temp;

        for (indexType i = (array_size - 1); i >= 0; i--) {
                for (indexType j = 1; j <= i; j++) {
                        if (array[j] < array[j - 1]) {
                                temp = array[j - 1];
                                array[j - 1] = array[j];
                                array[j] = temp;
                        }
                }
        }
}

template < class T >
/// sort list using quickSort
void quickSort (T array[], const int &leftarg, const int &rightarg) {
        if (leftarg < rightarg) {

                T pivotvalue = array[leftarg];
                int left = leftarg - 1;
                int right = rightarg + 1;

                for (;;) {

                        while (array[--right] > pivotvalue) {
                        };
                        while (array[++left] < pivotvalue) {
                        };

                        if (left >= right)
                                break;

                        T temp = array[right];
                        array[right] = array[left];
                        array[left] = temp;
                }

                int pivot = right;
                quickSort (array, leftarg, pivot);
                quickSort (array, pivot + 1, rightarg);
        }
}

template < class itemType, class indexType >
/*** sort list using shellSort
 * The indices left and right are inclusive.
 */
void shellSort (itemType array[], indexType left, indexType right) {
        static indexType i, j, h;
        static itemType v;

        for (h = 1; h <= (right - left) / 9; h = 3 * h + 1) {
        };
        for (; h > 0; h /= 3) {
                for (i = left + h; i <= right; ++i) {
                        for (j = i - h, v = array[i]; j >= left && array[j] > v; array[j + h] = array[j], j -= h) {
                        };
                        array[j + h] = v;
                }
        }
}

/// replace all occurces of a substring in a string
inline std::string replaceString (std::string subject, const std::string &search, const std::string &replacement) {
        size_t pos = 0;
        while ((pos = subject.find (search, pos)) != std::string::npos) {
                subject.replace (pos, search.length (), replacement);
                pos += replacement.length ();
        }
        return subject;
}

/// print a double value as bits
void printdoubleasbits (double double_value, bool add_newline = true);

/// calculate directory name for job splitted into parts
std::string splitDir (std::vector< int > tag_indices);

/// calculate file name of job splitted into parts
std::string splitFile (std::vector< int > tag_indices);

/// calculate tag for job splitted into parts
std::string splitTag (std::vector< int > tag_indices);


