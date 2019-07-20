/*! \file oaoptions.h
    \brief Contains options for the oa algorithm

        These options are compiled into the program. Making them a run-time option would be inefficient.
*/

#pragma once

#include <string>

/* max number of columns in algorithms */
#define MAXCOLS 60
#define MAXROWS 2048

/* add checks for nonsense input from the user */
#define OACHECK 1

/* add extra debugging checks */
//#define OADEBUG 1

/* check for overflow */
//#define OAOVERFLOW 1

/* use memory efficient versions of the algorithm */
/* skip initialization of the rootrow-level structure */
//#define OAMEM 1

/* safe level permutations (needed for intermediate lmc test) */
//#define SAFELPERM 1

/* frequency element cache using row, value pairs */
#define FREQELEM 1

/* in the extending column increase at most by current max plus 1 */
#define USE_SMALLSTEP

/* use symmetry blocks from first column */
#define SYMMBLOCKS 1

/* use special code for column (t+1) */
#define TPLUSCOLUMN 1

#define oacolSort flipSort
#define oacolSortName "flipSort"
//#define oacolSort bubbleSort
//#define oacolSortName "bubbleSort"

/**
 * Print the compile-time options to string.
 * @return String with compile time information
 */
std::string compile_information ();

/// Print version
std::string version ();

/// Print copyright statement
void print_copyright ();


/// Print copyright statement
void print_copyright_light ();

/// Print compile time options to output stream
void print_options (std::ostream &output_stream);
/// Print compile time options to stdout
void print_options ();

#ifdef OADEV
inline int oadevelop () { return 1; }
#else
inline int oadevelop () { return 0; }
#endif

/// Return true if zlib is available and compressed files can be read
int has_zlib();
