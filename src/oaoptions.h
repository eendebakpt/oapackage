/*! \file oaoptions.h
	\brief Contains options for the oa algorithm

	These options are compiled into the program. Making them a run-time option would be inefficient.
*/

#ifndef ALGO_OPTIONS_H
#define ALGO_OPTIONS_H

//#include <iostream>
#include <string>

// reduce options of programs to only allow classic algorithms
#ifdef OADEV
#else
//#define CLASSICCODE
#endif


/* max number of columns is algorithms */
#define MAXCOLS 60
#define MAXROWS 2048


/* clean memory management */
#define CLEAN_ARRAY_LINK 1

const int HACK = 0;

/* output analysis data */
//#define OAANALYZE 1
//#define OAANALYZE_DISCR 1

/* for developing purposes*/
//#define OAEXTRA 1

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

/* apply strength 1 check to speed up calculations */
//#define

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

/* use j-value check */
//#define JCHECK 1

//#define NOFINAL 1

/**
 * Print the compile-time options to string.
 * @return String with information
 */
std::string compile_information();

/// Print version
std::string version();

/// Print copyright statement
void print_copyright();

/// Print copyright statement
void print_copyright_old();

/// Print copyright statement
void print_copyright_light();

/// Print compile time options
void print_options(std::ostream &outx);
void print_options();

#ifdef OADEBUG
int globalHackOption(int, int = -1);
#endif

#ifdef OADEV
inline int oadevelop()
{
    return 1;
}
#else
inline int oadevelop()
{
    return 0;
}
#endif

#endif




