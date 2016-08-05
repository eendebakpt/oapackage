#include <sys/timeb.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ostream>
#include <iomanip>
#include <fstream>


#if defined(_MSC_VER)
#pragma warning(disable: 4018)
#pragma warning(disable: 4996)
#endif

#ifdef _WIN32
#else
#include <sys/time.h>
#include <sys/utsname.h>
#endif

using namespace std;

#include "printfheader.h"
#include "tools.h"
#include "mathtools.h"
#include "arraytools.h"


static nullStream staticNullStream;


/** @brief Return string describing the system
 */
string system_uname()
{
#ifdef _WIN32
    stringstream ss;
    ss << "Windows";
    return ss.str();
#else
    stringstream ss ;
    struct utsname  uts;
    uname ( &uts );

    ss << printfstring ( "uts.sysname: %s\n", uts.sysname );
    ss << printfstring ( "uts.nodename: %s\n", uts.nodename );
    ss << printfstring ( "uts.release: %s\n", uts.release );
    ss << printfstring ( "uts.version: %s\n", uts.version );
    ss << printfstring ( "uts.machine: %s\n", uts.machine );
    return ss.str();
#endif
}

/*!
  Shows an array of x by y elements
  \brief Shows array
  \param array Pointer to OA
  \param x Number of columns
  \param y Number of rows
  */
void show_array ( carray_t *array, const int x, const int y )
{
#ifdef FULLPACKAGE
    write_array_format ( stdout, array, y, x );
#else
    myprintf ( "show_array: not implemented...\n" );
#endif
}

/*!
  Shows an array of r times c
  \brief Shows array
  \param array Pointer to OA
  \param r Number of rows
  \param c Number of columns
  */
void print_array ( const array_t *array, const rowindex_t r, const colindex_t c )
{
#ifdef FULLPACKAGE
    write_array_format ( cout, array, r, c );
#else
    myprintf ( "print_array: not implemented...\n" );
#endif
}

void print_array ( const char *str, const array_t *array, const rowindex_t r, const colindex_t c )
{
    myprintf ( "%s", str );
    print_array ( array, r, c );
}

/*! @brief Print array, overloaded function */
void print_array ( const array_link &A )
{
#ifdef FULLPACKAGE
    write_array_format ( stdout, A.array, A.n_rows, A.n_columns );
#else
    myprintf ( "print_array: not implemented...\n" );
#endif
}


/*!
 * Counts the occurence of each element in array
 * @brief Counts elements
 * @param array Pointer to array where elements are counted in
 * @param nelements
 * @param maxval Maximum value that can occur
 * @param elements
 */
void countelements ( carray_t* array, const int nelements, const int maxval, int* elements )
{
    //NOTE: use ideas from http://codereview.stackexchange.com/questions/47566/optimization-for-histogram-computation-algorithm-in-c

    memset ( elements, 0, maxval * sizeof ( int ) );

    for ( int i=0; i<nelements; i++ )
        elements[array[i]]++;
//	  printf("countelements: done\n");
}



/*!
  Gives next combination for k elements out of n based on an algorithm from wikipedia.
  The generate is sorted.
  \brief Gives next combination
  \param comb Pointer to combination
  \param k Number of the current combination
  \param n Number of elements in combination
  */
int next_comb ( int *comb, int k, int n )
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


/// Alternative combination walking algorithm
/// Also safe for combination of 0 elements
int next_comb_s ( int *comb, int k, int n )
{
    if ( k==0 )
        return 0;

    int i;
    for ( i=0; i< ( k-1 ); i++ ) {
        if ( comb[i] < comb[i+1]-1 ) {
            comb[i]++;
            init_perm ( comb, i );
            return 1;
        }
    }
    if ( comb[i]< ( n-1 ) ) {
        comb[i]++;
        init_perm ( comb, i );
        return 1;
    }

    // no more perms, reset the data
    init_perm ( comb, k );

    return 0;

}



/** @brief Construct file string from a design
 *
 * @return std::string
 */
std::string oafilestring ( const arraydata_t *ad )
{
    return oafilestring ( ad->N, ad->ncols, ad->s );
}

/** @brief Construct file string from a design
 *
 * @return String
 */
string oafilestring ( rowindex_t rows, colindex_t cols, array_t *s )
{
    string fname = "";
    fname += itos ( rows );

    for ( int i = 0; i < cols; i++ ) {
        if ( i==0 )
            fname += '.';
        else
            fname += '-';
        fname += itos ( s[i] );
    }
    fname += ".oa";

    return fname;
}

#define XPFS
#ifndef XPFS
/**
 * @brief Function similar to printf returning C++ style string
 * @param message
 * @return
 */
std::string printfstring ( const char *message, ... )
{
    char buf[8*1024];

    va_list va;
    va_start ( va, message );
    vsprintf ( buf, message, va );
    va_end ( va );

    std::string str ( buf );
    return str;
}
#endif

/*!Returns the actual time with millisecond precission, can be used for measuring times of small functions
  \brief Returns the actual time with millisecond precission

  */
double get_time_ms()
{
    struct timeb	tb;
    ftime ( &tb );
    return ( double ) tb.time + ( ( double ) tb.millitm/1000.0f );
}

double get_time_ms ( double t0 )
{
    struct timeb	tb;
    ftime ( &tb );
    return ( double ) tb.time + ( ( double ) tb.millitm/1000.0f ) - t0;
}
const std::string whiteSpaces ( " \f\n\r\t\v" );


void trimRight ( std::string& str,
                 const std::string& trimChars = whiteSpaces )
{
    std::string::size_type pos = str.find_last_not_of ( trimChars );
    str.erase ( pos + 1 );
}


void trimLeft ( std::string& str, const std::string& trimChars = whiteSpaces )
{
    std::string::size_type pos = str.find_first_not_of ( trimChars );
    str.erase ( 0, pos );
}


void trim ( std::string& str, const std::string& trimChars )
{
    if ( trimChars.length() ==0 ) {
        trimRight ( str, whiteSpaces );
        trimLeft ( str, whiteSpaces );
    } else {
        trimRight ( str, trimChars );
        trimLeft ( str, trimChars );
    }
}


#ifdef OAANALYZE_DISCR
/*
 * This block of code contains functions and static variables for the analysis of the algorithm.
 *
 *
 */

#include <map>

#ifdef FULLPACKAGE
static int **discr_data1 = 0;
static int **discr_data2 = 0;
static int snr, snc;

typedef std::map<std::string, long> valueMap;
static valueMap analysis_values;

void analysis_init_values()
{
    analysis_values.clear();
    analysis_values.insert ( valueMap::value_type ( "root_level_perm", 0 ) );
}

void analysis_increase_counter ( const std::string &p )
{
    analysis_values[p]++;
}

void analysis_show_counter ( const std::string &p )
{
    printf ( "analysis: %s: %ld\n", p.c_str(), analysis_values[p] );
}



void init_discriminant ( int nr, int nc )
{
    if ( discr_data1==0 || snc!=nc || snr!=nr ) {
        cout << " init_discriminant" << printfstring ( ": %d %d", nr, nc ) << endl;
        snr = nr;
        snc=nc;
        //logstream(SYSTEM) << "analyse_discriminant: allocating memory " << endl;

        if ( discr_data1!=0 ) {
            free2d ( discr_data1 );
            free2d ( discr_data2 );
        }

        discr_data1 = malloc2d<int> ( nc,nr );
        discr_data2 = malloc2d<int> ( nc,nr );
        for ( int i=0; i<nc*nr; i++ )
            discr_data1[0][i]=0;
        for ( int i=0; i<nc*nr; i++ )
            discr_data2[0][i]=0;
    }
}

void analyse_discriminant ( int row, int col, lmc_t lmc, int nr, int nc )
{
    init_discriminant ( nr,nc );
    switch ( lmc ) {
    case LMC_LESS:
        discr_data1[col][row]++;
        break;
    case LMC_MORE:
        discr_data2[col][row]++;
        break;
    case LMC_EQUAL:
        // no not store data
        break;
    default
            :
        std::cout << "problem" << endl;
        break;
    }
}

void print_discriminant ( int nr, int nc )
{
    if ( discr_data1==0 || snc!=nc || snr!=nr ) {
        printf ( "ERROR: no init done %d %d %d %d!\n", nc, snc, nr, snr );
    }
    //cout << " print_discriminant" << printfstring(": %d %d", nr, nc) << endl;
    cout << "LESS" << endl;
    write_array_format<int> ( cout, discr_data1[0], snr, snc, 9 );
    cout << "MORE" << endl;
    write_array_format<int> ( cout, discr_data2[0], snr, snc, 9 );
    printf ( "MORE TOTALS %d %d\n", nr, nc );
    for ( int i=0; i<nc; i++ ) {
        long count=0;
        for ( int r=0; r<nr; r++ ) {
            count+=discr_data2[i][r];
        }
        printf ( "%9ld ", count );
    }
    printf ( "\n" );
}

void clear_discriminant ( int nr, int nc )
{
    //cout << " clear_discriminant" << printfstring(": %d %d", nr, nc) << endl;
    init_discriminant ( nr,nc );
    for ( int x=0; x<nr*nc; x++ ) {
        discr_data1[0][x]=0;
        discr_data2[0][x]=0;
    }
}

#endif

void analysis_print ( const int level, const char *message, ... )
{
#ifdef OAANALYZE_DISCR
    static int 	loglevel = ANALYSIS_DISCRIMINATOR;
    static FILE *fid = 0;
    va_list		va;

    if ( fid==0 ) {
        char buf[STRBUFSIZE];
        sprintf ( buf, "analysislog.txt" ); // , MPI::COMM_WORLD.Get_rank());
        fid = fopen ( buf, "w" );
    }

    va_start ( va, message );
    char str[STRBUFSIZE];
    vsprintf ( str, message, va );
    if ( loglevel==level ) {
        fprintf ( fid, "%s", str );
    }
    va_end ( va );
    fflush ( fid );
#else
    message=0;
    int x = level+1;
    x++; /* prevent compiler warnings */
#endif
}

#endif // FULLPACKAGE

/* this integer determines to level of logging */
static int streamloglvl = NORMAL;

bool checkloglevel ( int l )
{
    bool b;
    #pragma omp critical
    {
//#pragma omp atomic
        b = l<=streamloglvl;
    }
    return b;
}

int getloglevel()
{
    return streamloglvl;
}

void setloglevel ( int n )
{
    #pragma omp critical
    {
        streamloglvl = n;
    }
    log_print ( -n, "" );	// for log_print
}

#ifdef FULLPACKAGE
/** @brief Returns an output stream object
 *
 * Depending on the log level the stream is standard output or a null stream.
 *
 * @param level
 * @return
 */
ostream& logstream ( int level )
{
    if ( level <= streamloglvl )
        return std::cout;
    else {
        //nullStream ns; return ns;
        return staticNullStream;
    }
}
#endif

/*!
  The function is used to print messages with a different priority. The first time, a threshold is set to
  determine which messages are actually shown and which aren't. The initialisation is done by calling the function
  with the negative value of the threshold and no further arguments. The macro setloglevel can be used with the normal
  level for easier use.
  The function can just be used as the printf function, only the first argument is always the priority level.
  The rest is a standard format string with arguments.
  \brief Print a message with a set priority
  \param level The priority of the message given
  \param message Format string for the message, see printf
  \return
  */
int log_print ( const int level, const char *message, ... )
{
    int result = 0;

#ifdef FULLPACKAGE
    // TODO: convert loglevel to more generic streamloglvl...
#ifdef SWIGPYTHON
    static int 	_loglevel = ( int ) SYSTEM;
#else
    static int 	_loglevel = ( int ) QUIET;
#endif

    if ( level < 0 ) {		//level < 0 means set level, any message appended is printed as well
        #pragma omp critical(logprint)
        {
            va_list		va;
            va_start ( va, message );
            int mlevel=-level;
            _loglevel = mlevel;
            vprintf ( message, va );
            va_end ( va );
        }
    } else {
        if ( level <= _loglevel )	{ //if printing level is high enough, the message is shown
            #pragma omp critical(logprint)
            {
                va_list		va;
                va_start ( va, message );
                vprintf ( message, va );
                va_end ( va );
            }
        }
    }
    result = ( level <= _loglevel );
#endif // FULLPACKAGE

    return result;
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
