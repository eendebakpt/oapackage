#include <Eigen/Core>
#include <Eigen/SVD>

#include <iostream>
#include <sstream>
#include <string>

#include "arraytools.h"
#include "oaoptions.h"
#include "printfheader.h"
#include "tools.h"
#include "version.h"

namespace nauty {
#include "nauty.h"
}

/// Information about compile time options
std::string compile_information () {
        std::stringstream ss;
        print_options (ss);

        return ss.str ();
}

/** Return version of program */
std::string version () {
        std::string v = __version__;
        return v;
}

/** @brief Print copyright notice
 */
void print_copyright () {
        myprintf ("Orthogonal Arrays %s\n", version ().c_str ());
        myprintf ("For more details see the files README.txt and LICENSE.txt\n");
}
/** @brief Print brief copyright notice
 */
void print_copyright_light () { myprintf ("Orthogonal Array package %s\n", version ().c_str ()); }

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>

// http://sourceforge.net/p/predef/wiki/Compilers/
#if defined(__GNUC__)
#if defined(__GNUC_PATCHLEVEL__)
#define __GNUC_VERSION__ (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
#define __GNUC_VERSION__ (__GNUC__ * 10000 + __GNUC_MINOR__ * 100)
#endif
#else
#define __GNUC_VERSION__ "none"
#endif

#ifndef __INTEL_COMPILER
#define __INTEL_COMPILER "none"
#endif

/**
 * Print the compile-time options to a string.
 */
std::string print_options_string () {
        std::string tabsep = "  ";

        std::stringstream outx;

        outx << "Orthogonal Array Package " << version () << std::endl;

        outx << "Compile date: " << __DATE__ << " " << __TIME__ << std::endl;

        outx << tabsep << "void * type: sizeof(void *) " << sizeof (void *) << std::endl;
        outx << tabsep << "array_t type: sizeof(array_t) " << sizeof (array_t) << std::endl;
        outx << tabsep << "integer types: sizeof(short int) " << sizeof (short int) << std::endl;
        outx << tabsep << "integer types: sizeof(unsigned long int) " << sizeof (unsigned long int) << ","
             << " sizeof(int) " << sizeof (int) << std::endl;
        outx << tabsep << "integer types: sizeof(long long) " << sizeof (long long) << std::endl;
        outx << tabsep << "floating point type: sizeof(float) " << sizeof (float) << ", sizeof(double) "
             << sizeof (double) << ", sizeof(long double) " << sizeof (long double) << std::endl;
        outx << tabsep << "Eigen version: " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "."
             << EIGEN_MINOR_VERSION << std::endl;

        Eigen::MatrixXd mymatrix (1, 1);
        Eigen::FullPivLU< Eigen::MatrixXd > lu_decomp (mymatrix);

        outx << tabsep << "eigen: JacobiSVD threshold " << lu_decomp.threshold () << std::endl;

// http://sourceforge.net/p/predef/wiki/Compilers/
#ifdef WIN32
#else
        outx << tabsep << "Compiler: __VERSION__ " << __VERSION__ << std::endl;
#endif
        outx << tabsep << "Compiler: __GNUC_VERSION__ " << __GNUC_VERSION__ << std::endl;
        outx << tabsep << "Compiler: __INTEL_COMPILER " << __INTEL_COMPILER << std::endl;

#ifdef __OPTIMIZE__
        outx << tabsep << "Optimization: __OPTIMIZE__" << std::endl;
#endif
        outx << tabsep << "Compile time options: "; // << std::endl;
        std::string sep = ", ";

#ifdef USEZLIB
        outx << "USEZLIB" << sep;
#endif

#ifdef OACHECK
        outx << "OACHECK" << sep;
#endif

#ifdef OADEBUG
        outx << "OADEBUG" << sep;
#endif

#ifdef OAOVERFLOW
        outx << "OAOVERFLOW" << sep;
#endif

#ifdef OAMEM
        outx << "OAMEM" << sep;
#endif

#ifdef SAFELPERM
        outx << "SAFELPERM" << sep;
#endif

#ifdef FREQELEM
        outx << "FREQELEM" << sep;
#endif

#ifdef USE_ROW_SYMMETRY
        outx << "USE_ROW_SYMMETRY" << sep;
#endif

#ifdef USE_SMALLSTEP
        outx << "USE_SMALLSTEP" << sep;
#endif

#ifdef TPLUSCOLUMN
        outx << "USE_TPLUSCOLUMN" << sep;
#endif

#ifdef SYMMBLOCKS
        outx << "SYMMBLOCKS" << sep;
#endif

#ifdef JCHECK
        outx << "JCHECK" << sep;
#endif
#ifdef OADEV
        outx << "OADEV" << sep;
#endif
#ifdef SWIGCODE
        outx << "SWIG" << sep;
#endif
        outx << std::endl;

        outx << tabsep << "columns sorting method: " << oacolSortName << std::endl;

        outx << tabsep << "nauty: NAUTYVERSION " << NAUTYVERSION << std::endl;
        outx << tabsep << "nauty: SIZEOF_LONG " << SIZEOF_LONG << std::endl;

        if (SIZEOF_LONG != sizeof (long)) {
                outx << tabsep << "!! ERROR: sizeof(long) does not correspond to compile time size of long"
                     << std::endl;
        }
        const std::string s = outx.str ();
        return s;
}

void print_options () {
#ifdef RPACKAGE
#else
        std::string s = print_options_string ();
        myprintf ("%s", s.c_str ());
#endif
}

int has_zlib(){
#ifdef USEZLIB
	return true;
#else
	return false;
#endif
}

/**
 * Print the compile-time options to output stream.
 * @param out
 */
void print_options (std::ostream &out) {
        std::string s = print_options_string ();
        out << s;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
