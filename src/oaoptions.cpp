#include <Eigen/Core>
#include <Eigen/SVD>

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include "oaoptions.h"
#include "arraytools.h"

/// Information about compile time options
std::string compile_information()
{
	std::stringstream ss;
	print_options(ss);

	return ss.str();
}


/** Return version of program */
std::string version()
{
    std::string v = "2.0.0";
    return v;
}

/** @brief Print copyright notice
 */
void print_copyright_old()
{
	printf("Orthogonal Arrays %s: Copyright TNO Science & Industry (2010), Copyright Pieter Eendebak (2011-2015)\n", version().c_str());
	printf("For more details see the files README.txt and LICENSE.txt\n");
}

/** @brief Print copyright notice
 */
void print_copyright()
{
	printf("Orthogonal Arrays %s\n", version().c_str());
	printf("For more details see the files README.txt and LICENSE.txt\n");
}
/** @brief Print brief copyright notice
 */
void print_copyright_light()
{
	printf("Orthogonal Array package %s\n", version().c_str());
}

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>

// http://sourceforge.net/p/predef/wiki/Compilers/
#if defined(__GNUC__)
# if defined(__GNUC_PATCHLEVEL__)
#  define __GNUC_VERSION__ (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100 \
                            + __GNUC_PATCHLEVEL__)
# else
#  define __GNUC_VERSION__ (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100)
# endif
#else
#define __GNUC_VERSION__ "none"
#endif

#ifndef __INTEL_COMPILER
#define __INTEL_COMPILER "none"
#endif

void print_options()
{
 print_options(std::cout); 
}

/**
 * Print the compile-time options to output stream.
 * @param out 
 */
void print_options(std::ostream &outx)
{
  
  	outx << "Orthogonal Array Package " << version() << std::endl;
	
	outx << "Compile date: " << __DATE__ << " " << __TIME__ << std::endl;
	//outx << "SVN version: " << svn_version << std::endl;

	outx << "array_t type: sizeof(array_t) " << sizeof(array_t) << std::endl;
	//outx << "SVN version: " << svn_version << std::endl;

	outx << "Eigen version: " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;
	
	    Eigen::MatrixXd mymatrix(1,1);
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(mymatrix);
	
	outx << "eigen: JacobiSVD threshold " << lu_decomp.threshold () << std::endl;
	
	// http://sourceforge.net/p/predef/wiki/Compilers/
#ifdef WIN32
#else
	outx << "Compiler: __VERSION__ " << __VERSION__ << std::endl;
#endif
	outx << "Compiler: __GNUC_VERSION__ " << __GNUC_VERSION__ << std::endl;
	outx << "Compiler: __INTEL_COMPILER " << __INTEL_COMPILER << std::endl;
	
	
#ifdef __OPTIMIZE__
	outx << "Optimization: __OPTIMIZE__" << std::endl;	
#endif
	outx << "Compile time options: "; // << std::endl;
	std::string sep = ", ";

	#ifdef OAANALYZE
	outx << "OAANALYZE" << sep;
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
	
//	if (INIT_COL_PREVIOUS==1)
	//	out << "INIT_COL_PREVIOUS" << sep;
	
	#ifdef SYMMBLOCKS
	outx << "SYMMBLOCKS" << sep;
	#endif

	#ifdef JCHECK
	outx << "JCHECK" << sep;
	#endif
	#ifdef OADEV
	outx << "OADEV" << sep;
	#endif
	outx << std::endl;

	outx << "columns sorting method: " << oacolSortName << std::endl;

	
}


#ifdef OADEBUG
int hopts[10]={0,0,0,0,0,0,0,0,0,0};
int globalHackOption(int i, int val) {
  if (i>=0 && i<10) {
    if (val!=-1)
      hopts[i]=val;
    return hopts[i];
  }
  printf("globalHackOption: error: invalid index\n");
  return -1;
}
#endif
