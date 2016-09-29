/** \file oastreaming.cpp

 C++ program: oastreaming

 oastreaming can calculate series of non-isomorphic orthogonal arrays

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See COPYING file that comes with this distribution
*/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <list>

#ifdef _WIN32
#else
#include <stdbool.h>
#include <unistd.h>
#endif

#include "anyoption.h"
#include "mathtools.h"
#include "strength.h"
#include "extend.h"
#include "tools.h"
#include "lmc.h"


inline void print_progress ( int csol, int ninput, int nextensions, double Tstart, colindex_t col )
{
	time_t seconds;
	struct tm *tminfo;

	time ( &seconds );
	tminfo = localtime ( &seconds );
	log_print ( QUIET, "Extending column %d of array %i/%i ~ %4.1f%% (total collected %d), time %.2f s, %s", col+1, csol, ( int ) ninput, 100.0 * ( ( float ) csol ) / ( float ) ninput, ( int ) nextensions, ( get_time_ms()-Tstart ), asctime ( tminfo ) );

}


/// parse the command line options for the program
AnyOption * parseOptions ( int argc, char* argv[], algorithm_t &algorithm )
{
	AnyOption *opt =  new AnyOption;

	/* parse command line options */
	opt->setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
	opt->setOption ( "loglevel", 'l' );
	opt->setFlag ( "x", 'x' );   /* special debugging flag */
	opt->setOption ( "restart", 'r' );
	opt->setOption ( "oaconfig", 'c' ); /* file that specifies the design */
	opt->setOption ( "output", 'o' ); /* prefix for output files */
	opt->setFlag ( "generate", 'g' );	/* only generate extensions, do not perform LMC check */
	opt->setFlag ( "coptions" );	/* only generate extensions, do not perform LMC check */
	opt->setOption ( "format", 'f' ); /* format to write to */
	opt->setOption ( "mode", 'm' ); /* algorithm method */
	opt->setOption ( "maxk", 'K' ); /* max number of columns to extend to */
	opt->setOption ( "rowsymmetry", 'R' ); /* max number of columns to extend to */


	opt->setOption ( "initcolprev", 'I' ); /* algorithm method */

	opt->addUsage ( "Orthonal Arrays: extend orthogonal arrays" );
	opt->addUsage ( "Usage: oastreaming [OPTIONS]" );

	opt->addUsage ( "" );
	opt->addUsage ( " -h  --help  			Prints this help " );
	opt->addUsage ( " --coptions			Show compile time options used " );
	opt->addUsage ( " -r [FILE]  --restart [FILE]	Restart with results file " );
	opt->addUsage ( " -l [LEVEL] --loglevel [LEVEL]	Set loglevel to number " );
	opt->addUsage ( " -c [FILE]  --oaconfig [FILE]	Set file with design specification" );
	opt->addUsage ( " -g  --generate [FILE]		Only generate arrays, do not perform LMC check" );
	opt->addUsage ( "  --rowsymmetry [VALUE]		Use row symmetry in generation" );
	opt->addUsage ( " -o [STR]  --output [FILE]	Set prefix for output (default: result) " );
	opt->addUsage ( " -f [FORMAT]			Output format (default: TEXT, or BINARY) " );
	opt->addUsage ( " --initcolprev [INTEGER]	Initialization method of new column (default: 1)" );
	opt->addUsage ( " --maxk [INTEGER] Maximum number of columns to exten to (default: extracted from config file) " );

#ifdef CLASSICCODE
#else
	std::string ss = printfstring ( " -m [MODE]			Algorithm (" ) + algorithm_t_list() + ")\n" ;
	opt->addUsage ( ss.c_str() );
#endif
	//opt->printUsage();
	opt->addUsage ( "" );
	opt->addUsage ( "Example: ./oastreaming -r inputarrays.oa -l 2" );

	opt->processCommandArgs ( argc, argv );

	if ( opt->getValue ( "mode" ) != NULL || opt->getValue ( 'm' ) != NULL ) {
		int vv= atoi ( opt->getValue ( 'm' ) );	//set custom loglevel
		algorithm = ( algorithm_t ) vv;
	} else
		algorithm = MODE_AUTOSELECT;

	return opt;

}

/*!
  For restarting an extension, the reading of all arays and putting them into the memmory is handled by init_restart. It only needs
  the file descriptor of the restart file, the characteristic numbers of the design and a "pointer" to the list, where the arrays
  need to be stored.
  \brief Initialises from a previous solutions file
  \param fname File name of file with solutions
  \param p Characteristic numbers of design
  \param solutions List where solutions are going to be put into
  */
int init_restart ( const char *fname, colindex_t &cols, arraylist_t &solutions )
{
	int narrays = readarrayfile ( fname, &solutions, 1, &cols );
	log_print ( NORMAL, "init_restart: number of arrays: %i\n",narrays );

	return 0;
}



/**
 * @brief Main function for oaextendmpi and oaextendsingle
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] )
{

	double Tstart = get_time_ms();

	algorithm_t algorithm=MODE_INVALID;
	AnyOption *opt = parseOptions ( argc, argv, algorithm );
	print_copyright_old();

	if ( opt->getFlag ( "help" ) || ( opt->getValue ( "coptions" ) != NULL ) ) {
		if ( opt->getFlag ( "help" ) )
			opt->printUsage();
		if ( opt->getValue ( "coptions" ) != NULL )
			print_options ( cout );

		delete opt;
		return 0;
	}


	int loglevel = opt->getIntValue ( 'l', NORMAL );
	setloglevel ( loglevel );

	int userowsymm = opt->getIntValue ( "rowsymmetry", 1 );
	int maxk = opt->getIntValue ( "maxk", 100000 );
	int initcolprev = opt->getIntValue ( "initcolprev", 1 );

	bool restart = false;
	if ( opt->getValue ( "restart" ) != NULL || opt->getValue ( 'r' ) != NULL ) {
		restart = true;
	} else {
		printf ( "oastreaming: please enter an input file with the -r option\n" );
		exit ( 1 );
	}
	std::string restartfile = opt->getValue ( 'r' );

	const char *oaconfigfile = opt->getStringValue ( 'c', "oaconfig.txt" );
	const char *resultprefix = opt->getStringValue ( 'o', "result" );

	arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString ( opt->getStringValue ( 'f', "B" ) );

	myassert ( mode!=ATEXT, "error: cannot use text file in streaming mode" );

	OAextend oaextend;
	oaextend.setAlgorithm ( algorithm );
	oaextend.extendarraymode=OAextend::STOREARRAY;


	//J5_45
	int xx = opt->getFlag ( 'x' );
	if ( xx )
		oaextend.j5structure=J5_45;

	if ( userowsymm==0 ) {
		oaextend.use_row_symmetry=userowsymm;
		printf ( "use row symmetry -> %d\n", oaextend.use_row_symmetry );
	}
	if ( opt->getFlag ( 'g' ) ) {
		std::cout << "only generating arrays (no LMC check)" << endl;
		oaextend.checkarrays=0;
	}


	logstream ( QUIET ) << "#time start: "<< currenttime() << std::endl;

	arraydata_t *ad;
	ad = readConfigFile ( oaconfigfile );
	if ( ad==0 )
		return 1;
	//ad->show(3);

	log_print ( NORMAL, "Using design file: %s (runs %d, strength %d)\n", oaconfigfile, ad->N, ad->strength );


	if ( oaextend.getAlgorithm() == MODE_AUTOSELECT ) {
		//oaextend.algmode = OAextend::getPreferredAlgorithm(*ad);
		oaextend.setAlgorithmAuto ( ad );
	}

	if ( initcolprev==0 ) {
		log_print ( NORMAL, "setting oaextend.init_column_previous to %d\n", INITCOLUMN_ZERO );
		oaextend.init_column_previous = INITCOLUMN_ZERO;
	}


	log_print ( SYSTEM, "using algorithm %d (%s)\n", oaextend.getAlgorithm(), oaextend.getAlgorithmName().c_str() );
	if ( log_print ( NORMAL,"" ) ) {
		std::cout << oaextend.__repr__();
	}


	/*-----------MAIN EXTENSION LOOP-------------*/

	arrayfile_t inputarrays ( restartfile.c_str() );

	// check that oaconfig agrees with loaded arrays
	if ( ad->N != inputarrays.nrows ) {
		printf ( "error: oaconfig does not agree with loaded arrays! (specified N %d)\n", ad->N );
		exit ( 1 );
	}

	colindex_t col_start = inputarrays.ncols;
	int current_col=col_start;
	printf ( "## oastreaming: extension of %d arrays from %d to %d columns\n", inputarrays.narrays, col_start+1, col_start+2 );
	time_t seconds;

	fflush ( stdout );
	arraydata_t *adcol = new arraydata_t ( ad, current_col+1 );

	string fname = resultprefix;
	fname += "-streaming";
	fname += "-" + oafilestring ( adcol );
	logstream ( NORMAL ) << "oaextend: streaming mode: create file " << fname << std::endl;
	int nb = arrayfile_t::arrayNbits ( *ad );
	oaextend.storefile.createfile ( fname, adcol->N, adcol->ncols, -1, mode, nb );


	log_print ( SYSTEM, "Starting with column %d (total time: %.2f [s])\n", current_col + 1, get_time_ms()-Tstart );
	int nr_extensions = 0;

	int csol = 0;
	arraylist_t extensionsdummy;
	for ( int i=0; i<inputarrays.narrays; i++ ) {
		array_link a ( inputarrays.nrows, inputarrays.ncols, i );
		inputarrays.read_array ( a );

		long nextensions = oaextend.storefile.narraycounter;
		print_progress ( csol, inputarrays.narrays, nextensions, Tstart, current_col );
		nr_extensions += extend_array ( a.array, adcol, current_col,extensionsdummy, oaextend );

		csol++;	/* increase current solution */
	}

	long nextensions = oaextend.storefile.narraycounter;
	log_print ( SYSTEM, "Done with column %i, total of %i solutions (time %.2f s))\n", current_col+1, ( int ) nextensions, get_time_ms()-Tstart );

	oaextend.storefile.closefile();

	/* report time */
	time ( &seconds );
	struct tm *tminfo;
	tminfo = localtime ( &seconds );
	log_print ( SYSTEM, "TIME: %.2f s, %s", get_time_ms()-Tstart, asctime ( tminfo ) );
	logstream ( QUIET ) << "#time end: "<< currenttime() << std::endl;
	logstream ( QUIET ) << "#time total: " << printfstring ( "%.1f", get_time_ms()-Tstart ) << " [s]" << std::endl;

	delete adcol;
	delete ad;
	delete opt;

	return 0;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
