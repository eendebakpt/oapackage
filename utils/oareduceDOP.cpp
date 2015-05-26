/** \file oareduceDOP.cpp

 C++ program: oareduceDOP

 oareduceDOP: reduce arrays to delete-one-factor canonical form

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"


/**
 * @brief Read in files with arrays and join them into a single file
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] )
{
   AnyOption opt;

   /* parse command line options */
   opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
   opt.setOption ( "output", 'o' );
   opt.setFlag ( "extend", 'e' );
   opt.setFlag ( "prune", 'p' );
   opt.setOption ( "verbose", 'v' );
   opt.setOption ( "dolmc", 'd' );
   opt.setOption ( "format", 'f' );
   opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */
   opt.setOption ( "strength", 's' );

   opt.addUsage ( "Orthonal Array: oareduceDOP: testing platform" );
   opt.addUsage ( "Usage: oareduceDOP [OPTIONS] [FILES]" );
   opt.addUsage ( "" );
   opt.addUsage ( " -h --help  			Prints this help " );
   opt.addUsage ( " -f [FORMAT]					Output format (default: TEXT, or BINARY) " );
   opt.addUsage ( " -o [FILE] --output [FILE]	Output prefix (default: standard output) " );
   opt.addUsage ( " -s [STRENGTH]  --strength [STRENGTH] Set strength to use in checking" );
   opt.addUsage ( " -e --extend			Extend set of arrays before checking " );
   opt.addUsage ( " -p --prune			Prune" );
   opt.addUsage ( " -d --dolmc			Do LMC?" );
   opt.addUsage ( " " );
   opt.addUsage ( "Reduce a set of arrays to delete-one-factor projection form." );

   opt.processCommandArgs ( argc, argv );

   print_copyright();
   setloglevel ( NORMAL );

   /* parse options */
   if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() ==0 ) {
      opt.printUsage();
      exit ( 0 );
   }

   logstream ( QUIET ) << "#time start: "<< currenttime() << std::endl;

   const char *oaconfigfile = "oaconfig.txt";
   const char *resultprefix;

   if ( opt.getValue ( "oaconfig" ) !=NULL )
      oaconfigfile = opt.getValue ( 'c' );

   arrayfile::arrayfilemode_t mode = arrayfile::ATEXT;


   if ( opt.getValue ( "format" ) !=0 ) {
      std::string format = opt.getValue ( 'f' );
      mode = arrayfile_t::parseModeString ( format );
   }

   const char *outputprefix = "test.oa";
   if ( opt.getValue ( "output" ) !=0 )
      outputprefix = opt.getValue ( 'o' );

   int verbose = opt.getIntValue ( 'v', 2 );
   //double Afinal = opt.getDoubleValue('A', 0);
   int kfinal = opt.getIntValue ( 'k', 7 );
   int strength = opt.getIntValue ( 's', 2 );
   int dolmc = opt.getIntValue ( 'd', 1 );

   int doextend = opt.getFlag ( 'e' );
   int dopruning = opt.getFlag ( 'p' );


   /* read in the arrays */
   arraylist_t *arraylist = new arraylist_t;
   std::cout << "oareduceDOP: reading " << opt.getArgc() << " file(s)" << endl;
   if ( verbose>=3 ) {
      printf ( "oareduceDOP: strength %d\n", strength );
   }
   for ( int i = 0 ; i < opt.getArgc() ; i++ ) {
      int n = readarrayfile ( opt.getArgv ( i ), arraylist );
      cout << "file " <<  opt.getArgv ( i ) << ":   read " << n << " array(s)" << endl;
   }

   arraylist_t earrays;

   // extend arrays

   arraydata_t *ad;
   ad = readConfigFile ( oaconfigfile );

   arraylist_t extensions;
   arraylist_t *arraylist2 = arraylist;

   if ( doextend && arraylist->size() >0 ) {
      array_link al =arraylist->at ( 0 );

      OAextend oaextend;
      oaextend.setAlgorithm ( MODE_ORIGINAL );
      int ecol=al.n_columns;
      oaextend.init_column_previous=0;
      oaextend.checkarrays=0;

      extend_arraylist ( *arraylist, *ad, oaextend, ecol, extensions );
      if ( verbose )
         printf ( "  extended %d arrays to %d arrays\n", ( int ) arraylist->size(), ( int ) extensions.size() );
      arraylist2 = &extensions;

   }

   reduceArraysGWLP ( arraylist2, earrays, verbose, dopruning, strength, dolmc );

   // write arrays to disk
   std::string outfile = outputprefix;
   printf ( "  writing %ld arrays to %s\n", earrays.size(), outfile.c_str() );
   writearrayfile ( outfile.c_str(), &earrays, arrayfile::ABINARY );

   delete ad;
   delete arraylist;

   logstream ( QUIET ) << "#time end: "<< currenttime() << std::endl;

   return 0;
}
