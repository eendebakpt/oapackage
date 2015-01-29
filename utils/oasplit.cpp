/** \file oasplit.cpp

 C++ program: oasplit

 oasplit can split a file with arrays into multiple files. This is usefull when the total number of arrays exceeds the memory
 capacity of the system

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "arraytools.h"
#include "anyoption.h"
#include "tools.h"


const int NSPLIT_DEFAULT = 8;	/* default number of files to split into */

/**
 * @brief Read in a file with arrays and split it into several files
 * @param argc 
 * @param argv[] 
 * @return 
 */
int main(int argc, char* argv[])
{
	AnyOption opt;

	/* parse command line options */
	opt.setFlag(  "help", 'h' );   /* a flag (takes no argument), supporting long and short form */ 
	opt.setOption("input", 'i');
	opt.setOption("output", 'o');
	opt.setOption("number", 'n');
	opt.setOption("verbose", 'v');
	opt.setOption("nwritemax", 0);
	opt.setOption("format", 'f');

	opt.addUsage( "oasplit: Split an array file into several files" );
	opt.addUsage( "Usage: oasplit -i [FILE] [OPTIONS]" );
	opt.addUsage( "" );
	opt.addUsage( " -h  --help  		Prints this help " );
	opt.addUsage( " -i [FILE]  --input [FILE]		Input file " );
	opt.addUsage( " -n [NUMBER]  --number [NUMBER]	Number of files to split into (default: 8) " );
	opt.addUsage( " --nwritemax [NUMBER]  					Max number of files to actucally write (default: 1000) " );
	opt.addUsage( " --nb [NUMBER]  					Number of bits for binary file (1 or 8, default: 8) " );
	opt.addUsage( " -o [STR]  --output [FILE]		Output prefix (default: split) " );
	opt.addUsage( " -f [FORMAT]					Output format (default: TEXT, or BINARY; B) " );
	opt.addUsage( " -v [INTEGER]					Verbose (default: 2) " );
	opt.processCommandArgs(argc, argv);

	print_copyright();

	/* parse options */
	if( opt.getFlag( "help" ) || opt.getFlag( 'h' ) ) {
		opt.printUsage();
		exit(0);
	}

	if(opt.getValue("input") == NULL && opt.getValue('i') == 0) {
		printf("No input file specified, use oasplit -h to get help\n");
		exit(1);
	}

	const char *inputfile;
	inputfile = opt.getValue('i');

		const char *outputprefix = opt.getStringValue('o', "split");
	int nsplit = opt.getIntValue('n', NSPLIT_DEFAULT);
	int nmax = opt.getIntValue("nwritemax", 1000);
	int nb = opt.getIntValue("nb", 8);
	int verbose = opt.getIntValue("verbose", 2);

	std::string format = opt.getStringValue('f', "BINARY");
	arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString(format);
	
	if(mode==ABINARY_DIFFZERO) {
	  nb=1;
	}
	/* open the file containint the arrays */
	int rows, cols, narrays;
	arrayfile_t *afile = new arrayfile_t(inputfile);
	rows = afile->nrows;
	cols = afile->ncols;
	narrays = afile->narrays;

	if(! afile->isopen() ) {
		printf("oasplit: problem opening file %s\n", inputfile);
		exit(1);
	}
	const int nwrite = std::min(nmax, nsplit);
	if (verbose)
	  printf("oasplit: input %s, outputprefix %s, nsplit %d, nwrite %d\n", inputfile, outputprefix, nsplit, nwrite);

	/* prepare the output files */
	string *outputfiles = new string[nsplit];
	int *nsp = new int[nsplit];
	arrayfile_t ** outfid = 0;

	if (nmax>1000 ) {
	  	logstream(NORMAL) << printfstring("oasplit: warning nsplit might be larger then max number of file descriptors available...\n");		
	}
	
	
	outfid = new arrayfile_t* [nwrite];
	for(int i = 0; i < nwrite; i++) {
		outputfiles[i] = outputprefix;
		outputfiles[i] += "-split-" + itos(i) + ".oa";

		if (verbose>=2) {
		  logstream(NORMAL)  << "oasplit: outputfile " << i << ": " << outputfiles[i] << endl;
		}
		
		int nars = floor( (double)narrays/nsplit);
		if ((narrays%nsplit)>i)
			nars++;
		logstream(DEBUG) << printfstring("  creating arrayfile %d (%d arrays)\n", i, nars);
		//int nb = 8; //arrayfile_t::arrayNbits(arraylist->at(0));
		outfid[i] = create_arrayfile(outputfiles[i].c_str(), rows, cols, nars, mode, nb);
		if(! outfid[i]->isopen() ) {
		printf("oasplit: error opening %s, aborting program\n", outputfiles[i].c_str());
		  exit(1); 
		}
		
		nsp[i] = 0;
	}

	array_link al(rows, cols, 0);
	int target;
	for(int i = 0; i < narrays; i++)
	{
		if(i%5000==0 && verbose) 
			printf("  scanning array %d/%d\n", i, narrays);

		afile->read_array(al);
		
		target = (i%nsplit);
		if(target>=nmax)
		  continue;
		nsp[target]++;
		outfid[target]->append_array(al);
	}

	for(int i = 0; i < nwrite; i++)
	{
		//finish_arrayfile(outfid[i]);
		outfid[i]->finisharrayfile();
		delete outfid[i];
	}

	/* clean up memory */
	delete [] outfid;
	delete [] outputfiles;
	delete [] nsp;
	afile->closefile();
	delete afile;

	return 0;
}
