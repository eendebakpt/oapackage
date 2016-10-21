/** \file oaextend.cpp

 C++ program: oaextendmpi / oaextendsingle

 oaextendmpi and oaextendsingle can calculate series of non-isomorphic orthogonal arrays

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

 Copyright: See COPYING file that comes with this distribution
*/

#ifdef OAEXTEND_MULTICORE
// this must be included before stdio.h
#include <mpi.h>
#endif

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

#ifdef OAEXTEND_MULTICORE

#include "mpi.h"
#include "mpitools.h"

#else
/** number of master processor */
const int MASTER = 0;

/* OAEXTEND_SINGLECORE */
int extend_slave_code ( const int this_rank, OAextend const &oaextend ) {
    return 0;
}
#endif

inline void print_progress ( int csol, arraylist_t &solutions, arraylist_t &extensions, double Tstart, colindex_t col ) {
    time_t seconds;
    struct tm *tminfo;

    time ( &seconds );
    tminfo = localtime ( &seconds );
    log_print ( QUIET, "Extending column %d of array %i/%i ~ %4.1f%% (total collected %d), time %.2f s, %s", col+1, csol, ( int ) solutions.size(), 100.0 * ( ( float ) csol ) / ( float ) solutions.size(), ( int ) extensions.size(), ( get_time_ms()-Tstart ), asctime ( tminfo ) );

}


/// parse the command line options for the program
AnyOption * parseOptions ( int argc, char* argv[], algorithm_t &algorithm ) {
    AnyOption *opt =  new AnyOption;

    /* parse command line options */
    opt->setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt->setOption ( "loglevel", 'l' );
    opt->setOption ( "sort", 's' );
    opt->setFlag ( "x", 'x' );   /* special debugging flag */
    opt->setOption ( "restart", 'r' );
    opt->setOption ( "oaconfig", 'c' ); /* file that specifies the design */
    opt->setOption ( "output", 'o' ); /* prefix for output files */
    opt->setFlag ( "generate", 'g' );	/* only generate extensions, do not perform LMC check */
    opt->setFlag ( "coptions" );	/* print compile time options */
    opt->setOption ( "format", 'f' ); /* format to write to */
    opt->setOption ( "mode", 'm' ); /* algorithm method */
    opt->setOption ( "maxk", 'K' ); /* max number of columns to extend to */
    opt->setOption ( "rowsymmetry", 'R' ); /* max number of columns to extend to */

    opt->setFlag ( "streaming" ); /* operate in streaming mode */

    opt->setOption ( "initcolprev", 'I' ); /* algorithm method */

    opt->addUsage ( "Orthonal Arrays: extend orthogonal arrays" );
#ifdef OAEXTEND_SINGLECORE
    opt->addUsage ( "Usage: oaextendmpi [OPTIONS]" );
#else
    opt->addUsage ( "Usage: oaextendsingle [OPTIONS]" );
#endif

    opt->addUsage ( "" );
    opt->addUsage ( " -h  --help  			Prints this help " );
    opt->addUsage ( " --coptions			Show compile time options used " );
    opt->addUsage ( " -r [FILE]  --restart [FILE]	Restart with results file " );
    opt->addUsage ( " -l [LEVEL] --loglevel [LEVEL]	Set loglevel to number " );
    opt->addUsage ( " -s [INTEGER] --sort [INTEGER]	Sort input and output arrays (default: 1) " );
    opt->addUsage ( " -c [FILE]  --oaconfig [FILE]	Set file with design specification" );
    opt->addUsage ( " -g  --generate [FILE]		Only generate arrays, do not perform LMC check" );
    opt->addUsage ( "  --rowsymmetry [VALUE]		Use row symmetry in generation" );
    opt->addUsage ( " -o [STR]  --output [FILE]	Set prefix for output (default: result) " );
    opt->addUsage ( " -f [FORMAT]			Output format (default: TEXT, or BINARY, D, Z). Format D is binary difference, format Z is binary with zero-difference " );
    opt->addUsage ( " --initcolprev [INTEGER]	Initialization method of new column (default: 1)" );
    opt->addUsage ( " --maxk [INTEGER] Maximum number of columns to exten to (default: extracted from config file) " );
    opt->addUsage ( " --streaming			Operate in streaming mode. Generated arrays will be written to disk immediately. " );

#ifdef CLASSICCODE
    algorithm_t_list();
#else
    std::string ss = printfstring ( " -m [MODE]			Algorithm (" ) + algorithm_t_list() + ")\n" ;
    opt->addUsage ( ss.c_str() );
#endif
    //opt->printUsage();
    opt->addUsage ( "" );
    opt->addUsage ( "Example: ./oaextendsingle -l 2" );

    opt->processCommandArgs ( argc, argv );

    if ( opt->getValue ( "mode" ) != NULL || opt->getValue ( 'm' ) != NULL ) {
        int vv= atoi ( opt->getValue ( 'm' ) );	//set custom loglevel
        algorithm = ( algorithm_t ) vv;
    } else {
        algorithm = MODE_AUTOSELECT;
    }
//  printf("Setting algorithm to %d (val %s -> %d )\n", algorithm, opt->getValue("mode"), algorithm);

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
int init_restart ( const char *fname, colindex_t &cols, arraylist_t &solutions ) {
//	array_link        *sols;
//	int index;

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
int main ( int argc, char* argv[] ) {
#ifdef OAEXTEND_SINGLECORE
    const int n_processors = 1;
    const int this_rank = 0;
#else
    MPI::Init ( argc, argv );

    const int n_processors = MPI::COMM_WORLD.Get_size();
    const int this_rank = MPI::COMM_WORLD.Get_rank();
#endif

    //printf("MPI: this_rank %d\n", this_rank);

    /*----------SET STARTING ARRAY(S)-----------*/
    if ( this_rank != MASTER ) {
#ifdef OAEXTEND_MULTICORE
        slave_print ( QUIET, "M: running core %d/%d\n", this_rank, n_processors );
#endif
        algorithm_t algorithm = MODE_INVALID;
        AnyOption *opt = parseOptions ( argc, argv, algorithm );
        OAextend oaextend;
        oaextend.setAlgorithm ( algorithm );

#ifdef OAEXTEND_MULTICORE
        log_print ( NORMAL, "slave: receiving algorithm int\n" );
        algorithm = ( algorithm_t ) receive_int_slave();
        oaextend.setAlgorithm ( algorithm );
        //printf("slave %d: receive algorithm %d\n", this_rank, algorithm);
#endif

        //assert(algorithm==MODE_ORIGINAL);
        //oaextend.algmode = algorithm;
        extend_slave_code ( this_rank, oaextend );

        delete opt;
    } else {
        double Tstart = get_time_ms();
        int nr_extensions;
        arraylist_t		solutions, extensions;

        algorithm_t algorithm=MODE_INVALID;
        AnyOption *opt = parseOptions ( argc, argv, algorithm );
        print_copyright_old();

        int loglevel = opt->getIntValue ( 'l', NORMAL );
        setloglevel ( loglevel );
        //printf("set log level to %d\n", loglevel);

        int dosort = opt->getIntValue ( 's', 1 );
        int userowsymm = opt->getIntValue ( "rowsymmetry", 1 );
        int maxk = opt->getIntValue ( "maxk", 100000 );
        int initcolprev = opt->getIntValue ( "initcolprev", 1 );

        const bool streaming = opt->getFlag ( "streaming" );

        bool restart = false;
        if ( opt->getValue ( "restart" ) != NULL || opt->getValue ( 'r' ) != NULL ) {
            restart = true;
        }

        const char *oaconfigfile = opt->getStringValue ( 'c', "oaconfig.txt" );
        const char *resultprefix = opt->getStringValue ( 'o', "result" );

        arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString ( opt->getStringValue ( 'f', "T" ) );

        OAextend oaextend;
        oaextend.setAlgorithm ( algorithm );

        if ( streaming ) {
            logstream ( SYSTEM ) << "operating in streaming mode, sorting of arrays will not work " << std::endl;
            oaextend.extendarraymode=OAextend::STOREARRAY;
            //oaextend.storefile.createfile();
        }
        //oaextend.init_column_previous=1;

        //J5_45
        int xx = opt->getFlag ( 'x' );
        if ( xx ) {
            oaextend.j5structure=J5_45;
        }

        if ( userowsymm==0 ) {
            oaextend.use_row_symmetry=userowsymm;
            printf ( "use row symmetry -> %d\n", oaextend.use_row_symmetry );
        }
        if ( opt->getFlag ( 'g' ) ) {
            std::cout << "only generating arrays (no LMC check)" << endl;
            oaextend.checkarrays=0;
        }

        if ( opt->getFlag ( "help" ) || ( opt->getValue ( "coptions" ) != NULL ) ) {
            if ( opt->getFlag ( "help" ) ) {
                opt->printUsage();
            }
            if ( opt->getValue ( "coptions" ) != NULL ) {
                print_options ( cout );
            }
        } else {
            logstream ( QUIET ) << "#time start: "<< currenttime() << std::endl;

            arraydata_t *ad;
            ad = readConfigFile ( oaconfigfile );
            ad->lmc_overflow_check();

            if ( ad==0 ) {
                return 1;
            }

            log_print ( NORMAL, "Using design file: %s (runs %d, strength %d)\n", oaconfigfile, ad->N, ad->strength );

            if ( oaextend.getAlgorithm() == MODE_AUTOSELECT ) {
                oaextend.setAlgorithmAuto ( ad );
            }

#ifdef OAEXTEND_SINGLECORE
            if ( initcolprev==0 ) {
                log_print ( NORMAL, "setting oaextend.init_column_previous to %d\n", INITCOLUMN_ZERO );
                oaextend.init_column_previous = INITCOLUMN_ZERO;
            }
#endif

#ifdef OAEXTEND_MULTICORE
            int alg=oaextend.getAlgorithm();
            for ( int i=0; i<n_processors; i++ ) {
                if ( i==MASTER ) {
                    continue;
                }
                log_print ( NORMAL, "MASTER: sending algorithm %d to slave %d\n", alg, i );
                send_int ( alg, i );
            }
#endif
            colindex_t col_start;

            log_print ( SYSTEM, "using algorithm %d (%s)\n", oaextend.getAlgorithm(), oaextend.getAlgorithmName().c_str() );
            if ( log_print ( NORMAL,"" ) ) {
                cout << oaextend.__repr__();
            }

            if ( restart ) { //start from result file
                int initres =init_restart ( opt->getValue ( 'r' ), col_start, solutions );
                if ( initres == 1 ) {	//check if restarting went OK
                    logstream ( SYSTEM ) << "Problem with restart from " << opt->getValue ( 'r' ) << "" << endl;

#ifdef OAEXTEND_MPI
                    MPI_Abort ( MPI_COMM_WORLD, 1 );
#endif
                    exit ( 1 );
                }

                if ( dosort ) {
                    double Ttmp = get_time_ms();
                    sort ( solutions.begin(), solutions.end() ); //solutions.sort();

                    log_print ( QUIET, "   sorting of initial solutions: %.3f [s]\n", get_time_ms()-Ttmp );
                }

                // check that oaconfig agrees with loaded arrays
                if ( solutions.size() >0 ) {
                    if ( ad->N != solutions[0].n_rows ) {
                        printf ( "Problem: oaconfig does not agree with loaded arrays!\n" );
                        // free_sols(solutions); ??
                        solutions.clear();
                    }
                }

            } else {
                //starting with root
                if ( check_divisibility ( ad ) == false )	{	//Divisibility test
                    log_print ( SYSTEM, "ERROR: Failed divisibility test!\n" );

#ifdef OAEXTEND_MPI
                    MPI_Abort ( MPI_COMM_WORLD, 1 );
#endif
                    exit ( 1 );
                }
                create_root ( ad, solutions );
                col_start = ad->strength;
            }

            /*-----------MAIN EXTENSION LOOP-------------*/

            log_print ( SYSTEM, "M: running with %d procs\n", n_processors );

            maxk=std::min ( maxk, ad->ncols );

            //oaextend.info();  ad->show();

            time_t seconds;
            for ( colindex_t current_col = col_start; current_col < maxk; current_col++ ) {
                fflush ( stdout );
                arraydata_t *adcol = new arraydata_t ( ad, current_col+1 );

                if ( streaming ) {

                    string fname = resultprefix;
                    fname += "-streaming";
                    fname += "-" + oafilestring ( ad );
                    logstream ( NORMAL ) << "oaextend: streaming mode: create file " << fname << std::endl;
                    int nb = arrayfile_t::arrayNbits ( *ad );

                    //oaextend.storefile.setVerbose(2);
                    oaextend.storefile.createfile ( fname, adcol->N, ad->ncols, -1, ABINARY, nb );
                }


                //printf("hack: adcol current_col %d\n", current_col); adcol->show();

                //time(&seconds); tminfo = localtime(&seconds);
                //log_print(SYSTEM, "TIME: %5i\t%s\n", seconds - t0, asctime(tminfo));
                log_print ( SYSTEM, "Starting with column %d (%d, total time: %.2f [s])\n", current_col + 1, ( int ) solutions.size(),  get_time_ms()-Tstart );
                nr_extensions = 0;
                arraylist_t::iterator	cur_extension;

                int csol = 0;
                for ( cur_extension = solutions.begin(); cur_extension != solutions.end(); cur_extension++ ) {
                    print_progress ( csol, solutions, extensions, Tstart, current_col );
                    logstream ( NORMAL ) << cur_extension[0];

                    //if (cur_extension-solutions.begin()>2) break;

                    if ( n_processors==1 ) {
                        nr_extensions += extend_array ( cur_extension->array, adcol, current_col,extensions, oaextend );

                    } else {
#ifdef OAEXTEND_MPI
                        double Ttmp = get_time_ms();
                        int slave = collect_solutions_single ( extensions, adcol );
                        log_print ( DEBUG, "   time: %.2f, collect time %.2f\n", ( get_time_ms()-Tstart ), ( get_time_ms()-Ttmp ) );

                        /* OPTIMIZE: send multiple arrays to slave 1 once step */
                        extend_array_mpi ( cur_extension, 1, adcol, current_col, slave );
#endif
                    }

#ifdef OAEXTEND_MPI
                    if ( ( csol+1 ) % nsolsync ==0 ) {
                        collect_solutions_all ( extensions, adcol );
                    }
#endif
                    // OPTIMIZE: periodically write part of the solutions to disk
                    csol++;	/* increase current solution */
                }
#ifdef OAEXTEND_MPI
                collect_solutions_all ( extensions, adcol );
#endif

                if ( checkloglevel ( NORMAL ) ) {
                    csol=0;
                    for ( cur_extension = extensions.begin(); cur_extension != extensions.end(); cur_extension++ ) {
                        log_print ( DEBUG, "%i: -----\n", ++csol );
                        logstream ( DEBUG ) << cur_extension[0];
                    }
                }
                if ( dosort ) {
                    log_print ( DEBUG, "Sorting solutions, necessary for multi-processor code\n" );
                    double Ttmp = get_time_ms();
                    sort ( extensions.begin(), extensions.end() );
                    log_print ( DEBUG, "   sorting time: %.3f [s]\n", get_time_ms()-Ttmp );
                }

                save_arrays ( extensions, adcol, extensions.size(), n_processors, resultprefix, mode );
                solutions.swap ( extensions );	//swap old and newly found solutions
                extensions.clear();		//clear old to free up space for new extensions

                log_print ( SYSTEM, "Done with column %i, total of %i solutions (time %.2f s))\n", current_col+1, ( int ) solutions.size(), get_time_ms()-Tstart );

                delete adcol;

#ifdef OAANALYZE_DISCR
                logstream ( NORMAL ) << "Discriminant: " << endl;
                print_discriminant ( ad->N,current_col+1 );
#endif

                //log_print(NORMAL, "final: ");
                //print_fracs(NORMAL);
                clear_fracs();
            }
            /*------------------------*/
            time ( &seconds );
            struct tm *tminfo;
            tminfo = localtime ( &seconds );
            log_print ( SYSTEM, "TIME: %.2f s, %s", get_time_ms()-Tstart, asctime ( tminfo ) );
            logstream ( QUIET ) << "#time end: "<< currenttime() << std::endl;
            logstream ( QUIET ) << "#time total: " << printfstring ( "%.1f", get_time_ms()-Tstart ) << " [s]" << std::endl;

            solutions.clear();
            delete ad;

#ifdef OAEXTEND_MPI
            stop_slaves();
#endif
        }

        delete opt;
    } /* end of master code */

#ifdef OAEXTEND_MPI
    MPI_Finalize();
#endif

    //lmc_stats();

    return 0;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
