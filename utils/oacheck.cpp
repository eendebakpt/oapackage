
/** \file oacheck.cpp

 C++ program: oacheck

 oacheck allows to check wether a list of arrays in LMC form or not.

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2008

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arraytools.h"
#include "extend.h"
#include "lmc.h"
#include "oaoptions.h"
#include "tools.h"

#define stringify(name) #name

const int ncheckopts = 14;
enum checkmode_t {
        MODE_CHECK,
        MODE_REDUCE,
        MODE_REDUCERANDOM,
        MODE_REDUCETRAIN,
        MODE_REDUCETRAINRANDOM,
        MODE_CHECKJ4,
        MODE_REDUCEJ4,
        MODE_REDUCEJ4RANDOM,
        MODE_CHECKJ5X,
        MODE_REDUCEJ5X,
        MODE_HADAMARD,
        MODE_NONE,
        MODE_CHECK_SYMMETRY,
        MODE_CHECKJ5XFAST
};

/// return mode for oacheck in std::string
string modeString (checkmode_t checkmode) {
        std::string str = "invalid";
        switch (checkmode) {
        case MODE_CHECK:
		case MODE_REDUCE:
		case MODE_REDUCETRAIN:
		case MODE_REDUCETRAINRANDOM:
		case MODE_REDUCERANDOM:
		case MODE_HADAMARD:
		case MODE_NONE:
		case MODE_CHECKJ4:
		case MODE_CHECK_SYMMETRY:
		case MODE_REDUCEJ4:
		case MODE_CHECKJ5X:
		case MODE_CHECKJ5XFAST:
		case MODE_REDUCEJ5X:
		case MODE_REDUCEJ4RANDOM:
			str = stringify (checkmode);
            break;
        default:
                printf ("error: no such mode %d\n", checkmode);
                break;
        }
        return str;
}

/**
 * @brief Read in a list of array and check for each of the arrays whether they are in LMC form or not
 * @param argc
 * @param argv[]
 * @return
 */
int main (int argc, char *argv[]) {
        char *fname;
        int correct = 0, wrong = 0;

        /* parse command line options */
        AnyOption opt;
        opt.setFlag ("help", 'h'); /* a flag (takes no argument), supporting long and short form */
        opt.setOption ("loglevel", 'l');
        opt.setOption ("oaconfig", 'c');
        opt.setOption ("check");
        opt.setOption ("strength", 's');
        opt.setOption ("prune", 'p');
        opt.setOption ("output", 'o');
        opt.setOption ("mode", 'm');

        opt.addUsage ("oacheck: Check arrays for LMC form or reduce to LMC form");
        opt.addUsage ("Usage: oacheck [OPTIONS] [ARRAYFILE]");
        opt.addUsage ("");
        opt.addUsage (" -h  --help  			Prints this help");
        opt.addUsage (" -l [LEVEL]  --loglevel [LEVEL]	Set loglevel to number");
        opt.addUsage (" -s [STRENGTH]  --strength [STRENGTH] Set strength to use in checking");
        std::string ss = " --check [MODE]			Select mode: ";
        ss +=
            printfstring ("%d (default: check ), %d (reduce to LMC form), %d (apply random transformation and reduce)",
                          MODE_CHECK, MODE_REDUCE, MODE_REDUCERANDOM); //
        for (int i = 3; i < ncheckopts; ++i)
                ss += printfstring (", %d (%s)", static_cast< checkmode_t > (i), modeString ((checkmode_t)i).c_str ());
        opt.addUsage (ss.c_str ());
        opt.addUsage (" -o [OUTPUTFILE] Output file for LMC reduction");

        opt.processCommandArgs (argc, argv);

        algorithm_t algmethod = (algorithm_t)opt.getIntValue ("mode", MODE_AUTOSELECT);

        int loglevel = NORMAL;
        if (opt.getValue ("loglevel") != NULL || opt.getValue ('l') != NULL)
                loglevel = atoi (opt.getValue ('l')); // set custom loglevel
        setloglevel (loglevel);

        if (checkloglevel (QUIET)) {
                print_copyright ();
        }

        if (opt.getFlag ("help") || opt.getFlag ('h') || opt.getArgc () == 0) {
                opt.printUsage ();
                return 0;
        }

        checkmode_t mode = MODE_CHECK;
        if (opt.getValue ("check") != NULL)
                mode = (checkmode_t)atoi (opt.getValue ("check")); // set custom loglevel
        if (mode >= ncheckopts)
                mode = MODE_CHECK;

        logstream (QUIET) << "#time start: " << currenttime () << std::endl;

        double t = get_time_ms ();
        t = 1e3 * (t - floor (t));
        srand (t);

        int prune = opt.getIntValue ('p', 0);

        fname = opt.getArgv (0);
        setloglevel (loglevel);

        int strength = opt.getIntValue ('s', 2);
        if (strength < 1 || strength > 10) {
                printf ("Strength specfied (t=%d) is invalid\n", strength);
                exit (1);
        } else {
                log_print (NORMAL, "Using strength %d to increase speed of checking arrays\n", strength);
        }

        arraylist_t outputlist;
        char *outputfile = 0;
        bool writeoutput = false;
        if (opt.getValue ("output") != NULL) {
                writeoutput = true;
                outputfile = opt.getValue ('o');
        }

        arrayfile_t *afile = new arrayfile_t (fname);
        if (!afile->isopen ()) {
                printf ("Problem opening %s\n", fname);
                exit (1);
        }

        logstream (NORMAL) << "Check mode: " << mode << " (" << modeString (mode) << ")" << endl;

        /* start checking */
        double Tstart = get_time_ms (), dt;

        int index;
        array_link al (afile->nrows, afile->ncols, -1);
        for (int i = 0; i < afile->narrays; i++) {

                index = afile->read_array (al);
                array_t *array = al.array;

                arraydata_t arrayclass = arraylink2arraydata (al, 0, strength);
                arrayclass.lmc_overflow_check ();

                OAextend oaextend;

                lmc_t result = LMC_NONSENSE;

                LMCreduction_t *reduction = new LMCreduction_t (&arrayclass);
                LMCreduction_t *randtest = new LMCreduction_t (&arrayclass);
				array_link test_arraylink = al.clone();
				array_t *testarray = test_arraylink.array;
                randtest->transformation->randomize ();
                reduction->setArray ( al); 
                /* variables needed within the switch statement */
                switch (mode) {
                case MODE_CHECK:
				{
					/* LMC test with reduction code */
					reduction->mode = OA_TEST;
					result = LMCcheck(array, arrayclass, oaextend, *reduction);
					break;
				}
                case MODE_CHECKJ4: {
                        /* LMC test with special code */
                        reduction->mode = OA_TEST;
                        reduction->init_state = COPY;
                        array_link al (array, arrayclass.N, arrayclass.ncols, -10);
                        oaextend.setAlgorithm (MODE_J4, &arrayclass);
                        result = LMCcheck (al, arrayclass, oaextend, *reduction);
                        break;
                }
                case MODE_REDUCEJ4: {
                        /* LMC test with special code */
                        reduction->mode = OA_REDUCE;
                        array_link al (array, arrayclass.N, arrayclass.ncols, -10);
                        reduction->setArray (al);
                        result = LMCcheckj4 (al, arrayclass, *reduction, oaextend);
                        break;
                }
                case MODE_CHECK_SYMMETRY: {
					myprintf("MODE_CHECK_SYMMETRY not supported any more");
					exit(1);
                }
                case MODE_CHECKJ5X: {
                        oaextend.setAlgorithm (MODE_J5ORDERX, &arrayclass);

                        /* LMC test with special code */
                        reduction->mode = OA_TEST;
                        reduction->init_state = COPY;

                        result = LMCcheck (al, arrayclass, oaextend, *reduction);
                        break;
                }
                case MODE_CHECKJ5XFAST: {
                        oaextend.setAlgorithm (MODE_J5ORDERX, &arrayclass);

                        /* LMC test with special code */
                        reduction->mode = OA_TEST;
                        reduction->init_state = COPY;
                        if (1) {
                                array_link al (array, arrayclass.N, arrayclass.ncols, al.INDEX_NONE);
                                reduction->updateSDpointer (al);
                        }
                        result = LMCcheck (al, arrayclass, oaextend, *reduction);
                        break;
                }
                case MODE_REDUCEJ5X: {
                        OAextend oaextend;
                        oaextend.setAlgorithm (MODE_J5ORDERX, &arrayclass);

                        /* LMC test with special code */
                        printf ("oacheck: WARNING: MODE_CHECKJ5X: untested code, not complete\n");
                        reduction->mode = OA_REDUCE;
                        array_link al (array, arrayclass.N, arrayclass.ncols, -10);
                        result = LMCcheck (al, arrayclass, oaextend, *reduction);
                        break;
                }

                case MODE_REDUCE:
                        /* LMC reduction */
				{
					reduction->mode = OA_REDUCE;
					copy_array(array, reduction->array, arrayclass.N, arrayclass.ncols);
					dyndata_t dynd = dyndata_t(arrayclass.N);

					result = LMCreduction(array, array, &arrayclass, &dynd, reduction, oaextend);
					break;
				}
				case MODE_REDUCERANDOM: 
					{
						/* random LMC reduction */
						oaextend.setAlgorithm(MODE_ORIGINAL, &arrayclass);
						reduction->mode = OA_REDUCE;
						randtest->transformation->apply(array, testarray);
						copy_array(testarray, reduction->array, arrayclass.N, arrayclass.ncols);
						reduction->init_state = INIT;
						array_link alp = arrayclass.create_root(arrayclass.ncols, 1000);
						reduction->setArray(alp);
						dyndata_t dynd = dyndata_t(arrayclass.N);

						result = LMCreduction(testarray, testarray, &arrayclass, &dynd, reduction, oaextend);

						if (log_print(NORMAL, "")) {
							reduction->transformation->show();
						}
					} break;
				
                case MODE_REDUCETRAIN:
                        /* LMC reduction with train*/
                        reduction->mode = OA_REDUCE;
                        result = LMCreduction_train (al, &arrayclass, reduction, oaextend);
                        break;
                case MODE_REDUCETRAINRANDOM:
                        /* random LMC reduction with train*/
                        reduction->mode = OA_REDUCE;

                        randtest->transformation->apply (array, testarray);
                        result = LMCreduction_train (test_arraylink, &arrayclass, reduction, oaextend);
                        break;
                case MODE_HADAMARD:
						myprintf("MODE_HADAMARD not supported any more\n");
						exit(1);
                default:
                        result = LMC_NONSENSE;
                        std::cout << "function " << __FUNCTION__ << "line " << __LINE__ << "Unknown mode" << std::endl;
                        exit (1);
                        break;
                }

                bool aequal;
                switch (mode) {
                case MODE_CHECK:
                        if (result == LMC_LESS) {
                                ++wrong;
                                log_print (NORMAL, "Found array nr %i/%i NOT in lmc form:\n", i, afile->narrays);
                        } else {
                                ++correct;
                                log_print (NORMAL, "Found array nr %i/%i in lmc form.\n", i, afile->narrays);
                        }

                        break;
                case MODE_CHECKJ4:
                case MODE_CHECKJ5X:
                case MODE_CHECKJ5XFAST:
                case MODE_CHECK_SYMMETRY:
                        if (result == LMC_LESS) {
                                ++wrong;
                                log_print (NORMAL, "Found array nr %i/%i NOT in minimal form:\n", i, afile->narrays);
                        } else {
                                ++correct;
                                log_print (NORMAL, "Found array nr %i/%i in minimal form.\n", i, afile->narrays);
                        }

                        break;
                case MODE_REDUCE:
                case MODE_REDUCEJ4:
                case MODE_REDUCEJ5X:
                case MODE_REDUCETRAIN:
                        aequal = std::equal (array, array + arrayclass.N * arrayclass.ncols, reduction->array);

                        if ((!aequal) || reduction->state == REDUCTION_CHANGED) {
                                ++wrong;
                                log_print (QUIET, "Reduced array %i/%i to lmc form.\n", i, afile->narrays);
                                if (checkloglevel (NORMAL)) {
                                        log_print (QUIET, "Original:\n");

                                        print_array (array, afile->nrows, afile->ncols);
                                        print_array ("Reduction:\n", reduction->array, afile->nrows, afile->ncols);
                                        printf ("---------------\n");
                                }
                                if (checkloglevel (DEBUG)) {
                                        cout << "Transformation: " << endl;
                                        reduction->transformation->show (cout);

                                        rowindex_t row;
                                        colindex_t column;
                                        array_diff (array, reduction->array, arrayclass.N, arrayclass.ncols, row, column);
                                        cout << "Difference at: ";
                                        printf (" row %d, col %d\n", row, column);
                                }
                        } else {
                                ++correct;
                                log_print (QUIET, "Array %i/%i was already in lmc form:\n", i, afile->narrays);
                        }
                        break;
                case MODE_REDUCETRAINRANDOM:
                case MODE_REDUCERANDOM:
                case MODE_REDUCEJ4RANDOM:
                        aequal = std::equal (array, array + arrayclass.N * arrayclass.ncols, reduction->array);

                        if (!aequal) {
                                ++wrong;
                                log_print (SYSTEM,
                                           "Array %i/%i: reduced random transformation not equal to original (BUG!)\n",
                                           i, afile->narrays);
                                if (checkloglevel (NORMAL)) {
                                        print_array ("Original:\n", array, afile->nrows, afile->ncols);
                                        print_array ("Randomized:\n", testarray, afile->nrows, afile->ncols);
                                        print_array ("Reduction:\n", reduction->array, afile->nrows, afile->ncols);
                                        printf ("---------------\n");
                                }
                        } else {
                                ++correct;
                                log_print (QUIET, "Array %i/%i: reduced random transformation to original\n", i,
                                           afile->narrays);
                        }

                        break;
                case MODE_HADAMARD:
                        break;
                default:
                        std::cout << "function " << __FUNCTION__ << "line " << __LINE__ << "unknown mode" << std::endl;
                        break;
                }

                array_link al = array_link ((carray_t *)reduction->array, afile->nrows, afile->ncols, -1);
                if (result == LMC_MORE || !prune) {
                        outputlist.push_back (al);
                }
                delete randtest;
                delete reduction;
                destroy_array (testarray);
        }

        afile->closefile ();
        delete afile;

        if (writeoutput) {
                logstream (NORMAL) << "Writing output (" << outputlist.size () << " arrays)" << endl;
                writearrayfile (outputfile, &outputlist);
        }

        logstream (QUIET) << "#time end: " << currenttime () << std::endl;

        return 0;
}
