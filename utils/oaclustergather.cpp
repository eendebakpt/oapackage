/** \file oaclustergather.cpp

 C++ program: oaclustergather

 Dedicated tool to gather data on cluster

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2015

 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <iostream>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arrayproperties.h"
#include "arraytools.h"
#include "tools.h"

#include "timsort.hpp"

#include "evenodd.h"
#include "pareto.h"

// convert to split vector
inline std::vector< int > tovec (std::vector< int > ii, int i) {
        ii.push_back (i);
        return ii;
}



/** read number of arrays and number of Pareto arrays from text file
 *
 * The format of the text file is a number of lines of the form:
 * [ncolumns]: [narrays] [npareto]
 */
bool readNumbersFile (const char *numbersfile, std::vector< long > &na, std::vector< long > &npareto, int kmin = -1,
                      int kmax = -1) {
        FILE *fid = fopen (numbersfile, "rt");
        if (fid == 0) {
                // printf ( "error: read numbers file %s: error\n" );
                return 0;
        }
        long np, n, k;

        std::fill (na.begin (), na.end (), -1);
        std::fill (npareto.begin (), npareto.end (), -1);

        while (1) {
                int r = fscanf (fid, "k %ld: %ld %ld\n", &k, &n, &np);
                if (r < 3)
                        break;
                if (k >= (int)na.size ())
                        na.resize (k + 1);
                if (k >= (int)npareto.size ())
                        npareto.resize (k + 1);
                na[k] = n;
                npareto[k] = np;
        }

        fclose (fid);

        return 1;
}

/// write number of arrays and number of Pareto arrays to text file
void writeNumbersFile (const char *numbersfile, const std::vector< long > &na, const std::vector< long > &npareto,
                       int kmin = -1, int kmax = -1) {

        if (kmin < 0)
                kmin = 0;
        if (kmax < 0)
                kmax = na.size () + 1;

        FILE *fid = fopen (numbersfile, "wt");
        for (int k = kmin; k <= kmax; k++) {
                fprintf (fid, "k %d: %ld %ld\n", k, na[k], npareto[k]);
        }

        fclose (fid);
}

std::vector< int > getSplit (AnyOption &opt) {
        int nsplit0 = opt.getIntValue ("nsplit0", -1);
        int nsplit1 = opt.getIntValue ("nsplit1", -1);
        int nsplit2 = opt.getIntValue ("nsplit2", -1);
        int nsplit3 = opt.getIntValue ("nsplit3", -1);
        int nsplit4 = opt.getIntValue ("nsplit4", -1);

        assert (nsplit0 != -1); // prevent legacy code from running

        std::vector< int > nsplit;
        nsplit.push_back (nsplit0);
        nsplit.push_back (nsplit1);
        nsplit.push_back (nsplit2);
        nsplit.push_back (nsplit3);
        nsplit.push_back (nsplit4);

        return nsplit;
}

std::vector< int > getLevels (AnyOption &opt) {
        int split0 = opt.getIntValue ("split0", -1);
        int split1 = opt.getIntValue ("split1", -1);
        int split2 = opt.getIntValue ("split2", -1);
        int split3 = opt.getIntValue ("split3", -1);
        int split4 = opt.getIntValue ("split4", -1);

        std::vector< int > lvls;
        if (split0 >= 0) {
                lvls.push_back (split0);
                if (split1 >= 0) {
                        lvls.push_back (split1);
                        if (split2 >= 0) {
                                lvls.push_back (split2);
                                if (split3 >= 0) {
                                        lvls.push_back (split3);
                                }
                        }
                }
        }
        return lvls;
}

void paretoInfo (const array_link alx) {
        std::vector< int > j5 = alx.Jcharacteristics (5);
        int j5max = vectormax (j5, 0);

        int v1 = (j5max == alx.n_rows);
        int v2 = 1 - v1;

        int N = alx.n_rows;
        int rank = array2xf (alx).rank ();
        std::vector< int > F4 = alx.Fvalues (4);
        std::vector< double > gwlp = alx.GWLP ();
        printf ("pareto data: %d ; ", rank);
        printf (" %d ", (int)(N * N * gwlp[4]));
        printf (" ; ");
        display_vector (F4);
        printf (" ; ");
        printf (" %d ; %d", v1, v2);
        printf ("\n");
}

const std::string filesep = "/";

/**
 * @brief Gather subset of arrays using Pareto optimality
 *
 *
 * @param argc
 * @param argv[]
 * @return Returns zero if we have a good run, otherwise an errorcode
 */
int main (int argc, char *argv[]) {
        /* parse command line options */
        AnyOption opt;
        opt.setFlag ("help", 'h'); /* a flag (takes no argument), supporting long and short form */
        opt.setFlag ("coptions");  /* print compile time options */
        opt.setOption ("output", 'o');
        opt.setOption ("verbose", 'v');
        opt.setOption ("basedir", 'b');
        opt.setOption ("method");

        opt.setOption ("split0");
        opt.setOption ("split1");
        opt.setOption ("split2");
        opt.setOption ("split3");
        opt.setOption ("split4");
        opt.setOption ("numbersfile");
        opt.setOption ("cleanrun");
        opt.setOption ("config", 'c');
        opt.setOption ("nsplit0");
        opt.setOption ("paretomethod");
        opt.setOption ("nsplit1");
        opt.setOption ("nsplit2");
        opt.setOption ("nsplit3");
        opt.setOption ("nsplit4");
        opt.setOption ("nparetodiff"); // allow difference in pareto files for .oa files and numbers file
        opt.setOption ("kmin");
        opt.setOption ("kmax");
        opt.setOption ("format", 'f');
        opt.setOption ("debug");
        // for debugging
        opt.setOption ("hm");
        opt.setOption ("lm");
        opt.setOption ("dindex");

        opt.addUsage ("Orthonal Array Cluster Gather: special tool");
        opt.addUsage ("Usage: oaclustergather [OPTIONS] ");
        opt.addUsage ("");
        opt.addUsage (" -h --help  			Prints this help ");
        opt.addUsage (" -v --verbose  			Verbose level (default: 1) ");
        opt.addUsage (" -b --basedir [DIR]  			Base calculation dir ");
        opt.addUsage (" -c --config [CONFIGFILE]  		Config file to use ");
        opt.addUsage (" -f [FORMAT]					Output format (TEXT, or default:BINARY) ");
        opt.addUsage (
            " --method [METHOD]				Default: 0 (Pareto). Other options: 1 (J-statistics)");
        opt.addUsage (" --numbersfile [FILENAME] 	Output name of number of arrays ");
        opt.addUsage (
            " --cleanrun [INTEGER]		If set to 1 abort when not all files are found. If set to zero "
            "generate partial results (default: 1)");
        opt.addUsage (" --nsplit0 [NUMBER]		Number of split files at level 0");
        opt.addUsage (" --nsplit1 [NUMBER]		Number of split files at level 1");
        opt.addUsage (" --nsplit2 [NUMBER]		Number of split files at level 2");
        opt.addUsage (" --nsplit3 [NUMBER]		Number of split files at level 3");
        opt.addUsage (" --split[x] [NUMBER]		Split index at level [x]");
        opt.addUsage (" --kmin [NUMBER] --kmax [NUMBER]	Min and max number of columns (inclusive) to process");
        opt.addUsage (" --paretomethod [0, 1]		If 1 add J5 to Pareto criterium");
        opt.addUsage (" -o [FILE] --output [FILE]	Output prefix for filtered arrays (default: no output) ");
        opt.addUsage ("");
        opt.addUsage ("Returns zero output code on succesfull run.");
        opt.processCommandArgs (argc, argv);

        print_copyright ();

        /* parse options */
        if (opt.getFlag ("help") || opt.getFlag ('h') || (opt.getValue ("coptions") != 0)) {
                opt.printUsage ();
                if (opt.getValue ("coptions") != 0) {
                        print_options (cout);
                }
                exit (0);
        }

        const char *outputprefix = 0;
        if (opt.getValue ("output") != 0)
                outputprefix = opt.getValue ('o');

        enum pareto_method { EVENODDPARETO, J5STATS, CONFERENCESTATS };

        const std::string basedir = opt.getStringValue ("basedir", "");
        const char *configfile = opt.getStringValue ("config", "oaconfig.txt");
        int verbose = opt.getIntValue ('v', 1);
        pareto_method method = (pareto_method)opt.getIntValue ("method", 0);
        int debug = opt.getIntValue ("debug", 0);
        const int dindex = opt.getIntValue ("dindex", 24);

        int needcleanrun = opt.getIntValue ("cleanrun", 1);
        int paretomethod = opt.getIntValue ("paretomethod", 0);
        int allowparetodiff = opt.getIntValue ("nparetodiff", 0);
        int kmin = opt.getIntValue ("kmin", 9);
        int kmax = opt.getIntValue ("kmax", 24);
        const char *numbersfile = opt.getStringValue ("numbersfile", 0);

        arraydata_t *adata = readConfigFile (configfile);
        assert (adata != 0);

        std::vector< int > nsplit = getSplit (opt);
        std::vector< int > lvls = getLevels (opt);

        if (verbose) {
                std::string splittag = splitTag (lvls);

                printf ("oaclustergather: basedir %s, kmin %d, kmax %d, verbose %d, method %d\n", basedir.c_str (),
                        kmin, kmax, verbose, method);
#ifdef DOOPENMP
                printf ("oaclustergather: openmp: num threads %d, max num threads %d\n", omp_get_num_threads (),
                        omp_get_max_threads ());
#endif
        }
        if (verbose) {
                std::cout << "#time start: " << currenttime () << std::endl;
        }
        double time0 = get_time_ms ();
        if (verbose) {
                std::cout << "oaclustergather: levels " << splitTag (lvls) << std::endl;
        }

        /* Loop over all subdirectories for all number of columns */
        int cleanrun = 1; /// indicates whether all necessary files have been found

        const char *methodtag = 0;
        if (method == J5STATS) {

                methodtag = "jstats";
                const int jjval = 5;

                Jcounter jc (adata->N, 5);

                int level = lvls.size ();
                std::string splittag = splitTag (lvls);

                // loop over all subsections
                for (int jj = 0; jj < nsplit[level]; jj++) {
                        std::string subdir = splitDir (tovec (lvls, jj));

                        if (verbose >= 1) {
                                if (verbose >= 2)
                                        printf ("\n");
                                printf ("#### oaclustergather: block %d (time %.1f [s])\n", jj,
                                        get_time_ms () - time0);
                                fflush (0);
                        }

                        std::string nfilesub0 = "" + (methodtag + ("-" + splitTag (tovec (lvls, jj)) + ".txt"));
                        std::string nfilesub = basedir + filesep + subdir + filesep + nfilesub0;

                        Jcounter jck = readStatisticsFile (nfilesub.c_str (), verbose >= 2);
                        bool b = jck.isOpen ();

                        if (verbose >= 2) {
                                if (b) {
                                        printf ("   --> read numbers file %s\n", nfilesub.c_str ());
                                } else {
                                        printf ("   --> could not read read numbers file %s\n", nfilesub.c_str ());
                                }
                        }
                        if (b) {
                                jc += jck;
                        }
                        // loop over all columns
                        for (int k = kmin; k <= kmax; k++) {
                                arraydata_t adata0 (adata, k);

                                std::string subfile0 = splitFile (tovec (lvls, jj)) +
                                                       printfstring ("-extend-%s.oa", adata0.idstr ().c_str ());
                                const std::string afile = basedir + filesep + subdir + filesep + subfile0;

                                bool existfile = file_exists (afile.c_str ());
                                if (jck.hasColumn (k) && (!existfile)) {
                                        if (verbose >= 2) {
                                                printfd ("statistics file has data and no file, continuing\n");
                                        }
                                        continue;
                                }

                                if (jck.hasColumn (k) && (existfile)) {
                                        if (verbose >= 2) {
                                                printfd ("statistics file has data and file exists, raising error\n");
                                        }
                                }

                                // get numbers of arrays from array file
                                int nnarrays = nArrays (afile.c_str ());
                                if ((!b) && (nnarrays < 0)) {
                                        printf ("no numbers file and no array file (afile %s)\n", afile.c_str ());
                                        printf (" numbers file is %s\n", nfilesub.c_str ());
                                        cleanrun = 0;
                                        if (needcleanrun) {

                                                exit (1);
                                        } else {
                                                continue;
                                        }
                                }
                                // na[k] += nnarrays;

                                if (verbose >= 2) {
                                        printfd ("calculate statistics on file %s\n", afile.c_str ());
                                }
                                Jcounter jcounter = calculateJstatistics (afile.c_str (), jjval, verbose >= 2);
                                if (!jcounter.validData ()) {
                                        printfd ("could not read statistics on file %s\n", afile.c_str ());
                                        printfd ("b %d, nnarrays: %d\n", b, nnarrays);
                                        exit (1);
                                }
                                if (verbose >= 2)
                                        jcounter.show ();

                                if (jcounter.N >= 0) {
                                        jc += jcounter;
                                }
                        }
                }

                fflush (0);

                if (verbose) {
                        printf ("totals (cleanrun %d):\n", cleanrun);
                        jc.showcompact ();
                }

                if (numbersfile && cleanrun) {
                        if (verbose)
                                printf ("writing numbers file %s\n", numbersfile);
                        writeStatisticsFile (numbersfile, jc, 0);
                }

                if (verbose >= 1) {
                        printf ("  total number of arrays: %ld, %.1f Marrays/hour\n", jc.narrays (),
                                double(3600. / 1e6) * double(jc.narrays ()) / (get_time_ms () - time0));
                }
        }

        if (method == EVENODDPARETO) {
                methodtag = "pareto";

                {
                        pareto_cb_cache paretofunction = calculateArrayParetoJ5Cache< array_link >;
                        if (paretomethod)
                                paretofunction = calculateArrayParetoJ5Cache< array_link >;
                        assert (paretomethod == 1); // other methods not implemented at this moment...
                }
                pareto_cb paretofunction = calculateArrayParetoJ5< array_link >;
                assert (paretomethod == 1); // other methods not implemented at this moment...

                arrayfile::arrayfilemode_t arrayfilemode =
                    arrayfile_t::parseModeString (opt.getStringValue ('f', "TEXT"));

                std::vector< long > na (kmax + 1);      /// number of arrays
                std::vector< long > npareto (kmax + 1); /// total number of pareto arrays

                Combinations::initialize_number_combinations (30);

                // loop over all columns
                for (int k = kmin; k <= kmax; k++) {
                        int cleanrunK = 1; /// indicates whether all necessary files for k columns have been found

                        if (verbose >= 2)
                                printf ("\n");
                        if (verbose >= 2)
                                printf ("#### oaclustergather: %d columns (time %.1f [s])\n", k,
                                        get_time_ms () - time0);
                        Pareto< mvalue_t< long >, array_link > pset;

                        arraydata_t adata0 (adata, k);

                        int level = lvls.size ();

                        std::string splittag = splitTag (lvls);

                        if (verbose) {
                                printf (" \n## oaclustergather: %d columns, gathering results for stage %d: split %s "
                                        "(time %.1f [s])\n",
                                        k, level, splittag.c_str (), get_time_ms () - time0);
                                fflush (0);
                        }
                        // loop over all subsections
                        for (int jj = 0; jj < nsplit[level]; jj++) {
                                if (debug) {
                                        int hm = opt.getIntValue ("hm", 569);
                                        int lm = opt.getIntValue ("lm", 568);
                                        if (jj < lm)
                                                continue;
                                        if (jj > hm)
                                                continue;
                                }
                                std::string subdir = splitDir (tovec (lvls, jj));
                                std::string nfilesub0 = "numbers-" + splitTag (tovec (lvls, jj)) + ".txt";
                                std::string nfilesub = basedir + filesep + subdir + filesep + nfilesub0;
                                std::vector< long > nasub (kmax + 1);
                                std::vector< long > nparetosub (kmax + 1);
                                bool b = readNumbersFile (nfilesub.c_str (), nasub, nparetosub, kmin, kmax);

                                if (verbose >= 2) {
                                        if (b) {
                                                printf ("   --> read numbers file %s\n", nfilesub.c_str ());
                                        } else {
                                                printf ("   --> could not read read numbers file %s\n",
                                                        nfilesub.c_str ());
                                        }
                                }

                                std::string subfile0 = splitFile (tovec (lvls, jj)) +
                                                       printfstring ("-extend-%s.oa", adata0.idstr ().c_str ());
                                std::string subfilepareto0 = splitFile (tovec (lvls, jj)) +
                                                             printfstring ("-pareto-%s.oa", adata0.idstr ().c_str ());
                                const std::string afile = basedir + filesep + subdir + filesep + subfile0;
                                std::string psourcefile = afile;
                                const std::string parfile = basedir + filesep + subdir + filesep + subfilepareto0;

                                if (b) {
                                        if (verbose >= 3)
                                                printf ("  --> nasub[%d]=%ld\n", k, nasub[k]);

                                        if (nasub[k] < 0) {
                                                long nnarrays = nArrays (afile.c_str ());
                                                if (verbose >= 3)
                                                        printfd ("   --> adding %ld arrays\n", nnarrays);
                                                // nasub[k]=nnarrays;
                                                na[k] += nnarrays;
                                        } else {
                                                // get number of arrays from numbers file
                                                na[k] += nasub[k];

                                                assert (nasub[k] >= 0);
                                        }
                                } else {
                                        // get numbers of arrays from array file
                                        int nnarrays = nArrays (afile.c_str ());
                                        if (nnarrays < 0) {
                                                printf ("no numbers file and no array file (afile %s)\n",
                                                        afile.c_str ());
                                                cleanrun = 0;
                                                cleanrunK = 0;
                                                if (needcleanrun) {

                                                        exit (1);
                                                } else {
                                                        continue;
                                                }
                                        }
                                        na[k] += nnarrays;
                                }

                                if (verbose >= 2)
                                        printf ("  --> check pareto file %s\n", subfilepareto0.c_str ());
                                bool paretofile = 0;
                                if (oa_file_exists (parfile)) {
                                        if (verbose >= 2) {
                                                printf ("  switching to Pareto file %s->%s!\n", subfile0.c_str (),
                                                        subfilepareto0.c_str ());
                                        }
                                        psourcefile = parfile;
                                        paretofile = 1;
                                }

                                int nn = nArrays (psourcefile.c_str ());
                                if (verbose >= 2)
                                        printf ("   ### file %s: %d arrays\n", base_name (afile).c_str (), nn);

                                if (nn >= 0)
                                        if (paretofile && b) {
                                                if (verbose >= 3) {
                                                        printf ("  --> both numbers file and pareto file, checking "
                                                                "whether numbers are equal\n");
                                                }
                                                // check
                                                if (nparetosub[k] != nn) {
                                                        printfd (" \n### warning: column %d: nparetosub[%d] %d, nn %d "
                                                                 "(number of arrays in .oa file)\n\n",
                                                                 jj, k, nparetosub[k], nn);
                                                        if (!allowparetodiff) {
                                                                cleanrun = 0;
                                                                cleanrunK = 0;
                                                        }
                                                }
                                        } else {
                                        }
                                else {
                                        // no source of Pareto files....
                                        if (cleanrun || verbose >= 1) { // only report error if we are in a clean run
                                                fprintf (stdout, "   error: file %s\n", psourcefile.c_str ());
                                        }
                                        cleanrun = 0;
                                        cleanrunK = 0;
                                        if (needcleanrun) {
                                                printfd ("found an error, aborting the program...\n");
                                                exit (1);
                                        } else {
                                                continue;
                                        }
                                }

                                if (debug) {
                                        int apos = arrayInFile (exampleArray (dindex), psourcefile.c_str ());

                                        if (apos >= 0) {
                                                myprintf ("found array in source file! setting debug to 10...\n");
                                                arraylist_t ll = readarrayfile (psourcefile.c_str ());
                                                for (size_t ij = 0; ij < ll.size (); ij++) {
                                                        array_link alx = ll[ij];
                                                        printf ("array %d: ", (int)ij);
                                                        paretoInfo (alx);
                                                }
                                                pset.show (2);
                                                pset.verbose = 3;
                                                debug = 10;
                                                // exit(0);
                                        }
                                }

                                long naread = 0; // number of arrays read
                                {
                                        // blocked read of arrays
                                        const long narraymax = 50000; // max number of arrays to read in a block
                                        arrayfile_t afile (psourcefile.c_str (), 0);

                                        if (!afile.isopen ()) {
                                                if (verbose) {
                                                        myprintf ("oaclustergather: problem with file %s\n",
                                                                  afile.filename.c_str ());
                                                }
                                                return 0;
                                        }

                                        long narrays = afile.narrays;
                                        if (narrays < 0)
                                                narrays = arrayfile_t::NARRAYS_MAX;

                                        int loop = 0;
                                        while (true) {
                                                long n = std::min (narraymax, narrays - naread);
                                                arraylist_t arraylist = afile.readarrays (n);

                                                gfx::timsort (arraylist.begin (),
                                                              arraylist.end ()); // sorting the arrays makes the rank
                                                                                 // calculations with subrank re-use
                                                                                 // more efficient
                                                if (verbose >= 2)
                                                        printf ("oaclustergather: read arrays in block %d: %d arrays "
                                                                "(%ld/%ld)\n",
                                                                loop, (int)arraylist.size (), naread, narrays);
                                                if (arraylist.size () <= 0)
                                                        break;
                                                addArraysToPareto (pset, paretofunction, arraylist, jj, verbose);
                                                naread += arraylist.size ();
                                                loop++;
                                        }
                                }
                                if (debug == 10) {
                                        pset.show (2);

                                        std::vector< array_link > llx = pset.allindices ();
                                        arraylist_t ll (llx.begin (), llx.end ());
                                        int jx = arrayInList (exampleArray (dindex), ll);
                                        myprintf ("after merge...jj %d, index in pareto list %d\n", jj, jx);
                                        if (jx < 0) {
                                                exit (0);
                                        }
                                }

                                if (verbose >= 3 || ((jj % 60 == 0 || (jj == nsplit[level] - 1)) && verbose >= 2)) {
                                        printf ("oaclustergather: file %d/%d, %ld arrays: %d Pareto values, %d Pareto "
                                                "elements\n",
                                                jj, nsplit[level], naread, pset.number (), pset.numberindices ());
                                        fflush (0);
                                }
                        }

                        if (verbose) {
                                printf ("  final pareto set %d cols (%d files): ", k, nsplit[level]);
                                if (verbose >= 2)
                                        pset.show (2);
                                else
                                        pset.show (2);
                        }
                        // write pareto set to disk

                        arraylist_t pp = pset.allindicesdeque ();
                        if (debug) {
                                printf ("---\n");
                                arrayInList (exampleArray (dindex), pp);
                                std::vector< array_link > llx = pset.allindices ();
                                arraylist_t ll (llx.begin (), llx.end ());
                                int jx = arrayInList (exampleArray (dindex), ll);
                        }
                        npareto[k] = pset.numberindices ();

                        if ((cleanrunK || (!needcleanrun)) && (outputprefix != 0)) {
                                std::string cdir = splitDir (lvls); // printfstring ( "sp0-split-%d/", split0 );
                                std::string xxx = splitFile (lvls);
                                if (lvls.size () > 0) {
                                } else {
                                        xxx = "results-j5evenodd";
                                }
                                std::string pfile0 =
                                    printfstring ("%s-%s-%s.oa", xxx.c_str (), outputprefix, adata0.idstr ().c_str ());
                                if (!cleanrunK) {
                                        pfile0 = replaceString (pfile0, ".oa", "-partial.oa");
                                }
                                std::string pfile = basedir + "/" + cdir + pfile0;

                                if (verbose) {
                                        printf ("  writing pareto file %s (%ld/%ld arrays)\n", pfile0.c_str (),
                                                (long)npareto[k], (long)na[k]);
                                        if (verbose >= 2)
                                                printfd ("pfile %s\n", pfile.c_str ());
                                }
                                writearrayfile (pfile.c_str (), &pp, arrayfilemode, adata->N, k);
                                if (debug) {
                                }
                        }

                        fflush (0);
                }

                if (verbose)
                        printf ("totals (cleanrun %d):\n", cleanrun);
                for (int k = kmin; k <= kmax; k++) {
                        if (verbose) {
                                printf ("%d columns: %ld arrays, %ld pareto\n", k, na[k], npareto[k]);
                        }
                }

                if (numbersfile && cleanrun) {
                        if (verbose)
                                printf ("writing numbers file %s\n", numbersfile);
                        writeNumbersFile (numbersfile, na, npareto, kmin, kmax);
                }

                const long natotal = std::accumulate (na.begin (), na.end (), (long)0);

                if (verbose >= 1) {
                        printf ("  total number of arrays: %ld, %.1f Marrays/hour\n", natotal,
                                double(3600. / 1e6) * double(natotal) / (get_time_ms () - time0));
                }
        }

        if (verbose) {
                std::cout << "#time end: " << currenttime () << std::endl;
                std::cout << "#time total: " << printfstring ("%.1f", get_time_ms () - time0) << " [s]" << std::endl;
                fflush (0);
        }

        /* free allocated structures */
        delete adata;

        if (cleanrun)
                return 0;
        else
                return 1;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 8; ;
