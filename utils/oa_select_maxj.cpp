/** \file oa_select_maxj.cpp

 C++ program: oa_select_maxj

 Select designs from a file with maxj value and write to outputfile
 
 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "anyoption.h"
#include "arraytools.h"
#include "tools.h"

bool strings_ends_with(const string& a, const string& b) {
    if (b.size() > a.size()) return false;
    return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void write_map(const char *filename, const std::map<long,long> &m) {
    map<long, long>::const_iterator it;

    FILE *fid = fopen(filename,"wt");
    fprintf(fid, "{");
    for (it = m.begin(); it != m.end(); it++) {
        fprintf(fid, " \"%ld\": %ld", it->first, it->second);
        if (std::next(it)!=m.end())
            fprintf(fid, ",");
    }
    fprintf(fid, "}\n");
    fclose(fid);
}

int main (int argc, char *argv[]) {
        AnyOption opt;

        /* parse command line options */
        opt.setFlag ("help", 'h'); /* a flag (takes no argument), supporting long and short form */
        opt.setOption ("input", 'i');
        opt.setOption ("output", 'o');
        opt.setOption ("verbose", 'v');
        opt.setOption ("format", 'f');

        opt.addUsage ("oa_select_maxj: Select designs with max(J5)");
        opt.addUsage ("Usage: oa_select_maxj -i [FILE] -o file [OPTIONS]");
        opt.addUsage ("");
        opt.addUsage (" -h  --help  		Prints this help ");
        opt.addUsage (" -i [FILE]  --input [FILE]		Input file ");
        opt.addUsage (" -o [STR]  --output [FILE]		Output file ");
        opt.addUsage (" -f [FORMAT]					Output format (default: TEXT, or BINARY; B, or DIFF; D, or DIFFZERO; Z) ");
        opt.addUsage (" -v [INTEGER]					Verbose (default: 2) ");
        opt.processCommandArgs (argc, argv);

        print_copyright ();

        /* parse options */
        if (opt.getFlag ("help") || opt.getFlag ('h')) {
                opt.printUsage ();
                exit (0);
        }

        if (opt.getValue ("input") == NULL && opt.getValue ('i') == 0) {
                printf ("No input file specified, use oa_select_maxj -h to get help\n");
                exit (1);
        }

        const char *inputfile;
        inputfile = opt.getValue ('i');

        const char *outputfile = opt.getValue ('o');
        
        if (!strings_ends_with(outputfile, ".oa")) {
            throw_runtime_exception("outputfile should end with .oa");   
        }
        
        std::string outputjson = std::string(outputfile);
        replace(outputjson, std::string(".oa"), std::string(".json"));
        std::map<long,long> maxj_values;

        int verbose = opt.getIntValue ("verbose", 2);

        std::string format = opt.getStringValue ('f', "BINARY");
        arrayfile::arrayfilemode_t mode = arrayfile_t::parseModeString (format);

        int nb = 8;
        if (mode == ABINARY_DIFFZERO) {
                nb = 1;
        }
        /* open the file containing the arrays */
        arrayfile_t *afile = new arrayfile_t (inputfile);
        const int rows = afile->nrows;
        const int cols = afile->ncols;
        const int narrays = afile->narrays;

        if (!afile->isopen ()) {
                printf ("oa_select_maxj: problem opening file %s\n", inputfile);
                exit (1);
        }
        if (verbose)
                printf ("oa_select_maxj: input %s, output %s\n", inputfile, outputfile);

        /* prepare the output files */
        arrayfile_t *outfid = 0;

        outfid = new arrayfile_t(outputfile, rows, cols, -1, mode, nb);
        
        jstruct_t *js = 0;
        array_link al (rows, cols, 0);
        
        array_link active_al5;
        int active_value = -1;
        int active_maxj = - 1;
        for (int i = 0; i < narrays; i++) {
                if (i % 500000 == 0 && verbose)
                        printf ("  oa_select_maxj: scanning array %d/%d\n", i, narrays);

           
                afile->read_array (al);

                array_link al5 = al.selectFirstColumns(5);
                
                if (js==0) {
                 // initial initialization
                    js = new jstruct_t(al5, 5);
                    std::vector<int> v= js->Fval();
                    for(std::vector<int>::iterator it = v.begin(); it != v.end(); ++it) {
                        maxj_values[(*it)]= 0;
                    }
                }
                
                if (al5==active_al5) {
                    // nothing to do
        }
                else {
                    js->calcj5(al5);
                
                    active_al5= al5;
                    active_maxj = js->maxJ ();
                    active_value = js->calculateF()[0];
                }
                
                maxj_values[active_maxj]++;
                
                if (active_value)
                    outfid->append_array (al);
        }

        write_map(outputjson.c_str(), maxj_values);
        
        outfid->finisharrayfile ();

        /* clean up memory */
        afile->closefile ();
        outfid->closefile();
        delete afile;
        delete outfid;

        return 0;
}
