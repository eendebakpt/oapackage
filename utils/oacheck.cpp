
/** \file oacheck.cpp

 C++ program: oacheck

 oacheck allows to check wether a list of arrays in LMC form or not.

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2008

 Copyright: See LICENSE.txt file that comes with this distribution
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "oaoptions.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"
#include "arraytools.h"
#include "lmc.h"



#define stringify( name ) # name

#ifdef CLASSICCODE
const int ncheckopts=5;
#else
const int ncheckopts=14;
#endif
enum checkmode_t {MODE_CHECK, MODE_REDUCE, MODE_REDUCERANDOM, MODE_REDUCETRAIN, MODE_REDUCETRAINRANDOM, MODE_CHECKJ4, MODE_REDUCEJ4, MODE_REDUCEJ4RANDOM, MODE_CHECKJ5X, MODE_REDUCEJ5X, MODE_HADAMARD, MODE_NONE, MODE_CHECK_SYMMETRY, MODE_CHECKJ5XFAST};

/// return mode for oacheck in std::string
string modeString ( checkmode_t m ) {
    std::string str = "invalid";
    switch ( m ) {
    case MODE_CHECK:
        str =  stringify ( MODE_CHECK );
        break;
    case MODE_REDUCE:
        str =  stringify ( MODE_REDUCE );
        break;
    case MODE_REDUCETRAIN:
        str =  stringify ( MODE_REDUCETRAIN );
        break;
    case MODE_REDUCETRAINRANDOM:
        str =  stringify ( MODE_REDUCETRAINRANDOM );
        break;
    case MODE_REDUCERANDOM:
        str =  stringify ( MODE_REDUCERANDOM );
        break;
    case MODE_HADAMARD:
        str =  stringify ( MODE_HADAMARD );
        break;
    case MODE_NONE:
        str =  stringify ( MODE_NONE );
        break;
    case MODE_CHECKJ4:
        str =  stringify ( MODE_CHECKJ4 );
        break;
    case MODE_CHECK_SYMMETRY:
        str =  stringify ( MODE_CHECK_SYMMETRY );
        break;

    case MODE_REDUCEJ4:
        str =  stringify ( MODE_REDUCEJ4 );
        break;
    case MODE_CHECKJ5X:
        str =  stringify ( MODE_CHECKJ5X );
        break;
    case MODE_CHECKJ5XFAST:
        str =  stringify ( MODE_CHECKJ5XFAST );
        break;
    case MODE_REDUCEJ5X:
        str =  stringify ( MODE_REDUCEJ5X );
        break;
    case MODE_REDUCEJ4RANDOM:
        str =  stringify ( MODE_REDUCEJ4RANDOM );
        break;
    default:
        printf("error: no such mode %d\n", m);
        break;
    }
    return str;
}


void hadamardcheck(int i, array_t *array, const char *fname, const arraydata_t &ad, dyndata_t &dynd, LMCreduction_t *reduction, double Tstart, double &dt)
{
    int verbose=2;
    OAextend oaextend;

    /* variables needed within the switch statement */
    arraylist_t hlist;
    arraylist_t * classlist ;
    // it;
    lmc_t result;
    std::vector<int> indices;
    /* Hadamard */

    /* read list of arrays for comparing */
    classlist = new arraylist_t;
    readarrayfile ( fname, classlist );

    // with the current array make all Hadamard transformations
    reduction->mode = OA_REDUCE;
    dt = (get_time_ms() - Tstart);
    logstream ( NORMAL ) << "  time: " << printfstring("%.1f [s]", dt) << endl;

    indices.resize(ad.ncols);
    for(int ij=0; ij<ad.ncols; ij++) indices[ij]=ij;

    int nn=ad.ncols;
    int hcolmax=ad.ncols; // hcolmax=1;
    for ( colindex_t hcol=0; hcol<hcolmax; hcol++ ) {	// 4->
        logstream ( NORMAL ) << "array " << i << ": applying Hadamard transformation on column " << hcol << endl;

        array_t *cpy = clone_array ( array, ad.N, ad.ncols );
        apply_hadamard ( &ad, cpy, hcol );
        reduction->transformation->reset();
        result = LMCreduction_train ( cpy, &ad, &dynd, reduction, oaextend );
//        result = LMCreduce ( cpy, &ad, &dynd, reduction );

        int dbgval=-1;
        //array_link
        //reduction->transformation->apply(cpy);
        if (0) {
            printf("!! check reduction transformation\n");
            array_link tmparray(ad.N, ad.ncols, -1);
            array_link resultarray(reduction->array, ad.N, ad.ncols, -1);

            //reduction->transformation->inverse().apply(array, tmparray.array);
            reduction->transformation->apply(cpy, tmparray.array);
            if (verbose>=3) {
                tmparray.showarray();
                resultarray.showarray();
            }
            printf("equal: %d\n", tmparray==resultarray);
            dbgval = (tmparray==resultarray);

            return;//HACK

        }

        // LMCreduction_t reduction ( &ad );
        //reduction.mode=LMC_REDUCE;
        //lmc_t ret = LMCcheck ( alf, ad, oaextend, reduction );


        array_link link ( reduction->array, ad.N, ad.ncols, hcol );
        hlist.push_back ( link );
        destroy_array ( cpy );
        dt = (get_time_ms() - Tstart);
        logstream ( NORMAL ) << "  subtime: " << printfstring("%.1f [s]", dt) << endl;

        /* find index of this array in array list */
        arraylist_t::iterator it = find ( classlist->begin(),  classlist->end(), hlist[hcol] );
        indices[hcol] = (it-classlist->begin());

        if ( indices[hcol]==(int)classlist->size()) {
            cout << "Problem with col " << hcol << "!" << endl;
        }
        else
            logstream(NORMAL) << "   subsub: " << i << ", hcol " << hcol << ": " << indices[hcol] << printfstring(", dbgval %d", dbgval) << endl;
    }


    logstream(NORMAL) << "indices: ";
    for(int ii=0; ii<ad.ncols; ii++) {
        logstream(NORMAL) << indices[ii] << " ";
    }
    logstream(NORMAL) << endl;

    /* sort the resulting list of arrays, the first element is then a LMC representative */
    // sort ( hlist.begin(), hlist.end() );


    indexsort hlistsort(hlist );
    int findex = hlistsort.indices[0];

    // hlistsort.show(); printf("\n");

    if (0) {
        array_link a1 = classlist->at(101) ;
        array_link a2 =classlist->at(129);
        int vv = (a1<a2);
        printf("a1 < a2? %d\n", vv);

        printf("classlist[101]<classlist[129] %d\n", classlist->at(101) <  classlist->at(129) );
        printf("classlist[129]<classlist[101] %d\n",  classlist->at(129) <  classlist->at(101) );
        printf("classlist[121]<classlist[129] %d\n", classlist->at(121) <  classlist->at(129) );
        printf("classlist[121]<classlist[101] %d\n",  classlist->at(121) <  classlist->at(101) );

        for(int jj=0; jj<ad.ncols; jj++) {
            int val = indices[hlistsort.indices[jj] ] ;
            printf("jj %d: pos %d val %d\n", jj, hlistsort.indices[jj], val);
        }
    }

    if ( checkloglevel ( NORMAL ) ) {
        cout << "Found representative for Hadamard orbit:" << endl;
        //cout << hlist[0];
        cout << " hcol: " << hlist[findex].index << printfstring(", findex %d", findex)<< endl;
    }

    arraylist_t::iterator it = find ( classlist->begin(),  classlist->end(), hlist[findex] );
    cout << "Index of " << i << "  is " << ( it-classlist->begin() ) << "/" << classlist->size() << " (hcol " << findex << ")" << endl;

    hlist.clear();
}


/**
 * @brief Read in a list of array and check for each of the arrays whether they are in LMC form or not
 * @param argc
 * @param argv[]
 * @return
 */
int main ( int argc, char* argv[] ) {
    char	*fname;
    int correct = 0, wrong = 0;



    /* parse command line options */
    AnyOption opt;
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "loglevel", 'l' );
    opt.setOption ( "oaconfig", 'c' );
    opt.setOption ( "check" );
    opt.setOption ( "strength", 's' );
    opt.setOption ( "prune", 'p' );
    opt.setOption ( "output", 'o' );
    opt.setOption ( "mode", 'm' );

    opt.addUsage ( "oacheck: Check arrays for LMC form or reduce to LMC form" );
    opt.addUsage ( "Usage: oacheck [OPTIONS] [ARRAYFILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h  --help  			Prints this help" );
    opt.addUsage ( " -l [LEVEL]  --loglevel [LEVEL]	Set loglevel to number" );
    opt.addUsage ( " -s [STRENGTH]  --strength [STRENGTH] Set strength to use in checking" );
    std::string ss = " --check [MODE]			Select mode: ";
    ss += printfstring("%d (default: check ), %d (reduce to LMC form), %d (apply random transformation and reduce)", MODE_CHECK, MODE_REDUCE, MODE_REDUCERANDOM); //
    for(int i=3; i<ncheckopts; ++i)
        ss += printfstring(", %d (%s)", static_cast<checkmode_t>(i), modeString( (checkmode_t)i).c_str() );
    opt.addUsage ( ss.c_str() );
    opt.addUsage ( " -o [OUTPUTFILE] Output file for LMC reduction" );


    opt.processCommandArgs ( argc, argv );

    algorithm_t algmethod =  (algorithm_t)opt.getIntValue("mode", MODE_AUTOSELECT);

    int loglevel = NORMAL;
    //int verbose = loglevel;
    if ( opt.getValue ( "loglevel" ) != NULL || opt.getValue ( 'l' ) != NULL )
        loglevel = atoi ( opt.getValue ( 'l' ) );	//set custom loglevel
    setloglevel ( loglevel );

    if (checkloglevel(QUIET)) {
        print_copyright();
    }

    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() ==0 ) {
        opt.printUsage();
        return 0;
    }

    checkmode_t mode = MODE_CHECK;
    if ( opt.getValue ( "check" ) != NULL )
        mode = ( checkmode_t ) atoi ( opt.getValue ( "check" ) );	//set custom loglevel
    if (mode>=ncheckopts) mode = MODE_CHECK;

    //logstream(QUIET) << "oacheck: mode  " << modeString(mode) << std::endl;
    logstream(QUIET) << "#time start: "<< currenttime() << std::endl;

    double t=get_time_ms();
    t= 1e3*(t-floor(t));
    srand(t);

    int prune = opt.getIntValue('p', 0);

    fname = opt.getArgv ( 0 );
    setloglevel ( loglevel );

    int strength = opt.getIntValue('s', 2);
    if (strength<1 || strength>10) {
        printf("Strength specfied (t=%d) is invalid\n", strength );
        exit(1);
    }
    else {
        log_print( NORMAL, "Using strength %d to increase speed of checking arrays\n", strength );
    }

    arraylist_t outputlist;
    char *outputfile = 0;
    bool writeoutput = false;
    if(opt.getValue("output")!=NULL) {
        writeoutput = true;
        outputfile = opt.getValue('o');
    }

    arrayfile_t *afile = new arrayfile_t ( fname );
    if ( ! afile->isopen() ) {
        printf ( "Problem opening %s\n", fname );
        exit ( 1 );
    }

    logstream ( NORMAL ) << "Check mode: " << mode << " (" << modeString ( mode ) << ")" << endl;


#ifdef OAANALYZE_DISCR
    clear_discriminant(afile->nrows,afile->ncols);
    analysis_init_values();
#endif

    /* start checking */
    double Tstart = get_time_ms(), dt;

    int index;
    array_link al(afile->nrows, afile->ncols, -1);
    for ( int i = 0; i < afile->narrays; i++ ) {

        index = afile->read_array(al);
        array_t *array=al.array;

        arraydata_t ad = arraylink2arraydata(al, 0,strength);
	ad.lmc_overflow_check();
	
        dyndata_t dynd = dyndata_t ( ad.N );
        OAextend oaextend;
//		oaextend.setAlgorithm(algmethod, &ad);

        lmc_t result = LMC_NONSENSE;

        LMCreduction_t *reduction = new LMCreduction_t ( &ad );
        LMCreduction_t *randtest = new LMCreduction_t ( &ad );
        array_t *testarray = clone_array ( array, ad.N, ad.ncols );
        randtest->transformation->randomize(); // HACK, set back!
        //randtest->transformation->randomizecolperm();
        reduction->setArray(al);	// FIXME: for reduction of arrays not in root-form we need to add the special initialization

        /* variables needed within the switch statement */
        switch ( mode ) {
        case MODE_CHECK:
            /* LMC test with reduction code */
            reduction->mode = OA_TEST;
            //result = LMCreduceFull ( array, array, &ad, &dynd, reduction );
            result = LMCreduce ( array, array, &ad, &dynd, reduction, oaextend );
            break;
        case MODE_CHECKJ4:
        {
            /* LMC test with special code */
            reduction->mode = OA_TEST;
            reduction->init_state = COPY;
            array_link al(array, ad.N, ad.ncols, -10);
            oaextend.setAlgorithm(MODE_J4, &ad);
            //reduction->show(0);
            result = LMCcheck ( al, ad, oaextend, *reduction);
            break;
        }
        case MODE_REDUCEJ4:
        {
            /* LMC test with special code */
//					printf("WARNING: not complete\n");
            reduction->mode = OA_REDUCE;
            array_link al(array, ad.N, ad.ncols, -10);
            reduction->setArray(al);
            //copy_array(al.array, reduction->array, ad.N, ad.ncols);
            result = LMCcheckj4 ( al, ad, *reduction, oaextend );
            break;
        }
        case MODE_CHECK_SYMMETRY:
        {
            const int verbose=0;
            //OAextend oaextend;

       {
                oaextend.setAlgorithm(MODE_LMC_SYMMETRY, &ad);

                oaextend.j5structure=J5_45;
		if (verbose>=2) {
                printf("oacheck : MODE_CHECK_SYMMETRY: array %d...\n", i );
		}

                if (verbose>=2) {
                    printf("oacheck : MODE_CHECK_SYMMETRY...\n" );
                    // al.showarray();
                }
                /* LMC test with special code */
                reduction->reset();
                reduction->mode = OA_TEST;
                reduction->init_state=INIT;
                array_link al(array, ad.N, ad.ncols, -10);
                reduction->setArray(al);
                result = LMCcheck ( al, ad, oaextend, *reduction );
		 if (verbose) {
                    printf("oacheck : MODE_CHECK_SYMMETRY: result %d\n", (int) result );
		 }
            }

            arraydata_t adata(ad);
    //mydebug2(al, adata, oaextend, 0);
    
            if (0) {
                if (verbose)
                    printf("oacheck : MODE_CHECK_SYMMETRY: checking algorithm using direct calculation...\n" );
                     oaextend.setAlgorithm(MODE_J5ORDERX, &ad);
                //    oaextend.j5structure=J5_45;

                int dverbose=0;

                OAextend x = oaextend;
                //arraydata_t adx(ad);
                //x.setAlgorithm(MODE_J5ORDERX, &adx);
                LMCreduction_t reductionsub = calculateSymmetryGroups( al.deleteColumn(-1), ad,  oaextend, 0);

                if (dverbose) {
                    printf("LMCcheck: MODE_LMC_SYMMETRY: testing: strength %d, al ", ad.strength);
                    al.show();
                    reductionsub.symms.showColperms(1);
                }
                reduction->clearSymmetries();
                reduction->setArray(al);

//		        rn = LMCcheckXX(al, adata, oaextend, reduction, reductionsub, dverbose) ;

//		copy_array(al.array, reduction.array, adata.N, adata.ncols); // hack?
                reduction->setArray(al);
                lmc_t result3 = LMCcheckSymmetryMethod(al, ad, oaextend, *reduction, reductionsub, dverbose) ;
                if(result!=result3) {
                    printf("oacheck: MODE_CHECK_SYMMETRY: %d (~ %d)\n ", result, result3);
                }
              //  result=result3;

            }

            if (0)
            {
                oaextend.setAlgorithm(MODE_J5ORDERX, &ad);
		        LMCreduction_t *reduction2 = new LMCreduction_t ( &ad );
	      reduction2->setArray(al);
                reduction2->mode = OA_TEST;
				//    reduction2->symms.showSymmetries();
				 //   exit(0);

                lmc_t result2 = LMCcheck ( al, ad, oaextend, *reduction2 );
                if (result!=result2) {
                    printfd("ERROR: oacheck: compare with original: MODE_CHECK_SYMMETRY: %d (should be %d)\n ", result, result2);
                    reduction2->transformation->show();
		    
		    //reduction->c
		    
		    LMCreduction_t reductionsub= calculateSymmetryGroups( al.deleteColumn(-1),  ad,  oaextend, 0);
		    reductionsub.symms.show(2);
		    reductionsub.symms.showSymmetries(2);
		    
		    
		    delete reduction2;
                }
            }
            break;
        }
        case MODE_CHECKJ5X:
        {
            oaextend.setAlgorithm(MODE_J5ORDERX, &ad);

            /* LMC test with special code */
            reduction->mode = OA_TEST;
            reduction->init_state=COPY;
	    
            result = LMCcheck ( al, ad, oaextend, *reduction );
            break;
        }
        case MODE_CHECKJ5XFAST:
        {
            //OAextend oaextend;
            oaextend.setAlgorithm(MODE_J5ORDERX, &ad);

            /* LMC test with special code */
            reduction->mode = OA_TEST;
            reduction->init_state=COPY;
	    if (1) {
	      // TODO: make special MODE for this
            array_link al(array, ad.N, ad.ncols, -10);
	    reduction->updateSDpointer(al);
	   // reduction->sd->show(2);
	    }
            result = LMCcheck ( al, ad, oaextend, *reduction );
            break;
        }
        case MODE_REDUCEJ5X:
        {
            OAextend oaextend;
            oaextend.setAlgorithm(MODE_J5ORDERX, &ad);

            /* LMC test with special code */
            printf("oacheck: WARNING: MODE_CHECKJ5X: untested code, not complete\n");
            reduction->mode = OA_REDUCE;
            array_link al(array, ad.N, ad.ncols, -10);
            result = LMCcheck ( al, ad, oaextend, *reduction );
            break;
        }
        case MODE_REDUCEJ4RANDOM:
        {
            /* LMC test with special code */
            printf("WARNING: untested code\n");
            reduction->mode = OA_REDUCE;

            randtest->transformation->apply ( array, testarray );

            copy_array(testarray, reduction->array, ad.N, ad.ncols);
            check_root_update(testarray, ad, reduction->array);

            array_link al(testarray, ad.N, ad.ncols, -10);
            result = LMCcheckj4 ( al, ad, *reduction , oaextend);
            break;
        }
        case MODE_REDUCE:
            /* LMC reduction */
            reduction->mode = OA_REDUCE;
            copy_array ( array, reduction->array, ad.N, ad.ncols );

            //cout << "input: " << endl;
            //show_array(array, afile->ncols, afile->nrows );
            //cout << "input: " << endl;
            //show_array(reduction->array, afile->ncols, afile->nrows );

            result = LMCreduce ( array, array, &ad, &dynd, reduction, oaextend );
            break;
        case MODE_REDUCERANDOM:
        {
            /* random LMC reduction */
            oaextend.setAlgorithm(MODE_ORIGINAL, &ad);
            reduction->mode = OA_REDUCE;
            randtest->transformation->apply ( array, testarray );
            copy_array ( testarray, reduction->array, ad.N, ad.ncols );
            reduction->init_state=INIT;
            array_link alp = rootPlus(ad);
            reduction->setArray(alp);

            //oaextend.info();
            //reduction->show();
            result = LMCreduce ( testarray, testarray, &ad, &dynd, reduction, oaextend);
            //printf("oacheck: debug: result %d\n", result);

            if(log_print(NORMAL, "")) {
                reduction->transformation->show();
            }
        }
        break;
        case MODE_REDUCETRAIN:
            /* LMC reduction with train*/
            reduction->mode = OA_REDUCE;

            result = LMCreduction_train ( array, &ad, &dynd, reduction , oaextend);
            break;
        case MODE_REDUCETRAINRANDOM:
            /* random LMC reduction with train*/
            reduction->mode = OA_REDUCE;

            randtest->transformation->apply ( array, testarray );
            result = LMCreduction_train ( testarray, &ad, &dynd, reduction , oaextend);
            break;
        case MODE_HADAMARD:
            /* Hadamard */
            oaextend.setAlgorithm(MODE_J4, &ad);
            printf("MODE_HADAMRD: oaextend alg: %s\n", oaextend.getAlgorithmName().c_str() );

            if (1)
                hadamardcheck(i,array, fname, ad, dynd, reduction, Tstart, dt);

            break;
        default:
            result = LMC_NONSENSE;
            std::cout << "function " << __FUNCTION__ << "line " << __LINE__ << "Unknown mode" << std::endl;
            exit ( 1 );
            break;
        }


        bool aequal;
        switch ( mode ) {
        case MODE_CHECK:
            if ( result == LMC_LESS ) {
                ++wrong;
                log_print (NORMAL, "Found array nr %i/%i NOT in lmc form:\n", i, afile->narrays );
                //show_array(array, afile->ncols, afile->nrows );
                //printf ( "---------------\n" );
            } else {
                ++correct;
                log_print (NORMAL, "Found array nr %i/%i in lmc form.\n", i, afile->narrays );
            }

            break;
        case MODE_CHECKJ4:
        case MODE_CHECKJ5X:
        case MODE_CHECKJ5XFAST:
        case MODE_CHECK_SYMMETRY:
            if ( result == LMC_LESS ) {
                ++wrong;
                log_print (NORMAL, "Found array nr %i/%i NOT in minimal form:\n", i, afile->narrays );
                //show_array(array, afile->ncols, afile->nrows );
                //printf ( "---------------\n" );
            } else {
                ++correct;
                log_print (NORMAL, "Found array nr %i/%i in minimal form.\n", i, afile->narrays );
            }

            break;
        case MODE_REDUCE:
        case MODE_REDUCEJ4:
        case MODE_REDUCEJ5X:
        case MODE_REDUCETRAIN:
            aequal = std::equal(array, array+ad.N*ad.ncols, reduction->array);

            if ( ( !aequal ) || reduction->state==REDUCTION_CHANGED ) {
                ++wrong;
                log_print (QUIET, "Reduced array %i/%i to lmc form.\n", i, afile->narrays );
                if (checkloglevel(NORMAL)) {
                    log_print (QUIET, "Original:\n" );

                    //print_array(array, afile->ncols, afile->nrows );
                    print_array(array, afile->nrows, afile->ncols );
                    //print_array("Reduction:\n", reduction->array, afile->ncols, afile->nrows );
                    print_array("Reduction:\n", reduction->array, afile->nrows, afile->ncols );
                    printf ( "---------------\n" );
                }
                if (checkloglevel(DEBUG)) {
                    cout << "Transformation: " << endl;
                    reduction->transformation->show(cout);

                    rowindex_t r;
                    colindex_t c;
                    array_diff(array, reduction->array, ad.N, ad.ncols, r, c);
                    cout << "Difference at: " ;
                    printf(" row %d, col %d\n", r, c);
                }
            } else {
                ++correct;
                log_print (QUIET, "Array %i/%i was already in lmc form:\n", i, afile->narrays );
            }
            break;
        case MODE_REDUCETRAINRANDOM:
        case MODE_REDUCERANDOM:
        case MODE_REDUCEJ4RANDOM:
            aequal = std::equal(array, array+ad.N*ad.ncols, reduction->array);

            if ( !aequal ) {
                ++wrong;
                log_print (SYSTEM, "Array %i/%i: reduced random transformation not equal to original (BUG!)\n", i, afile->narrays );
                if (checkloglevel(NORMAL)) {
                    print_array("Original:\n", array, afile->nrows, afile->ncols );
                    print_array("Randomized:\n", testarray, afile->nrows, afile->ncols );
                    //printf("Reduced:\n"); reduction->transformation->print_transformed(testarray);
                    print_array("Reduction:\n",reduction->array, afile->nrows, afile->ncols );
                    printf ( "---------------\n" );
                }
            } else {
                ++correct;
                log_print (QUIET, "Array %i/%i: reduced random transformation to original\n", i, afile->narrays );
            }

            break;
        case MODE_HADAMARD:
            break;
        default:
            std::cout << "function " << __FUNCTION__ << "line " << __LINE__ << "unknown mode" << std::endl;
            break;
        }

        array_link al = array_link((carray_t *)reduction->array, afile->nrows, afile->ncols, -1);
        if (result==LMC_MORE || !prune) {
            outputlist.push_back(al);
        }
        delete randtest;
        delete reduction;
        destroy_array ( testarray );
    }

#ifdef OAANALYZE_DISCR
    logstream(NORMAL) << "Discriminant: " << endl;
    print_discriminant(afile->nrows,afile->ncols);
    analysis_show_counter("root_level_perm");

#endif

    //  delete [] s;
    //destroy_array ( array );
    //close_arrayfile ( afile );
    afile->closefile();
    delete afile;

    if (writeoutput) {
        logstream(NORMAL) << "Writing output (" << outputlist.size() << " arrays)" << endl;
        writearrayfile(outputfile, &outputlist);
    }



    logstream(QUIET) << "#time end: "<< currenttime() << std::endl;

    return 0;
}
