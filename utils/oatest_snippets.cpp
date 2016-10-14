
    if ( 1 ) {
        array_link  al = exampleArray ( 2, 1 );
        al=al.randomrowperm();
        //array_link  al = exampleArray(22, 1);
        array_link out = sortrows ( al );
        al.showarray();
        out.showarray();
        exit ( 0 );
        al.show();
        if ( 1 ) {
            printf ( "## 3\n" );
            std::vector<int> f3 = al.FvaluesConference ( 3 );
            display_vector ( f3 );
            printf ( "\n" );
        }
        const int N = al.n_rows;
        jstructconference_t js ( N, 4 );
        jstructconference_t js2 ( al, 4 );
        js2.showdata ( 2 );

        printf ( "## 4\n" );
        std::vector<int> f4 = al.FvaluesConference ( 4 );
        std::vector<int> j4 = js.Jvalues();

        for ( int i=0; i< ( int ) js2.jvalues.size();  i++ ) {
            printf ( "i %d: %d -> %d\n", i, js2.jvalues[i], js2.jvalue2index.find ( js2.jvalues[i] )->second );
        }

        //printf("j4.size %d\n" , (int) j4.size() );
        display_vector ( j4 );
        printf ( "\n" );
        display_vector ( f4 );
        printf ( "\n" );
        exit ( 0 );
    }
    print_copyright();
    //cout << system_uname();
    setloglevel ( NORMAL );

    /* parse options */
    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
        opt.printUsage();
        exit ( 0 );
    }

    setloglevel ( SYSTEM );
    {
        array_link al = exampleArray ( 23, 1 );
        {
            jstructconference_t js ( al.n_rows, 4 );
            js.show();
            js.showdata();
        }
        jstructconference_t js ( al, 4 );
        js.show();
        js.showdata();


        printf ( "size %d\n", ( int ) js.values.size() );

        return 0;
    }

    array_link al = exampleArray ( 2 );
    al.showarray();
    Eigen::MatrixXd m1 = arraylink2eigen ( array2xf ( al ) );
    //std::cout << (m1) << std::endl; exit(0);

    Eigen::MatrixXd m2 = arraylink2eigen ( array2xf2 ( al ) );

    std::cout << ( m1-m2 ) << std::endl;

    for ( int i=0; i<1000000; i++ ) {
        //Eigen::MatrixXd m = array2xfeigen ( al);
        //Eigen::MatrixXd m = arraylink2eigen(array2xf ( al) );
        Eigen::MatrixXd m = arraylink2eigen ( array2xf2 ( al ) );
    }

    exit ( 0 );

    if ( 1 ) {
        arraylist_t ll = readarrayfile ( "dummy-24-4.oa" );
        array_link al=ll[0];
        al.show();

        int N=al.n_rows;
        conference_t ctype ( N, N );


        array_link al2 = ctype.create_root();
        array_link al3 = ctype.create_root_three();

        //arraylist_t lst; lst.push_back(al3); writearrayfile("test.oa", lst, arrayfile::ABINARY);

        int extcol=6;
        cperm_list ee= generateConferenceExtensions ( al2, ctype, extcol, 0, 0, 1 );
        printf ( "generated %d\n", ( int ) ee.size() );

        ee= generateConferenceExtensions ( al3, ctype, extcol, 0, 0, 1 );
        printf ( "generated %d\n", ( int ) ee.size() );

        DconferenceFilter filter ( al, 0,1,0 );

        filter.filterList ( ee,1 );

        exit ( 0 );
    }

    if ( 1 ) {

// FIXME: make unittest of this one
        //arraylist_t ll = readarrayfile ( "/home/eendebakpt/oatmp/conf/dconferencej1j3-36-8.oa" );
        //arraylist_t ll = readarrayfile ( "/home/eendebakpt/oatmp/conf/dconferencej1j3-28-4.oa" );
        arraylist_t ll = readarrayfile ( "/home/eendebakpt/oatmp/conf/dconferencej1j3-32-8.oa" );

        //ll[52]=ll[51];
        //ll[51]=ll[52];

        //ll[0].Fvalues(4);
        //ll[0].Fvalues(5);

        array_link alx = ll[0];
        alx.showproperties();
        int N = alx.n_rows;

        int fi=51;
        printf ( "first:\n" );
        ll[fi].showarray();
        printf ( "next :\n" );
        ll[fi+1].showarray();
        printf ( "  first diff index: %d\n", ll[fi].firstColumnDifference ( ll[fi+1] ) );


        conference_t ct ( N, 2*N );
        CandidateGeneratorDouble cgenerator ( array_link() , ct );
        cgenerator.verbose=1;
        cgenerator.showCandidates();

        double t0;
        if ( 0 ) {
            t0=get_time_ms();

            for ( size_t i=0; i<ll.size(); i++ ) {
                const array_link &al = ll[i];
                std::vector<cperm> cc1 = generateDoubleConferenceExtensionsInflate ( al, ct, 1, 1, 1 );
            }
            printf ( "   dtt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );
        }

        t0=get_time_ms();
        std::vector<cperm> cc1;
        for ( size_t i=0; i<ll.size(); i++ ) {
            printf ( "--- i %d -------\n", ( int ) i );
            const array_link &al = ll[i];
            int nc1=-1;
            if ( 1 ) {
                printfd ( "   _________________________\n" );
                cc1 = generateDoubleConferenceExtensionsInflate ( al, ct, 1, 1, 1 );
                nc1=cc1.size();
                printfd ( "   _________________________\n" );
            }
            //cperm tmp= cc1[0]; printf("size cc1: %d\n", (int)tmp.size() );

            cgenerator.verbose=2;
            std::vector<cperm> cc2 = cgenerator.generateDoubleConfCandidates ( al );
            //printf("size cc2: %d\n", cc2[0].size());



            int nc2=cc2.size();

            printf ( "%d: number of candidates: %d/%d\n", ( int ) i, nc1, nc2 );
            cgenerator.showCandidates();

            if ( i==52 && 0 ) {
                printf ( "cc1: \n" );
                showCandidates ( cc1 );
                printf ( "cc2: \n" );
                showCandidates ( cc2 );
                break;
            }
        }
        printf ( "   dtt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );

        exit ( 0 );
    }

    {
        int filterip=1;
        int filtersymm=1;
        int filterj3 = 1;
        array_link al = exampleArray ( r, 1 );
        al.show();

        al=al.selectFirstColumns ( xx );
        int kstart=xx-1;
        //al=al.selectFirstColumns(4);		int kstart=3;
        array_link als = al.selectFirstColumns ( kstart );

        //als.show();

        printf ( "find extensions of:\n" );
        al.showarray();
        al.row_symmetry_group().show();

        int N = al.n_rows;

        conference_t				ct = conference_t ( N, N );
        ct.j1zero=1;
        ct.j3zero=1;
        ct.itype=CONFERENCE_RESTRICTED_ISOMORPHISM;
        ct.ctype=conference_t::DCONFERENCE;


        //for(int z=0; z<10; z++)
        std::vector<cperm> cci = generateDoubleConferenceExtensionsInflate ( al, ct, verbose, 1, 1 );


        //printf ( "no symm:\n" );
        //cc = generateDoubleConferenceExtensions ( als, ct, verbose, 0, filterip, filterj3, 0 );

        if ( ix )
            exit ( 0 );
        printf ( "## full array (with symm):\n" );

        t0=get_time_ms();
        std::vector<cperm> cc3 = generateDoubleConferenceExtensions ( al, ct, verbose, 0, filterip, filterj3, 1 );
        for ( size_t i=0; i<cc3.size(); i++ ) {
            printf ( "  %d: ", ( int ) i );
            print_cperm ( cc3[i] );
            printf ( "\n" );
        }
        printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );
        exit ( 0 );

        printf ( "## full array (no symm):\n" );
        t0=get_time_ms();
        cc3 = generateDoubleConferenceExtensions ( al, ct, verbose, 0, filterip, filterj3, 0 );
        printf ( "   dt %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );


        //	printf ( "extend_conference: extended array %d/%d to %d arrays\n", ( int ) i, ( int ) lst.size(), nn );
        exit ( 0 );
    }



    {


        //		long imax = std::numeric_limits<long>::max(); printf("max for long %ld\n" , imax);
        //		 imax = std::numeric_limits<int>::max();
        //	printf("max for int %ld\n" , imax);
        //	exit(0);

        arraylist_t lst= readarrayfile ( input );
        printf ( "read %d arrays\n", ( int ) lst.size() );
        std::sort ( lst.begin(), lst.end() );



        Jcounter jcounter = calculateJstatistics ( input, jj, 1 );
        //Jcounter jcounter2 = calculateJstatistics ( input, jj, 1 ); jcounter += jcounter2;

        printf ( "--- results ---\n" );
        jcounter.show();
        jcounter.showPerformance();

        writeStatisticsFile ( "numbers-J.txt", jcounter, 1 );

        const char *numbersfile = "numbers-J.txt";
        Jcounter jc = readStatisticsFile ( numbersfile, verbose );
        jc.show();
        exit ( 0 );


        int jj=0;
        if ( xx ) {
            Pareto<mvalue_t<long>,array_link> pset;
            for ( int i=0; i<niter; i++ )
                addArraysToPareto ( pset, calculateArrayParetoJ5Cache<array_link>, lst, jj, verbose );
        }

        exit ( 0 );
    }




    /*

    if ( 1 ) {

    	arraydata_t ad;

    	conference_t ctype ( 8, 3 );
    	ctype.itype=CONFERENCE_RESTRICTED_ISOMORPHISM;
    	ctype.ctype=conference_t::DCONFERENCE;
    	arraylist_t lst= readarrayfile ( "test.oa" );
    	int verbose=1;

    	for ( int i=0; i< ( int ) lst.size() ; i++ ) {
    		array_transformation_t t = reduceOAnauty ( lst[i]+1, 2 );

    		array_link A = t.apply ( lst[i]+1 ) + ( -1 );
    		printf ( "array %d\n", i );
    		lst[i].showarray();
    		printf ( "array %d reduced\n", i );
    		A.showarray();
    	}

    	arraylist_t lst2 = addConstant ( lst, 0 );

    	arraylist_t outlist = selectConferenceIsomorpismClasses ( lst2, verbose, ctype.itype );
    	outlist = addConstant ( outlist, 0 );
    	writearrayfile ( "test2.oa", outlist );
    	exit ( 0 );
    }
    */

    if ( 1 ) {

        array_link al = exampleArray ( 5 );
        array_link alx = al;
        alx.randomperm();

        array_transformation_t t1 = reduceOAnauty ( al, 1 );
        //t1.show();	return 0;

        array_link alr1 = t1.apply ( al );


        array_transformation_t t2 = reduceOAnauty ( alx, 1 );
        array_link alr2 = t2.apply ( alx );


        printf ( "reduced:\n" );
        alr1.showarray();
        printf ( "random reduced:\n" );
        alr2.showarray();

        if ( alr1 != alr2 )
            printf ( "error: reductions unequal!\n" );

        return 0;
    }
    if ( 0 ) {

        arraylist_t ll = readarrayfile ( "/home/eendebakpt/tmp/sp0-split-10/sp0-split-10-pareto-64.2-2-2-2-2-2-2-2-2-2.oa" );
        array_link al=ll[0];

        int r0= ( al ).rank();
        int r=array2xf ( al ).rank();
        printf ( "rank: %d %d\n",  r0,r );

        arraydata_t arrayclass=arraylink2arraydata ( al, 1 );

        OAextend oaextend=OAextend();
        oaextend.checkarrays=0;
        oaextend.setAlgorithm ( MODE_J5ORDERXFAST, &arrayclass );
        setloglevel ( NORMAL );

        arrayclass.show();
        oaextend.info();

        printf ( "extend!\n" );
        al.show();

        int current_col=al.n_columns;
        arraylist_t extensions;
        int nr_extensions = extend_array ( al.array, &arrayclass, current_col,extensions, oaextend );

        //arraylist_t ww=extend_array(al, arrayclass, oaextend);

        return 0;
    }


    {
        arraylist_t lst = readarrayfile ( input );
        arraylist_t lstgood  = selectConferenceIsomorpismClasses ( lst, verbose );

        return 0;
    }


