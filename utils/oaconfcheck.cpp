/** \file oatest.cpp

C++ program: oatest

oatest: tool for testing new algorithms

Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2014

Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <map>

#include "arraytools.h"
#include "arrayproperties.h"
#include "anyoption.h"
#include "tools.h"
#include "extend.h"

#include "evenodd.h"
//<<<<<<< HEAD
//#include "oadevelop.h"
//=======
//>>>>>>> upstream/master
#include "lmc.h"

#include "conference.h"

/* Define functions that might be in the mathtools.h  file */
template <class numtype>
/**
 * Initialiaze a permutation
 * @param perm
 * @param len
 */
void init_signperm ( std::vector<numtype> &signperm)
{
    for ( size_t i=0; i<signperm.size(); i++ )
        signperm[i]=1;
}

/* Start the real functions*/
/*
rowlevel_permutation: find the row-level permutation of a column
al: array link, pointer
rowperm: current row permutation, pointer, rowsort_t
colperm: current column permutation, pointer
rowsignperm: current rowsign permutation, pointer
column: the column to find the permutation, int
n_rows: number of rows, n_rows, index
*/
//template<class NumType, class NumTypeIn>
void rowlevel_permutation ( const array_link &al, rowsort_t *rowperm, std::vector<int> &colperm, std::vector<int> &rowsignperm, std::vector<int> &colsignperm, const rowindex_t n_rows, int column ) {

    //std::vector<NumType> rowsignone ( nrows );
    //init_signperm ( rowsignone );
    //see what level-permutation we need to get a column with all ones and one zero
    // We don't copy the array but use the current modifications given by the rowsignperm, colperm, and rowperm
    for ( rowindex_t r=0; r < n_rows; r++ ) {
        if ( ( (rowsignperm[rowperm[r].val]*colsignperm[colperm[column]]) * al.atfast( rowperm[r].val, colperm[column] )) < 0 )
            rowsignperm[ rowperm[r].val ] = -1; // assign new values to the rowsignperm vector
    }

}

int LMC0check_firstrow ( const array_link &al, std::vector<int> &colperm,  int ci ) {
    printf ( "not implemented...\n" );
    return LMC_LESS;

}

/* calc_rowsort: Helper function to sort rows of a conference array
al: link to the array, pointer
colperm: current column permutation of the array, pointer
n_rows: number of rows of the array, int
n_cols: number of columns of the array, int
*/
indexsort calc_rowsort(const array_link &al, rowsort_t *rowperm, std::vector<int> &colperm, std::vector<int> &rowsignperm, std::vector<int> &colsignperm, const rowindex_t n_rows, const colindex_t n_cols)
{
    // Find
    std::vector<mvalue_t<int> > rr;
    for ( int i=0; i < n_rows; i++ ) {
        mvalue_t<int> m;
        for ( int k=0; k < n_cols; k++ )
            // We transform the elements (0,1,-1) to (0,1,2)
            // To perform the sort. We use tha transformations of the array
            m.v.push_back ( ( ( (colsignperm[colperm[k]]*rowsignperm[rowperm[i].val])*al.at( rowperm[i].val, colperm[k] ) )+3) % 3 );
        rr.push_back ( m );
    }
    indexsort is ( rr );
    return rr;
}

/*
LMC0_sortrows: Compute the new row sort for the array
\brief Calculate the new order and assign it to the existing
al: link to the array, pointer
colperm: current column permutation of the array, pointer
n_rows: number of rows of the array, int
n_cols: number of columns of the array, int
*/
void LMC0_sortrows ( const array_link &al, rowsort_t *rowperm, std::vector<int> &colperm, std::vector<int> &rowsignperm, std::vector<int> &colsignperm, const rowindex_t n_rows, const colindex_t n_cols )
{
    indexsort aa = calc_rowsort(al, rowperm, colperm, rowsignperm, colsignperm, n_rows, n_cols);
    // Assign the new sorting of the rows
    for (rowindex_t j = 0; j < n_rows; j++){
        rowperm[j].val = aa.indices[j];
    }

}

lmc_t lmc0_compare_columns ( const array_link &al, rowsort_t *rowperm, std::vector<int> &colperm, int column, std::vector<int> &rowsignperm, std::vector<int> &colsignperm ) {
    const int nrows=al.n_rows;
    // Compare columns with respect to the current column permutation and row ordering
    // We also include the current rowsign permutation
    //printf("compare_columns: source %d, target %d\n", column, colperm[column]);
    /* Compare the two columns according to the LMC0 ordering*/
    for ( int r=0; r<nrows; r++ ) {
    // Auxiliar variable to store the current value of the array
        int value_cdesign_trans = (colsignperm[colperm[column]]*rowsignperm[rowperm[r].val]) * al.atfast ( rowperm[r].val, colperm[column] );
        if ( ((al.atfast ( r, column )+3) % 3) <  ( (value_cdesign_trans+3) % 3) ) // Transform the elements from (0, 1, -1) to (0, 1, 2)
            return LMC_MORE;
        if ( ((al.atfast ( r, column )+3) % 3) >  ( (value_cdesign_trans+3) % 3) )
            return LMC_LESS;
    }
    return LMC_EQUAL;
}

lmc_t LMC0_columns ( const array_link &al, rowsort_t *rowperm, std::vector<int> &colperm, int column, std::vector<int> &rowsignperm, std::vector<int> &colsignperm, const int ncols, const int nrows, int verbose=0 ) {
    //const int ncols = al.n_columns;
    //const int nrows = al.n_rows;
    lmc_t r = LMC_NONSENSE;

    if ( verbose ) {
        printf ( "LMC0_columns: column %d, colperm ", column );
        display_vector ( colperm );
        printf ( "\n" );
    }
    for ( int c=column; c<ncols ; c++ ) {
        // select first column. Swap columns
        int col = colperm[c];
        colperm[c]=colperm[column];
        colperm[column]=col;
        /* i. Apply the correct column level permutation to make the element X(1,k) equal to 1*/
        // Current value in the X(1,k) position
        int current_sign_col = colsignperm[colperm[column]]; // save current sign
        int current_val_firstrow = (rowsignperm[rowperm[0].val]*current_sign_col)*(al.at(rowperm[0].val, colperm[column]));
        // Assign correct column sign permutation
        colsignperm[ colperm[column] ] = colsignperm[ colperm[column] ] * current_val_firstrow;

        /* ii. Sort rows using the ordering 0, 1, -1 */
        LMC0_sortrows ( al, rowperm, colperm, rowsignperm, colsignperm, nrows, ncols );

        // compare the current pair of columns
        r = lmc0_compare_columns ( al, rowperm, colperm, column, rowsignperm, colsignperm );

        //printf("LMC0_columns: column comparison at %d: %d\n", column, r);
        if ( r==LMC_EQUAL ) {
            // columns are equal, so go one column deeper
            // keep column permutation, and column sign permutation
            if ( verbose>=2 )
                printf ( "LMC0_columns: go to col %d\n", column+1 );
            r = LMC0_columns ( al, rowperm, colperm, column+1, rowsignperm, colsignperm );
        }
        if ( r==LMC_LESS ) {
            // we already know the array is not in minimal form
            break;
        }

        /* result is LMC_MORE, continue with the calculation */
        // restore column sign permutation for the first row
        colsignperm[ colperm[column] ] = current_sign_col;
        // swap back column. Restore column permutation
        colperm[c]=col;
        colperm[column]=colperm[col];
        // restore column sign permutation for the first row
        colsignperm[ colperm[column] ] = current_sign_col;
    }
    return r;
}



/** Algorithm: check an array X for LMC0 form
 *
 * The tag [B] indicates branching
 *
 *
 *
 * 1. Select the first column [B]
 * 2. Apply row level permutations such that the first column only contains ones (and a single zero)
 * 3. Sort the rows of the design
 * 4. Select one of two possible sign permutations for the first row [B]
 * 5. Select the next column k [B]
 * 6. With the selected column:
 *    i. Apply the correct column level permutation to make the element X(1, k) equal to 1
 *    ii. Sort the rows using the ordering 0, 1, -1
 *    iii. Compare the column to the original column k of the design. If less, then abort with LMC_LESS, if more, continue to 5, if equal then go to 5 with k+1
 * 7. Return LMC_MORE
 *
 *
 * When selecting columns we use a std::vector to keep track of selected indices
 * When sorting the rows to do not sort the rows, but use a rowsort_t structure to keep track of the order of rows
 **/

/**
 * Dummy algorithm: only allow column permutations
 *
 *
 */
lmc_t LMC0check ( const array_link &al ) {
    /*0. Initialize data */
    // & implies a pointer, thus al is a pointer
    // the * is used to access to the variable sotred at the pointer
    lmc_t result = LMC_MORE;
    // Get size of the array
    const int ncols = al.n_columns; // number of columns, rowindex_t
    const int nrows = al.n_rows; // number of rows, colindex_t

    // Initialize the column permutation
    std::vector<int> colperm ( ncols ); // create pointer for the column permutations
    init_perm ( colperm );
    //initialize column level permutation
    std::vector<int> colsignperm ( ncols );
    init_signperm ( colsignperm );
    //initialize row-level permutation
    std::vector<int> rowsignperm ( nrows );
    init_signperm ( rowsignperm );
    // Initialize row order
    dyndata_t rowperm_data = dyndata_t ( nrows );//rowperm_data.show();
    rowsort_t *rowsort = rowperm_data.rowsort;
    for (rowindex_t i = 0; i < nrows; i++){
        rowsort[i].val = i;
    } //print_rowsort(rowsort, nrows);

    for (int sel_col = 0; sel_col < ncols; sel_col++){
        /*1. Select the first (sel_col) column */
        // make this column the first in the array
        int first_col = colperm[ 0 ]; // copy current first column index
        colperm[ 0 ] = colperm[ sel_col ]; // assign to be the first one
        colperm[ sel_col ] = first_col; // finish the swap

        /*2. Find row-level permutation such that the first column only contains ones */
        rowlevel_permutation ( al, rowsort, colperm, rowsignperm, colsignperm, nrows, 0 ); //print_perm(rowsignperm);
        /* 3. Find permutation to sort the array*/
        LMC0_sortrows( al, rowsort, colperm, rowsignperm, colsignperm, nrows, ncols ); //print_rowsort(rowsort, nrows);
        int value_rowsign_firstrow = rowsignperm[ rowsort[0].val ];
        /* 4. Select one of two possible sign permutations for the first row */
        for (int r_sign = 0; r_sign < 2; r_sign++){

            // replace the element k by the possible sign permutation
            rowsignperm[ rowsort[0].val ] = 2*r_sign - 1;

            /* 5. Select the next column */
            lmc_t result = LMC0_columns (al, rowperm, colperm, 0, rowsignperm, colsignperm, ncols, nrows );
            if ( result==LMC_LESS ) {
            /* FINISH TEST */
            // we already know the array is not in minimal form
                break; // FINISH TEST

            } // end if
            /* RESTORE  values for colsignperm, colperm, and sort rows*/
        } // end for
        // restore value of the first_row
        rowsignperm[ rowsort[0].val ] = value_rowsign_firstrow;

        // swap back column. Restore column permutation
        colperm[ sel_col ] = colperm[ 0 ];
        colperm[ 0 ] = first_col;

    }
    // call the recursive part starting at column 0
    //lmc_t r = LMC0_columns ( al, colperm, 0 );

    return result;
}



int main ( int argc, char* argv[] ) {
    AnyOption opt;
    /* parse command line options */
    opt.setFlag ( "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt.setOption ( "output", 'o' );
    opt.setOption ( "input", 'I' );
    opt.setOption ( "rand", 'r' );
    opt.setOption ( "verbose", 'v' );
    opt.setOption ( "ii", 'i' );
    opt.setOption ( "jj" );
    opt.setOption ( "rows" );
    opt.setOption ( "cols" );
    opt.setOption ( "oaconfig", 'c' ); /* file that specifies the design */

    opt.addUsage ( "Orthonal Array: oaconfcheck: testing platform" );
    opt.addUsage ( "Usage: oaconfcheck [OPTIONS] [FILE]" );
    opt.addUsage ( "" );
    opt.addUsage ( " -h --help  			Prints this help " );
    opt.processCommandArgs ( argc, argv );


    double t0=get_time_ms(), dt=0;
    int randvalseed = opt.getIntValue ( 'r', 1 );
    int ix = opt.getIntValue ( 'i', 1 );
    int jj = opt.getIntValue ( "jj", 5 );

    int xx = opt.getIntValue ( 'x', 0 );
    int niter  = opt.getIntValue ( "niter", 10 );
    int verbose  = opt.getIntValue ( "verbose", 1 );

    char *input = opt.getValue ( 'I' );
    if ( input==0 )
        input="test.oa";

    srand ( randvalseed );
    if ( randvalseed==-1 ) {
        randvalseed=time ( NULL );
        printf ( "random seed %d\n", randvalseed );
        srand ( randvalseed );
    }


    print_copyright();

    /* parse options */
    if ( opt.getFlag ( "help" ) || opt.getFlag ( 'h' ) || opt.getArgc() <0 ) {
        opt.printUsage();
        exit ( 0 );
    }

    setloglevel ( SYSTEM );

    /* test on known set of arrays */
    /*
    printf("### test code on known array\n");
    array_link al = exampleArray ( 2 );
    lmc_t r=LMC0check ( al );
    printf ( "minimal form check result: %d (should be %d)\n", r, LMC_MORE );

    for ( int i=0; i<3; i++ ) {
        printf ( "round %d: apply random transformation to the array\n", i );
        array_link altmp =  al.randomcolperm();
        lmc_t r=LMC0check ( altmp );

        if ( altmp==al ) {
            printf ( "  minimal form check: %d (should be %d)\n", r, LMC_MORE );
        } else        printf ( "  minimal form check: %d (should be %d)\n", r, LMC_LESS );

    }
    */

    printf("### read from file\n");

    arraylist_t ll= readarrayfile ( input );

    for ( size_t i=0; i<ll.size(); i++ ) {
        array_link al = ll[i];

        //al = al.randomrowperm();
        //al = al.randomcolperm();
        //cout << al; print array
        lmc_t r = LMC0check ( al );
        printf ( "array %d: result %d\n", (int) i, ( int ) r );
    }


    return 0;

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;

