#include <iostream>
#include <math.h>
#include <algorithm>

#include "nonroot.h"
#include "lmc.h"
#include "extend.h"


// IDEA: fold non_root and non_root_j4 into each other: only difference: TPLUSONECOLUMN?

inline void cpy_dyndata_rowsort ( const dyndata_t *src, dyndata_t *dest )
{
#ifdef OADEBUG
    if ( src->N!=dest->N ) {
        myprintf ( "error: source and target dyndata_t have unequal number of rows!\n" );
        exit ( 0 );
    }
#endif

    memcpy ( dest->rowsort, src->rowsort, sizeof ( rowsort_t ) *dest->N );
}

inline void cpy_dyndata_rowsortl ( const dyndata_t *src, dyndata_t *dest )
{
#ifdef OADEBUG
    if ( src->N!=dest->N ) {
        myprintf ( "error: source and target dyndata_t have unequal number of rows!\n" );
        exit ( 0 );
    }
    if ( src->rowsortl==0 ) {
        myprintf ( "cpy_dyndata_rowsortl: ERROR: source %ld\n", ( long ) src->rowsortl );
        exit ( 0 );
    }
    if ( dest->rowsortl==0 ) {
        myprintf ( "cpy_dyndata_rowsortl: ERROR: target %ld\n", ( long ) dest->rowsortl );
        exit ( 0 );
        return;
    }
#endif

    memcpy ( dest->rowsortl, src->rowsortl, sizeof ( rowsortlight_t ) *dest->N );
}

void debug_check_rowsort_overflow ( const arraydata_t *ad, const rowsort_t *rowsort, const dyndata_t* dd, rowindex_t cur_row )
{
#ifdef OAOVERFLOW
    assert ( rowsort[cur_row].val >= 0 );
#endif
#ifdef OADEBUG
    /* the rowsort value should not exceed the maximum int value! */
    if ( std::numeric_limits<rowsort_value_t>::max() /ad->s[dd->col]<2*rowsort[cur_row].val )
        log_print ( SYSTEM, " LMC_check_col: possible integer overflow detected! s[col]=%ld, rowsort[cur_row].val=%ld\n", ad->s[dd->col], ( long ) rowsort[cur_row].val );
#endif
}


/* internal functions */
//inline lmc_t LMC_check_col_tplusplus(const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd);
/// Perform LMC check on a single column (t+1), also apply a level transformation
inline lmc_t LMC_check_col_tplus ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd );
/// Perform LMC check on a single column, also apply a level transformation
lmc_t LMC_check_col ( const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd );
/// Perform LMC check on a single column, special case for 2-level arrays
lmc_t LMC_check_col_ft_2level ( const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose=0 );
//lmc_t __attribute__ ((noinline))   LMC_check_col_ft ( const array_t *originalcol, carray_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose=0 );

/// compare 2 columns
inline lmc_t LMC_check_col ( const array_t *original, const array_t *array, const arraydata_t *ad, const dyndata_t *dd );


/// Perform LMC check on a single column, perform full sorting
inline lmc_t LMC_check_col_complete ( const array_t *original, const array_t *array, const arraydata_t *ad, const dyndata_t *dd );
/// Perform LMC check on a single column, a buffer is used for the rowsort value
inline lmc_t LMC_check_col_buffer ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, rowsort_value_t *rsbuffer );

inline lmc_t LMC_check_col_less ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd );

// note: this one passes the complete array, with column permutation data!
inline lmc_t LMC_check_col_j5order ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd );

/**
* @brief Perform LMC check on a single column (t+1), also apply a level transformation
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col_tplus ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd )
{
    assert ( ad->s[0]==2 );
    lmc_t ret = LMC_EQUAL;
    int cur_row=0, rowp;
    const int oaindex = ad->oaindex;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
    prefetch ( rowsort );
#endif

    /* we check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        /* sort rows according to new information */
        for ( int k = 0; k < oaindex; k++ ) {
            cur_row = j*oaindex+k; // OPTIMIZE: use cur_row++

            rowp = rowsort[cur_row].r;
            rowsort[cur_row].val = lperm[array[rowp]];


        }

        // this is the special condition for column (t+1)
        // we only need to check the ordering for the first N/oaindex rows. the remaining rows can be brought to normal form if the LMC check is satisfied
        // see the appendix A in the article by Deng and Tang
        if ( j<1 ) {
            // OPTIMIZE: combine sorting and comparing
            // by checking whether we have flipped any elements
            // other sort functions: insertionSort, bubbleSort, bubbleSort2, quickSort, shellSort, std::sort, std::__insertion_sort

            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );

            for ( int k = 0; k < oaindex; k++ ) {
                cur_row = j*oaindex+k;
                rowp = rowsort[cur_row].r;

#ifdef SAFELPERM
                if ( original[cur_row ] < safe_lperm<array_t> ( array[rowp], lperm, ad->s[dd->col] ) ) {
                    ret = LMC_MORE;
                    break;
                } else if ( original[cur_row] > safe_lperm<array_t> ( array[rowp], lperm, ad->s[dd->col] ) ) { //the permuted array is lex. less than original
                    ret = LMC_LESS;
                    break;
                }
#else
                if ( original[cur_row ] < lperm[array[rowp]] ) {
                    ret = LMC_MORE;	// first check for LMC_MORE, this happens much more often then LMC_LESS
                    break;
                } else if ( original[cur_row] > lperm[array[rowp]] ) { //the permuted array is lex. less than original
                    ret = LMC_LESS;
                    //myprintf("ret: %d, cur_row %d= j(%d)*oaindex(%d)+k(%d);", ret, cur_row,j,oaindex,k);
                    //myprintf(" val %d    val %d\n",  original[cur_row ] ,  lperm[array[rowp]] );
                    break;
                }
#endif
            }
            if ( ret!=LMC_EQUAL ) {
                break;
            }
        }
    }

#ifdef OAANALYZE_DISCR
    analyse_discriminant ( cur_row, dd->col, ret, ad->N, ad->ncols );
#endif
    return ret;
}

/**
* @brief Perform LMC check on a single column, also apply a level transformation
*
* Immediately return once LMC_MORE is detected
*
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col_less ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd )
{
    //myprintfd("LMC_check_col: start: lperm "); print_perm(lperm, ad->s[0]);
    lmc_t ret = LMC_EQUAL;
    int cur_row=0, rowp;
    const int oaindex = ad->oaindex;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;

    /* we check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        if ( oaindex>1 ) {
            /* sort rows according to new information */
            for ( int k = 0; k < oaindex; k++ ) {
                cur_row = j*oaindex+k; // OPTIMIZE: use cur_row++
                rowp = rowsort[cur_row].r;

                debug_check_rowsort_overflow ( ad, rowsort, dd, cur_row );

                /* OPTIMIZE: the multiplication can be taken out of this function */
                rowsort[cur_row].val = ad->s[dd->col]*rowsort[cur_row].val+lperm[array[rowp]];
            }

            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        }

        if ( ret==LMC_EQUAL ) {
            for ( int k = 0; k < oaindex; k++ ) {
                cur_row = j*oaindex+k;
                rowp = rowsort[cur_row].r;

                if ( original[cur_row ] < lperm[array[rowp]] ) {
                    ret = LMC_MORE;
                    break;
                } else if ( original[cur_row] > lperm[array[rowp]] ) { //the permuted array is lex. less than original
                    ret = LMC_LESS;
                    //printfd("LMC_check_col: LMC_LESS: lperm "); print_perm(lperm, ad->s[0]);
                    break;
                }
            }
        }
        if ( ret==LMC_MORE )
            break;
    }

    //myprintf("LMC_check_col: at col %d, ret %d (LMC_EQUAL 1, LMC_LESS 0)\n", dd->col,ret); printf("sorted column: \n");
    //print_col_sorted(array, lperm, rowsort, ad->N);

#ifdef OAANALYZE_DISCR
    analyse_discriminant ( cur_row, dd->col, ret, ad->N, ad->ncols );
#endif
    return ret;
}

/**
* @brief Perform LMC check on a single column, also apply a level transformation, no sorting
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col_nosort ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd )
{
    lmc_t ret = LMC_EQUAL;
    int rowp;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
    prefetch ( rowsort );
#endif

    /* we check in blocks of oaindex */
    for ( int cur_row = 0; cur_row < ( nrows ); cur_row++ ) {
        rowp = rowsort[cur_row].r;
        array_t v = lperm[array[rowp]];
        if ( original[cur_row ] < v ) {
            ret = LMC_MORE;
            break;
        } else if ( original[cur_row] > v ) { //the permuted array is lex. less than original
            ret = LMC_LESS;
            break;
        }
    }

    return ret;
}
/**
* @brief Perform LMC check on a single column
*
* Also apply a level transformation. A buffer is used for the rowsort value
*
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col_buffer ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, rowsort_value_t *rsbuffer )
{
    //myprintf("LMC_check_col: start\n");
    lmc_t ret = LMC_EQUAL;
    int cur_row, rowp;
    const int oaindex = ad->oaindex;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
    prefetch ( rowsort );
#endif

    /* we check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        if ( oaindex>1 ) {
            /* sort rows according to new information */
            for ( int k = 0; k < oaindex; k++ ) {
                cur_row = j*oaindex+k;

                rowp = rowsort[cur_row].r;

#ifdef OAOVERFLOW
                assert ( rowsort[cur_row].val>0 );
                /* the rowsort value should not exceed the maximum int value! */
                if ( std::numeric_limits<rowsort_value_t>::max() /ad->s[dd->col]<2*rowsort[cur_row].val )
                    log_print ( SYSTEM, " LMC_check_col: integer overflow detected! s[col]=%ld, rowsort[cur_row].val=%ld\n", ad->s[dd->col], ( long ) rowsort[cur_row].val );
#endif /* OAOVERFLOW */

                /* OPTIMIZE: the multiplication can be taken out of this function */
                rowsort[cur_row].val = rsbuffer[cur_row]+lperm[array[rowp]];
            }
            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        }

        for ( int k = 0; k < oaindex; k++ ) {
            cur_row = j*oaindex+k;
            rowp = rowsort[cur_row].r;

            if ( original[cur_row ] < lperm[array[rowp]] ) {
                ret = LMC_MORE;
                break;
            } else if ( original[cur_row] > lperm[array[rowp]] ) { //the permuted array is lex. less than original
                ret = LMC_LESS;
                break;
            }
        }
        if ( ret!=LMC_EQUAL )
            break;
    }
    return ret;
}

// check range of elements (inclusive)
inline lmc_t checkLMChelperSorted ( int ix, int iy, const array_t *arraycol,  levelperm_t lperm, const arraydata_t *ad )
{
    //myprintf("checkLMChelper: %d to %d (inclusive)\n", ix, iy);
    array_t prevval=0; //lperm[0];

    for ( int cur_row = ix; cur_row <= iy; cur_row++ ) {
        array_t cval=lperm[arraycol[cur_row]];
        if ( prevval==1 && cval==0 )
            return LMC_LESS;

        prevval=cval;
    }

    return LMC_EQUAL;
}

// check range of elements (inclusive)
inline lmc_t checkLMChelperSorted ( int ix, int iy, const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad,  rowsortlight_t *rowsort )
{
    int v1=0, v2=0;

    array_t prevval=0;
    for ( int cur_row = ix; cur_row <= iy; cur_row++ ) {
        int rowp = rowsort[cur_row];
        array_t cval=originalcol[cur_row ];	// current value

        // count check
        v1 +=originalcol[cur_row ];
        v2 += lperm[arraycol[rowp]];

        // sort check
        if ( prevval==1 && cval==0 ) {
            if ( 0 ) {
                myprintf ( "checkLMChelperSorted line %d: %d to %d: ", __LINE__, ix, iy );
                for ( int v=ix; v<=iy; v++ ) {
                    int rowp = rowsort[v];
                    myprintf ( "%d ", lperm[arraycol[rowp]] );
                }
                myprintf ( "\n" );
            }
            return LMC_LESS;
        }
        prevval=cval;
    }
    // if (lperm[0]) v2=(1+iy-ix)-v2;

    if ( v1==v2 )
        return LMC_EQUAL;
    if ( v1<v2 )
        return LMC_MORE;
    return LMC_LESS;
}

// check range of elements (inclusive)
inline lmc_t checkLMChelper ( int ix, int iy, const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad,  rowsortlight_t *rowsortl )
{
    //myprintf("checkLMChelper: %d to %d (inclusive)\n", ix, iy);
    int v1=0, v2=0;

#ifdef OADEBUG
    myassert ( rowsortl!=0, "checkLMChelper: need rowsortl structure\n" );
#endif
    for ( int cur_row = ix; cur_row <= iy; cur_row++ ) {
        int rowp = rowsortl[cur_row];

        v1 +=originalcol[cur_row ];
        v2 += lperm[arraycol[rowp]];
        // v2 += arraycol[rowp];
    }

    // if (lperm[0]) v2=(1+iy-ix)-v2;

    if ( v1==v2 )
        return LMC_EQUAL;
    if ( v1<v2 )
        return LMC_MORE;
    return LMC_LESS;
}

// check range of elements (inclusive)
inline lmc_t checkLMChelper ( int ix, int iy, const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad,  rowsort_t *rowsort )
{
    //myprintf("checkLMChelper: %d to %d (inclusive)\n", ix, iy);
    int v1=0, v2=0;

    for ( int cur_row = ix; cur_row <= iy; cur_row++ ) {
        int rowp = rowsort[cur_row].r;

        v1 +=originalcol[cur_row ];
        v2 += lperm[arraycol[rowp]];
        // v2 += arraycol[rowp2];
    }

    // if (lperm[0]) v2=(1+iy-ix)-v2;

    if ( v1==v2 )
        return LMC_EQUAL;
    if ( v1<v2 )
        return LMC_MORE;
    return LMC_LESS;
}

// sort for 2-level arrays (low and high are inclusive)
inline void sortzerooneR ( rowsortlight_t *rs, int low, int high, const array_t *arraycol )
{
    while ( low < high ) {
        if ( arraycol[rs[low]] == 1 ) {
            low ++;
            continue;
        }
        while ( arraycol[rs[high]] == 0 ) {
            high --;
        }
        if ( low < high ) {
            //swap arr[low], arr[high]
            std::swap ( rs[low], rs[high] );
        }
    }
}
// sort for 2-level arrays (low and high are inclusive)
inline void sortzeroone ( rowsortlight_t *rs, int low, int high, const array_t *arraycol )
{
    while ( low < high ) {
        if ( arraycol[rs[low]] == 0 ) {
            low ++;
            continue;
        }
        while ( arraycol[rs[high]] == 1 ) {
            high --;
        }
        if ( low < high ) {
            //swap arr[low], arr[high]
            std::swap ( rs[low], rs[high] );
        }
    }
}

// sort for 2-level arrays (low and high are inclusive)
inline void sortzerooneR ( rowsort_t *rs, int low, int high, const array_t *arraycol )
{
    while ( low < high ) {
        if ( arraycol[rs[low].r] == 1 ) {
            low ++;
            continue;
        }
        while ( arraycol[rs[high].r] == 0 ) {
            high --;
        }
        if ( low < high ) {
            //swap arr[low], arr[high]
            std::swap ( rs[low].r, rs[high].r );
        }
    }
}
// sort for 2-level arrays (low and high are inclusive)
inline void sortzeroone ( rowsort_t *rs, int low, int high, carray_t *arraycol )
{
    while ( low < high ) {
        if ( arraycol[rs[low].r] == 0 ) {
            low ++;
            continue;
        }
        while ( arraycol[rs[high].r] == 1 ) {
            high --;
        }
        if ( low < high ) {
            //swap arr[low], arr[high]
            std::swap ( rs[low].r, rs[high].r );
        }
    }
}

/// check column for row symmetry exchanges
lmc_t LMC_check_col_rowsymm ( const array_t *arraycol,  const arraydata_t *ad, const symmdata &sd, int col, int dverbose )
{
    lmc_t ret = LMC_EQUAL;
    const int scol = col-1;

    int nb = sd.ft.atfast ( sd.ft.n_rows-1, scol ); // myprintf("LMC_check_col_ft: nb %d\n",nb);
    array_t * sdp=sd.ft.array+scol*sd.ft.n_rows;

    array_t lperm[2] = { ( array_t ) 0,1};

    /* we check in blocks determined by the ft */
    for ( int j = 0; j < nb; j++ ) {
        int x1 = sdp[2*j];
        int x2 = sdp[2*j+1];

        ret = checkLMChelperSorted ( x1, x2-1, arraycol, lperm, ad );

        if ( ret!=LMC_EQUAL ) {
            return ret;
        }
    }
    return ret;
}

/// Perform LMC check on a single column, special case for 2-level arrays
lmc_t LMC_check_col_ft_2level_rowsymm ( const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose )
{
    //myassert(sd!=0, "LMC_check_col_ft");

    lmc_t ret = LMC_EQUAL;
    const int scol = dd->col-1;
    rowsortlight_t *rowsort = dd->rowsortl;

    int nb = sd.ft.atfast ( sd.ft.n_rows-1, scol ); // myprintf("LMC_check_col_ft: nb %d\n",nb);
    array_t * sdp=sd.ft.array+scol*sd.ft.n_rows;

    /* we check in blocks determined by the ft */
    for ( int j = 0; j < nb; j++ ) {
        int x1 = sdp[2*j];
        int x2 = sdp[2*j+1];

        ret = checkLMChelperSorted ( x1, x2-1, originalcol, arraycol, lperm, ad, rowsort );

        if ( ret!=LMC_EQUAL ) {
            return ret;
        }
    }


    // new method
    if ( lperm[0]==0 ) {
        for ( int j = 0; j < nb; j++ ) {
            int x1 = sdp[2*j];
            int x2 = sdp[2*j+1];
            sortzeroone ( rowsort, x1, x2-1, arraycol );
        }
    } else {
        for ( int j = 0; j < nb; j++ ) {
            int x1 = sdp[2*j];
            int x2 = sdp[2*j+1];
            if ( x1==x2 && 0 ) {
                myprintf ( "\n### \n" );
                sd.ft.showarray();
                myprintf ( "calling sortzerooneR: scol %d, j %d, x1 %d, x2 %d\n", scol, j, x1, x2 );
            }
            sortzerooneR ( rowsort, x1, x2-1, arraycol );
        }
    }

    return ret;

}


/** check 2-level column in fast mode
 *
 * We assume the sd pointer has been set and that only comparison is needed. The rowsort_t structure is not updated.
 * We also assume the column is sorted with respect to the previous columns, so trivial cases are not detected!
 */
lmc_t LMC_check_col_ft_2level ( const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose )
{
#ifdef OADEBUG
    if ( dd->rowsortl==0 ) {
        printfd ( "OADEBUG: error: dd->rowsortl==%ld\n", ( long ) dd->rowsortl ) ;
        lperm[-1]=-1; // trick to get feedback from valgrind
    }
    myassert ( dd->rowsortl!=0, "LMC_check_col_ft_2level: need rowsortl structure\n" );
#endif

    //myassert(sd!=0, "LMC_check_col_ft");

    lmc_t ret = LMC_EQUAL;
    const int scol = dd->col-1;
    rowsortlight_t *rowsort = dd->rowsortl;

    int nb = sd.ft.atfast ( sd.ft.n_rows-1, scol ); // myprintf("LMC_check_col_ft: nb %d\n",nb);
    array_t * sdp=sd.ft.array+scol*sd.ft.n_rows;

    /* we check in blocks determined by the ft */
    for ( int j = 0; j < nb; j++ ) {
        int x1 = sdp[2*j];
        int x2 = sdp[2*j+1];

        ret = checkLMChelper ( x1, x2-1, originalcol, arraycol, lperm, ad, rowsort );

        if ( ret!=LMC_EQUAL ) {
            return ret;
        }
    }


    // new method
    if ( lperm[0]==0 ) {
        for ( int j = 0; j < nb; j++ ) {
            int x1 = sdp[2*j];
            int x2 = sdp[2*j+1];
            sortzeroone ( rowsort, x1, x2-1, arraycol );
        }
    } else {
        for ( int j = 0; j < nb; j++ ) {
            int x1 = sdp[2*j];
            int x2 = sdp[2*j+1];
            if ( x1==x2 && 0 ) {
                myprintf ( "\n### \n" );
                sd.ft.showarray();
                myprintf ( "calling sortzerooneR: scol %d, j %d, x1 %d, x2 %d\n", scol, j, x1, x2 );
            }
            sortzerooneR ( rowsort, x1, x2-1, arraycol );
        }
    }

    return ret;
}

lmc_t LMC_check_col_ft_testing ( const array_t *originalcol, carray_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose )
{
    //myassert(sd!=0, "LMC_check_col_ft");

    lmc_t ret = LMC_EQUAL;
    const int nrows = ad->N;
    const int scol = dd->col-1;
    rowsort_t *rowsort = dd->rowsort;

    int nb = sd.ft.atfast ( sd.ft.n_rows-1, scol ); // myprintf("LMC_check_col_ft: nb %d\n",nb);
    array_t * sdp=sd.ft.array+scol*sd.ft.n_rows;

    /* we check in blocks determined by the ft */
    for ( int j = 0; j < nb; j++ ) {
        //int x1 = sd.ft.atfast(2*j, scol); int x2 = sd.ft.atfast(2*j+1, scol);
        int x1 = sdp[2*j];
        int x2 = sdp[2*j+1];
        //if (col>4) 	    myprintf("  calling checkLMChelper: %d to %d (inclusive)\n", x1, x2-1);

        ret = checkLMChelper ( x1, x2-1, originalcol, arraycol, lperm, ad, rowsort );
        if ( ret!=LMC_EQUAL ) {
            return ret;
        }
    }

    // we have LMC_EQUAL: prepare data for deeper levels
    if ( 0 ) {
        const rowsort_value_t sval = ad->s[dd->col];

        for ( int cur_row=0; cur_row<ad->N; cur_row++ ) {
            int rowp = rowsort[cur_row].r;
            rowsort[cur_row].val = lperm[arraycol[rowp]];
            //rowsort[cur_row].val = sval*rowsort[cur_row].val+lperm[arraycol[rowp]];
        }
    }

    if ( 0 ) {
        for ( int j = 0; j < nb; j++ ) {
            int x1 = sdp[2*j];
            int x2 = sdp[2*j+1];
            if ( ( x2-x1 ) <=1 )
                continue;
            //int x1 = sd.ft.atfast(2*j, scol);  int x2 = sd.ft.atfast(2*j+1, scol);
            for ( int cur_row=x1; cur_row<x2; cur_row++ ) {
                rowsort[cur_row].val = lperm[arraycol[rowsort[cur_row].r]];
            }
            oacolSort ( rowsort, x1, x2-1 );
        }
    } else {
        // new method

        if ( lperm[0]==0 ) {
            for ( int j = 0; j < nb; j++ ) {
                int x1 = sdp[2*j];
                int x2 = sdp[2*j+1];
                sortzeroone ( rowsort, x1, x2-1, arraycol );
            }
        } else {
            for ( int j = 0; j < nb; j++ ) {
                int x1 = sdp[2*j];
                int x2 = sdp[2*j+1];
                if ( x1==x2 && 0 ) {
                    myprintf ( "\n### \n" );
                    sd.ft.showarray();
                    myprintf ( "calling sortzerooneR: scol %d, j %d, x1 %d, x2 %d\n", scol, j, x1, x2 );
                }
                sortzerooneR ( rowsort, x1, x2-1, arraycol );
            }
        }
    }

    return ret;
}



inline lmc_t LMC_check_col_ft_X ( const array_t *originalcol, carray_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd, const symmdata *sd, int dverbose )
{
    const int newalg=1;

    lmc_t ret = LMC_EQUAL;
    int cur_row, rowp;
    const int oaindex = ad->oaindex;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;
    const rowsort_value_t sval = ad->s[dd->col];

//lmc_t rete=LMC_EQUAL; int fh=0;

    /* we check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        /* sort rows according to new information */
        rowsort_value_t prevvalue=rowsort[j*oaindex].val;
        int previndex=0;
        for ( int k = 0; k < oaindex; k++ ) {
            cur_row = j*oaindex+k; // OPTIMIZE: use cur_row++
            rowp = rowsort[cur_row].r;

            //if (dverbose) {
//		myprintf("##\nk %d: prevvalue %d, rowsort[cur_row].val %d\n", k, prevvalue, rowsort[cur_row].val);
//	      }

            if ( rowsort[cur_row].val!=prevvalue ) {
                if ( newalg ) {
                    //  if ( ( (k-1)-previndex)>-2) {
                    lmc_t v = checkLMChelper ( j*oaindex+previndex, j*oaindex+k-1, originalcol, arraycol, lperm, ad, rowsort );
                    if ( v!=LMC_EQUAL ) {
                        ret=v;
                        return ret;
                    }
                    //  }
                } else {
                    oacolSort ( rowsort+ ( j*oaindex ), previndex, k-1 );

                    for ( int kk = previndex; kk < k; kk++ ) {
                        int cur_row2 = j*oaindex+kk;
                        int rowp2 = rowsort[cur_row2].r;

                        if ( originalcol[cur_row2 ] != lperm[arraycol[rowp2]] ) {

                            if ( originalcol[cur_row2 ] < lperm[arraycol[rowp2]] ) {
                                if ( dverbose )  {
                                    myprintf ( "early abort MORE: j %d, k %d, previndex %d, kk %d\n", j, k, previndex, kk );
                                    myprintf ( "  %d %d (arraycol[rowp2] %d %d)\n", originalcol[cur_row2 ] , lperm[arraycol[rowp2]], arraycol[rowp2], rowp2 );
                                }
                                // continue;
                                //if (!fh) {   rete = LMC_MORE;  fh=1;  }
                                ret=LMC_MORE;

                                return ret;
                                break;
                            } else if ( originalcol[cur_row2] > lperm[arraycol[rowp2]] ) {
                                //continue;
                                if ( dverbose ) {
                                    myprintf ( "early abort: LESS j %d, k %d, previndex %d, kk %d\n", j, k, previndex, kk );
                                }
                                // if (!fh) { fh=1;  rete = LMC_LESS;  }
                                ret = LMC_LESS;
                                return ret;
                                break;
                            }
                        }

                    }
                }
                previndex=k;
                prevvalue=rowsort[cur_row].val;
            }

            // if (dverbose)
            //   myprintf("k: %d, rowsort.val %d->%d (val array %d)\n", k, rowsort[cur_row].val, sval*rowsort[cur_row].val+lperm[arraycol[rowp]], lperm[arraycol[rowp]] );
            if ( !newalg ) {
                rowsort[cur_row].val = sval*rowsort[cur_row].val+lperm[arraycol[rowp]];
                if ( ret!=LMC_EQUAL )
                    break;
            }
        }
        // if (ret!=LMC_EQUAL) break;

        if ( newalg ) {
            lmc_t v = checkLMChelper ( j*oaindex+previndex, j*oaindex+oaindex-1, originalcol, arraycol, lperm, ad, rowsort );
            if ( v!=LMC_EQUAL ) {
                ret=v;
                return ret;
            }
            //oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        } else {
            //oacolSort ( rowsort+ ( j*oaindex ), previndex, oaindex-1 );
            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );

            //myprintf("previndex %d, oaindex %d\n", previndex, oaindex);
            int k;
            for ( k = previndex; k < oaindex; k++ ) {
                cur_row = j*oaindex+k;
                rowp = rowsort[cur_row].r;

                // OPTIMIZE: use != method for other functions as well
                if ( originalcol[cur_row ] != lperm[arraycol[rowp]] ) {
                    if ( originalcol[cur_row ] < lperm[arraycol[rowp]] ) {
                        ret = LMC_MORE;
                        break;
                    } else if ( originalcol[cur_row] > lperm[arraycol[rowp]] ) { //the permuted array is lex. less than original
                        ret = LMC_LESS;
                        break;
                    }
                }
            }
            if ( ret!=LMC_EQUAL ) {
                if ( dverbose ) {
                    myprintf ( "abort: j %d, k %d, row %d, ret %d\n", j, k, j*oaindex+k, ret );
                }
                break;
            }
        }
    }

    if ( newalg ) {
        // we have LMC_EQUAL: prepare data for deeper levels
        // TODO: this can be blocked (so we do not need the rowsort values any more!!
        for ( int i=0; i<ad->N; i++ ) {
            cur_row=i;
            int rowp = rowsort[i].r;
            rowsort[i].val = sval*rowsort[i].val+lperm[arraycol[rowp]];
        }
        for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        }
    }

    return ret;
}

// OPTIMIZE: LMC_check_col: if no row symmetries, we can forget about the rowsort

/**
* @brief Perform LMC check on a single column, also apply a level transformation
*
* Immediately return once LMC_LESS or LMC_MORE is detected
*
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col ( const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd )
{
    //myprintf("LMC_check_col: start\n");
    lmc_t ret = LMC_EQUAL;
    int cur_row, rowp;
#ifdef OAANALYZE_DISCR
    cur_row=0;
#endif
    const int oaindex = ad->oaindex;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
    prefetch ( rowsort );
#endif

#ifdef OADEBUG
    for ( int xx=0; xx<nrows; xx++ ) {
        int rowp = rowsort[xx].r;
        if ( rowp>=nrows ) {
            printfd ( "error: rowp %d >= nrows %d\n", rowp, nrows );
            exit ( 0 );
        }
    }
#endif

    const rowsort_value_t sval = ad->s[dd->col];

    /* we check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        if ( oaindex>1 ) {
            /* sort rows according to new information */
            for ( int k = 0; k < oaindex; k++ ) {
                cur_row = j*oaindex+k; // OPTIMIZE: use cur_row++
                rowp = rowsort[cur_row].r;

#ifdef OADEBUG
                debug_check_rowsort_overflow ( ad, rowsort, dd, cur_row );
#endif

                /* OPTIMIZE: the multiplication can be taken out of this function */
#ifdef SAFELPERM
                rowsort[cur_row].val = ad->s[dd->col]*rowsort[cur_row].val+safe_lperm<array_t> ( array[rowp], lperm, ad->s[dd->col] );
#else
                rowsort[cur_row].val = sval*rowsort[cur_row].val+lperm[arraycol[rowp]];
#endif
            }

            // OPTIMIZE: combine sorting and comparing
            // other sort functions: insertionSort, bubbleSort2, quickSort, shellSort, std::sort, std::__insertion_sort
            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        }

#ifdef OADEBUG
        for ( int xx=0; xx<nrows; xx++ ) {
            int rowp = rowsort[xx].r;
            if ( rowp>=nrows ) {
                printfd ( "error: row %d: rowp %d >= nrows %d\n", xx, rowp, nrows );
                exit ( 0 );
            }
        }
#endif

        for ( int k = 0; k < oaindex; k++ ) {
            cur_row = j*oaindex+k;
            rowp = rowsort[cur_row].r;

#ifdef SAFELPERM
            if ( original[cur_row ] < safe_lperm<array_t> ( array[rowp], lperm, ad->s[dd->col] ) ) {
                ret = LMC_MORE;	// first check for LMC_MORE, this happens much more often
                break;
            } else if ( original[cur_row] > safe_lperm<array_t> ( array[rowp], lperm, ad->s[dd->col] ) ) { //the permuted array is lex. less than original
                ret = LMC_LESS;
                break;
            }
#else
            if ( originalcol[cur_row ] < lperm[arraycol[rowp]] ) {
                ret = LMC_MORE;
                break;
            } else if ( originalcol[cur_row] > lperm[arraycol[rowp]] ) { //the permuted array is lex. less than original
                ret = LMC_LESS;
                break;
            }
#endif
        }

        if ( ret!=LMC_EQUAL )
            break;
    }

    //printf("LMC_check_col: at col %d, ret %d (LMC_EQUAL 1, LMC_LESS 0)\n", dd->col,ret); printf("sorted column: \n");
    //print_col_sorted(array, lperm, rowsort, ad->N);

#ifdef OAANALYZE_DISCR
    analyse_discriminant ( cur_row, dd->col, ret, ad->N, ad->ncols );
#endif
    return ret;
}




/**
* @brief Perform LMC check on a single column
* @param original Pointer to original OA
* @param array Array
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
lmc_t LMC_check_col ( const array_t *original, const array_t *array, const arraydata_t *ad, const dyndata_t *dd )
{
    //myprintf("LMC_check_col: start\n");
    lmc_t ret = LMC_EQUAL;
    int cur_row, rowp;
    int oaindex = ad->oaindex;
    int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;

    /* we check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        if ( oaindex>1 ) {
            /* sort rows according to new information */
            for ( int k = 0; k < oaindex; k++ ) {
                cur_row = j*oaindex+k;

                rowp = rowsort[cur_row].r;
#ifdef OADEBUG
                assert ( rowsort[cur_row].val>0 );
                /* check for integer overflow */
                if ( std::numeric_limits<rowsort_value_t>::max() /ad->s[dd->col]<2*rowsort[cur_row].val )
                    log_print ( SYSTEM, " LMC_check_col: integer overflow detected! s[col]=%ld, rowsort[cur_row].val=%ld\n", ad->s[dd->col], ( long ) rowsort[cur_row].val );
#endif
                /* OPTIMIZE: the multiplication can be taken out of the function */
                rowsort[cur_row].val = ad->s[dd->col]*rowsort[cur_row].val+array[rowp];
            }

            // OPTIMIZE: select best sort algorithm
            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        }

        for ( int k = 0; k < oaindex; k++ ) {
            cur_row = j*oaindex+k;
            rowp = rowsort[cur_row].r;

            if ( original[cur_row] > array[rowp] ) { //the permuted array is lex. less than original
                ret = LMC_LESS;
                break;
            } else if ( original[cur_row ] < array[rowp] ) {
                ret = LMC_MORE;
                break;
            }
        }
        if ( ret!=LMC_EQUAL )
            break;
    }
    return ret;
}


/**
* @brief Perform LMC check on a single column, perform full sorting
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col_complete ( const array_t *original, carray_t *array, const arraydata_t *ad, const dyndata_t *dd )
{
    lmc_t ret = LMC_EQUAL;
    const int oaindex = ad->oaindex;
    const int nrows = ad->N;
    rowsort_t *rowsort = dd->rowsort;

    /* we sort check in blocks of oaindex */
    for ( int j = 0; j < ( nrows/oaindex ); j++ ) {
        if ( oaindex>1 ) {
            /* sort rows according to new information */
            for ( int k = 0; k < oaindex; k++ ) {
                int cur_row = j*oaindex+k;
                int rowp = rowsort[cur_row].r;

                /* OPTIMIZE: the multiplication can be taken out of this function */
                rowsort[cur_row].val = ad->s[dd->col]*rowsort[cur_row].val+array[rowp];
            }
            oacolSort ( rowsort+ ( j*oaindex ), 0, oaindex-1 );
        }
    }

    for ( rowindex_t cur_row = 0; cur_row < nrows; cur_row++ ) {
        int rowp = rowsort[cur_row].r;

        if ( original[cur_row ] < array[rowp] ) {
            ret = LMC_MORE;
            break;
        } else if ( original[cur_row] > array[rowp] ) {
            ret = LMC_LESS;
            break;
        }
        if ( ret!=LMC_EQUAL )
            break;
    }

    return ret;
}

/** Check column with ordering based on j5 */
inline lmc_t LMC_check_col_j5order ( const array_t *original, const array_t *array, levelperm_t lperm, const arraydata_t *ad, const dyndata_t *dd )
{
//   myprintf("LMC_check_col_j5order: check code path...\n");

    // OPTIMIZE: pre-calculate J5 values!
    // OPTIMIZE: use first m rows to calculate j5?
    assert ( ad->order==ORDER_J5 );
    //myprintf("LMC_check_col_j5order!\n");

    const int cpoffset = ad->N*dd->colperm[dd->col];
    if ( dd->col<4 ) { // XXX
        lmc_t ret = LMC_check_col ( original+dd->col*+ad->N, array+cpoffset, lperm, ad, dd );
        return ret;
    }

    if ( log_print ( DEBUG,"" ) ) {
        log_print ( NORMAL, "LMC_check_col_j5order col %d, colperm: ", dd->col );
        print_perm ( dd->colperm, ad->ncols );
    }

    //reduction->array+dyndata->col*+ad->N, original+cpoffset


    //int jbase = abs(predictJ(original, ad->N, 0));
    //int jcol = abs(predictJrowsort(array, ad->N, dd->rowsort, lperm));
    colperm_t pp = new_perm_init<colindex_t> ( 5 );
    if ( 0 ) {
        myprintf ( "original:\n" );
        print_array ( original, ad->N, ad->ncols );
    }
    pp[4]=dd->col;

    int jbase= abs ( jvaluefast ( original, ad->N, 5, pp ) );
    for ( int x=0; x<5; x++ )
        pp[x]=dd->colperm[x];
    pp[4]=dd->colperm[dd->col];
    int jcol= abs ( jvaluefast ( array, ad->N, 5, pp ) );	// NOTE: tricky lperms have not been done!
    if ( 1 && dd->col>4 ) {
        //myprintf("array:\n");
        //print_array(array, ad->N, ad->ncols);
        if ( log_print ( DEBUG,"" ) ) {
            myprintf ( "  xxx col %d, jbase %d, jcol %d: colperm: ", dd->col, jbase, jcol );
            print_perm ( pp, 5 );
        }

    }
    delete_perm ( pp );


    lmc_t ret;
    if ( jbase==jcol ) {
        ret = LMC_check_col ( original+dd->col*+ad->N, array+cpoffset, lperm, ad, dd );
        return ret;
    } else {
        if ( log_print ( DEBUG,"" ) ) {
            myprintf ( "  j5order: col %d, jbase %d, jcol %d: colperm: ", dd->col, jbase, jcol );
            print_perm ( pp, 5 );
        }

        if ( jbase ORDER_J5_SMALLER jcol )
            ret = LMC_MORE;
        else
            ret = LMC_LESS;
        return ret;
    }
}

/* Main functions */

lmc_t LMCreduce_non_root_j4 ( const array_t * original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic )
{
    //logstream ( DEBUG+1 ) << printfstring ( "LMCreduce_non_root: col %d, reduction->mode %d, perm ",  dyndata->col, reduction->mode );
    // print_perm ( dyndata->colperm, ad->ncols );


    /* static allocation of data structures for non-root stage of LMC check */
    levelperm_t *lperm_p = 0;
    colperm_t *colperm_p = 0;
    colperm_t *localcolperm_p = 0;

    /* static init of dyndata */
    dyndata_t ** dynd_p;
    int dynd_p_nelem;

    /* static allocation of buffer for the column check */
    array_t *colbuffer = 0;

    tmpStatic.init_nonroot_stage ( lperm_p, colperm_p, localcolperm_p, dynd_p, dynd_p_nelem, colbuffer, ad );

    const int remsize = ad->ncols-dyndata->col;	/* number of columns remaining */
    const int nlevels = ad->s[dyndata->col];	/* number of levels at this column */

    levelperm_t lperm = lperm_p[dyndata->col];
    colperm_t colpermloop = localcolperm_p[dyndata->col];
    dyndata_t *dyndatacpy = dynd_p[dyndata->col];
    const int cg = ad->get_col_group ( dyndata->col );	/* current column group */
    /* number of columns remaining in this group */
    const int ncolsremgroup = ad->colgroupsize[cg]+ad->colgroupindex[cg]-dyndata->col;
    //myprintf("LMCreduce_non_root:   ad->colgroupsize[cg] %d, ad->colgroupindex[cg] %d, \n", ad->colgroupsize[cg], ad->colgroupindex[cg]);

#ifdef OADEBUG
    myassertdebug ( dyndata->col!=ad->ncols, "LMC_non_root this code should not be reached! dyndata->col!=ad->ncols\n" );
#endif

    // OPTIMIZE: check whether changing order of this loop makes a difference
    lmc_t ret = LMC_MORE;
    init_perm<colindex_t> ( colpermloop, remsize );
    const int col = dyndata->col;
    for ( int i=0; i<ncolsremgroup; i++ ) {
        std::swap ( colpermloop[0], colpermloop[i] );
        //	 myprintf("\n\n## colperm loop i=%d (col %d):\n", i, dyndata->col);

        // keep track of applied column permutation
        copy_perm ( dyndata->colperm, dyndatacpy->colperm, ad->ncols );
        perform_inv_perm<colindex_t> ( dyndata->colperm+col, dyndatacpy->colperm+col, remsize, colpermloop );

        const int nlevelperms = init_perm_n<array_t, int> ( lperm, nlevels );
        const int cpoffset = ad->N*dyndatacpy->colperm[dyndata->col];


        /* loop over all level permutations */
        for ( int j=0; j<nlevelperms; j++ ) {
            /* copy dyndata rowsort data */
            // OPTIMIZE: eliminate this copy by moving it into check_col, saves up to 15% for large run sizes
            cpy_dyndata_rowsort ( dyndata, dyndatacpy );
            dyndatacpy->col = dyndata->col;	// TODO: needed?

            /* LMC_check_col performs level permutations, updates the sorting structure and compares with the original array on blocks of oaindex */
            // prefetch(reduction->array);

            if ( reduction->mode>=OA_REDUCE ) {
                /* for reduce check perform full column perm!!!! */ /* NOTE: why? */
#ifdef SAFELPERM
                safe_perform_level_perm<array_t> ( original+cpoffset, colbuffer, ad->N, lperm, ( int ) ad->s[dyndata->col] );
#else
                perform_level_perm ( original+cpoffset, colbuffer, ad->N, lperm );
#endif
                //ret = LMC_check_col_complete ( reduction->array+dyndata->col*ad->N, colbuffer, ad, dyndatacpy );
                ret = LMC_check_col_less ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );

            } else {

                if ( ad->order==ORDER_LEX ) {
                    // TODO: clean this up
#ifdef TPLUSCOLUMN
                    //ret =  LMC_check_col ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );
                    ret =  LMC_check_col_tplus ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );
#else
                    ret =  LMC_check_col ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );
#endif
                } else {
                    //myprintf("j5 order\n");
                    ret =  LMC_check_col_j5order ( reduction->array, original, lperm, ad, dyndatacpy );
                }

            }

            if ( ret==LMC_LESS ) {
                reduction->lastcol=col;
            }

            if ( reduction->mode>=OA_REDUCE && ret==LMC_LESS ) {
                reduction->updateLastCol ( col );
                reduction->updateFromLoop ( *ad, dyndatacpy, lperm_p, original );
                if ( ( reduction->nred%150000 ) ==0 )
                    log_print ( DEBUG-1, "found reduction: column %d, number is %ld\n", dyndata->col, reduction->nred );


                ret=LMC_EQUAL;	/* since the best array so far is updated, we need to continue */
            }
            if ( ret==LMC_EQUAL ) {
                // this column could not decide, go one level deeper

                dyndatacpy->col += 1;
                if ( dyndatacpy->col==ad->ncols ) {
                    /* this was the last column */
                    ret = LMC_MORE;
                } else {
                    if ( reduction->maxdepth==dyndata->col )
                        ret = LMC_MORE;
                    else
                        /* pass the old cpy as an argument, since it will not be changed */
                        ret = LMCreduce_non_root ( original, ad, dyndatacpy, reduction, oaextend, tmpStatic );
                }
            }
            if ( ret==LMC_LESS ) {
                break;
            }
            // else: LMC_MORE or LMC_EQUAL, continue with loop
            next_perm ( lperm, nlevels );
        }

        if ( ret==LMC_LESS )
            break;

        std::swap ( colpermloop[0], colpermloop[i] ); // swap back
    }

    return ret;
}


/** @perform The non-root stage of the LMC check or LMC reduction.
*
* @param original
* @param ad
* @param dyndata
* @param reduction
* @return
*/
lmc_t LMCreduce_non_root ( const array_t * original, const arraydata_t* ad, dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, const LMC_static_struct_t &tmpStatic )
{
    /* static allocation of data structures for non-root stage of LMC check */
    levelperm_t *lperm_p = 0;
    colperm_t *colperm_p = 0;
    colperm_t *localcolperm_p = 0;

    /* static init of dyndata */
    dyndata_t ** dynd_p;
    int dynd_p_nelem;

    /* static allocation of buffer for the column check */
    array_t *colbuffer = 0;

    tmpStatic.init_nonroot_stage ( lperm_p, colperm_p, localcolperm_p, dynd_p, dynd_p_nelem, colbuffer, ad );

    const int remsize = ad->ncols-dyndata->col;	/* number of columns remaining */
    const int nlevels = ad->s[dyndata->col];	/* number of levels at this column */

    levelperm_t lperm = lperm_p[dyndata->col];
    colperm_t colpermloop = localcolperm_p[dyndata->col];
    dyndata_t *dyndatacpy = dynd_p[dyndata->col];
    const int cg = ad->get_col_group ( dyndata->col );	/* current column group */
    /* number of columns remaining in this group */
    const int ncolsremgroup = ad->colgroupsize[cg]+ad->colgroupindex[cg]-dyndata->col;

#ifdef OADEBUG
    myassertdebug ( dyndata->col!=ad->ncols, "LMC_non_root this code should not be reached!  dyndata->col!=ad->ncols\n" );
#endif

    // OPTIMIZE: use buffer for rowsort values to take multiplication out of LMC_check_col
    // OPTIMIZE: check whether changing order of this loop makes a difference
    lmc_t ret = LMC_MORE;
    init_perm<colindex_t> ( colpermloop, remsize );	// intialize for entire loop
    const int col = dyndata->col;
    for ( int i=0; i<ncolsremgroup; i++ ) {
        std::swap ( colpermloop[0], colpermloop[i] );	// swap 2 columns, we swap them back at the end of the loop so we do not need to initialize each time
        //	 printf("\n\n## colperm loop i=%d (col %d):\n", i, dyndata->col);

        // keep track of applied column permutation
        copy_perm ( dyndata->colperm, dyndatacpy->colperm, ad->ncols );
        perform_inv_perm<colindex_t> ( dyndata->colperm+col, dyndatacpy->colperm+col, remsize, colpermloop );

        // OPTIMIZE: this yields too many level perms, we can reduce it by pre-selecting the first row
        // OPTIMIZE: even better: if there are no row symmetries any more, we can just select the best permutation
        const int nlevelperms = init_perm_n<array_t, int> ( lperm, nlevels );
        const int cpoffset = ad->N*dyndatacpy->colperm[dyndata->col];

        /* loop over all level permutations */
        for ( int j=0; j<nlevelperms; j++ ) {
            /* copy dyndata rowsort data */
            // OPTIMIZE: eliminate this copy by moving it into check_col, saves up to 15% for large run sizes
            cpy_dyndata_rowsort ( dyndata, dyndatacpy );
            dyndatacpy->col = dyndata->col;  // TODO: needed?

            /* LMC_check_col performs level permutations, updates the sorting structure and compares with the original array on blocks of oaindex */
            // prefetch(reduction->array);

            if ( reduction->mode>=OA_REDUCE ) {
                /* for reduce check perform full column perm!!!! */ /* NOTE: why? */
                // OPTIMIZE: in case of LMC_MORE we do not need the full column
                // NOTE: VERY important for random LMC reductions!

                if ( 1 ) {
                    ret = LMC_check_col_less ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );

                    if(ret==LMC_LESS) {
                        //printfd("less! "); dyndatacpy->show();
                    }
                } else {
#ifdef SAFELPERM
                    safe_perform_level_perm<array_t> ( original+cpoffset, colbuffer, ad->N, lperm, ( int ) ad->s[dyndata->col] );
#else
                    perform_level_perm ( original+cpoffset, colbuffer, ad->N, lperm );
#endif
                    //printfd("...\n");

                    ret = LMC_check_col_complete ( reduction->array+dyndata->col*ad->N, colbuffer, ad, dyndatacpy );
                }
            } else {
                if ( ad->order==ORDER_LEX ) {

                    if ( ( oaextend.getAlgorithm() ==MODE_LMC_2LEVEL || oaextend.getAlgorithm() ==MODE_J5ORDERXFAST ) && reduction->sd!=0 ) { // || oaextend.getAlgorithm()==MODE_J4) {
                        myassert ( reduction->sd!=0, "LMC_check_col_ft" );
                        ret =  LMC_check_col_ft_2level ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy, * ( reduction->sd ) );
//                     			  cpy_dyndata_rowsort ( dyndata, dyndatacpy );
                        //    ret =  LMC_check_col ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );
                    } else {
                        ret =  LMC_check_col ( reduction->array+dyndata->col*+ad->N, original+cpoffset, lperm, ad, dyndatacpy );
                    }
                    dyndatacpy->col = dyndata->col;  // TODO: needed?

                } else 	{
                    ret =  LMC_check_col_j5order ( reduction->array, original, lperm, ad, dyndatacpy );
                }
            }

            //myprintf("LMCreduce_non_root: column permutation %d, level permutation %d, ret %d\n", i,j, ret);

            if ( ret==LMC_LESS ) {
                //logstream ( DEBUG ) << printfstring ( "  LMC_LESS: column %d\n", col );
                reduction->lastcol=col;
            }

            if ( reduction->mode>=OA_REDUCE && ret==LMC_LESS ) {
                reduction->lastcol=col;
                reduction->updateFromLoop ( *ad, dyndatacpy, lperm_p, original );
                if ( ( reduction->nred%150000 ) ==0 )
                    log_print ( DEBUG-1, "found reduction: column %d, number is %ld\n", dyndata->col, reduction->nred );


                ret=LMC_EQUAL;	/* since the best array so far is updated, we need to continue */

                reduction->mincol=std::min ( reduction->mincol, col );

                if ( 0 && reduction->mode==OA_REDUCE_PARTIAL ) {
                    myprintf ( "debug: reduction->mode %d\n", reduction->mode );
                    if ( col<reduction->targetcol ) {
                        myprintf ( "  LMC_REDUCE_PARTIAL: abort loop!\n" );
                        ret=LMC_LESS;
                    }
                }
            }
            if ( ret==LMC_EQUAL ) {
                //myprintf("equal colgroup: "); print_perm( dyndatacpy->colperm, dyndatacpy->col+1 ); // @PTE
                reduction->symms.storeSymmetryPermutation ( dyndatacpy );

                //reduction->storeColumnPermutation(dyndatacpy->colperm, dyndatacpy->col+1 );
                // this column could not decide, go one level deeper

                dyndatacpy->col += 1;
                if ( dyndatacpy->col==ad->ncols ) {
                    /* this was the last column */
                    ret = LMC_MORE;
                } else {
                    if ( reduction->maxdepth==dyndata->col )
                        ret = LMC_MORE;
                    else {
                        /* pass the old cpy as an argument, since it will not be changed */
                        ret = LMCreduce_non_root ( original, ad, dyndatacpy, reduction, oaextend, tmpStatic );
                    }
                }
            }
            //}
            if ( ret==LMC_LESS ) {
                //logstream(DEBUG) << printfstring("  LMC_LESS: column %d\n", col);
                break; // break from levelperms loop
            }
            // else: LMC_MORE or LMC_EQUAL, continue with loop
            next_perm ( lperm, nlevels );
        }

        if ( ret==LMC_LESS ) {
            break; // break from loop over available columns
        }


        std::swap ( colpermloop[0], colpermloop[i] );	// swap 2 columns back
    }

    return ret;
}


/** @perform The non-root stage of the LMC check or LMC reduction.
*
* @param original
* @param ad
* @param dyndata
* @param reduction
* @return
*/
lmc_t LMCreduce_non_root_2level ( const array_t * original, const arraydata_t* ad, dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, const LMC_static_struct_t &tmpStatic )
{
    const int dverbose=0;

    if ( dverbose )
        myprintf ( "LMCreduce_non_root_2level: level %d\n", dyndata->col );

    /* static allocation of data structures for non-root stage of LMC check */
    levelperm_t *lperm_p = 0;
    colperm_t *colperm_p = 0;
    colperm_t *localcolperm_p = 0;

    /* static init of dyndata */
    dyndata_t ** dynd_p;
    int dynd_p_nelem;

    /* static allocation of buffer for the column check */
    array_t *colbuffer = 0;

    tmpStatic.init_nonroot_stage ( lperm_p, colperm_p, localcolperm_p, dynd_p, dynd_p_nelem, colbuffer, ad );

#ifdef OADEBUG
    if ( dynd_p==0 ) {
        printfd ( "dyndata->col %d\n", dyndata->col );
        tmpStatic.ad->show();
        printfd ( "error: dynd_p not initialized!  tmpStatic.dyndata_p %ld \n", ( long ) tmpStatic.dyndata_p );
    }

#endif

    const int remsize = ad->ncols-dyndata->col;	/* number of columns remaining */
    const int nlevels = 2; // ad->s[dyndata->col];	/* number of levels at this column */

    levelperm_t lperm = lperm_p[dyndata->col];
    //colperm_t colpermloopXX = localcolperm_p[dyndata->col];
    dyndata_t *dyndatacpy = dynd_p[dyndata->col];
    dyndatacpy->allocrowsortl();
    dyndatacpy->col = dyndata->col;  // TODO: needed?

    const int cg = ad->get_col_group ( dyndata->col );	/* current column group */
    /* number of columns remaining in this group */
    const int ncolsremgroup = ad->colgroupsize[cg]+ad->colgroupindex[cg]-dyndata->col;
    //myprintf("LMCreduce_non_root:   ad->colgroupsize[cg] %d, ad->colgroupindex[cg] %d, \n", ad->colgroupsize[cg], ad->colgroupindex[cg]);

    //debug_print_transformations(reduction, original, *array, ad, dyndata);

    lmc_t ret = LMC_MORE;
    //init_perm<colindex_t> ( colpermloopXX, remsize );	// intialize for entire loop
    const int col = dyndata->col;
    copy_perm ( dyndata->colperm, dyndatacpy->colperm, ad->ncols );

    const int nlevelperms = 2;	//special case for 2-level arrays
    for ( int i=0; i<ncolsremgroup; i++ ) {
        //if (i==0) printfd(" col %d: ncolsremgroup %d\n", dyndata->col, ncolsremgroup);

        // keep track of applied column permutation
        std::swap ( dyndatacpy->colperm[col], dyndatacpy->colperm[col+i] );

        lperm[0]=0;
        lperm[1]=1; // special init for 2-level arrays
        const int cpoffset = ad->N*dyndatacpy->colperm[dyndata->col];
        const array_t *reductionarrayoffset = reduction->array+dyndata->col*ad->N;

        //bool specialf=false;
        /* loop over all level permutations */
        for ( int j=0; j<nlevelperms; j++ ) {
            /* copy dyndata rowsort data */
            cpy_dyndata_rowsortl ( dyndata, dyndatacpy );  //NOTE: not for col_ft!?

            dyndatacpy->col = dyndata->col;  // TODO: needed?

            /* LMC_check_col performs level permutations, updates the sorting structure and compares with the original array on blocks of oaindex */

            if ( reduction->mode>=OA_REDUCE ) {
                myprintf ( "LMCreduce_non_root_2level: mode OA_REDUCE not implemented\n" );
                return LMC_NONSENSE;
            } else {
                if ( ad->order==ORDER_LEX ) {

#ifdef OADEBUG
                    myassert ( reduction->sd!=0, "LMC_check_col_ft: reduction->sd is not set" );
#endif
                    if ( dyndatacpy->col==ad->ncols-1 ) {
                        ret =  LMC_check_col_ft_2level_rowsymm ( reductionarrayoffset, original+cpoffset, lperm, ad, dyndatacpy, * ( reduction->sd ) , 0 );
                        //  myprintf("LMCreduce_non_root_2level: special check for final column: col %d (ncols %d) return %d!!!!\n", dyndatacpy->col, ad->ncols, ret);

                    } else {

//#define OADEBUGL
                        ret =  LMC_check_col_ft_2level ( reductionarrayoffset, original+cpoffset, lperm, ad, dyndatacpy, * ( reduction->sd ) , 0 );
                    }
                    dyndatacpy->col = dyndata->col;  // TODO: needed?



                } else 	{
                    ret =  LMC_check_col_j5order ( reduction->array, original, lperm, ad, dyndatacpy );
                }
            }

            //myprintf("LMCreduce_non_root: column permutation %d, level permutation %d, ret %d\n", i,j, ret);

            if ( ret==LMC_LESS ) {
                reduction->lastcol=col;
            }

            if ( ret==LMC_EQUAL ) {
                //myprintf("equal colgroup: "); print_perm( dyndatacpy->colperm, dyndatacpy->col+1 ); // @PTE
                reduction->symms.storeSymmetryPermutation ( dyndatacpy );

                // this column could not decide, go one level deeper
                dyndatacpy->col += 1;
                if ( dyndatacpy->col==ad->ncols ) {
                    /* this was the last column */
                    ret = LMC_MORE;
                } else {
                    if ( reduction->maxdepth==dyndata->col )
                        ret = LMC_MORE;
                    else {
                        /* pass the old cpy as an argument, since it will not be changed */
                        ret = LMCreduce_non_root_2level ( original, ad, dyndatacpy, reduction, oaextend, tmpStatic );
                    }
                }
                /*
                // special case for check of 2-level arrays: if the first permutation is LMC_EQUAL, then the second permutation can only be LMC_MORE
                if(ret==LMC_MORE && 0) {
                myprintf("early brake: col %d\n", dyndata->col);
                	 dyndatacpy->rowsortl2rowsort();
                	 reduction->updateTransformation( *ad, dyndatacpy, lperm_p, original );
                	 array_link tmp(original, ad->N, ad->ncols);
                	 array_link ww = reduction->transformation->apply(tmp);
                	 ww.selectFirstColumns(dyndata->col+1).showarray();
                	}
                	if (specialf) {
                	 array_link ww(original, ad->N, ad->ncols);
                	 myprintf("at col %d: level perm leads to same value\n", dyndata->col);
                	 ww.selectFirstColumns(dyndata->col+2).showarray();

                	}
                	 if(ret==LMC_MORE) {
                	  specialf=1;
                	 }
                           //  break;
                             */
            }
            //}
            if ( ret==LMC_LESS ) {
                break; // break from levelperms loop
            }
            // else: LMC_MORE or LMC_EQUAL, continue with loop
            next_perm_twoperm ( lperm, nlevels );
        }

        if ( ret==LMC_LESS ) {
            break; // break from loop over available columns
        }

        std::swap ( dyndatacpy->colperm[col], dyndatacpy->colperm[col+i] );
    }

    return ret;
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
