#include <algorithm>
#include <iostream>
#include <math.h>

#include "extend.h"
#include "lmc.h"
#include "nonroot.h"

inline void cpy_dyndata_rowsort (const dyndata_t *src, dyndata_t *dest) {
#ifdef OADEBUG
        if (src->N != dest->N) {
                myprintf ("error: source and target dyndata_t have unequal number of rows!\n");
                throw_runtime_exception ("error: source and target dyndata_t have unequal number of rows!");
        }
#endif

        memcpy (dest->rowsort, src->rowsort, sizeof (rowsort_t) * dest->N);
}

/// copy rowsortl variable of dyndata_t
inline void cpy_dyndata_rowsortl (const dyndata_t *src, dyndata_t *dest) {
#ifdef OADEBUG
        if (src->N != dest->N) {
                myprintf ("error: source and target dyndata_t have unequal number of rows!\n");
                throw_runtime_exception("");
        }
        if (src->rowsortl == 0) {
                myprintf ("cpy_dyndata_rowsortl: ERROR: source %ld\n", (long)src->rowsortl);
                throw_runtime_exception("");
        }
        if (dest->rowsortl == 0) {
                myprintf ("cpy_dyndata_rowsortl: ERROR: target %ld\n", (long)dest->rowsortl);
                throw_runtime_exception("");
        }
#endif

        memcpy (dest->rowsortl, src->rowsortl, sizeof (rowsortlight_t) * dest->N);
}

void debug_check_rowsort_overflow (const arraydata_t *ad, const rowsort_t *rowsort, const dyndata_t *dd,
                                   rowindex_t cur_row) {
#ifdef OADEBUG
        /* the rowsort value should not exceed the maximum int value! */
        if (std::numeric_limits< rowsort_value_t >::max () / ad->s[dd->col] < 2 * rowsort[cur_row].val)
                log_print (
                    SYSTEM,
                    " LMC_check_col: possible integer overflow detected! s[col]=%ld, rowsort[cur_row].val=%ld\n",
                    ad->s[dd->col], (long)rowsort[cur_row].val);
#endif
}

/// Perform LMC check on a single column (t+1), also apply a level transformation
inline lmc_t LMC_check_col_tplus (const array_t *original, const array_t *array, levelperm_t lperm,
                                  const arraydata_t *ad, const dyndata_t *dd);
/// Perform LMC check on a single column, also apply a level transformation
lmc_t LMC_check_col (const array_t *originalcol, const array_t *arraycol, levelperm_t lperm, const arraydata_t *ad,
                     const dyndata_t *dd);
/// Perform LMC check on a single column, special case for 2-level arrays
lmc_t LMC_check_col_ft_2level (const array_t *originalcol, const array_t *arraycol, levelperm_t lperm,
                               const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose = 0);

/**
* @brief Perform LMC check on a single column (t+1), also apply a level transformation
* @param original Pointer to original OA
* @param array Array
* @param lperm Level permutation to be applied before comparison
* @param ad Static array data
* @param dd Dynamic array data
* @return
*/
inline lmc_t LMC_check_col_tplus (const array_t *original, const array_t *array, levelperm_t lperm,
                                  const arraydata_t *ad, const dyndata_t *dd) {
        assert (ad->s[0] == 2);
        lmc_t ret = LMC_EQUAL;
        int cur_row = 0, rowp;
        const int oaindex = ad->oaindex;
        const int nrows = ad->N;
        rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
        prefetch (rowsort);
#endif

        /* we check in blocks of oaindex */
        for (int j = 0; j < (nrows / oaindex); j++) {
                /* sort rows according to new information */
                for (int k = 0; k < oaindex; k++) {
                        cur_row = j * oaindex + k; // OPTIMIZE: use cur_row++

                        rowp = rowsort[cur_row].r;
                        rowsort[cur_row].val = lperm[array[rowp]];
                }

                // this is the special condition for column (t+1)
                // we only need to check the ordering for the first N/oaindex rows. the remaining rows can be brought
                // to normal form if the LMC check is satisfied
                // see the appendix A in the article by Deng and Tang
                if (j < 1) {
                        oacolSort (rowsort + (j * oaindex), 0, oaindex - 1);

                        for (int k = 0; k < oaindex; k++) {
                                cur_row = j * oaindex + k;
                                rowp = rowsort[cur_row].r;

#ifdef SAFELPERM
                                if (original[cur_row] < safe_lperm< array_t > (array[rowp], lperm, ad->s[dd->col])) {
                                        ret = LMC_MORE;
                                        break;
                                } else if (original[cur_row] >
                                           safe_lperm< array_t > (
                                               array[rowp], lperm,
                                               ad->s[dd->col])) { // the permuted array is lex. less than original
                                        ret = LMC_LESS;
                                        break;
                                }
#else
                                if (original[cur_row] < lperm[array[rowp]]) {
                                        ret = LMC_MORE; // first check for LMC_MORE, this happens much more often then
                                                        // LMC_LESS
                                        break;
                                } else if (original[cur_row] >
                                           lperm[array[rowp]]) { // the permuted array is lex. less than original
                                        ret = LMC_LESS;
                                        break;
                                }
#endif
                        }
                        if (ret != LMC_EQUAL) {
                                break;
                        }
                }
        }


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
inline lmc_t LMC_check_col_less (const array_t *original, const array_t *array, levelperm_t lperm,
                                 const arraydata_t *ad, const dyndata_t *dd) {
        lmc_t ret = LMC_EQUAL;
        int cur_row = 0, rowp;
        const int oaindex = ad->oaindex;
        const int nrows = ad->N;
        rowsort_t *rowsort = dd->rowsort;

        /* we check in blocks of oaindex */
        for (int j = 0; j < (nrows / oaindex); j++) {
                if (oaindex > 1) {
                        /* sort rows according to new information */
                        for (int k = 0; k < oaindex; k++) {
                                cur_row = j * oaindex + k; // OPTIMIZE: use cur_row++
                                rowp = rowsort[cur_row].r;

                                debug_check_rowsort_overflow (ad, rowsort, dd, cur_row);

                                rowsort[cur_row].val = ad->s[dd->col] * rowsort[cur_row].val + lperm[array[rowp]];
                        }

                        oacolSort (rowsort + (j * oaindex), 0, oaindex - 1);
                }

                if (ret == LMC_EQUAL) {
                        for (int k = 0; k < oaindex; k++) {
                                cur_row = j * oaindex + k;
                                rowp = rowsort[cur_row].r;

                                if (original[cur_row] < lperm[array[rowp]]) {
                                        ret = LMC_MORE;
                                        break;
                                } else if (original[cur_row] >
                                           lperm[array[rowp]]) { // the permuted array is lex. less than original
                                        ret = LMC_LESS;
                                        break;
                                }
                        }
                }
                if (ret == LMC_MORE)
                        break;
        }


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
inline lmc_t LMC_check_col_nosort (const array_t *original, const array_t *array, levelperm_t lperm,
                                   const arraydata_t *ad, const dyndata_t *dd) {
        lmc_t ret = LMC_EQUAL;
        int rowp;
        const int nrows = ad->N;
        rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
        prefetch (rowsort);
#endif

        /* we check in blocks of oaindex */
        for (int cur_row = 0; cur_row < (nrows); cur_row++) {
                rowp = rowsort[cur_row].r;
                array_t v = lperm[array[rowp]];
                if (original[cur_row] < v) {
                        ret = LMC_MORE;
                        break;
                } else if (original[cur_row] > v) { // the permuted array is lex. less than original
                        ret = LMC_LESS;
                        break;
                }
        }

        return ret;
}

// check range of elements (inclusive)
inline lmc_t checkLMChelperSorted (int ix, int iy, const array_t *originalcol, const array_t *arraycol,
                                   levelperm_t lperm, const arraydata_t *ad, rowsortlight_t *rowsort) {
        int v1 = 0, v2 = 0;

        array_t prevval = 0;
        for (int cur_row = ix; cur_row <= iy; cur_row++) {
                int rowp = rowsort[cur_row];
                array_t cval = originalcol[cur_row]; // current value

                // count check
                v1 += originalcol[cur_row];
                v2 += lperm[arraycol[rowp]];

                // sort check
                if (prevval == 1 && cval == 0) {
                        return LMC_LESS;
                }
                prevval = cval;
        }

        if (v1 == v2)
                return LMC_EQUAL;
        if (v1 < v2)
                return LMC_MORE;
        return LMC_LESS;
}

// check range of elements (inclusive)
inline lmc_t checkLMChelper (int ix, int iy, const array_t *originalcol, const array_t *arraycol, levelperm_t lperm,
                             const arraydata_t *ad, rowsortlight_t *rowsortl) {
        int v1 = 0, v2 = 0;

        for (int cur_row = ix; cur_row <= iy; cur_row++) {
                int rowp = rowsortl[cur_row];

                v1 += originalcol[cur_row];
                v2 += lperm[arraycol[rowp]];
        }

        if (v1 == v2)
                return LMC_EQUAL;
        if (v1 < v2)
                return LMC_MORE;
        return LMC_LESS;
}

// sort for 2-level arrays (low and high are inclusive)
inline void sortzerooneR (rowsortlight_t *rs, int low, int high, const array_t *arraycol) {
        while (low < high) {
                if (arraycol[rs[low]] == 1) {
                        low++;
                        continue;
                }
                while (arraycol[rs[high]] == 0) {
                        high--;
                }
                if (low < high) {
                        std::swap (rs[low], rs[high]);
                }
        }
}
// sort for 2-level arrays (low and high are inclusive)
inline void sortzeroone (rowsortlight_t *rs, int low, int high, const array_t *arraycol) {
        while (low < high) {
                if (arraycol[rs[low]] == 0) {
                        low++;
                        continue;
                }
                while (arraycol[rs[high]] == 1) {
                        high--;
                }
                if (low < high) {
                        std::swap (rs[low], rs[high]);
                }
        }
}

// sort for 2-level arrays (low and high are inclusive)
inline void sortzerooneR (rowsort_t *rs, int low, int high, const array_t *arraycol) {
        while (low < high) {
                if (arraycol[rs[low].r] == 1) {
                        low++;
                        continue;
                }
                while (arraycol[rs[high].r] == 0) {
                        high--;
                }
                if (low < high) {
                        std::swap (rs[low].r, rs[high].r);
                }
        }
}
// sort for 2-level arrays (low and high are inclusive)
inline void sortzeroone (rowsort_t *rs, int low, int high, carray_t *arraycol) {
        while (low < high) {
                if (arraycol[rs[low].r] == 0) {
                        low++;
                        continue;
                }
                while (arraycol[rs[high].r] == 1) {
                        high--;
                }
                if (low < high) {
                        std::swap (rs[low].r, rs[high].r);
                }
        }
}

/// sort a column using rowsort structure and symmetry data
void _sort_array_column(rowsortlight_t *rowsort, int nb, const levelperm_t lperm, array_t* symmdata_column_pointer, const array_t * arraycol)
{
	if (lperm[0] == 0) {
		for (int j = 0; j < nb; j++) {
			int x1 = symmdata_column_pointer[2 * j];
			int x2 = symmdata_column_pointer[2 * j + 1];
			sortzeroone(rowsort, x1, x2 - 1, arraycol);
		}
	}
	else {
		for (int j = 0; j < nb; j++) {
			int x1 = symmdata_column_pointer[2 * j];
			int x2 = symmdata_column_pointer[2 * j + 1];
			sortzerooneR(rowsort, x1, x2 - 1, arraycol);
		}
	}
}

/// Perform LMC check on a single column, special case for 2-level arrays
lmc_t LMC_check_col_ft_2level_rowsymm (const array_t *originalcol, const array_t *arraycol, levelperm_t lperm,
                                       const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose) {
        lmc_t ret = LMC_EQUAL;
        const int scol = dd->col - 1;
        rowsortlight_t *rowsort = dd->rowsortl;

        int nb = sd.ft.atfast (sd.ft.n_rows - 1, scol);
        array_t *symmdata_column_pointer = sd.ft.array + scol * sd.ft.n_rows;

        /* we check in blocks determined by the ft */
        for (int j = 0; j < nb; j++) {
                int x1 = symmdata_column_pointer[2 * j];
                int x2 = symmdata_column_pointer[2 * j + 1];

                ret = checkLMChelperSorted (x1, x2 - 1, originalcol, arraycol, lperm, ad, rowsort);

                if (ret != LMC_EQUAL) {
                        return ret;
                }
        }

		_sort_array_column(rowsort, nb, lperm, symmdata_column_pointer, arraycol);

		return ret;
}

/** check 2-level column in fast mode
 *
 * We assume the sd pointer has been set and that only comparison is needed. The rowsort_t structure is not updated.
 * We also assume the column is sorted with respect to the previous columns, so trivial cases are not detected!
 */
lmc_t LMC_check_col_ft_2level (const array_t *originalcol, const array_t *arraycol, levelperm_t lperm,
                               const arraydata_t *ad, const dyndata_t *dd, const symmdata &sd, int dverbose) {
        myassert (dd->rowsortl != 0, "LMC_check_col_ft_2level: need rowsortl structure\n");

        lmc_t ret = LMC_EQUAL;
        const int scol = dd->col - 1;
        rowsortlight_t *rowsort = dd->rowsortl;

        int nb = sd.ft.atfast (sd.ft.n_rows - 1, scol); 
        array_t *symmdata_column_pointer = sd.ft.array + scol * sd.ft.n_rows;

        /* we check in blocks determined by the ft */
        for (int j = 0; j < nb; j++) {
                int x1 = symmdata_column_pointer[2 * j];
                int x2 = symmdata_column_pointer[2 * j + 1];

                ret = checkLMChelper (x1, x2 - 1, originalcol, arraycol, lperm, ad, rowsort);

                if (ret != LMC_EQUAL) {
                        return ret;
                }
        }

		_sort_array_column(rowsort, nb, lperm, symmdata_column_pointer, arraycol);

        return ret;
}

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
inline lmc_t LMC_check_col (const array_t *originalcol, const array_t *arraycol, levelperm_t lperm,
                            const arraydata_t *ad, const dyndata_t *dd) {
        lmc_t ret = LMC_EQUAL;
        int cur_row, rowp;

        const int oaindex = ad->oaindex;
        const int nrows = ad->N;
        rowsort_t *rowsort = dd->rowsort;
#ifdef DOPREFETCH
        prefetch (rowsort);
#endif

        const rowsort_value_t sval = ad->s[dd->col];

        /* we check in blocks of oaindex */
        for (int j = 0; j < (nrows / oaindex); j++) {
                if (oaindex > 1) {
                        /* sort rows according to new information */
                        for (int k = 0; k < oaindex; k++) {
                                cur_row = j * oaindex + k; // OPTIMIZE: use cur_row++
                                rowp = rowsort[cur_row].r;

#ifdef OADEBUG
                                debug_check_rowsort_overflow (ad, rowsort, dd, cur_row);
#endif

/* OPTIMIZE: the multiplication can be taken out of this function */
#ifdef SAFELPERM
                                rowsort[cur_row].val = ad->s[dd->col] * rowsort[cur_row].val +
                                                       safe_lperm< array_t > (array[rowp], lperm, ad->s[dd->col]);
#else
                                rowsort[cur_row].val = sval * rowsort[cur_row].val + lperm[arraycol[rowp]];
#endif
                        }

                        // OPTIMIZE: combine sorting and comparing
                        oacolSort (rowsort + (j * oaindex), 0, oaindex - 1);
                }


                for (int k = 0; k < oaindex; k++) {
                        cur_row = j * oaindex + k;
                        rowp = rowsort[cur_row].r;

#ifdef SAFELPERM
                        if (original[cur_row] < safe_lperm< array_t > (array[rowp], lperm, ad->s[dd->col])) {
                                ret = LMC_MORE; // first check for LMC_MORE, this happens much more often
                                break;
                        } else if (original[cur_row] >
                                   safe_lperm< array_t > (
                                       array[rowp], lperm,
                                       ad->s[dd->col])) { // the permuted array is lex. less than original
                                ret = LMC_LESS;
                                break;
                        }
#else
                        if (originalcol[cur_row] < lperm[arraycol[rowp]]) {
                                ret = LMC_MORE;
                                break;
                        } else if (originalcol[cur_row] >
                                   lperm[arraycol[rowp]]) { // the permuted array is lex. less than original
                                ret = LMC_LESS;
                                break;
                        }
#endif
                }

                if (ret != LMC_EQUAL)
                        break;
        }

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
lmc_t LMC_check_col (const array_t *original, const array_t *array, const arraydata_t *ad, const dyndata_t *dd) {
        lmc_t ret = LMC_EQUAL;
        int cur_row, rowp;
        int oaindex = ad->oaindex;
        int nrows = ad->N;
        rowsort_t *rowsort = dd->rowsort;

        /* we check in blocks of oaindex */
        for (int j = 0; j < (nrows / oaindex); j++) {
                if (oaindex > 1) {
                        /* sort rows according to new information */
                        for (int k = 0; k < oaindex; k++) {
                                cur_row = j * oaindex + k;

                                rowp = rowsort[cur_row].r;
#ifdef OADEBUG
                                assert (rowsort[cur_row].val > 0);
                                /* check for integer overflow */
                                if (std::numeric_limits< rowsort_value_t >::max () / ad->s[dd->col] <
                                    2 * rowsort[cur_row].val)
                                        log_print (SYSTEM, " LMC_check_col: integer overflow detected! s[col]=%ld, "
                                                           "rowsort[cur_row].val=%ld\n",
                                                   ad->s[dd->col], (long)rowsort[cur_row].val);
#endif
                                rowsort[cur_row].val = ad->s[dd->col] * rowsort[cur_row].val + array[rowp];
                        }

                        oacolSort (rowsort + (j * oaindex), 0, oaindex - 1);
                }

                for (int k = 0; k < oaindex; k++) {
                        cur_row = j * oaindex + k;
                        rowp = rowsort[cur_row].r;

                        if (original[cur_row] > array[rowp]) { // the permuted array is lex. less than original
                                ret = LMC_LESS;
                                break;
                        } else if (original[cur_row] < array[rowp]) {
                                ret = LMC_MORE;
                                break;
                        }
                }
                if (ret != LMC_EQUAL)
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
inline lmc_t LMC_check_col_complete (const array_t *original, carray_t *array, const arraydata_t *ad,
                                     const dyndata_t *dd) {
        lmc_t ret = LMC_EQUAL;
        const int oaindex = ad->oaindex;
        const int nrows = ad->N;
        rowsort_t *rowsort = dd->rowsort;

        /* we sort check in blocks of oaindex */
        for (int j = 0; j < (nrows / oaindex); j++) {
                if (oaindex > 1) {
                        /* sort rows according to new information */
                        for (int k = 0; k < oaindex; k++) {
                                int cur_row = j * oaindex + k;
                                int rowp = rowsort[cur_row].r;

                                rowsort[cur_row].val = ad->s[dd->col] * rowsort[cur_row].val + array[rowp];
                        }
                        oacolSort (rowsort + (j * oaindex), 0, oaindex - 1);
                }
        }

        for (rowindex_t cur_row = 0; cur_row < nrows; cur_row++) {
                int rowp = rowsort[cur_row].r;

                if (original[cur_row] < array[rowp]) {
                        ret = LMC_MORE;
                        break;
                } else if (original[cur_row] > array[rowp]) {
                        ret = LMC_LESS;
                        break;
                }
                if (ret != LMC_EQUAL)
                        break;
        }

        return ret;
}

/** Check column with ordering based on j5 */
inline lmc_t LMC_check_col_j5order (const array_t *original, const array_t *array, levelperm_t lperm,
                                    const arraydata_t *ad, const dyndata_t *dd) {

        assert (ad->order == ORDER_J5);

        const int cpoffset = ad->N * dd->colperm[dd->col];
        if (dd->col < 4) { 
                lmc_t ret = LMC_check_col (original + dd->col * +ad->N, array + cpoffset, lperm, ad, dd);
                return ret;
        }

        if (log_print (DEBUG, "")) {
                log_print (NORMAL, "LMC_check_col_j5order col %d, colperm: ", dd->col);
                print_perm (dd->colperm, ad->ncols);
        }

        colperm_t pp = new_perm_init< colindex_t > (5);

        pp[4] = dd->col;

        int jbase = abs (jvaluefast (original, ad->N, 5, pp));
        for (int x = 0; x < 5; x++)
                pp[x] = dd->colperm[x];
        pp[4] = dd->colperm[dd->col];
        int jcol = abs (jvaluefast (array, ad->N, 5, pp)); 
        if ( dd->col > 4) {
                if (log_print (DEBUG, "")) {
                        myprintf ("  xxx col %d, jbase %d, jcol %d: colperm: ", dd->col, jbase, jcol);
                        print_perm (pp, 5);
                }
        }
        delete_perm (pp);

        lmc_t ret;
        if (jbase == jcol) {
                ret = LMC_check_col (original + dd->col * +ad->N, array + cpoffset, lperm, ad, dd);
                return ret;
        } else {
                if (log_print (DEBUG, "")) {
                        myprintf ("  j5order: col %d, jbase %d, jcol %d: colperm: ", dd->col, jbase, jcol);
                        print_perm (pp, 5);
                }

                if (jbase ORDER_J5_SMALLER jcol)
                        ret = LMC_MORE;
                else
                        ret = LMC_LESS;
                return ret;
        }
}

lmc_t LMCreduce_non_root_j4 (const array_t *original, const arraydata_t *ad, const dyndata_t *dyndata,
                             LMCreduction_t *reduction, const OAextend &oaextend, LMCreduction_helper_t &helper_structure) {

        /* static allocation of data structures for non-root stage of LMC check */
        levelperm_t *lperm_p = 0;
        colperm_t *colperm_p = 0;
        colperm_t *localcolperm_p = 0;

        /* static init of dyndata */
        dyndata_t **dynd_p;
        int dynd_p_nelem;

        /* static allocation of buffer for the column check */
        array_t *colbuffer = 0;

        helper_structure.init_nonroot_stage (lperm_p, colperm_p, localcolperm_p, dynd_p, dynd_p_nelem, colbuffer, ad);

        const int number_remaining = ad->ncols - dyndata->col; 
        const int number_levels = ad->s[dyndata->col];      

        levelperm_t lperm = lperm_p[dyndata->col];
        colperm_t colpermloop = localcolperm_p[dyndata->col];
        dyndata_t *dyndatacpy = dynd_p[dyndata->col];
        const int current_column_group = ad->get_col_group (dyndata->col); 
        /* number of columns remaining in this group */
        const int ncolsremgroup = ad->colgroupsize[current_column_group] + ad->colgroupindex[current_column_group] - dyndata->col;


        lmc_t ret = LMC_MORE;
        init_perm< colindex_t > (colpermloop, number_remaining);
        const int col = dyndata->col;
        for (int i = 0; i < ncolsremgroup; i++) {
                std::swap (colpermloop[0], colpermloop[i]);

                // keep track of applied column permutation
                copy_perm (dyndata->colperm, dyndatacpy->colperm, ad->ncols);
                perform_inv_perm< colindex_t > (dyndata->colperm + col, dyndatacpy->colperm + col, number_remaining,
                                                colpermloop);

                const int nlevelperms = init_perm_n< array_t, int > (lperm, number_levels);
                const int cpoffset = ad->N * dyndatacpy->colperm[dyndata->col];

                /* loop over all level permutations */
                for (int j = 0; j < nlevelperms; j++) {
                        /* copy dyndata rowsort data */
                        cpy_dyndata_rowsort (dyndata, dyndatacpy);
                        dyndatacpy->col = dyndata->col; // TODO: needed?

                        /* LMC_check_col performs level permutations, updates the sorting structure and compares with
                         * the original array on blocks of oaindex */

                        if (reduction->mode >= OA_REDUCE) {
                        /* for reduce check perform full column perm!!!! */ /* NOTE: why? */
#ifdef SAFELPERM
                                safe_perform_level_perm< array_t > (original + cpoffset, colbuffer, ad->N, lperm,
                                                                    (int)ad->s[dyndata->col]);
#else
                                perform_level_perm (original + cpoffset, colbuffer, ad->N, lperm);
#endif
                                ret = LMC_check_col_less (reduction->array + dyndata->col * +ad->N,
                                                          original + cpoffset, lperm, ad, dyndatacpy);

                        } else {

                                if (ad->order == ORDER_LEX) {
#ifdef TPLUSCOLUMN
                                        ret = LMC_check_col_tplus (reduction->array + dyndata->col * +ad->N,
                                                                   original + cpoffset, lperm, ad, dyndatacpy);
#else
                                        ret = LMC_check_col (reduction->array + dyndata->col * +ad->N,
                                                             original + cpoffset, lperm, ad, dyndatacpy);
#endif
                                } else {
                                        ret =
                                            LMC_check_col_j5order (reduction->array, original, lperm, ad, dyndatacpy);
                                }
                        }

                        if (ret == LMC_LESS) {
                                reduction->lastcol = col;
                        }

                        if (reduction->mode >= OA_REDUCE && ret == LMC_LESS) {
                                reduction->updateLastCol (col);
                                reduction->updateFromLoop (*ad, dyndatacpy, lperm_p, original);
                                if ((reduction->nred % 150000) == 0)
                                        log_print (DEBUG - 1, "found reduction: column %d, number is %ld\n",
                                                   dyndata->col, reduction->nred);

                                ret = LMC_EQUAL; /* since the best array so far is updated, we need to continue */
                        }
                        if (ret == LMC_EQUAL) {
                                // this column could not decide, go one level deeper

                                dyndatacpy->col += 1;
                                if (dyndatacpy->col == ad->ncols) {
                                        /* this was the last column */
                                        ret = LMC_MORE;
                                } else {
                                        if (reduction->maxdepth == dyndata->col)
                                                ret = LMC_MORE;
                                        else
                                                /* pass the old cpy as an argument, since it will not be changed */
                                                ret = LMCreduce_non_root (original, ad, dyndatacpy, reduction,
                                                                          oaextend, helper_structure);
                                }
                        }
                        if (ret == LMC_LESS) {
                                break;
                        }
                        // else: LMC_MORE or LMC_EQUAL, continue with loop
                        next_perm (lperm, number_levels);
                }

                if (ret == LMC_LESS)
                        break;

                std::swap (colpermloop[0], colpermloop[i]); // swap back
        }

        return ret;
}

/** @perform The non-root stage of the LMC check or LMC reduction.
*
* @param original Pointer to array
* @param arrayclass Specification of array class
* @param dyndata
* @param reduction Transformation generated
* @param oaextend Options for the algorithm
* @param tmpStatic Static part of the algorithm
* @return
*/
lmc_t LMCreduce_non_root (const array_t *original, const arraydata_t *arrayclass, dyndata_t *dyndata,
                          LMCreduction_t *reduction, const OAextend &oaextend, const LMCreduction_helper_t &helper_structure) {
        /* static allocation of data structures for non-root stage of LMC check */
        levelperm_t *lperm_p = 0;
        colperm_t *colperm_p = 0;
        colperm_t *localcolperm_p = 0;

        /* static init of dyndata */
        dyndata_t **dynd_p;
        int dynd_p_nelem;

        /* static allocation of buffer for the column check */
        array_t *colbuffer = 0;

        helper_structure.init_nonroot_stage (lperm_p, colperm_p, localcolperm_p, dynd_p, dynd_p_nelem, colbuffer, arrayclass);

        const int number_remaining = arrayclass->ncols - dyndata->col; 
        const int nlevels = arrayclass->s[dyndata->col];      /* number of levels at this column */

        levelperm_t lperm = lperm_p[dyndata->col];
        colperm_t colpermloop = localcolperm_p[dyndata->col];
        dyndata_t *dyndatacpy = dynd_p[dyndata->col];
        const int current_column_group = arrayclass->get_col_group (dyndata->col);
        /* number of columns remaining in this group */
        const int ncolsremgroup = arrayclass->colgroupsize[current_column_group] + arrayclass->colgroupindex[current_column_group] - dyndata->col;

        myassert (dyndata->col != arrayclass->ncols,
                       "LMC_non_root this code should not be reached!  dyndata->col!=ad->ncols\n");

        lmc_t ret = LMC_MORE;
        init_perm< colindex_t > (colpermloop, number_remaining); // intialize for entire loop
        const int col = dyndata->col;
        for (int i = 0; i < ncolsremgroup; i++) {
                std::swap (colpermloop[0], colpermloop[i]); // swap 2 columns, we swap them back at the end of the loop
                                                            // so we do not need to initialize each time

                // keep track of applied column permutation
                copy_perm (dyndata->colperm, dyndatacpy->colperm, arrayclass->ncols);
                perform_inv_perm< colindex_t > (dyndata->colperm + col, dyndatacpy->colperm + col, number_remaining,
                                                colpermloop);

                const int nlevelperms = init_perm_n< array_t, int > (lperm, nlevels);
                const int cpoffset = arrayclass->N * dyndatacpy->colperm[dyndata->col];

                /* loop over all level permutations */
                for (int j = 0; j < nlevelperms; j++) {
                        /* copy dyndata rowsort data */
                        cpy_dyndata_rowsort (dyndata, dyndatacpy);
                        dyndatacpy->col = dyndata->col; 

                        /* LMC_check_col performs level permutations, updates the sorting structure and compares with
                         * the original array on blocks of oaindex */

                        if (reduction->mode >= OA_REDUCE) {
                                if (1) {
                                        ret = LMC_check_col_less (reduction->array + dyndata->col * +arrayclass->N,
                                                                  original + cpoffset, lperm, arrayclass, dyndatacpy);

                                } else {
#ifdef SAFELPERM
                                        safe_perform_level_perm< array_t > (original + cpoffset, colbuffer, ad->N,
                                                                            lperm, (int)ad->s[dyndata->col]);
#else
                                        perform_level_perm (original + cpoffset, colbuffer, arrayclass->N, lperm);
#endif
                                        ret = LMC_check_col_complete (reduction->array + dyndata->col * arrayclass->N,
                                                                      colbuffer, arrayclass, dyndatacpy);
                                }
                        } else {
                                if (arrayclass->order == ORDER_LEX) {

                                        if ((oaextend.getAlgorithm () == MODE_LMC_2LEVEL ||
                                             oaextend.getAlgorithm () == MODE_J5ORDER_2LEVEL) &&
                                            reduction->sd != 0) { 
                                                myassert (reduction->sd != 0, "LMC_check_col_ft");
                                                ret = LMC_check_col_ft_2level (
                                                    reduction->array + dyndata->col * +arrayclass->N, original + cpoffset,
                                                    lperm, arrayclass, dyndatacpy, *(reduction->sd)); 
                                        } else {
                                                ret = LMC_check_col (reduction->array + dyndata->col * +arrayclass->N,
                                                                     original + cpoffset, lperm, arrayclass, dyndatacpy);
                                        }
                                        dyndatacpy->col = dyndata->col; 

                                } else {
                                        ret =
                                            LMC_check_col_j5order (reduction->array, original, lperm, arrayclass, dyndatacpy);
                                }
                        }

                        if (ret == LMC_LESS) {
                                reduction->lastcol = col;
                        }

                        if (reduction->mode >= OA_REDUCE && ret == LMC_LESS) {
                                reduction->lastcol = col;
                                reduction->updateFromLoop (*arrayclass, dyndatacpy, lperm_p, original);
                                if ((reduction->nred % 150000) == 0)
                                        log_print (DEBUG - 1, "found reduction: column %d, number is %ld\n",
                                                   dyndata->col, reduction->nred);

                                ret = LMC_EQUAL; /* since the best array so far is updated, we need to continue */

                                reduction->mincol = std::min (reduction->mincol, col);
                        }
                        if (ret == LMC_EQUAL) {
                                // this column could not decide, go one level deeper
                                dyndatacpy->col += 1;
                                if (dyndatacpy->col == arrayclass->ncols) {
                                        /* this was the last column */
                                        ret = LMC_MORE;
                                } else {
                                        if (reduction->maxdepth == dyndata->col)
                                                ret = LMC_MORE;
                                        else {
                                                /* pass the old cpy as an argument, since it will not be changed */
                                                ret = LMCreduce_non_root (original, arrayclass, dyndatacpy, reduction,
                                                                          oaextend, helper_structure);
                                        }
                                }
                        }
                        if (ret == LMC_LESS) {
                                break; // break from levelperms loop
                        }
                        // else: LMC_MORE or LMC_EQUAL, continue with loop
                        next_perm (lperm, nlevels);
                }

                if (ret == LMC_LESS) {
                        break; // break from loop over available columns
                }

                std::swap (colpermloop[0], colpermloop[i]); // swap 2 columns back
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
lmc_t LMCreduce_non_root_2level (const array_t *original, const arraydata_t *ad, dyndata_t *dyndata,
                                 LMCreduction_t *reduction, const OAextend &oaextend,
                                 const LMCreduction_helper_t &tmpStatic) {
        const int dverbose = 0;

        if (dverbose)
                myprintf ("LMCreduce_non_root_2level: level %d\n", dyndata->col);

        /* static allocation of data structures for non-root stage of LMC check */
        levelperm_t *lperm_p = 0;
        colperm_t *colperm_p = 0;
        colperm_t *localcolperm_p = 0;

        /* static init of dyndata */
        dyndata_t **dynd_p;
        int dynd_p_nelem;

        /* static allocation of buffer for the column check */
        array_t *colbuffer = 0;

        tmpStatic.init_nonroot_stage (lperm_p, colperm_p, localcolperm_p, dynd_p, dynd_p_nelem, colbuffer, ad);

        if (dynd_p == 0) {
                throw_runtime_exception ("error: dynd_p not initialized!");
        }

        const int remsize = ad->ncols - dyndata->col; /* number of columns remaining */
        const int nlevels = 2; 		/* number of levels at this column */

        levelperm_t lperm = lperm_p[dyndata->col];
        dyndata_t *dyndatacpy = dynd_p[dyndata->col];
        dyndatacpy->allocate_rowsortl ();
        dyndatacpy->col = dyndata->col;

        const int current_column_group = ad->get_col_group (dyndata->col); 
        /* number of columns remaining in this group */
        const int ncolsremgroup = ad->colgroupsize[current_column_group] + ad->colgroupindex[current_column_group] - dyndata->col;


        lmc_t ret = LMC_MORE;
        const int col = dyndata->col;
        copy_perm (dyndata->colperm, dyndatacpy->colperm, ad->ncols);

        const int nlevelperms = 2; // special case for 2-level arrays
        for (int i = 0; i < ncolsremgroup; i++) {
                // keep track of applied column permutation
                std::swap (dyndatacpy->colperm[col], dyndatacpy->colperm[col + i]);

                lperm[0] = 0;
                lperm[1] = 1; // special init for 2-level arrays
                const int cpoffset = ad->N * dyndatacpy->colperm[dyndata->col];
                const array_t *reductionarrayoffset = reduction->array + dyndata->col * ad->N;

                /* loop over all level permutations */
                for (int j = 0; j < nlevelperms; j++) {
                        /* copy dyndata rowsort data */
                        cpy_dyndata_rowsortl (dyndata, dyndatacpy); 

                        dyndatacpy->col = dyndata->col; 

                        /* LMC_check_col performs level permutations, updates the sorting structure and compares with
                         * the original array on blocks of oaindex */

                        if (reduction->mode >= OA_REDUCE) {
                                myprintf ("LMCreduce_non_root_2level: mode OA_REDUCE not implemented\n");
                                return LMC_NONSENSE;
                        } else {
                                if (ad->order == ORDER_LEX) {

                                        if (dyndatacpy->col == ad->ncols - 1) {
                                                ret = LMC_check_col_ft_2level_rowsymm (
                                                    reductionarrayoffset, original + cpoffset, lperm, ad, dyndatacpy,
                                                    *(reduction->sd), 0);
                                        } else {
                                                ret = LMC_check_col_ft_2level (reductionarrayoffset,
                                                                               original + cpoffset, lperm, ad,
                                                                               dyndatacpy, *(reduction->sd), 0);
                                        }
                                        dyndatacpy->col = dyndata->col; 

                                } else {
                                        ret =
                                            LMC_check_col_j5order (reduction->array, original, lperm, ad, dyndatacpy);
                                }
                        }

                        if (ret == LMC_LESS) {
                                reduction->lastcol = col;
                        }

                        if (ret == LMC_EQUAL) {
                                // this column could not decide, go one level deeper
                                dyndatacpy->col += 1;
                                if (dyndatacpy->col == ad->ncols) {
                                        /* this was the last column */
                                        ret = LMC_MORE;
                                } else {
                                        if (reduction->maxdepth == dyndata->col)
                                                ret = LMC_MORE;
                                        else {
                                                /* pass the old cpy as an argument, since it will not be changed */
                                                ret = LMCreduce_non_root_2level (original, ad, dyndatacpy, reduction,
                                                                                 oaextend, tmpStatic);
                                        }
                                }
                                
                        }
                        if (ret == LMC_LESS) {
                                break; // break from levelperms loop
                        }
                        // else: LMC_MORE or LMC_EQUAL, continue with loop
                        next_perm_twoperm (lperm, nlevels);
                }

                if (ret == LMC_LESS) {
                        break; // break from loop over available columns
                }

                std::swap (dyndatacpy->colperm[col], dyndatacpy->colperm[col + i]);
        }

        return ret;
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
