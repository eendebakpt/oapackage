#ifdef OAEXTEND_MULTICORE
#include <mpi.h>
#endif
#include <algorithm>
#include <errno.h>
#include <list>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <math.h>
#else
#include <stdbool.h>
#include <unistd.h>
#endif

#include "arrayproperties.h"
#include "arraytools.h"
#include "extend.h"
#include "lmc.h"
#include "mathtools.h"
#include "strength.h"
#include "tools.h"

using namespace std;

/*!
 * Counts the occurence of each element in array
 * @brief Counts elements
 * @param array Pointer to array where elements are counted in
 * @param nelements
 * @param maxval Maximum value that can occur
 * @param elements
 */
void countelements (carray_t *array, const int nelements, const int maxval, int *elements) {

        memset (elements, 0, maxval * sizeof (int));

        for (int i = 0; i < nelements; i++)
                elements[array[i]]++;
}

std::vector< int > dextend_t::filterArrays (const array_link &al, const arraylist_t &earrays, arraylist_t &earraysout,
                                            std::vector< std::vector< double > > &edata, int verbose) {
        dextend_t &dextend = *this;
        int nn = dextend.filter.size ();

        std::vector< int > ctype (nn);
        int ngc = 0;
        for (int i = 0; i < nn; i++) {
                ctype[i] = dextend.filter[i] * (dextend.lmctype[i] == LMC_MORE);
        }

        int ngoodcombined = std::accumulate (ctype.begin (), ctype.end (), 0);

        array_link tmparray (al.n_rows, al.n_columns + 1, -1);
        std::copy (al.array, al.array + al.n_columns * al.n_rows, tmparray.array);

	cprintf (verbose >= 3, "dextend_t::filterArrays: filtering fraction %.3f: %d/%d\n", double(ngoodcombined) / nn, ngoodcombined, nn);
        flush_stdout ();
        const int lastcolx = al.n_columns;

        for (size_t idx = 0; idx < ctype.size (); idx++) {
                if (ctype[idx]) {
                        tmparray.setcolumn (lastcolx, earrays[idx]);
                        earraysout.push_back (tmparray);
                }
        }

        return ctype;
}

void dextend_t::DefficiencyFilter (double Dfinal, int k, int kfinal, double Lmax, int verbose) {
        int kn = k + 1;
        dextend_t &dextend = *this;

        double Cfinal = Dvalue2Cvalue (Dfinal, kfinal);
        double Cthr = Cfinal;
        double CthrMulti = Cfinal / pow (Lmax, kfinal - kn);

        double Lmaxmulti = pow (Lmax, kfinal - kn);

        int nn = dextend.Deff.size ();
        for (int ii = 0; ii < (int)nn; ii++) {

                double Ci = Dvalue2Cvalue (dextend.Deff[ii], kn);

                int chk = 1;
                switch (dextend.filtermode) {
                case DFILTER_NONE:
                        chk = 1;
                        break;
                case DFILTER_BASIC:
                        chk = Ci >= Cfinal;
                        break;
                case DFILTER_MULTI:
                        // chk= Ci >= Cfinalmulti;
                        chk = Lmaxmulti * Ci >= Cfinal;
                        break;
                }
                dextend.filter[ii] = chk;
        }
}

/*!  \brief Branche Information for extension algorithm
*
* This struct is used to keep track of the branches of different possibilities. If more than one option is available,
* this location should be stored at st[count] and count should be increased. For short, this could be done by
* st[count++]. If a column if filled, the highest element in use, gives the last branching point. Struct is used to
* represent a stack, since no random access is needed.
*/
struct split {
	///	Pointer to the stack
	rowindex_t *st; // alternative for a stack object
					///	Counter for number of elements on the stack and the number of options at the last spot, is obsolete?
	rowindex_t count;
	/// size of stack
	rowindex_t stacksize;

	/// number of valid positions
	array_t *nvalid;
	array_t *cvalidpos;
	///  valid positions
	array_t **valid;

	/**
	* @brief Constructor for stack structure
	* @param stacksize
	* @param maxs Maximum value in array
	*/
	split(rowindex_t stacksize, int maxs = 26) {
		log_print(DEBUG + 1, "split constructor: size %d\n", stacksize);
		this->stacksize = stacksize;
		this->st = new rowindex_t[stacksize];
		this->nvalid = new array_t[stacksize];
		this->cvalidpos = new array_t[stacksize];
		this->valid = malloc2d< array_t >(stacksize, maxs + 1);
		this->count = 0; // initialize stack
	}

	~split() {
		delete[] this->st;
		delete[] this->nvalid;
		delete[] this->cvalidpos;
		free2d(this->valid); //, this->stacksize);
	}

	inline rowindex_t position() { return this->count - 1; };

	/// show information about stack
	void show() {
		myprintf("stack: %d elements\n", this->count);
		for (int i = 0; i < this->count; i++) {
			myprintf("row %d (pos %d): ", this->st[i], this->cvalidpos[i]);
			print_perm(this->valid[i], this->nvalid[i]);
		}
	}
};

void OAextend::updateArraydata (arraydata_t *arrayclass) const {
        if (arrayclass == 0)
                return;

        switch (this->algmode) {
        case MODE_LMC_SYMMETRY:
        case MODE_ORIGINAL:
        case MODE_LMC_2LEVEL:
                arrayclass->order = ORDER_LEX;
                break;
        case MODE_J4:
                arrayclass->order = ORDER_LEX;
                break;
        case MODE_J5ORDER:
                arrayclass->order = ORDER_J5;
                break;
        case MODE_J5ORDERX:
        case MODE_J5ORDER_2LEVEL:
                arrayclass->order = ORDER_J5;
                arrayclass->order = ORDER_LEX;
                break;
        case MODE_AUTOSELECT:
        default:
                printfd ("error: OAextend::updateArraydata: no such mode %d for algorithm_t\n", this->algmode);
                break;
        }
}

std::string OAextend::__repr__ () const {
        std::string split = "\n";
        std::string s = printfstring ("OAextend: checkarrays %d", checkarrays) + split;
        s += printfstring ("OAextend: singleExtendTime %.1f [s], nLMC %d", singleExtendTime, nLMC) + split;
        s += printfstring ("OAextend: init_column_previous %d", init_column_previous) + split;
        s += printfstring ("OAextend: algorithm %s", this->getAlgorithmName ().c_str ()) + split;
        if (this->algmode == MODE_J5ORDERX || this->algmode == MODE_J5ORDER_2LEVEL) {
                s += printfstring ("OAextend: special: j5structure %d", this->j5structure) + split;
        }
        return s;
}

void OAextend::setAlgorithmAuto (arraydata_t *ad) {
        if (ad == 0) {
                myprintf ("setAlgorithmAuto: zero pointer\n");
                return;
        }
        algorithm_t x = OAextend::getPreferredAlgorithm (*ad);

        this->setAlgorithm (x, ad);
}

void OAextend::info(int verbose) const {
	std::cout << __repr__();

	if (verbose > 1) {
		myprintf("OAextend: use_row_symmetry: %d\n", this->use_row_symmetry);
	}
	std::cout << std::endl;
}

void OAextend::setAlgorithm (algorithm_t algorithm, arraydata_t *ad) {
        this->algmode = algorithm;

        switch (this->algmode) {
        case MODE_J5ORDER:
        case MODE_J5ORDERX:
        case MODE_J5ORDER_2LEVEL:
                init_column_previous = INITCOLUMN_J5;
                j5structure = J5_45;
                break;
        case MODE_LMC_SYMMETRY:
        case MODE_LMC_2LEVEL:
        case MODE_ORIGINAL:
        case MODE_J4:
                init_column_previous = INITCOLUMN_PREVIOUS;
                break;
        case MODE_AUTOSELECT:
                init_column_previous = INITCOLUMN_PREVIOUS;
                break;
        default:
                init_column_previous = INITCOLUMN_PREVIOUS;
                myprintf ("OAextend::setAlgorithm: error: algorithm (%d) unknown\n", algorithm);
                break;
        }

        this->updateArraydata (ad);
}


/// Find the range of the elements that are allowed
void get_range (array_t *array, extendpos *p, extend_data_t *es, int use_row_symmetry);

/**
 * @brief Returns the value of (part of) a row
 * @param array
 * @param start_idx
 * @param end_idx
 * @param n_rows
 * @param index Value index for each of the columns of the array
 * @return
 */
static inline array_t row_rank_partial (carray_t *array, const colindex_t start_idx, const colindex_t end_idx,
                                        const rowindex_t row, const rowindex_t n_rows, const int *index) {
        register int i, sum = 0, j = row;
        j += n_rows * start_idx;
        for (i = start_idx; i < end_idx; i++) {
                sum += index[i] * array[j];
                j += n_rows;
        }
        return sum;
}

/** @brief Find row symmetry structure of an array
 *
 * It is assumed the array is already in LMC form
 *
 * @param array
 * @param ad
 * @param col
 * @param gidx
 * @param gstart
 * @param gsize
 */
rowindex_t find_row_symm (carray_t *array, const arraydata_t *ad, const colindex_t col, rowindex_t *&gidx,
                          rowindex_t *&gstart, rowindex_t *&gsize) {
        int nsg = 0; // number of symmetry groups
        array_t *rowvalues;

        vindex_t *vindex = new_valueindex< vindex_t > (ad->s, col);

        rowvalues = new array_t[ad->N];
        for (int i = 0; i < ad->N; i++) {
                // NOTE: max value of rowvalues is of the order s^ncols
                rowvalues[i] = row_rank_partial (array, 0, col, i, ad->N, vindex);
        }

        nsg = symm_group_index_plain (rowvalues, ad->N, gidx, gstart, gsize);

        delete[] rowvalues;
        delete[] vindex;
        return nsg;
}

int compare_array_block (carray_p A, carray_p B, rowindex_t N, rowindex_t rs1, rowindex_t rs2, rowindex_t nrows,
                         colindex_t cs, colindex_t ncols) {
        for (int y = 0; y < ncols; y++) {
                int coloffset = (cs + y) * N;
                for (int x = 0; x < nrows; x++) {
                        if (A[x + rs1 + coloffset] < B[x + rs2 + coloffset]) {
                                return -1;
                        }
                        if (A[x + rs1 + coloffset] > B[x + rs2 + coloffset]) {
                                return 1;
                        }
                }
        }
        return 0;
}

/**
 * @brief Check whether a block exchange (from the first column blocks) leads to an LMC less array
 * @param array
 * @param N
 * @param blocksize
 * @param bidx1
 * @param bidx2
 * @param ncols
 * @return
 */
int check_block_exchange (carray_p array, rowindex_t N, int blocksize, int bidx1, int bidx2, colindex_t cs,
                          colindex_t clast) {
        rowindex_t rs1 = bidx1 * blocksize;
        rowindex_t rs2 = bidx2 * blocksize;

        colindex_t ncols = clast - cs + 1;
        if (compare_array_block (array, array, N, rs1, rs2, blocksize, cs, ncols) > 0) {
                logstream (DEBUG + 1) << "check_block_exchange: found exchange"
                                      << printfstring (" block %d, %d", bidx1, bidx2) << endl;
                return 1;
        }
        return 0;
}

/**
 * @brief Return next element from the stack
 * @param stack
 * @return
 */
inline array_t stack_next_element (const split *stack) {
        return stack->valid[stack->count - 1][stack->cvalidpos[stack->count - 1]];
}

/**
 * We check for all possible extensions. There are several criteria that can be used to check whether an element is
 * valid or not.
 * These are: strength check, strength 1 check (i.e. the number of elements), range check (due to LMC test)
 *
 * @param es
 * @param array
 * @param p
 * @param frequencies
 * @param esdyn
 * @param stack
 * @return
 */
int check_branch (extend_data_t *es, carray_t *array, extendpos *p, split *stack, array_t &first,
                  int use_row_symmetry) {
        // number of options left according to number of elements
        int npos = 0;

        const int posmax = p->ad->N / p->ad->s[p->col];
        const int coloffset = p->ad->N * p->col;

        array_t vmin;
        if (use_row_symmetry) {
                vmin = array[coloffset + p->row] >= es->range_low ? array[coloffset + p->row] : es->range_low;
        } else {
                vmin = (array[coloffset + p->row] >= 0) ? array[coloffset + p->row] : 0;
        }

#ifdef USE_SMALLSTEP
        const array_t vmax = std::min< array_t > (p->ad->s[p->col] - 1, es->range_high);
#else
        const array_t vmax = p->ad->s[p->col] - 1;
#endif

        /* loop over all positions specified by range */
        for (int i = vmin; i <= vmax; i++) {
#ifdef COUNTELEMENTCHECK
                /* strength 1 check */
                if (es->elements[i] < posmax) {
#else
                if (1) {
#endif
                        p->value = i; 
                        /* strength t check */
                        if (valid_element (es, p, array)) {
                                stack->valid[stack->count][npos] = p->value;
                                npos++;
                        } else {
                                // printf("reject on strength: %d value %d (row %d col %d)\n", range, p->value, p->
                                // row, p->col);
                        }
                } else {
                        // printf("reject on strength 1 check: value %d (row %d col %d)\n", p->value, p-> row, p->col);
                }
        }

        /* assign first valid element to variable */
        first = stack->valid[stack->count][0];

        /* if multiple positions then create a branch */
        if (npos > 1) {
                stack->st[stack->count] = p->row;
                stack->nvalid[stack->count] = npos;
                stack->cvalidpos[stack->count] = 0;
                stack->count++;
        }

        return npos;
}

/** Special version for 2-level array
 *
 * We check for all possible extensions. There are several criteria that can be used to check whether an element is
 * valid or not.
 * These are: strength check, strength 1 check (i.e. the number of elements), range check (due to LMC test)
 *
 * @param es
 * @param array
 * @param p
 * @param frequencies
 * @param esdyn
 * @param stack
 * @return
 */
int check_branch_2level (extend_data_t *es, carray_t *array_colstart, extendpos *p, split *stack, array_t &first,
                         int use_row_symmetry) {
        // number of options left according to number of elements
        int npos = 0;

        array_t vmin;
        if (use_row_symmetry) {
                vmin = array_colstart[p->row] >= es->range_low ? array_colstart[p->row] : es->range_low;
        } else {
                vmin = (array_colstart[p->row] >= 0) ? array_colstart[p->row] : 0;
        }

        const array_t vmax = 1;

        /* loop over all positions specified by range */
        for (int i = vmin; i <= vmax; i++) {
#ifdef COUNTELEMENTCHECK
                /* strength 1 check */
                if (es->elements[i] < posmax) {
#else
                if (1) {
#endif
                        p->value = i; 
                        /* strength t check */
                        if (valid_element_2level (es, p)) {
                                stack->valid[stack->count][npos] = p->value;
                                npos++;
                        } else {
                                // printf("reject on strength: %d value %d (row %d col %d)\n", range, p->value, p->
                                // row, p->col);
                        }
                } else {
                        // printf("reject on strength 1 check: value %d (row %d col %d)\n", p->value, p-> row, p->col);
                }
        }

        /* assign first valid element to variable */
        first = stack->valid[stack->count][0];

        /* if multiple positions then create a branch */
        if (npos > 1) {
                stack->st[stack->count] = p->row;
                stack->nvalid[stack->count] = npos;
                stack->cvalidpos[stack->count] = 0;
                stack->count++;
        }

        return npos;
}

/** @brief Find the range of the elements that are allowed
 * A lower bound is given by the symmetry structure of the array, an upper bound is given by the fact
 * that the elements cannot increase by more than 1 with respect to the current maximum
 *
 * @param array
 * @param p
 * @param esdyn
 */
void get_range (array_t *array, extendpos *p, extend_data_t *es, int use_row_symmetry) {
        array_t *array_col = array + p->ad->N * p->col;

        if (use_row_symmetry) {
                rowindex_t symmidx = es->gidx[p->row];
                array_t minvalue = 0;
                for (rowindex_t x = es->gstart[symmidx]; x < p->row; x++) {
                        minvalue = max (minvalue, array_col[x]);
                }
                es->range_low = minvalue;
        } else {
                es->range_low = 0;
        }

        /* determine maximum value */
        int maxval = -1;
        if (p->row < es->oaindextmin) {
                for (int j = 0; j < p->row; j++) {
                        array_t v = array_col[j];
                        maxval = (maxval < v) ? v : maxval;
                }
        } else {
                maxval = es->adata->s[es->extcolumn] - 1;
        }
        es->range_high = maxval + 1;
}

/*!
  For every column, some standard operations need to be performed:
  - first element is always zero
  - count its value in the number of elements
  - set remaining values in array to unused (= -1)
  - set parameters to start at second element
  \brief Prepare columnn for generating extensions
  \param array Pointer to OA
  \param p Characteristic numbers of OA
  \param col_offset Offset for first element in column, for faster calculations
  \param stack
  */
int *init_column (array_t *array, extendpos *p, int *col_offset, split *&stack) {
        log_print (DEBUG, "init_column: p->col = %d\n", p->col);
        int *elements;

        *col_offset = p->ad->N * p->col; // faster position calculation
#ifdef COUNTELEMENTCHECK

        elements = (int *)calloc (p->ad->s[p->col], sizeof (int));

        elements[0] = 1; // first element in column is always 0
#else
        elements = 0;
#endif
        array[*col_offset] = 0;                                                  // count it and set the value
        memset (&array[*col_offset + 1], -1, (p->ad->N - 1) * sizeof (array_t)); // set remainder of row to -1: unused
        p->row = 1;                                                              // continue with next element

        /* define the stack for the extensions tree */
        stack = new split (p->ad->N);

        return elements;
}

/** @brief Copy the frequency count table
*/
inline void copy_freq_table(strength_freq_table source, strength_freq_table target, int ftsize) {
	memcpy((void *)target[0], (void *)source[0], ftsize * sizeof(int));
}

/**
 * @brief Initialize extension with column
 * @param array
 * @param p
 * @param col_offset
 * @param stack
 * @param es
 * @param esdyn
 */
void init_column_full (array_t *array, extendpos *p, int &col_offset, split *&stack, extend_data_t *es) {
#ifdef COUNTELEMENTCHECK
        es->elements = init_column (array, p, &col_offset, stack);
#else
        init_column (array, p, &col_offset, stack);
#endif

        /* set cache structures */
        recount_frequencies (es->freqtable, es, p->col, 0, -1, array);

        for (int j = 1; j < p->row; j++) {
                p->row = j;
                add_element_freqtable (es, p->row - 1, array, es->freqtable);
                copy_freq_table (es->freqtable, es->freqtable_cache[j], es->freqtablesize);
        }
}

/**
 * @brief Initialize extension with previous column
 * @param array
 * @param p
 * @param col_offset
 * @param stack
 * @param es
 * @param esdyn
 */
void init_column_previous (array_t *array, extendpos *p, int &col_offset, split *&stack, extend_data_t *es,
                           const OAextend &oaextend) {
        log_print (DEBUG, "init_column_previous: p->col = %d\n", p->col);
        const rowindex_t N = p->ad->N;

        col_offset = N * p->col; // faster position calculation
#ifdef COUNTELEMENTCHECK
        es->elements = (int *)calloc (p->ad->s[p->col], sizeof (int));
#endif
        /* copy previous column */
        memcpy (&array[col_offset], &array[N * (p->col - 1)], N * sizeof (array_t));

        /* define the stack for the extensions tree */
        stack = new split (N);

        recount_frequencies (es->freqtable, (extend_data_t *)es, p->col, 0, -1, array); /* init freq to zero */

        for (int j = 1; j < N; j++) {
                p->row = j;
                get_range (array, p, es, oaextend.use_row_symmetry);
#ifdef COUNTELEMENTCHECK
                countelements (&array[col_offset], p->row, p->ad->s[p->col], es->elements);
#endif
                add_element_freqtable (es, p->row - 1, array, es->freqtable);

                array_t firstpos;
                int npos = check_branch (es, array, p, stack, firstpos, oaextend.use_row_symmetry);

                if (npos > 1) {
                        /* a branch was created, update the data tables */
                        copy_freq_table (es->freqtable, es->freqtable_cache[p->row], es->freqtablesize);
                }

                if (npos == 0) {
                        /* the current branch is empty, this happend for strength > 1 since the current column and
                           * the previous one cannot be equal */
                        firstpos = -1;
                }
                if (firstpos != array[col_offset + p->row]) {
                        /* trace back */
                        // if npos==0, then trace back and add 1 to stack pointer
                        for (int z = p->row; z < N; z++) {
                                array[col_offset + z] = -1;
                        }

                        if (npos > 1) {
                                stack->count--;
                        }
                        break;
                }
        }

        if (stack->count == 0) {
                /* since all elements of the column are set the extend loop thinks it will return to a braching point,
                 * since there is no branching point set we make the last element -1 */
                log_print (DEBUG, "init_column_previous: there is a unique extention to the array??\n");
                log_print (DEBUG, "stack->count %d, p->row %d (N %d)\n", stack->count, p->row, p->ad->N);
                array[col_offset + N - 1] = -1;
        }

        /* we set this for designs with strength 1 (here columns can be repeated) */
        array[col_offset + N - 1] = -1;
}

/** @brief Print progress on current column
 */
double progress_column (array_t *A, extendpos *p) {
        cout << "column: ";
        for (int i = 0; i < p->ad->N; i++) {
                cout << A[p->ad->N * p->col + i] << " ";
        }
        cout << endl;
        return 0;
}

/** @brief Predict j4(1,2,3,k) using the theorem from Deng
* This works only for 2-level arrays. The 0 corresponds to a +
*
*/
inline int predictJ(const array_t *array, const int N, const int k) {
	int t = N / 4;
	int tt = t / 2;

	int number_of_plus_values = 0;
	for (int i = 0; i < tt; i++) {
		if (array[k * N + i] == 0) {
			number_of_plus_values++;
		}
	}
	for (int i = tt; i < t; i++) {
		if (array[k * N + i] == 1) {
			number_of_plus_values++;
		}
	}

	return 8 * number_of_plus_values - N;
}

/// perform Jcheck, return 1 if branch should be cut
int Jcheck (carray_t *array, const rowindex_t N, const int jmax, const extendpos *p) {
        if (p->row == N / 4 + 1) {
                int Jpred = predictJ (array, N, p->col);
                if (abs (Jpred) > abs (jmax)) {
                        logstream (NORMAL) << printfstring ("   cutting branch jpred %d, jmax %d!\n", Jpred, jmax);
                        return 1;
                }
        }
        return 0;
}

/// return string with current time
inline std::string printtime() {
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	return printfstring("%s", asctime(timeinfo));
}

/**
 *
 * @brief Check the stack of extensions
 * @param
 * @param p
 * @param array
 * @return
 */
bool return_stack (split *stack, extendpos *p, array_t *array, int col_offset) {
        int more_branches = true;
        if (stack->count > 0) {
                /* more branch points on stack */
                p->row = stack->st[stack->count - 1];
                memset (&array[col_offset + p->row + 1], -1, (p->ad->N - p->row - 1) * sizeof (array_t));
        } else {
                /* no more branch points: stop */
                more_branches = false;
        }

        return more_branches;
}

inline void showLoopProgress (array_t *array, const int col_offset, const rowindex_t N, const int node_rank = 0,
                              int nlmcarrays = -1) {
#ifdef _OPENMPX
#define SAFEOMP
#ifdef SAFEOMP
        static long _nloops[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        int tid = omp_get_thread_num ();
        _nloops[tid % 16]++;
        long nloops = _nloops[tid % 16]; 
#else
        static long nloops = 0;
        nloops++;
#endif
#else
        static long nloops = 0;
        nloops++;
#endif
        if (nloops % 10000 == 0) {
                fflush (stdout);

                if (nloops % (500 * 1000 * 1000) == 0) {
                        std::cout << "node [" << node_rank << "]: extend loop " << nloops / (1000 * 1000);
                        cout << "m, ";
                        print_perm (array + col_offset, N, 28);
                }
        }
}

typedef std::vector< array_link > extensioncol_list_t;


// This function is quite long (and complicated) since it performs the entire extension of arrays. The LMC test is
// performed by a function call to LMCreduce.
//
// For calculation efficiency several cache systems are used counting the number of occurences of elements and
// t-tuples.
//
int extend_array (const array_link &input_array, const arraydata_t *fullad, const colindex_t extensioncol,
                  arraylist_t &extensions, OAextend const &oaextend) {

        carray_t *origarray = input_array.array;
        const int start_number = extensions.size ();
        const colindex_t ncolsextension = extensioncol + 1;

        /* make copy of the array with one additional column */
        array_t *array = create_array (fullad->N, ncolsextension);
        copy_array (origarray, array, fullad->N, extensioncol);

#ifdef OACHECK
        if (fullad->strength < 1) {
                log_print (SYSTEM, " extend_array: error: function not defined for strength < 1\n");
                throw_runtime_exception("extend_array: strength should be >=1");
        }
#endif

        /* array data */
        arraydata_t *ad = new arraydata_t (fullad, ncolsextension);
        oaextend.updateArraydata (ad);

        extendpos *p = new extendpos (extensioncol, ad);

        /* extension data */
        extend_data_t *es = new extend_data_t (ad, p->col);

        rowindex_t nsg = find_row_symm (array, ad, ncolsextension - 1, es->gidx, es->gstart, es->gsize);
        log_print (DEBUG, "  number of row symmetry groups: %d\n", nsg); // print_perm(esdyn->gidx, p->N);

        /* the stack contains data about the branches found during the extension */
        split *stack = 0;

#ifdef FREQELEM
        /* set elem frequencies */
        es->init_frequencies (array);
#endif

        /* check whether we are in the same column group */
        int col_offset;
        if (oaextend.init_column_previous == INITCOLUMN_PREVIOUS && ad->s[p->col] == ad->s[p->col - 1]) {
                init_column_previous (array, p, col_offset, stack, es, oaextend);
        } else {
                if (oaextend.init_column_previous == INITCOLUMN_J5) {
                        if (extensioncol > 5) {
                                // printf("extensioncol %d: here\n", extensioncol);
                                init_column_previous (array, p, col_offset, stack, es, oaextend);
                        } else
                                init_column_full (array, p, col_offset, stack, es);

                } else {
                        init_column_full (array, p, col_offset, stack, es);
                }
        }

        /* set the frequencies table */
        recount_frequencies (es->freqtable, es, p->col, 0, p->row - 2,
                             array); // NOTE: is this not also done in init_column_xxx?

        log_print (DEBUG, "Starting extension\n");

        array_t *array_colstart = array + p->col * ad->N;
        bool more_branches = true;

        double extendTime = get_time_ms ();
        unsigned long narrays = 0, nlmcarrays = 0;

#ifdef COUNTELEMENTCHECK
        countelements (array_colstart, 0, p->ad->s[p->col], es->elements);
#endif

        const rowindex_t N = p->ad->N;
#ifdef OAEXTEND_MULTICORE
        const int node_rank = MPI::COMM_WORLD.Get_rank ();
#else
        const int node_rank = 0;
#endif

        LMCreduction_t reduction (ad);
        reduction.init_state = COPY;
        reduction.initStatic (); // needed to make the extend code thread safe

        do {
                showLoopProgress (array, col_offset, N, node_rank, nlmcarrays);

                if ((extensions.size () * ad->N * extensioncol) > (1024 * 1024 * 1024 / sizeof (array_t))) {
                        myprintf ("extend_array: ERROR: running out of memory! extensions.size() %d, ad->N %d\n",
                                  (int)extensions.size (), ad->N);
                        array_link al (origarray, ad->N, extensioncol);
                        myprintf ("row symmetry group of array: \n");
                        al.row_symmetry_group ().show ();
                        throw_runtime_exception("number of extensions too large");
                }

                if (p->row < N) { /*column is not yet full */

#ifdef SYMMBLOCKS
                        const int bsize = N / p->ad->s[0];

                        if (((p->row % bsize) == 0) && (p->row > (bsize * 2))) {
                                int bidx1 = (p->row / bsize) - 2;
                                int bidx2 = (p->row / bsize) - 1;

                                if (check_block_exchange (array, p->ad->N, bsize, bidx1, bidx2, p->ad->strength,
                                                          p->col) == 1) {
                                        more_branches = return_stack (stack, p, array, col_offset);

                                        if (more_branches)
                                                continue;
                                        else {
                                                logstream (DEBUG) << "cutbranch: no more branches" << endl;
                                                break;
                                        }
                                }
                        }
#endif

                        if (array_colstart[p->row] == -1) { // current entry is empty
                                // frequencies table from previous loop is clean, only add a single element
                                add_element_freqtable_col (es, p->row - 1, array_colstart, es->freqtable);

                                get_range (array, p, es, oaextend.use_row_symmetry);
#ifdef COUNTELEMENTCHECK
                                addelement (array_colstart[p->row - 1], es->elements);
#endif

                                array_t firstpos;
                                int npos = check_branch (es, array, p, stack, firstpos, oaextend.use_row_symmetry);

                                if (npos > 1)
                                        copy_freq_table (es->freqtable, es->freqtable_cache[p->row],
                                                         es->freqtablesize);

                                if (npos == 0) { /* reached dead end, goto branch point */
                                        more_branches = return_stack (stack, p, array, col_offset);
                                } else {
                                        array_colstart[p->row] = firstpos;
                                        p->row++;
                                }
                        } else { /* came back from branche */
                                /* copy frequencies from cache */
                                copy_freq_table (es->freqtable_cache[p->row], es->freqtable, es->freqtablesize);
#ifdef COUNTELEMENTCHECK
                                /* range and elements do not need to be copied, they are recalculated */
                                countelements (array_colstart, p->row, p->ad->s[p->col], es->elements);
#endif

                                stack->cvalidpos[stack->count - 1]++;

                                if (stack->cvalidpos[stack->count - 1] == stack->nvalid[stack->count - 1]) {
                                        array_colstart[p->row] = -1; /* should be done in return_stack ? */
                                        stack->count--;              /* remove element from stack */
                                        more_branches = return_stack (stack, p, array, col_offset);
                                } else {
                                        array_colstart[p->row] = stack_next_element (stack);
                                        p->row++;
                                }
                        }
                } else { // reached end of column
                        narrays++;

                        if ((narrays % oaextend.nLMC == 0) ||
                            ((get_time_ms () - extendTime) > oaextend.singleExtendTime)) {
                                extendTime = get_time_ms ();
                                if (log_print (QUIET, "")) {
                                        logstream (QUIET) << printfstring ("  OA extension: ") << narrays - 1
                                                          << " arrays checked, " << extensions.size ()
                                                          << " solutions so far";
                                        logstream (QUIET) << ", time " << printtime ();
                                        myprintf ("  OA extension: current row: ");
                                        print_perm (array + col_offset, N, 36);
                                }
                        }

                        lmc_t lmc;
                        reduction.reset ();
                        reduction.mode = OA_TEST;

                        if (oaextend.checkarrays == 0)
                                lmc = LMC_MORE;
                        else {
                                lmc = LMCcheck (array, *ad, oaextend, reduction);
                        }

                        if (lmc == LMC_MORE || lmc == LMC_EQUAL) {
                                if (checkloglevel (DEBUG)) {
                                        std::cout << "Found array:" << endl;
                                        print_array (array, ad->N, ad->ncols);
                                }
                                /* the extension found is LMC */
                                switch (oaextend.extendarraymode) {
                                case OAextend::APPENDFULL: {
                                        array_link tmp_extension (array, N, p->col + 1, nlmcarrays);
                                        extensions.push_back (tmp_extension);
                                } break;
                                case OAextend::APPENDEXTENSION: {
                                        array_link tmp_extension (array + p->col * N, N, 1, nlmcarrays);
                                        extensions.push_back (tmp_extension);
                                } break;
                                case OAextend::STOREARRAY: {
                                        array_link tmp_extension (array, N, p->col + 1, nlmcarrays);
                                        arrayfile_t *storefile =
                                            (arrayfile_t *)&oaextend.storefile; // trick to prevent const warnings
                                                                                /** Note:
                                                                                *
                                                                                * If you only want to store a selection of arrays, here one should insert
                                                                                * a check on the properties of the array.
                                                                                *
                                                                                */
                                        storefile->append_array (tmp_extension);
                                } break;
                                case OAextend::NONE: {

                                        // do nothing
                                } break;
                                default:
                                        myprintf ("extend_array: no such mode!\n");
                                        break;
                                }
                                nlmcarrays++;
                        }
                        more_branches = return_stack (stack, p, array, col_offset);

                        if (oaextend.check_maximal) {
                                // abort the algorithm
                                more_branches = false;
                        }
                }
        } while (more_branches); /* continue as long as more branches are left */

        log_print (QUIET, "Found %lu arrays in total, only %d in LMC form, total is %d\n", narrays,
                   (int)(extensions.size () - start_number), (int)extensions.size ());

        delete es;
        destroy_array (array);
        delete p->ad;
        delete p;
        delete stack;

        return narrays;
}

arraylist_t extend_array (const array_link &al, arraydata_t &fullad, OAextend const &oaextend) {
        arraylist_t alist;
        alist.push_back (al);
        return extend_arraylist (alist, fullad, oaextend);
}

arraylist_t extend_array (const array_link &al, arraydata_t &fullad) {
        OAextend oaextend;
        oaextend.setAlgorithm (MODE_ORIGINAL);
        arraylist_t alist;
        alist.push_back (al);
        return extend_arraylist (alist, fullad, oaextend);
}

arraylist_t extend_arraylist (const arraylist_t &alist, const arraydata_t &arraydata) {
        arraydata_t atmp (arraydata); // make a copy so arraydata can be const
        OAextend oaextend (atmp);
        return extend_arraylist (alist, atmp, oaextend);
}

arraylist_t extend_arraylist (const arraylist_t &alist, arraydata_t &fullad, OAextend const &oaextend) {
        colindex_t extensioncol = -1;

        arraylist_t ll;
        extend_arraylist (alist, fullad, oaextend, extensioncol, ll);
        return ll;
}

int extend_arraylist (const arraylist_t &alist, arraydata_t &fullad, OAextend const &oaextend,
                      colindex_t extensioncolx, arraylist_t &extensions) {
        colindex_t extensioncol;
        int n = 0;
        for (arraylist_t::const_iterator al = alist.begin (); al != alist.end (); al++) {
                if (extensioncolx == -1)
                        extensioncol = al->n_columns;
                else
                        extensioncol = extensioncolx;

                if (extensioncol != al->n_columns) {
                        myprintf ("problem: array has %d columns, extension to column %d is requested\n",
                                  al->n_columns, extensioncol);
                        continue;
                }
                n += extend_array (*al, &fullad, extensioncol, extensions, oaextend);
        }
        fflush (stdout);
        return n;
}

arraylist_t runExtendRoot (arraydata_t arrayclass, int max_number_columns, int verbose) {

        array_link al = arrayclass.create_root ();
        OAextend oaoptions;
        oaoptions.setAlgorithmAuto (&arrayclass);

        arraylist_t solutions;
        solutions.push_back (al);

        for (int ii = arrayclass.strength; ii < max_number_columns; ii++) {
                arraylist_t extensions;
                extend_arraylist (solutions, arrayclass, oaoptions, ii, extensions);
                if (verbose) {
                        myprintf ("runExtend: ncols %d: %ld arrays\n", ii + 1, (long)extensions.size ());
                }
                solutions = extensions;
        }
        return solutions;
}
