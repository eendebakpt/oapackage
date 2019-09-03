#include <algorithm>
#include <math.h>

#include "arraytools.h"
#include "extend.h"
#include "lmc.h"
#include "mathtools.h"
#include "oaoptions.h"
#include "tools.h"

#include "nonroot.h"

#ifndef DOOPENMP
// no OpenMP
// disable omp pragma warnings
#ifdef WIN32
#pragma warning(disable : 4068)
#elif defined(__APPLE__)
#pragma warning(disable : 4068)

bool operator!= (symmdataPointer const &ptr, int x) {
        // A friendly reminder to not test pointers against any values except 0 (NULL)
        myassert (!x, "invalid pointer in comparison");
        return ptr;
}

#else
// http://choorucode.com/2014/03/11/how-to-ignore-specific-warning-of-gcc/
#pragma GCC diagnostic ignored "-Wenum-compare"
#endif
#endif

#ifdef DOOPENMP
#include <omp.h>
#endif

using namespace std;

// worker pool of static objects

static object_pool< LMCreduction_helper_t > LMCreduction_pool;

LMCreduction_helper_t *acquire_LMCreduction_object () {
#ifdef NOREUSE
        LMCreduction_helper_t *pp = new LMCreduction_helper_t ();
        pp->id = 100 * omp_get_thread_num () + int(random () % 40);
        return pp;
#endif


        LMCreduction_helper_t *p = 0;
#pragma omp critical
        { p = LMCreduction_pool.New (); }
        return p;
}
void release_LMCreduction_object (LMCreduction_helper_t *p) {
#ifdef NOREUSE
        delete p;
        return;
#endif

#pragma omp critical
        {
                LMCreduction_pool.Delete (p);
        }
}

void clear_LMCreduction_pool () { LMCreduction_pool.reset (); }

LMCreduction_helper_t::~LMCreduction_helper_t () {
        this->freeall ();
}

LMCreduction_helper_t::LMCreduction_helper_t () {
        /* initialize static structure values to zero */
        this->LMC_non_root_init = 0;
        this->LMC_root_init = 0;
        this->LMC_root_rowperms_init = 0;
        this->rootrowperms_full = 0;

        this->ad = 0;
        this->colbuffer = 0;
        this->nrootrowperms = 0;
        this->rootrowperms = 0;
        this->dyndata_p = 0;
        this->current_trans = 0;
        this->colperm_p = 0;
        this->localcolperm_p = 0;
}

/// return name of the algorithm
std::string algnames (algorithm_t m) {
        std::string str;
        switch (m) {
        case MODE_ORIGINAL:
                str = stringify (MODE_ORIGINAL);
                break;
        case MODE_LMC_2LEVEL:
                str = stringify (MODE_LMC_DEBUG);
                break;
        case MODE_J4:
                str = stringify (MODE_J4);
                break;
        case MODE_LMC_SYMMETRY:
                str = stringify (MODE_LMC_SYMMETRY);
                break;
        case MODE_AUTOSELECT:
                str = stringify (MODE_AUTOSELECT);
                break;
        case MODE_J5ORDER:
                str = stringify (MODE_J5ORDER);
                break;
        case MODE_J5ORDERX:
                str = stringify (MODE_J5ORDERX);
                break;
        case MODE_J5ORDER_2LEVEL:
                str = stringify (MODE_J5ORDERXFAST);
                break;
        default:
        case MODE_INVALID:
                printfd ("MODE_INVALID!");
                str = stringify (MODE_INVALID);
                break;
        }
        return str;
}

template < class numtype, class objecttype >
/** Convert number to a permutation and apply to an array
* See also http://en.wikipedia.org/wiki/Permutation#Numbering_permutations
* @param permutation_index Index of permutation to apply
* @param array Array of objects to be permuted
* @param length Length of array
* @return
*/
void permutationLex (numtype permutation_index, objecttype *array, numtype length) {
        numtype fact = factorial (length - 1);
        objecttype tempj, temps;

        for (int j = 0; j < length - 1; j++) {
                tempj = (permutation_index / fact) % (length - j);
                temps = array[j + tempj];
                for (int i = j + tempj; i >= j + 1; i--) {
                        array[i] = array[i - 1]; // shift the chain right
                }
                array[j] = temps;
                fact = fact / (length - (j + 1));
        }
}

/**
* @brief From a given permutation index construct the corresponding level permutations
* @param permindex
* @param ad
* @param lperms
*/
void root_row_permutation_from_index (int permindex, const arraydata_t *ad, levelperm_t *lperms) {
        for (colindex_t i = 0; i < ad->strength; i++) {
                const int col = ad->strength - i - 1;
                int fac = factorial< int > (ad->s[col]);

                int idx = (permindex % fac);
                int rem = (permindex / fac);

                init_perm (lperms[col], ad->s[col]);

                permutationLex< array_t > (idx, lperms[col], ad->s[col]);
                permindex = rem;
        }
}

/**
* @brief Perform a row permutation on an array using a rowsort_t structure
* @param source
* @param target
* @param rs
* @param nrows
* @param ncols
*/
void perform_row_permutation (const carray_t *source, array_t *target, rowsort_t *rs, rowindex_t nrows,
                              colindex_t ncols) {
        int t;
        for (colindex_t i = 0; i < ncols; i++) {
                for (rowindex_t j = 0; j < nrows; j++) {
                        t = rs[j].r; // perm[j];
                        target[nrows * i + j] = source[nrows * i + t];
                }
        }
}


void debug_show_reduction (LMCreduction_t *reduction, carray_t *original, const arraydata_t *ad,
                           const char *str = "##\n") {
        myprintf ("####3 new reduction:\n");
        reduction->transformation->show ();
        print_array ("original:\n", original, ad->ncols, ad->N);
        reduction->transformation->apply (original, reduction->array);
        print_array ("update:\n", reduction->array, ad->ncols, ad->N);
}

void debug_print_transformations (LMCreduction_t *reduction, carray_t *original, carray_t *array,
                                  const arraydata_t *ad, const dyndata_t *dyndata) {
        myprintf ("### LMCreduce_non_root: col %d\n", dyndata->col);
        myprintf ("current array:\n");

        reduction->transformation->show ();
        myprintf ("orignal:\n");
        print_array (original, ad->ncols, ad->N);
        myprintf ("current:\n");
        reduction->transformation->print_transformed (original);
}

#ifdef OAOVERFLOW
int check_bounds_rowsort (arraydata_t *ad) {
        long long maxval = 1;

        for (int x = 0; x < ad->ncols; x++) {
                if (maxval > (std::numeric_limits< rowsort_value_t >::max () / ad->s[x])) {
                        logstream (SYSTEM) << "check_bounds_rowsort: possible overflow";
                }
                maxval *= ad->s[x];
        }
        return 0;
}
#endif



LMCreduction_t::LMCreduction_t (const arraydata_t *adp) {
        mode = OA_TEST;

        transformation = new array_transformation_t (adp);
        array = create_array (adp);
        sd = symmdataPointer ((symmdata *)0);
        maxdepth = -1;

        staticdata = 0;

        ncols = adp->ncols;
        nrows = adp->N;
        reset ();
#ifdef SWIGCODE
        init_state = COPY;
#else
        init_state = INIT_STATE_INVALID; // remove this later
#endif
}

LMCreduction_t::~LMCreduction_t () {
        free ();

        releaseStatic ();
}

void LMCreduction_t::free () {
        delete transformation;
        transformation = 0;
        destroy_array (array);
        array = 0;
}

#define MINCOLMAX 1000

LMCreduction_t::LMCreduction_t (const LMCreduction_t &at) {
        mode = at.mode;
        state = at.state;
        init_state = at.init_state;

        maxdepth = at.maxdepth;
        lastcol = at.lastcol;
        nred = at.nred;
        ncols = at.ncols;
        nrows = at.nrows;

        targetcol = at.targetcol;
        mincol = at.mincol;
        staticdata = 0;

        sd = symmdataPointer ((symmdata *)0);

        transformation = new array_transformation_t (*(at.transformation));
        array = create_array (transformation->ad);
}

LMCreduction_t &LMCreduction_t::operator= (const LMCreduction_t &at) 
{
        mode = at.mode;
        state = at.state;
        init_state = at.init_state;
        maxdepth = at.maxdepth;
        lastcol = at.lastcol;
        nred = at.nred;

        targetcol = at.targetcol;
        mincol = at.mincol;
        ncols = at.ncols;
        nrows = at.nrows;

        staticdata = at.staticdata; 

        free ();

        transformation = new array_transformation_t (*(at.transformation));
        array = create_array (transformation->ad);

        return *this;
}

void LMCreduction_t::reset () {
        lastcol = -1;
        nred = 0;
        state = REDUCTION_INITIAL;
        init_state = COPY;
        transformation->reset ();

        targetcol = 0;
        mincol = MINCOLMAX;

}

void LMCreduction_t::show(int verbose) const {
	myprintf("LMCreduction_t: mode %d, state %d (REDUCTION_INITIAL %d, REDUCTION_CHANGED %d), init_state "
		"%d, lastcol %d\n",
		this->mode, this->state, REDUCTION_INITIAL, REDUCTION_CHANGED, this->init_state,
		this->lastcol);
	if (verbose >= 1) {
		myprintf("LMCreduction_t: nred %ld\n", nred);
		print_array("array:\n", this->array, this->transformation->ad->N,
			this->transformation->ad->ncols);
	}
	if (verbose >= 2)
		this->transformation->show();
}

std::string LMCreduction_t::__repr__ () const {
	std::string ss = printfstring ("LMCreduction_t: mode %d, state %d (REDUCTION_INITIAL %d, "
					"REDUCTION_CHANGED %d), init_state %d, lastcol %d\n",
					this->mode, this->state, REDUCTION_INITIAL, REDUCTION_CHANGED,
					this->init_state, this->lastcol);
	return ss;
}
void LMCreduction_t::updateTransformation (const arraydata_t &ad, const dyndata_t &dyndatacpy, levelperm_t *lperm_p,
                                           const array_t *original) {

        nred++;
        this->state = REDUCTION_CHANGED;
        /* copy level permutations, note: expensive operation */
        for (colindex_t c = 0; c < ad.ncols; c++)
                copy_perm< array_t > (lperm_p[c], this->transformation->lperms[c], ad.s[c]);

        /* copy row permutation */
        dyndatacpy.getRowperm (this->transformation->rperm);

        /* copy column permutation */
        copy_perm (dyndatacpy.colperm, this->transformation->cperm, ad.ncols);

}

void LMCreduction_t::updateFromLoop (const arraydata_t &ad, const dyndata_t &dyndatacpy, levelperm_t *lperm_p,
                                     const array_t *original) {

        nred++;
        this->state = REDUCTION_CHANGED;
        /* copy level permutations, note: expensive operation */
        for (colindex_t c = 0; c < ad.ncols; c++)
                copy_perm< array_t > (lperm_p[c], this->transformation->lperms[c], ad.s[c]);

        /* copy row permutation */
        dyndatacpy.getRowperm (this->transformation->rperm);

        /* copy column permutation */
        copy_perm (dyndatacpy.colperm, this->transformation->cperm, ad.ncols);

        /* update array used for comparing */

        this->transformation->apply (original, this->array);

}

/** Apply a random transformation to an array (inplace) **/
void random_transformation (array_t *array, const arraydata_t *adp) {
        array_transformation_t *transformation = new array_transformation_t (adp);
        transformation->randomize ();

        array_t *cpy = clone_array (array, adp->N, adp->ncols);
        transformation->apply (cpy, array);
}

/**
* @brief Apply Hadamard transformation to orthogonal array
* @param source
* @param target
* @param hcol
*/
void apply_hadamard (const arraydata_t *ad, array_t *array, colindex_t hcol) {
        if ((ad->N - 1) != ad->ncols) {
                myprintf ("WARNING: array does not have size N x (N-1), Hadamard transformation makes no sense\n");
				throw_runtime_exception("apply_hadamard: input array does not have proper shape");
        }
        for (colindex_t z = 0; z < ad->ncols; z++) {
                if (ad->s[z] != 2) {
                        myprintf ("error: Hadamard transformation only defined for 2-level arrays\n");
                        return;
                }
        }

        carray_t idperm[2] = {0, 1};
        carray_t swapperm[2] = {1, 0};
        carray_t *perm = 0;

        for (rowindex_t r = 0; r < ad->N; r++) {
                /* detemine action on this row */
                if (array[r + ad->N * hcol] == 0)
                        perm = idperm;
                else
                        perm = swapperm;

                for (colindex_t c = 0; c < ad->ncols; c++) {
                        if (c == hcol)
                                continue;

                        array[c * ad->N + r] = perm[array[ad->N * c + r]];
                }
        }
}

/// Apply Hadamard transformation to orthogonal array
void apply_hadamard (array_link &al, colindex_t hcol) {
        arraydata_t adata = arraylink2arraydata (al);
        apply_hadamard (&adata, al.array, hcol);
}

void dyndata_t::reset () {
        init_perm (this->colperm, this->N);

        for (int i = 0; i < N; i++) {
                this->rowsort[i].r = i;
                this->rowsort[i].val = 0;
        }
}

void dyndata_t::getRowperm(rowpermtypelight &rp) const {
	rp.resize(this->N);
	if (this->rowsortl == 0) {
		for (int i = 0; i < this->N; i++)
			rp[i] = this->rowsort[i].r;
	}
	else {
		for (int i = 0; i < this->N; i++)
			rp[i] = this->rowsortl[i];
	}
}

void dyndata_t::getRowperm(rowperm_t &rperm) const {
	if (this->rowsortl == 0) {
		for (rowindex_t x = 0; x < this->N; x++)
			rperm[x] = this->rowsort[x].r;
	}
	else {
		for (rowindex_t x = 0; x < this->N; x++)
			rperm[x] = this->rowsortl[x];
	}
}

rowpermtypelight dyndata_t::getRowperm() const {
	rowpermtypelight rp(this->N);
	this->getRowperm(rp);
	return rp;
}

colpermtypelight dyndata_t::getColperm() const {
	colpermtypelight cp(this->colperm, this->col + 1);
	return cp;
}
void dyndata_t::getColperm(colpermtypelight &cp) const {
	cp.resize(this->col + 1);
	std::copy(this->colperm, this->colperm + this->col + 1, cp.data_pointer);
}

/// allocate structure to keep track of row sorting
rowsort_t * allocate_rowsort(int N) {
  return (rowsort_t *)malloc (sizeof (rowsort_t) * N);
}

void deallocate_rowsort(rowsort_t *& rowsort) {
  free (rowsort);
  rowsort = 0;
}

rowsorter_t::rowsorter_t(int number_of_rows) {
      this->number_of_rows = number_of_rows;
      this->rowsort = allocate_rowsort(number_of_rows);	
	  this->reset_rowsort();
  }  

void rowsorter_t::reset_rowsort()
{
	for (int i = 0; i < this->number_of_rows; i++) {
		this->rowsort[i].r = i;
		this->rowsort[i].val = 0;
	}
}

rowsorter_t::~rowsorter_t() {
 deallocate_rowsort(this->rowsort); 
}

dyndata_t::dyndata_t (int N_, int col_) {
        this->N = N_;
        this->col = col_;
        this->rowsort = allocate_rowsort(this->N);
        this->colperm = new_perm_init< colindex_t > (this->N); 
        this->rowsortl = 0;

        for (int i = 0; i < N; i++) {
                this->rowsort[i].r = i;
                this->rowsort[i].val = 0; 
        }
}

/**
* @brief Make a a deep copy of a dyndata_t structure
* @param dd
*/
dyndata_t::dyndata_t (const dyndata_t &dd) {
        this->rowsortl = 0;
        this->rowsort = 0;
        initdata (dd);
}
/**
* @brief Make a a deep copy of a dyndata_t structure
* @param dd
*/
dyndata_t::dyndata_t (const dyndata_t *dd) {
        this->rowsortl = 0;
        this->rowsort = 0;
        initdata (*dd);
}

dyndata_t::~dyndata_t () {
        deallocate_rowsort(this->rowsort);
        delete_perm (colperm);

        deleterowsortl ();
}

void dyndata_t::initdata (const dyndata_t &dd) {
        this->col = dd.col;
        this->N = dd.N;
        this->colperm = new_perm< colindex_t > (this->N);
        copy_perm (dd.colperm, this->colperm, N);

        if (dd.rowsort != 0) {
                this->rowsort = allocate_rowsort(this->N);
                memcpy (this->rowsort, dd.rowsort, N * sizeof (rowsort_t));
        }
        if (dd.rowsortl != 0) {
                this->rowsortl = new_perm< rowindex_t > (this->N); 
                copy_perm (dd.rowsortl, this->rowsortl, this->N);
        }
}

void dyndata_t::copydata (const dyndata_t &dd) {
        if (this->col == dd.col) {
                copy_perm (dd.colperm, this->colperm, N);
        } else {
                this->col = dd.col;
                delete_perm (colperm);
                this->colperm = new_perm< colindex_t > (this->N);
                copy_perm (dd.colperm, this->colperm, N);
        }

        if (this->N == dd.N) {
                if (dd.rowsort != 0) {
                        memcpy (this->rowsort, dd.rowsort, N * sizeof (rowsort_t));
                }
                if (dd.rowsortl != 0) {
                        copy_perm (dd.rowsortl, this->rowsortl, this->N);
                }
        } else {
                this->N = dd.N;
                if (dd.rowsort != 0) {
                        if (this->rowsort != 0) {
                                free (this->rowsort);
                                this->rowsort = 0;
                        }

                        this->rowsort = allocate_rowsort(this->N);
                        memcpy (this->rowsort, dd.rowsort, N * sizeof (rowsort_t));
                }
                if (dd.rowsortl != 0) {
                        deleterowsortl ();
                        this->rowsortl = new_perm< rowindex_t > (this->N); 
                        copy_perm (dd.rowsortl, this->rowsortl, this->N);
                }
        }
}

dyndata_t &dyndata_t::operator= (const dyndata_t &dd) {
        this->copydata (dd);
        return *this;
}

void dyndata_t::show () const {
        myprintf ("dyndata_t: N %d, col %d; ", this->N, this->col);
        myprintf ("  colperm: ");
        print_perm (this->colperm, this->col);
        myprintf ("  rowsort:\n");
        print_rowsort (this->rowsort, this->N);
}

/**
* @brief Convert column permutation of the root to a sorting of the rows (to make the array in root form again)
* @param tmprperm Current rowperm
* @param rperm New permutation
* @param cp Column permutation
*/
inline void cperm_lperms_to_rowsort (const rowperm_t tmprperm, levelperm_t *lperms, rowperm_t rperm,
                                     const colperm_t cp, const arraydata_t *ad) {
        array_t *mults = new array_t[ad->strength + 1];

        /// invert column permutation
        colperm_t cpi = invert_permutation (cp, ad->strength);

        int mult = 1; // mults[strength]=mult;
        for (int l = ad->strength - 1; l >= 0; l--) {
                mults[l] = mult;
                mult *= ad->s[cpi[l]];
        }

        // loop over all rows
        for (rowindex_t k = 0; k < ad->N; k++) {
                rowindex_t v = k;

                /* calculate the value of the permutation, modulo blocks of oaindex */
                rperm[k] = 0;
                v /= ad->oaindex;
                for (int l = ad->strength - 1; l >= 0; l--) {
                        int nl = cp[l];             // new column
                        int vmod = v % (ad->s[nl]); // value in this row
                        v /= ad->s[l];
                        rperm[k] += mults[nl] * lperms[l][vmod];
                }
                /* translate value into new position */
                rperm[k] *= ad->oaindex;
                rperm[k] += (k % ad->oaindex);
        }

        delete_perm (cpi);
        delete[] mults;
}

/**
* @brief Convert a set of level permutations to the corresponding row permutations (for an array in root form)
* @param rperm
* @param lperms
* @param ad
*/
inline void lperms_to_rowsort (rowperm_t rperm, levelperm_t *lperms, const arraydata_t *ad) {
        // loop over all rows
        for (rowindex_t k = 0; k < ad->N; k++) {
                int mult = 1;
                rowindex_t v = k;

                /* calculate the value of the permutation, modulo blocks of oaindex */
                rperm[k] = 0;
                v /= ad->oaindex;
                for (int l = ad->strength - 1; l >= 0; l--) {
                        int vmod = v % (ad->s[l]);
                        rperm[k] += mult * lperms[l][vmod];
                        v /= ad->s[l];
                        mult *= ad->s[l];
                }
                /* translate value into new position */
                rperm[k] *= ad->oaindex;
                rperm[k] += (k % ad->oaindex);
        }
}

/** @brief Helper function for create_root_permutations_index()
*/
void create_root_permutations_index_helper (rowperm_t *rperms, levelperm_t *lperms, const arraydata_t *ad, int level,
                                            int &permcounter) {
        if (level == ad->strength) {
                /* combine level permutations into root row permutations */

                lperms_to_rowsort (rperms[permcounter], lperms, ad);

                permcounter++;
        } else {
                int nlperms = init_perm_n< array_t, int > (lperms[level], ad->s[level]);

                // loop over permutations
                for (int i = 0; i < nlperms; i++) {
                        create_root_permutations_index_helper (rperms, lperms, ad, level + 1, permcounter);
                        next_perm (lperms[level], ad->s[level]);
                }
        }

        return;
}

/**
* @brief Create in index of all root level permutations
* @param ad Parameters of the array
* @param totalpermsr Returns the total number of permutations found
* @return
*/
rowperm_t *create_root_permutations_index (const arraydata_t *ad, int &totalpermsr) {
        /* OPTIMIZE: can this function be made simpler by noting that all permutations of N/oaindex occur ??*/

        totalpermsr = 1;
        for (int i = 0; i < ad->strength; i++)
                totalpermsr *= factorial< int > (ad->s[i]);

        log_print (DEBUG, "create_root_permutations_index: allocating rootrowperms table of size "
                          "%d*%d*sizeof(rowsort_t)=%ld=%.1f MB (sizeof(rowsort_t)=%d)\n",
                   totalpermsr, ad->N, (long)totalpermsr * (long)ad->N * (long)sizeof (rowsort_t),
                   (long)totalpermsr * (long)ad->N * ((double)sizeof (rowsort_t) / (1024 * 1024.0)),
                   sizeof (rowsort_t));

        rowperm_t *rperms = malloc2d< rowindex_t > (totalpermsr, ad->N);
        array_t **lperms = malloc2d_irr< array_t > (ad->strength, (ad->s));

        int totalperms2 = 0;
        create_root_permutations_index_helper (rperms, lperms, ad, 0, totalperms2);

        free2d (lperms);
        return rperms;
}

/**
* @brief Create in index of all root level permutations (both level and column permutations)
* @param ad Parameters of the array
* @param totalpermsr Returns the total number of permutations found
* @return
*/
rowperm_t *create_root_permutations_index_full (const arraydata_t *ad, int &totalpermsr) {
        if (!(ad->colgroupsize[0] >= ad->strength)) {
                if (log_print (NORMAL, "")) {
                        myprintf ("mixed array, not allocating\n");
                        myprintf ("ad->strength %d, ad->colgroupsize[0] %d\n", ad->strength, ad->colgroupsize[0]);
                        ad->show ();
                        ad->show_colgroups ();
                }
                return 0;
        }

        totalpermsr = static_cast< int > (pow (2.0, ad->strength) * factorial< int > (ad->strength));

        log_print (DEBUG, "create_root_permutations_index_full: allocating rootrowperms table of size "
                          "%d*%d*sizeof(rowsort_t)=%ld=%.1f MB (sizeof(rowsort_t)=%d)\n",
                   totalpermsr, ad->N, (long)totalpermsr * (long)ad->N * (long)sizeof (rowsort_t),
                   (long)totalpermsr * (long)ad->N * ((double)sizeof (rowsort_t) / (1024 * 1024.0)),
                   sizeof (rowsort_t));

        rowperm_t *rperms = malloc2d< rowindex_t > (totalpermsr, ad->N);
        array_t **lperms = malloc2d_irr< array_t > (ad->strength, (ad->s));

        rowperm_t tmprperm = new_perm< rowindex_t > (ad->N);

        colperm_t cp = new_perm_init< colindex_t > (ad->strength);
        int permcounter = 0;
        for (int i = 0; i < factorial< int > (ad->strength); i++) {
                // loop over all columns permutations
                for (int j = 0; j < pow (double(2), ad->strength); j++) {
                        // loop over all possible level permutations
                        root_row_permutation_from_index (j, ad, lperms); 

                        cperm_lperms_to_rowsort (tmprperm, lperms, rperms[permcounter], cp, ad);
                        permcounter++;
                }
                next_perm (cp, ad->strength);
        }
        delete_perm (cp);
        delete_perm (tmprperm);

        if (log_print (DEBUG + 1, "")) {
                myprintf ("create_root_permutations_index_j4: strength %d\n", ad->strength);
                for (int k = 0; k < permcounter; k++) {
                        myprintf ("perm %3d: ", k);
                        print_perm (rperms[k], ad->N);
                }
        }

        free2d (lperms); // we do not need to store the lperms
        return rperms;
}

template < class numtype >
/* create permutation that induces a certain combination selection at a starting point */
void create_perm_from_comb2 (numtype *perm, const numtype *comb, int k, int ncols) {
        int *tmpperm = new_perm_init< int > (ncols);
        for (int i = 0; i < k; i++) {
                perm[i] = tmpperm[comb[i]];
                tmpperm[comb[i]] = -1;
        }
        int j = 0;
        for (int i = k; i < ncols; i++) {
                while (tmpperm[j] < 0)
                        j++;
                perm[i] = tmpperm[j];
                j++;
        }
        delete_perm (tmpperm);
}

template <class numtype>
/* create permutation that induces a certain combination selection at a starting point. Note we assume the combination is ordered! */
void create_perm_from_comb ( numtype *perm, const numtype *comb, int k, int ncols, numtype startcol )
{
        // initialize
        init_perm (perm, ncols);

#ifdef OADEBUG
        for (int i = 0; i < (k - 1); i++)
                if (comb[i] >= comb[i + 1])
                        myprintf ("create_perm_from_comb: ERROR!\n");
#endif
        for (int i = 0; i < k; i++)
                std::swap (perm[i + startcol], perm[comb[i] + startcol]);
}

#ifdef OAMEM
inline void static_init_rp (rowperm_t &rp, rowindex_t N) {
        static rowindex_t Np;
        static rowperm_t rowp = 0;

        if (N != Np) {
                Np = N;
                if (rowp != 0)
                        delete_perm (rowp);
                rowp = new_perm< rowindex_t > (N);
        }

        rp = rowp;
}
#endif

void LMCreduction_helper_t::init_root_stage (levelperm_t *&lperm_p, colperm_t *&colperm_p, const arraydata_t *adp) {
        /* permutations buffer */
        lperm_p = this->current_trans->lperms;
        colperm_p = this->colperm_p;

        this->LMC_root_init = 1;
}

/**
* @brief Do the static initialization of structures necessary for LMC test and reduction
*
* The argument is used to check whether the static structures need to be updated or not
*
* @param adp
*/
void LMCreduction_helper_t::freeall () {

        /* clear old structures */
        if (this->current_trans != 0) {
                delete this->current_trans;
                this->current_trans = 0;
        }

        if (this->colperm_p != 0) {
                free2d< colindex_t > (this->colperm_p);
        }
        if (this->localcolperm_p != 0) {
                free2d< colindex_t > (this->localcolperm_p);
        }

        if (this->dyndata_p != 0) {
                if (this->ad != 0) {

                        for (int z = 0; z < this->ad->ncols; z++)
                                delete this->dyndata_p[z];
                        delete[] this->dyndata_p;
                }
        }

        if (this->rootrowperms != 0)
                free2d (this->rootrowperms);
        if (this->rootrowperms_full != 0)
                free2d (this->rootrowperms_full);

        if (this->ad != 0)
                delete this->ad;

        if (this->colbuffer != 0)
                free (this->colbuffer);

        /* make sure all variables are initialized again */
        this->LMC_non_root_init = 0;
        this->LMC_root_init = 0;
        this->LMC_root_rowperms_init = 0;
}

void LMCreduction_helper_t::init (const arraydata_t *adp) {
        log_print (DEBUG, "static_init_LMC: adp->ncols %d, adp->strength %d\n", adp->ncols, adp->strength);

        /* clear old structures */
        this->freeall ();

        /* allocate new structures */
        this->ad = new arraydata_t (*adp);
        this->current_trans = new array_transformation_t (adp);
        this->colperm_p = malloc2d< colindex_t > (this->ad->ncols, this->ad->ncols);
        this->localcolperm_p = malloc2d< colindex_t > (this->ad->ncols, this->ad->ncols);

/* level permutations of the root structure */
#ifdef OAMEM
/* for memory efficient version of the algorithm do not create a root_permutations index */
#else
        this->rootrowperms = create_root_permutations_index (this->ad, this->nrootrowperms);
        this->rootrowperms_full = create_root_permutations_index_full (this->ad, this->nrootrowperms_full);

#endif

        this->dyndata_p = new dyndata_t *[adp->ncols];
        for (int z = 0; z < adp->ncols; z++)
                this->dyndata_p[z] = new dyndata_t (adp->N);

        this->colbuffer = (array_t *)malloc (sizeof (array_t) * adp->N);

        /* make sure all variables are initialized again */
        this->LMC_non_root_init = 0;
        this->LMC_root_init = 0;
        this->LMC_root_rowperms_init = 0;
}

int LMCreduction_helper_t::needUpdate (const arraydata_t *adp) const {
        int update = 0;
        if (this->ad == 0)
                update = 1;
        else {
                if (!(*(this->ad) == *adp)) {
                        update = 1;
                }
        }
        return update;
}

int LMCreduction_helper_t::update (const arraydata_t *adp) {
        int update = this->needUpdate (adp);

        if (update) {
                this->init (adp);
        }
        return update;
}

void LMCreduction_helper_t::init_nonroot_stage (levelperm_t *&lperm_p, colperm_t *&colperm_p, colperm_t *&localcolperm_p,
                                              dyndata_t **&dynd_p, int &dynd_p_nelem, array_t *&colbuffer,
                                              const arraydata_t *adp) const {
        int updatevars = 1; // always update variables since there is a second static init function

        // we assume the static initialization has already been done

        // updatevars=1;
        if (updatevars != 0) {
                /* permutations buffer */
                lperm_p = this->current_trans->lperms;
                colperm_p = this->colperm_p;
                localcolperm_p = this->localcolperm_p;

                /* dynamic_data buffer */
                dynd_p = this->dyndata_p;
                dynd_p_nelem = this->ad->ncols;

                colbuffer = this->colbuffer;
        }
}

/**
 * @brief Returns the value of (part of) a row
 * @param array
 * @param start_idx
 * @param end_idx
 * @param n_rows
 * @param index Value index for each of the columns of the array
 * @return
 */
static inline int row_rank_partial (const carray_t *array, rowindex_t n_rows, const vindex_t *index,
                                    colindex_t start_idx, colindex_t end_idx, const colperm_t &colperm,
                                    rowindex_t row) {
        int sum = 0;
        const array_t *ar = array + row;
        for (colindex_t i = start_idx; i < end_idx; i++) {
                sum += index[i] * ar[n_rows * colperm[i]];
        }
        return sum;
}

/**
* @brief Sort an array based on the values of the root
* @param array
* @param ad
* @param dyndata Contains information about current transformation of the array
* @param rowsort
*/
inline void LMC_root_sort (const carray_t *array, const arraydata_t *ad, const dyndata_t *dyndata) {
        vindex_t *valueindex = new_valueindex< vindex_t > (ad->s, ad->strength); // OPTIMIZE: make static allocation?

        rowsort_t *rowsort = dyndata->rowsort;
        /* perform initial sort, after the initial sort we can sort using pre-calculated structures */
        for (rowindex_t j = 0; j < ad->N; j++) {
                rowsort[j].val =
                    row_rank_partial (array, ad->N, valueindex, (colindex_t)0, ad->strength, dyndata->colperm, j);
        }
        delete[] valueindex;

        std::stable_sort (rowsort, rowsort + ad->N);
}

/**
 * @brief Returns the value of (part of) a row
 * @param array
 * @param start_idx
 * @param end_idx
 * @param n_rows
 * @param index Value index for each of the columns of the array
 * @return
 */
static inline array_t row_rank_partial (carray_t *array, colindex_t start_idx, colindex_t end_idx, rowindex_t n_rows,
                                        const vindex_t *index) {
        array_t sum = 0;
        int j = 0;
        j += n_rows * start_idx;
        for (colindex_t i = start_idx; i < end_idx; i++) {
                sum += index[i] * array[j];
                j += n_rows;
        }
        return sum;
}

/**
* @brief Sort an array based on the values of the root
* @param array
* @param ad
* @param rowsort
*/
inline void LMC_root_sort (carray_t *array, const arraydata_t *ad, rowsort_t *rowsort) {
        int *valueindex = new_valueindex< int > (ad->s, ad->strength); // OPTIMIZE: make static allocation?

        /* perform initial sort, after the initial sort we can sort using pre-calculated structures */
        for (rowindex_t j = 0; j < ad->N; j++) {
                rowsort[j].val = row_rank_partial (array + rowsort[j].r, 0, ad->strength, ad->N, valueindex);
        }
        delete[] valueindex;

        std::sort (rowsort, rowsort + ad->N);
}

/** @brief Static initialization of level permutations
 */
inline void static_init_lperms (const arraydata_t *adp, levelperm_t *&lperm_p, LMCreduction_helper_t &tmpStatic) {
        /* no static update, we assume this has been done already */
        tmpStatic.update (adp);

        lperm_p = tmpStatic.current_trans->lperms;
}

/** @brief Static initialization of root row permutations
*/
inline void static_init_rootrowperms (const arraydata_t *adp, int &totalperms, rowperm_t *&rootrowperms,
                                      levelperm_t *&lperm_p, LMCreduction_helper_t &tmpStatic) {
        /* no static update, we assume this has been done already */

        totalperms = tmpStatic.nrootrowperms;
        rootrowperms = tmpStatic.rootrowperms;
        lperm_p = tmpStatic.current_trans->lperms;

        tmpStatic.LMC_root_rowperms_init = 1;
}

/// show array with rowsort and colperm
void show_array (carray_t *array, const int ncols, const int nrows, colperm_t colperm, rowsort_t *rs) {
        int count;
        for (int j = 0; j < nrows; j++) {
                count = j;
                for (int k = 0; k < ncols; k++) {
                        const char *s = (k < ncols - 1) ? " " : "\n";
                        int val = array[nrows * k + rs[j].r];
                        myprintf ("%3i%s", static_cast< int > (val), s);
                        count += nrows;
                }
        }
}

/**
* @brief Perform root level permutations in LMC check (with full root permutation group)
*
* This method assumes J4
*
* @param original
* @param array
* @param ad Static array information
* @param dyndata Dynamic array information
* @return
* @see LMC
*/
lmc_t LMCreduce_root_level_perm_full (carray_t const *original, const arraydata_t *ad, const dyndata_t *dyndata,
                                      LMCreduction_t *reduction, const OAextend &oaextend,
                                      LMCreduction_helper_t &tmpStatic) {
        lmc_t ret = LMC_EQUAL;
        rowsort_t *rowsort = dyndata->rowsort;

        /* perform initial sorting of the root */
        LMC_root_sort (original, ad, dyndata);

        int totalperms = 0;
        rowperm_t *rootrowperms = 0;
        levelperm_t *lperm_p = 0;
        tmpStatic.init_rootrowperms_full (totalperms, rootrowperms, lperm_p);
        assert (rootrowperms);

        /* perform fake level permutations */
        dyndata_t dyndatatmp = dyndata_t (dyndata);

        for (int l = 0; l < totalperms; l++) {
                // update sort structure
                for (rowindex_t k = 0; k < ad->N; k++) {
                        dyndatatmp.rowsort[rootrowperms[l][k]].r = rowsort[k].r;
                }

                /* pass to non-root stage */
                ret = LMCreduce_non_root_j4 (original, ad, &dyndatatmp, reduction, oaextend, tmpStatic);

                if (ret == LMC_LESS) {
                        break;
                }
        }

        return ret;
}

/**
* @brief Perform root level permutations in LMC check
* @param original
* @param array
* @param ad Static array information
* @param dyndata Dynamic array information
* @return
* @see LMC
*/
lmc_t LMCreduce_root_level_perm (array_t const *original, const arraydata_t *ad, const dyndata_t *dyndata,
                                 LMCreduction_t *reduction, const OAextend &oaextend, LMCreduction_helper_t &tmpStatic) {

        lmc_t ret = LMC_EQUAL;
        rowsort_t *rowsort = dyndata->rowsort;

        /* perform initial sorting of the root */
        LMC_root_sort (original, ad, dyndata);

        int totalperms = 0;
        rowperm_t *rootrowperms = 0; // pointer to root row permutations
        levelperm_t *lperm_p = 0;
        static_init_rootrowperms (ad, totalperms, rootrowperms, lperm_p, tmpStatic);

        /* perform fake level permutations */
        dyndata_t dyndatatmp = dyndata_t (dyndata);

        for (int l = 0; l < totalperms; l++) { // loop over root permutations (in levels)
                // update sort structure

                for (rowindex_t k = 0; k < ad->N; k++) {
                        // 	this assumes that after the LMC_root_sort the root is already in blocks
                        dyndatatmp.rowsort[rootrowperms[l][k]].r = rowsort[k].r;
                }

                // update level permutations structure
                if (reduction->mode >= OA_REDUCE) {
                        root_row_permutation_from_index (l, ad, lperm_p);
                }

                /* pass to non-root stage */
                if (oaextend.getAlgorithm () == MODE_LMC_2LEVEL ||
                    (oaextend.getAlgorithm () == MODE_J5ORDERX && reduction->sd != 0) ||
                    (oaextend.getAlgorithm () == MODE_J5ORDER_2LEVEL && reduction->sd != 0)) {
                        dyndatatmp.initrowsortl ();
                        ret = LMCreduce_non_root_2level (original, ad, &dyndatatmp, reduction, oaextend,
                                                         tmpStatic); 
                } else
                        ret = LMCreduce_non_root (original, ad, &dyndatatmp, reduction, oaextend,
                                                  tmpStatic); 

                if (ret == LMC_LESS) {
                        break;
                }
        }

        return ret;
}

#ifdef OAMEM
lmc_t LMCreduce_root_level_perm_ME (carray_t const *original, const arraydata_t *ad, const dyndata_t *dyndata,
                                    LMCreduction_t *reduction) {
        lmc_t ret = LMC_EQUAL;
        rowsort_t *rowsort = dyndata->rowsort;

        levelperm_t *lperm_p = 0;
        static_init_lperms (ad, lperm_p);

        logstream (DEBUG + 1) << "LMCreduce_root_level_perm_ME: col " << dyndata->col << endl;
        if (dyndata->col == 0) {
                /* perform initial sorting of the root */
                LMC_root_sort (original, ad, dyndata);

                /* initialize level permutations */
                for (int x = 0; x < ad->strength; x++) {
                        init_perm< array_t > (lperm_p[x], ad->s[x]);
                }
        }

        if (dyndata->col == ad->strength) {
                /* pass to non-root stage */

                dyndata_t dyndatatmp = dyndata_t (dyndata); // OPTIMIZE: make static allocation
                rowperm_t rp;
                static_init_rp (rp, ad->N);

                lperms_to_rowsort (rp, lperm_p, ad);
                for (rowindex_t k = 0; k < ad->N; k++)
                        dyndatatmp.rowsort[rp[k]].r = dyndata->rowsort[k].r;

                /* pass to non-root stage */
                ret = LMCreduce_non_root (original, ad, &dyndatatmp, reduction);

                return ret;
        }

        colindex_t col = dyndata->col;

        /* perform level permutations */
        bool np = true;
        unsigned long nperms = 0;
        while (np) {
                dyndata_t *dyndatatmp = new dyndata_t (dyndata);
                dyndatatmp->col++;
                ret = LMCreduce_root_level_perm_ME (original, ad, dyndatatmp, reduction);
                delete dyndatatmp;

                if (ret == LMC_LESS) {
                        break;
                }
                np = next_permutation (lperm_p[col], lperm_p[col] + ad->s[col]);

                if (col == 0 && (nperms % 10000 == 0 || nperms < 20)) {
                        logstream (DEBUG) << "LMCreduce_root_level_perm_ME: doing permutation " << nperms
                                          << " (of first col)" << endl;
                }
                nperms++;
        }
        return ret;
}

#endif

vector< int > LMCcheckLex (arraylist_t const &list, arraydata_t const &ad, int verbose) {
        vector< int > check (list.size ());

        for (unsigned int x = 0; x < list.size (); x++) {
                if (verbose)
                        myprintf ("LMCcheck: array %d\n", x);
                check[x] = LMCcheckLex (list[x], ad);
        }
        return check;
}


/** @brief Calculate J-characteristic for a column combination
*
* We assume the array has values 0 and 1
*/
int fastj3 (carray_t *array, rowindex_t N, const int J, const colindex_t *pp) {
        int jval = 0;

        for (rowindex_t r = 0; r < N; r++) {
                array_t tmp = 0;
                for (int i = 0; i < J; i++) {
                        tmp += array[r + pp[i] * N];
                }
                tmp %= 2;
                jval += tmp;
        }
        jval = 2 * jval - N;
        return (jval);
}

// forward declaration
lmc_t LMCreduction(array_t const *original, array_t const *array, const arraydata_t *ad, const dyndata_t *dyndata,
	LMCreduction_t *reduction, const OAextend &oaextend);

lmc_t LMCcheckj4 (array_link const &al, arraydata_t const &adin, LMCreduction_t &reduction, const OAextend &oaextend,
                  int jj) {
        LMCreduction_helper_t &tmpStatic = * (acquire_LMCreduction_object () );


        const int maxjj = 40;
        assert (jj < maxjj);

        if (reduction.init_state == INIT_STATE_INVALID) {
                myprintf ("LMCreduce: reduction.init_state is INIT_STATE_INVALID, please set it to INIT_STATE::COPY "
                          "or INIT!\n");
                fflush (stdout);
                std::cout.flush ();
                throw_runtime_exception("LMCreduce: reduction.init_state is INIT_STATE_INVALID");
        }
        if (reduction.init_state == COPY) {
                reduction.setArray (al);
        }

        assert (adin.ncolgroups == 1);
        assert (adin.s[0] == 2);
        assert (jj >= adin.strength);
        if (adin.strength != 3) {
                throw_runtime_exception ("LMCcheckj4: error: strength not equal to 3\n");
        }

        if (al.n_columns < jj) {
                // fall back to original method
                dyndata_t dyndata (adin.N);

                lmc_t ret = LMCreduction (al.array, al.array, &adin, &dyndata, &reduction, oaextend);
                return ret;
        }

        arraydata_t ad (adin);
        ad.complete_arraydata_splitn (jj); // split column group!

        tmpStatic.update (&ad);

        int nc = ncombs (ad.ncols, jj);
        colindex_t *comb = new_comb_init< colindex_t > (jj);
        int nrootcombs = ncombs (jj, ad.strength);
        colindex_t *combroot = new_comb_init< colindex_t > (ad.strength); //
        colindex_t combx[maxjj]; // final combination for root selection + extra columns
        colindex_t comby[maxjj];

        carray_t *array = al.array;

        lmc_t ret = LMC_EQUAL;

        /* loop over all possible column combinations for the first jj columns */
        int goodj = abs (jvaluefast (array, ad.N, jj, comb));

        /* initialize variables outside of loop */
        colperm_t perm = new_perm< colindex_t > (ad.ncols);
        colperm_t pp = new_perm_init< colindex_t > (jj);

        for (int i = 0; i < nc; i++) {
                int jval = abs (jvaluefast (array, ad.N, jj, comb));

                if (jval > goodj) {
                        // this combination will lead to LMC less
                        if (reduction.mode >= OA_REDUCE) {
                                ret = LMC_EQUAL;
                                goodj = jval; // update state

                        } else {
                                ret = LMC_LESS;
                                break;
                        }
                } else if (jval < goodj) {
                        // this combination can only lead to LMC more

                        next_combination (comb, jj, ad.ncols); // increase combination
                        continue;
                }


                init_comb (combroot, ad.strength, jj);
                for (int kk = 0; kk < nrootcombs; kk++) {
                        // prepare for all root combinations

                        // convert combination to full permutation
                        create_perm_from_comb (combx, combroot, ad.strength, jj, 0);
                        perform_inv_perm (comb, comby, jj, combx);

                        init_perm (perm, ad.ncols);
                        create_perm_from_comb2< colindex_t > (perm, comby, jj, ad.ncols);

                        dyndata_t dyndata (ad.N);
                        dyndata.col = ad.strength; // set at col after root
                        copy_perm (perm, dyndata.colperm, ad.ncols);

                        ret = LMCreduce_root_level_perm_full (array, &ad, &dyndata, &reduction, oaextend, tmpStatic);

                        if (reduction.mode >= OA_REDUCE) {
                                if (ret == LMC_LESS)
                                        reduction.state = REDUCTION_CHANGED;
                                ret = LMC_EQUAL;

                        } else {
                                if (ret == LMC_LESS) {
                                        break;
                                }
                        }
                        next_combination (combroot, ad.strength, jj); // increase combination
                }
                if (ret == LMC_LESS) {
                        break;
                }

                next_combination (comb, jj, ad.ncols); // increase combination
        }
        delete_perm (perm);
        delete_perm (pp);
        delete_comb (combroot);
        delete_comb (comb);

	release_LMCreduction_object(&tmpStatic);

        return ret;
}

/// create splits for a symmetry group
std::vector< int > symmetrygroup2splits (const symmetry_group &sg, int ncols, int verbose = 0) {
        int jj = sg.n;

        std::vector< int > splits;

        if (verbose >= 2) {
                myprintf ("symmetry indices: %d: ", sg.ngroups);
                display_vector (sg.gstart);
                myprintf ("\n");
        }
        for (int j = 0; j < sg.ngroups; j++) {
                splits.push_back (sg.gstart[j]);
        }
        if (ncols > jj) {
                splits.push_back (jj);
        }
        return splits;
}

void _calculate_j4_values(std::vector<int> &j4_values, carray_t *array, const int N, const colperm_t comb)
{
	colindex_t lc[4];
	init_perm(lc, 4);
	colindex_t lc2[4];
	init_perm(lc2, 4);
	for (size_t i = 0; i < 5; i++) {
		perform_inv_perm(comb, lc2, 4, lc);

		j4_values[i] = abs(jvaluefast(array, N, 4, lc2));
		next_comb(lc, 4, 5);
	}
	// we reverse the values (to keep ordering matched to the original orderings)
	std::reverse(j4_values.begin(), j4_values.end());
}
/// Perform check or reduction using ordering based on J5 and delete-one-factor J4 values
int jj45split (carray_t *array, rowindex_t N, int jj, const colperm_t comb, const arraydata_t &ad,
               const OAextend &oaextend, LMCreduction_helper_t &tmpStatic, LMCreduction_t &reduction, int verbose = 0) {
        myassert (jj == 5, "jj should be equal to 5");
        lmc_t ret;

        const int maxjj = 40;
        colindex_t comby[maxjj];
        colindex_t perm[maxjj];
        if (verbose)
                myprintf ("-- jj45split --\n");

        // allocate buffer to hold the values
        std::vector< int > j4_values (5);

        /* calculate the J4 values */
		_calculate_j4_values(j4_values, array, N, comb);

        indexsort s (5);
        s.sort (j4_values);
        std::vector< int > ordered_j4_values = s.sorted (j4_values);

        // calculate symmetry group of column permutations in the first jj columns
        symmetry_group sg (ordered_j4_values, true, verbose >= 3);
        if (verbose >= 2) {
                myprintf ("jj45split: values sorted ");
                display_vector (ordered_j4_values);
                myprintf ("\n  ");
                sg.show ();
        }

        // create the sorted column combination
        perform_inv_perm (comb, comby, jj, s.indices);
        if (verbose >= 2) {
                myprintf ("initial comb: ");
                print_perm (comb, jj);
                myprintf (" result comb: ");
                print_perm (comby, jj);
        }
        // make splits
        init_perm (perm, ad.ncols);
        create_perm_from_comb2< colindex_t > (perm, comby, jj, ad.ncols);

        if (verbose >= 2) {
                myprintf ("perm: ");
                print_perm (perm, ad.ncols);
        }

        dyndata_t dyndata (ad.N);
        dyndata.col = 0; 
        dyndata.setColperm (perm, ad.ncols);

        if (oaextend.getAlgorithm () == MODE_J5ORDERX || oaextend.getAlgorithm () == MODE_J5ORDER_2LEVEL) {
                dyndata.initrowsortl ();
        }

        arraydata_t adfix = ad;

        if (verbose >= 4) {
                adfix.show_colgroups ();
        }

        lmc_t ret1;

        std::vector< int > splits = symmetrygroup2splits (sg, ad.ncols, verbose);

        adfix.set_colgroups (splits);

        if (verbose >= 2) {
                adfix.show_colgroups ();
        }

        ret = LMCreduction (array, array, &adfix, &dyndata, &reduction, oaextend);

        if (verbose)
                myprintf ("  ret %d\n", ret);

        ret1 = ret;

        double val = ret1;
        return val;
}

/** Convert J5 and tuple of J4 values to a single number
 *
 * \param ww Vector with first element J5 and then the J4 values with the combination (2,3,4,5) first and (1,2,3,4) last.
 */
inline double jj452double (const double *ww) {
        // the maximum value of the J4 characteristics is N. we create a unique value out of the pair by multiplication
        double val = 0;
        for (size_t i = 0; i < 6; i++)
                val = val * 64 + ww[i];
        return val;
}

/// return true if target is in root form, otherwise return false
inline bool check_root_form(const array_t *array, const arraydata_t &arrayclass) {
	array_t *root = create_array(arrayclass.N, arrayclass.strength);
	create_root(root, &arrayclass);
	if (std::equal(array, array + arrayclass.N * arrayclass.strength, root)) {
		destroy_array(root);
		return true;
	}
	else {
		destroy_array(root);
		return false;
	}
}

/// return 0 if target is equal to original, otherwise return 1 and copy root initialization + 1
inline int check_root_update(carray_t *original, const arraydata_t &arrayclass, array_t *target) {
	int changed = 0;

	array_t *root = create_array(arrayclass.N, arrayclass.strength);
	create_root(root, &arrayclass);
	if (!std::equal(original, original + arrayclass.N * arrayclass.strength, root)) {
		copy_array(root, target, arrayclass.N, arrayclass.strength);
		if (arrayclass.ncols>arrayclass.strength) {
		  for (int j = 0; j < arrayclass.N; j++)
			  target[arrayclass.N * arrayclass.strength + j] = arrayclass.s[arrayclass.strength] + 100;
	}
		changed = 1;
	}
	destroy_array(root);

	return changed;
}

/// helper function to calculate J-values
inline int fastJupdateValue(rowindex_t N, carray_t *tmpval) {
	int jval = 0;
	for (rowindex_t r = 0; r < N; r++) {
		jval += tmpval[r] % 2;
	}
	jval = 2 * jval - N;
	return (jval);
}

/** return value of subarray based on J4-J5 ordering
 *
 * \param array Pointer to array
 * \param N Number of rows
 * \param comb Combination of 5 columns to use in calculation
 * \param dosort If True, then sort the tuple of J4 values before calculation of the value
 */
jj45_t jj45val (carray_t *array, rowindex_t N, const colperm_t comb, int dosort = 1) {
	const int jj = 5;
        double ww[6];

        array_t tmpval[MAXROWS];

        std::fill (tmpval, tmpval + N, 0);
        fastJupdate (array, N, jj, comb, tmpval);
        ww[0] = abs (fastJupdateValue (N, tmpval));

        for (size_t i = 0; i < 5; i++) {
                int ii = 5 - i - 1;
                fastJupdate (array, N, 1, comb + ii, tmpval);
                ww[i + 1] = abs (fastJupdateValue (N, tmpval));
                fastJupdate (array, N, 1, comb + ii, tmpval); 
        }

        if (dosort) {
                std::sort (ww + 1, ww + 6, std::greater< int > ());
        }
        double val = jj452double (ww);

        return val;
}


lmc_t LMCcheckj5 (array_link const &al, arraydata_t const &adin, LMCreduction_t &reduction, const OAextend &oaextend) {
        const int dverbose = 0;
        LMCreduction_helper_t &tmpStatic = reduction.getReferenceReductionHelper ();

        const int jj = 5;
        const int maxjj = 10;
        myassert (jj <= maxjj, "jj is larger than maximum allowed value");

        // perform basic integrity checks
        myassert (adin.s[0] == 2, "array should be 2-level array");
        myassert (jj >= adin.strength, "jj should be at least equal to the strength");
        myassert (adin.strength >= 3, "LMCcheckj5: strength should be >=3!");

        if (reduction.mode >= OA_REDUCE) {
                throw_runtime_exception ("error: reduction mode not implemented for J5 ordering!!!\n");
        }

        if (reduction.init_state == COPY) {
                reduction.setArray (al);
                reduction.init_state = COPY;
        }
        if (reduction.init_state == INIT_STATE_INVALID) {
                throw_runtime_exception ("LMCcheckj5: reduction.init_state is INVALID\n");
        }

        // NOTE: the reduction.array is reduced, so reduction.array is a good value to start from.
        // NOTE: for arrays not in root form we set update the reduction.array to make sure that the root is set (this
        // part is not checked!)
        bool rootform = check_root_form (reduction.array, adin);

        if (!rootform) {
                log_print (
                    SYSTEM,
                    "error: LMC test or LMC reduction for arrays not in root form needs special initialization !\n");
                print_array (reduction.array, adin.N, adin.ncols);
        }

        if (al.n_columns < jj) {
                // fall back to original method
                dyndata_t dyndata (adin.N);
                lmc_t ret = LMCreduction (al.array, al.array, &adin, &dyndata, &reduction, oaextend);
                return ret;
        }

        arraydata_t ad (adin);
        // split column group!
        // the first 5 columns are selected in this function (and if necessary we loop over all permutations). the
        // remaining columns are then processed with the original LMC algorithm
        ad.complete_arraydata_splitn (jj);

        tmpStatic.update (&ad);

        colindex_t *firstcolcomb = new_comb_init< colindex_t > (jj);
        int nrootcombs = ncombs (jj, ad.strength);
        colindex_t *combroot = new_comb_init< colindex_t > (ad.strength); //
        colindex_t combx[maxjj]; // final combination for root selection + extra columns
        colindex_t comby[maxjj];

        carray_t *array = al.array;

        lmc_t ret = LMC_EQUAL;

        /* loop over all possible column combinations for the first jj columns */
        int jbase = abs (jvaluefast (array, ad.N, jj, firstcolcomb));
        jj45_t j54_base = jj45val (array, ad.N, firstcolcomb, 0);

        /* initialize variables outside of loop */
        colperm_t perm = new_perm< colindex_t > (ad.ncols);
        colperm_t pp = new_perm_init< colindex_t > (jj);

        const int orig = oaextend.j5structure == J5_ORIGINAL;

        int ncolsfirst = adin.colgroupsize[0];
        int nc = ncombs (ncolsfirst, jj);     

        if (dverbose) {
                myprintf ("LMCcheckj5: selected ncolsfirst %d, nc %d, jbase %d, wbase %f\n", ncolsfirst, nc, jbase,
                          j54_base);
        }

        // loop over all possible combinations for the first jj columns
        for (int i = 0; i < nc; i++) {

                // determine the J-characteristic for the selected combination
                int j5val = abs (jvaluefast (array, ad.N, jj, firstcolcomb));

                if (dverbose >= 2) {
                        myprintf ("   LMCcheckj5 column comb %d/%d (jj=%d): j5val %d, jbase %d\n", i, nc, jj, j5val,
                                  jbase);
                }

                // check order based on J5
                if (j5val ORDER_J5_SMALLER jbase) {
                        ret = LMC_LESS;
                        reduction.updateLastCol (4);

                        if (dverbose) {
                                printfd ("LMCcheckj5: ret %d: j5val %d, jbase %d\n", ret, j5val, jbase);
                        }
                        break;
                } else if (j5val ORDER_J5_GREATER jbase) {
                        // this combination can only lead to LMC more
                        next_combination (firstcolcomb, jj, ad.ncols); // increase combination
                        continue;
                }

                if (dverbose >= 2) {
                        myprintf ("   LMCcheckj5 column comb %d/%d : going to jjsplit\n", i, nc);
                }

                const int dverbose = 0;
                // check order based on J4-J5 value
                jj45_t j54 = jj45val (array, ad.N, firstcolcomb);
                if (j54 ORDER_J45_SMALLER j54_base) {
                        ret = LMC_LESS;
                        if (dverbose)
                                printfd ("LMCcheckj5: ret %d: w %ld, wbase %ld\n", (int)ret, (long)j54,
                                            (long)j54_base);
                        reduction.updateLastCol (4);
                        break;
                } else if (j54 ORDER_J45_GREATER j54_base) {
                        // this combination can only lead to LMC more
                        next_combination (firstcolcomb, jj, ad.ncols); // increase combination
                        continue;
                }

                if (dverbose >= 2) {
                        myprintf ("   LMCcheckj5 column comb %d/%d : going to jjsplit\n", i, nc);
                }

                if (orig) {
                        // do classic check
                        init_comb (combroot, ad.strength, jj);
                        for (int kk = 0; kk < nrootcombs; kk++) {
                                // prepare for all root combinations

                                // convert combination to full permutation
                                create_perm_from_comb (combx, combroot, ad.strength, jj, 0);
                                perform_inv_perm (firstcolcomb, comby, jj, combx);
                                init_perm (perm, ad.ncols);
                                create_perm_from_comb2< colindex_t > (perm, comby, jj, ad.ncols);

                                dyndata_t dyndata (ad.N);
                                dyndata.col = ad.strength;
                                copy_perm (perm, dyndata.colperm, ad.ncols);

                                ret = LMCreduce_root_level_perm_full (array, &ad, &dyndata, &reduction, oaextend,
                                                                      tmpStatic);

                                if (ret == LMC_LESS) {
                                        printfd ("LMCcheckj5: note: code unchecked...\n");
                                        break;
                                }

                                next_combination (combroot, ad.strength, jj); 
                        }
                }
                int retx = ret;
                if (!orig) {

                        retx = jj45split (array, ad.N, jj, firstcolcomb, ad, oaextend, tmpStatic, reduction);

                }
                ret = (lmc_t) (int)retx;

                if (ret == LMC_LESS) {
                        break;
                }

                next_combination (firstcolcomb, jj, ncolsfirst); 
        }
        delete_perm (perm);
        delete_perm (pp);
        delete_comb (combroot);
        delete_comb (firstcolcomb);

        return ret;
}

lmc_t LMCcheckLex (array_link const &al, arraydata_t const &ad) {
        int N = al.n_rows;

        if (al.n_columns != ad.ncols) {
                myprintf ("error: number of columns of matrix inconsistent!\n");
                return LMC_NONSENSE;
        }

        dyndata_t dynd = dyndata_t (N);
        lmc_t result = LMC_NONSENSE;
        arraydata_t adx = arraydata_t (&ad, al.n_columns);
        LMCreduction_t *reduction = new LMCreduction_t (&adx);
        OAextend oaextend;
        result = LMCreduction (al.array, al.array, &ad, &dynd, reduction, oaextend);

        return result;
}

template < class numtype >
/// Convert selection of elements to extended permutation
void combadd2perm (const larray< numtype > &comb, int newidx, int n, larray< numtype > &target,
                   larray< numtype > &wtmp) {
        for (int i = 0; i < n; i++) {
                wtmp[i] = i;
        }
        int j = comb.size ();
        for (int i = 0; i < j; i++) {
                target[i] = comb[i];
                wtmp[comb[i]] = -1;
        }
        if (j < 0)
                myprintf ("combadd2perm: j<0! j %d\n", j);
        target[j] = newidx;
        wtmp[newidx] = -1;
        j++;

        for (int i = 0; i < n; i++) {
                if (wtmp[i] >= 0) {
                        target[j] = wtmp[i];
                        j++;
                }
        }
}

template < class numtype >
/// Convert selection of elements to extended permutation
std::vector< numtype > comb2perm(const std::vector< numtype > comb, int n) {
	std::vector< numtype > w(n);
	std::copy(comb.begin(), comb.end(), w.begin());
	numtype j = comb.size();
	for (int i = 0; i < n; i++) {
		bool c = std::find(comb.begin(), comb.end(), i) != comb.end();
		if (!c) {
			w[j] = i;
			j++;
		}
	}
	return w;
}

lmc_t LMCcheck(const array_link &al) {
	myassert(al.is_orthogonal_array(), "input array is not an orthogonal array");
	arraydata_t arrayclass = arraylink2arraydata(al);

	OAextend oaextend(arrayclass);
	LMCreduction_t reduction(&arrayclass);
	reduction.init_state = COPY;
	return LMCcheck(al, arrayclass, oaextend, reduction);
}

lmc_t LMCcheckOriginal (const array_link &al) {
        myassert (al.is2level (), "input array is not a 2-level array");
        arraydata_t ad = arraylink2arraydata (al);
        int strength = al.strength ();

        OAextend oaextend (ad);
        LMCreduction_t reduction (&ad);
        reduction.init_state = INIT;
        copy_array (al.array, reduction.array, ad.N, ad.ncols);

        int changed = check_root_update (al.array, ad, reduction.array);

        return LMCcheck (al.array, ad, oaextend, reduction);
}

lmc_t LMCcheck (const array_link &al, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction) {
        myassert (ad.N == al.n_rows, "LMCcheck: wrong number of rows");
        mycheck (ad.ncols <= al.n_columns, "LMCcheck: wrong number of columns al %d, adata %d", al.n_columns,
                 ad.ncols);

        return LMCcheck (al.array, ad, oaextend, reduction);
}

lmc_t LMCcheck (const array_t *array, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction) {
        lmc_t lmc = LMC_NONSENSE;
        dyndata_t dynd = dyndata_t (ad.N);

        if (reduction.init_state == INIT_STATE_INVALID) {
                throw_runtime_exception("LMCcheck: reduction.init_state is INVALID, please set it to COPY or INIT");
        }
        if (reduction.init_state == SETROOT) {
                log_print (DEBUG, "LMCcheck: reduction.init_state is SETROOT\n");

                array_link alp = ad.create_root (ad.ncols, 1000);
                reduction.setArray (alp);
                reduction.init_state = INIT;
        }
        if (reduction.init_state == COPY) {
                reduction.setArray (array, ad.N, ad.ncols);
        }

        switch (oaextend.getAlgorithm ()) {
        case MODE_ORIGINAL: {
                lmc = LMCreduction (array, array, &ad, &dynd, &reduction, oaextend);
        } break;
        case MODE_J4: {
                array_link al (array, ad.N, ad.ncols, -20);

                lmc = LMCcheckj4 (al, ad, reduction, oaextend);
        } break;
        case MODE_J5ORDER: {
                myprintf ("LMCcheck: algorithm MODE_J5ORDER: code path untested\n");
                lmc = LMCreduction (array, array, &ad, &dynd, &reduction, oaextend);
        } break;
        case MODE_LMC_2LEVEL: {
                array_link al (array, ad.N, ad.ncols, -20);
                reduction.updateSDpointer (al, true);
                lmc = LMCreduction (array, array, &ad, &dynd, &reduction, oaextend);
        } break;
        case MODE_LMC_SYMMETRY: {
			    myprintf("MODE_LMC_SYMMETRY not supported any more\n");
				throw_runtime_exception("MODE_LMC_SYMMETRY not supported any more");
        } break;
        case MODE_J5ORDERX: {
                array_link al (array, ad.N, ad.ncols, -20);
                copy_array (array, reduction.array, ad.N, ad.ncols);
                lmc = LMCcheckj5 (al, ad, reduction, oaextend);
        } break;
        case MODE_J5ORDER_2LEVEL: {
                array_link al (array, ad.N, ad.ncols, -20);
                copy_array (array, reduction.array, ad.N, ad.ncols);

                OAextend oaextend2 = oaextend;
                oaextend2.setAlgorithm (MODE_LMC_2LEVEL);
                reduction.updateSDpointer (al, false);

                lmc = LMCcheckj5 (al, ad, reduction, oaextend2);
        } break;

        default: {
                myprintf ("file %s: LMCcheck: line %d: error: no such mode %d\n", __FILE__, __LINE__,
                          oaextend.getAlgorithm ());
                myprintf ("  algorithm list: %s\n", algorithm_t_list ().c_str ());
        } break;
        }

        return lmc;
}

/**
* @brief Perform an LMC reduction in several steps.
* @param original
* @param array
* @param ad
* @param dyndata
* @param reduction
* @return
*/
lmc_t LMCreduction_train (const array_t *original, const arraydata_t *ad, const dyndata_t *dyndata,
                          LMCreduction_t *reduction, const OAextend &oaextend) {
        lmc_t ret = LMC_NONSENSE;
        int n;
        colindex_t *cols;

        switch (ad->ncols) {
        case 15:
                /* tight */
                n = std::max ((int)floor ((double)(ad->ncols - 7) / 2), 1);
                cols = new colindex_t[n];
                for (int z = 0; z < n; z++) {
                        cols[z] = ad->ncols - (n - z - 1);
                }
                break;
        default:
                /* default mode */
                n = std::max ((int)floor ((double)(ad->ncols - 7) / 2), 1);
                cols = new colindex_t[n];
                for (int z = 0; z < n; z++) {
                        cols[z] = ad->ncols - (n - z - 1) * 2;
                }
                break;
        }

        logstream (DEBUG) << "   cols selected: ";
        if (checkloglevel (DEBUG))
                print_perm (cols, n);

        array_t *atmp = create_array (ad);

        LMCreduction_t reductiontmp (ad);
        reduction->mode = OA_REDUCE;
        reductiontmp.mode = OA_REDUCE;

        copy_array (original, reduction->array, ad->N, ad->ncols);
        copy_array (original, reductiontmp.array, ad->N, ad->ncols);
        copy_array (original, atmp, ad->N, ad->ncols);

        const int verbose = 0;

        std::vector< array_link > aldbg (100);

        for (colindex_t z = 0; z < n; z++) {
                copy_array (reductiontmp.array, atmp, ad->N, ad->ncols);

                logstream (NORMAL) << printfstring ("   LMCreduction_train: stage %d/%d", z, n)
                                   << printfstring (" (current column %d, number of columns %d)", cols[z], cols[n - 1])
                                   << endl;

                reductiontmp.maxdepth = cols[z];
                reductiontmp.reset ();
                reductiontmp.transformation->reset ();
                reductiontmp.mode = OA_REDUCE;
                reductiontmp.nred = 0;

                ret = LMCreduction (atmp, atmp, ad, dyndata, &reductiontmp, oaextend);

                if (verbose >= 3) {
                        myprintf ("### LMCreduction_train: concatenate the transformations (z=%d, column %d)\n", z,
                                  cols[z]);
                        myprintf ("  accum -> \n");
                        reduction->transformation->show ();
                        myprintf ("  tmp -> \n");
                        reductiontmp.transformation->show ();
                }

                *(reduction->transformation) = (*(reductiontmp.transformation)) * (*(reduction->transformation));
        }
        /* final iteration */
        copy_array (reductiontmp.array, atmp, ad->N, ad->ncols);
        reduction->maxdepth = -1;
        reductiontmp.transformation->reset ();

        ret = LMCreduction (atmp, atmp, ad, dyndata, &reductiontmp, oaextend);

        if (verbose >= 2) {
                myprintf ("### LMCreduction_train: concatenate the transformations (final round)\n");
                myprintf ("  accum -> \n");
                reduction->transformation->show ();
                myprintf ("  tmp -> \n");
                reductiontmp.transformation->show ();
        }
        *(reduction->transformation) = (*(reductiontmp.transformation)) * (*(reduction->transformation));

        if (verbose >= 2) {
                myprintf ("  final -> \n");
                reduction->transformation->show ();
        }

        copy_array (reductiontmp.array, reduction->array, ad->N, ad->ncols);

        destroy_array (atmp);
        delete[] cols;
        return ret;
}

lmc_t LMCreduction_train(const array_link &al, const arraydata_t *ad, LMCreduction_t *reduction,
	const OAextend &oaextend) {
	dyndata_t dyndata(ad->N, 0);
	return LMCreduction_train(al.array, ad, &dyndata, reduction, oaextend);
}

/// full reduction, no root-trick
lmc_t LMCreduceFull (carray_t *original, const array_t *array, const arraydata_t *adx, const dyndata_t *dyndata,
                     LMCreduction_t *reduction, const OAextend &oaextend, LMCreduction_helper_t &tmpStatic) {
        arraydata_t *ad = new arraydata_t (*adx); 
        ad->oaindex = ad->N;                      // NOTE: this is to prevent processing on blocks in LMC_check_col

        if (dyndata->col == 0) {
                /* update information of static variables */
                tmpStatic.update (ad);
        }

        array_t *oarray = clone_array (original, ad->N, ad->ncols);
        dyndata_t dyndatatmp = dyndata_t (dyndata);
        if (log_print (NORMAL, "")) {
                myprintf ("LMCreduceFull:\n");
                print_array (oarray, ad->N, ad->ncols);
                dyndata->show ();
        }
        copy_array (original, reduction->array, ad->N, ad->ncols);
        lmc_t ret = LMCreduce_non_root (oarray, ad, &dyndatatmp, reduction, oaextend, tmpStatic);

        destroy_array (oarray);
        delete ad;
        return ret;
}

bool is_root_form(const array_link &array, int strength) {
	arraydata_t arrayclass = arraylink2arraydata(array, strength);

	array_link root_array = array.selectFirstColumns(strength);
	root_array.create_root(arrayclass);

	return array.selectFirstColumns(strength) == root_array;
}

/*!
  LMC performs an LMC test or LMC reduction on a given array.

  As input we require the array to be tested to be in root form.

  First we select the necessary column permutations. This is done in steps because we cannot interchange columns
    within a single column group.

  \brief LMC test
  \param original Pointer to the original OA
  \param array Pointer to the array to be checked
  \param ad Pointer to the OA data, these data are static
  \param dyndata Pointer to the dynamic OA data
  \return Value for LMC test
  */
lmc_t LMCreduction (const array_t *original, const array_t *array, const arraydata_t *ad, const dyndata_t *dyndata,
                 LMCreduction_t *reduction, const OAextend &oaextend) {

        LMCreduction_helper_t &tmpStatic = reduction->getReferenceReductionHelper();

        if (dyndata->col == 0) {
                /* update information of static variables */
                tmpStatic.update (ad);

                if (reduction->init_state == INIT_STATE_INVALID) {
                        throw_runtime_exception("LMCreduce: reduction.init_state is INIT_STATE_INVALID");
                }

                if (reduction->init_state == COPY) {
                        reduction->setArray (array, ad->N, ad->ncols);
                }

                bool rootform = check_root_form (reduction->array, *ad);

                if (!rootform) {
						if (checkloglevel(DEBUG) )
	                        printfd ("LMCreduce: WARNING: LMC test or LMC reduction for arrays not in root form needs "
                                 "special initialization! reduction->mode %d (OA_TEST %d, OA_REDUCE %d, OA_REDUCE_PARTIAL %d)\n",
                                 reduction->mode, OA_TEST, OA_REDUCE, OA_REDUCE_PARTIAL);

                        int changed = check_root_update (original, *ad, reduction->array);

                        if (checkloglevel (DEBUG) || 0) {
                                myprintf ("original:\n");
                                print_array (original, ad->N, ad->ncols);
                                myprintf ("reduction:\n");
                                print_array (reduction->array, ad->N, ad->ncols);
                        }
                }
        }

#ifdef OACHECK
        if (ad->ncols <= ad->strength) {
                myprintf ("LMC reduction not defined for arrays with number of columns (%d) less or equal to the "
                          "strength (%d)\n",
                          ad->ncols, ad->strength);
                return LMC_NONSENSE;
        }
#endif

        lmc_t ret = LMC_NONSENSE;

        levelperm_t *lperm_p = 0;
        colperm_t *colperm_p = 0;
        tmpStatic.init_root_stage (lperm_p, colperm_p, ad);
        log_print (DEBUG + 1, "LMC: entering root LMC phase (level %d) \n", dyndata->col);

        /* check whether the root is full */
        if (dyndata->col == ad->strength) {
                if (log_print (DEBUG, "")) {
                        log_print (DEBUG, "## LMC: root column selection complete \n");
                        myprintf ("   perm: ");
                        print_perm (dyndata->colperm, ad->ncols);
                }

#ifdef OAMEM
                dyndata_t *dd = new dyndata_t (dyndata);
                dd->col = 0;
                ret = LMCreduce_root_level_perm_ME (original, ad, dd, reduction);
                delete dd;
#else
                ret = LMCreduce_root_level_perm (original, ad, dyndata, reduction, oaextend, tmpStatic);
#endif

                return ret;
        }

        int col = dyndata->col; // current column
        int cg = ad->get_col_group (col);
        const int ncolsremroot = ad->strength - col;
        const int ncolsremgroup = ad->colgroupsize[cg] + ad->colgroupindex[cg] - col;
        int k = min (ncolsremroot, ncolsremgroup); /* number of columns to select in this stage */

        log_print (DEBUG + 1, "LMC: root stage: active column %d, selecting %d columns\n", col, k);
        /* select k columns from the entire group */

        int nc = ncombs (ncolsremgroup, k);
        colindex_t *comb = new_comb_init< colindex_t > (k);
        colindex_t *localcolperm = new_perm< colindex_t > (ad->ncols);

        array_t *cpy = create_array (ad);

        dyndata_t newdyndata = dyndata_t (dyndata);
        /* loop over all possible column combinations for the root */

        for (int i = 0; i < nc; i++) {
                create_perm_from_comb< colindex_t > (localcolperm, comb, k, ad->ncols, col);

                // loop over permutations of selected columns
                for (int z = 0; z < factorial< int > (k); z++) {

                        // set necessary components of dyndata (these might be changed inside LMCreduce)
                        newdyndata.copydata (*dyndata);

                        newdyndata.col = newdyndata.col + k;

                        perform_inv_perm (dyndata->colperm, colperm_p[col], ad->ncols, localcolperm);
                        copy_perm (colperm_p[col], newdyndata.colperm, ad->ncols);

                        if (log_print (DEBUG, "")) {
                                myprintf ("  colperm root selection ");
                                print_perm (newdyndata.colperm, ad->ncols);
                        }

                        ret = LMCreduction (original, array, ad, &newdyndata, reduction, oaextend);

                        if (ret == LMC_LESS) {
                                log_print (DEBUG + 1, "LMC at level %d: received LMC_less, doing break\n",
                                           dyndata->col);
                                break;
                        }

                        next_perm (localcolperm + col, k);
                }
                if (ret == LMC_LESS) {
                        break;
                }
                next_combination (comb, k, ncolsremgroup); // increase combination
        }

        destroy_array (cpy);
        delete_perm (localcolperm);
        delete_comb (comb);

        return ret;
}


void print_rowsort (rowsort_t *rowsort, int N) {
        myprintf ("row: ");
        for (int x = 0; x < N; x++)
                myprintf ("%2ld ", static_cast< long > (rowsort[x].r));
        myprintf ("\nval: ");
        for (int x = 0; x < N; x++)
                myprintf ("%2ld ", static_cast< long > (rowsort[x].val));
        myprintf ("\n");
}

void print_column_rowsort (const array_t *arraycol, rowsort_t *rowsort, int N) {
        myprintf ("column: \n");
        for (int x = 0; x < N; x++)
                myprintf ("%d: %d (value %ld)\n", x, arraycol[rowsort[x].r], (long)rowsort[x].val);
}

/// Print a column of an array with order determined by rowsort structure
void print_col_sorted (carray_t *array, levelperm_t lperm, rowsort_t *rs, int N) {
        myprintf ("do lperm: ");
        for (int i = 0; i < N; i++)
                myprintf ("%d ", lperm[array[rs[i].r]]);
        myprintf ("\n");
        myprintf ("no lperm: ");
        for (int i = 0; i < N; i++)
                myprintf ("%d ", array[rs[i].r]);
        myprintf ("\n");
}

/// Predict J-characteristic of an array based on rowsort of the root
inline int predictJrowsort (const array_t *array, const int N, const rowsort_t *rs, levelperm_t lperm) {
        int t = N / 4;
        int tt = t / 2;

        // x1 is the number of + entries in (s)
        int x1 = 0;
        for (int i = 0; i < tt; i++) {
                int ii = rs[i].r;
                if (lperm[array[ii]] == 0)
                        x1++;
        }
        for (int i = tt; i < t; i++) {
                int ii = rs[i].r;
                if (lperm[array[ii]] == 1)
                        x1++;
        }

        return 8 * x1 - N;
}

/// select the unique arrays in a list, the original list is sorted in place
void selectUniqueArrays (arraylist_t &input_arrays, arraylist_t &output_arrays, int verbose) {
        sort (input_arrays.begin (), input_arrays.end ()); 
        if (input_arrays.size () > 0) {
                std::vector< int > vv (input_arrays.size ());
                vv[0] = 1;

                for (size_t i = 0; i < input_arrays.size () - 1; i++) {
                        int e = input_arrays[i] == input_arrays[i + 1];
                        if (!e) {
                                vv[i + 1] = 1;
                        }
                }

                for (size_t i = 0; i < input_arrays.size (); i++) {
                        if (vv[i]) {
                                if (verbose >= 4) {
                                        myprintf ("  selecting array %d\n", (int)i);
                                }
                                output_arrays.push_back (input_arrays[i]);
                        }
                }
        }
}

array_link reduceLMCform (const array_link &array) {
		if (! (array.strength() >= 2) ) {
			throw_runtime_exception("strength should be at least 2");
		}
        int strength = 2; // assume strength is  2

        arraydata_t ad = arraylink2arraydata (array, 0, strength);
        LMCreduction_t *reduction = new LMCreduction_t (&ad);
        dyndata_t dynd = dyndata_t (ad.N);
        reduction->mode = OA_REDUCE;
        reduction->setArray (array);
        OAextend oaextend;
		const array_t *array_pointer = array.array;
		lmc_t result = LMCreduction (array_pointer, array_pointer, &ad, &dynd, reduction, oaextend);

        array_link nf = reduction->transformation->apply (array);

        return nf;
}

array_link reduceDOPform (const array_link &al, int verbose) {
        int dopruning = 0;
        int dolmc = 1;
        int strength = 2;
        arraylist_t lst;
        lst.push_back (al);
        arraylist_t reduced_arrays;
        if (verbose >= 2)
                myprintf ("reduceDOPform: calling reduceArraysGWLP\n");
        reduceArraysGWLP (lst, reduced_arrays, verbose, dopruning, strength, dolmc);
        return reduced_arrays[0];
}


/** Calculate mixed-level projection GWLP values from normal GWLP values. 
 * 
 * These are the normal projection values, with added the factor level of the removed column.
 *
 */
std::vector< GWLPvalue > mixedProjGWLP (const std::vector< GWLPvalue > dopgwp, const arraydata_t &ad,
                                        int verbose = 0) {
        std::vector< GWLPvalue > GWLPvalues;
        int nc = dopgwp.size ();
        for (int i = 0; i < nc; i++) {
                GWLPvalue w = dopgwp[i];
                std::vector< double > t;
                t.push_back (-ad.s[i]);

                t.insert (t.end (), w.values.begin (), w.values.end ());

                if (verbose >= 2) {
                        myprintf ("reduceArraysGWLP: mixed array: column %d: s[%d]=%d: \n   ", i, i, ad.s[i]);
                        display_vector< double > (t);
                        myprintf ("\n");
                }
                GWLPvalues.push_back (t);
        }

        return GWLPvalues;
}

std::vector< GWLPvalue > projectionDOFvalues (const array_link &array, int verbose ) {
	
	arraydata_t arrayclass=arraylink2arraydata(array);
	std::vector< GWLPvalue > projection_dof_values = projectionGWLPs(array);
	
	if (arrayclass.ismixed() ) {
		projection_dof_values =  mixedProjGWLP(projection_dof_values, arrayclass, verbose);
	}
	return projection_dof_values;
}

void helper_compare_dof_reductions(array_link &alf, array_link &lm, int verbose, lmc_t ret) {
if (alf == lm) {
	if (verbose >= 3)
		myprintf("  array unchanged (ret %d)\n", ret);
	if (verbose >= 3) {
		myprintf("input sort:\n");
		alf.showarray();
		myprintf("output:\n");
		lm.showarray();
	}
}
else {
	if (verbose >= 3) {
		myprintf("input sort:\n");
		alf.showarray();
		myprintf("output:\n");
		lm.showarray();
	}
	if (verbose >= 2)
		myprintf("  array changed\n");
}
}

bool _check_dof_order_minimal(std::vector<GWLPvalue> &dopgwp, int ncols, int verbose) 
{
	GWLPvalue x = *(min_element(dopgwp.begin(), dopgwp.begin() + ncols - 1));
	if (verbose >= 2) {
		myprintf("  delete-1 GWP sequence:        ");
		display_vector< GWLPvalue >(dopgwp);

		myprintf("\n");
		std::cout << "  pruning check: last element " << dopgwp[ncols - 1] << " minimal element " << x << std::endl;
	}
	if (dopgwp[ncols - 1] > x) {
		if (verbose >= 2) {
			myprintf("  reject based on dopgwp ordering\n");
		}
		return false;
	}
	return true;
}

array_transformation_t reductionDOP (const array_link &input_array, int verbose) {
        arraydata_t ad = arraylink2arraydata (input_array);

        std::vector< GWLPvalue > dopgwp = projectionGWLPs (input_array);
		std::vector< GWLPvalue > dofvalues = dopgwp;
		if (ad.ismixed ()) {
                dofvalues = mixedProjGWLP (dopgwp, ad, verbose);
        }

		indexsort is(dofvalues);
        std::vector< DOFvalue > sorted_dofvalues = is.sorted (dofvalues);

        if (verbose >= 2) {
                myprintf ("reductionDOP: delete-1-factor values sorted: ");
                display_vector< DOFvalue > (sorted_dofvalues);
                std::cout << endl;
        }

        symmetry_group sg (sorted_dofvalues, 0);
        if (verbose >= 2)
                sg.show ();

        array_transformation_t column_transformation (&ad);
        column_transformation.setcolperm (is.indices);

        array_link alf (input_array, is.indices);

        if (verbose) {
                myprintf ("is sorting: ");
                display_vector< int > (is.indices);
                myprintf ("\n");
        }
        // NOTE: since the array was re-ordered, the symmetry group has to be re-ordered
        ad.set_colgroups (sg);

        OAextend oaextend;
        oaextend.setAlgorithm (OAextend::getPreferredAlgorithm (ad));
        oaextend.setAlgorithm (MODE_ORIGINAL);

        LMCreduction_t reduction (&ad);

        if (verbose >= 3)
                ad.show (2);

        reduction.mode = OA_REDUCE;
        reduction.init_state = COPY;
		if (1)
        {
                reduction.init_state = INIT;
                reduction.setArray (alf);
                int changed = check_root_update (alf.array, ad, reduction.array);
                copy_array (alf.array, reduction.array, input_array.n_rows, input_array.n_columns); 
        }

        lmc_t ret = LMCcheck (alf, ad, oaextend, reduction);
        array_transformation_t array_transformation = (*(reduction.transformation)) * column_transformation;
        return array_transformation;
}

void reduceArraysGWLP (const arraylist_t &input_arrays, arraylist_t &output_arrays, int verbose, int dopruning, int strength,
                       int dolmc) {
        OAextend oaextend;
        arraylist_t xlist;
        int npruned = 0;

        if (verbose >= 2)
                myprintf ("reduceArraysGWLP: start\n");

        for (size_t i = 0; i < (size_t)input_arrays.size (); i++) {
                if (verbose >= 2)
                        myprintf ("reduceArrays: array %d/%d\n", (int)i, (int)input_arrays.size ());
                const array_link input_array = input_arrays.at (i);
                int ncols = input_array.n_columns;

                arraydata_t ad = arraylink2arraydata (input_array, 0, strength);

				std::vector< GWLPvalue > dopgwp = projectionGWLPs(input_array);

                if (dopruning) {
					if (! _check_dof_order_minimal(dopgwp, ncols, verbose) ) {
						npruned++;
						continue;
					}
                }

                if (verbose >= 2)
                        printfd ("reduceArraysGWLP:\n");

                std::vector< DOFvalue > dofvalues = dopgwp;
                if (ad.ismixed ()) {
                        dofvalues = mixedProjGWLP (dopgwp, ad, verbose);
                }

				indexsort is(dofvalues);
				std::vector< DOFvalue > sorted_dofvalues = is.sorted (dofvalues);

                if (verbose >= 2) {
                        myprintf ("  delete-1-factor values sorted: ");
                        display_vector< DOFvalue > (sorted_dofvalues);

                        myprintf ("\n");
                }

                symmetry_group sg (sorted_dofvalues, 0);
                if (verbose >= 2)
                        sg.show ();


                array_link alf (input_array, is.indices);
                ad.set_colgroups (sg);

                OAextend oaextend;
                oaextend.setAlgorithm (OAextend::getPreferredAlgorithm (ad));
                oaextend.setAlgorithm (MODE_ORIGINAL);

                LMCreduction_t reduction (&ad);

                array_link lm (input_array);

                reduction.mode = OA_REDUCE;
                reduction.init_state = COPY;
                {
                        reduction.init_state = INIT;
                        reduction.setArray (alf);

                        int changed = check_root_update (alf.array, ad, reduction.array);
                        if (changed && verbose >= 4) {
                                printfd ("reduceArraysGWLP: array changed %d.\n", changed);
                                reduction.getArray ().showarray ();
                        }
                }


                lmc_t ret = LMCcheck (alf, ad, oaextend, reduction);
                copy_array (reduction.array, lm.array, lm.n_rows, lm.n_columns);
                if (verbose >= 2)
                        printfd ("reduceArraysGWLP: ret %d (LMC_MORE %d, strength %d)\n", ret, LMC_MORE, ad.strength);

				helper_compare_dof_reductions(alf, lm, verbose, ret);
			
                xlist.push_back (lm);
        }
        if (verbose) {
                myprintf ("reduceArraysGWLP: input size %d, output size %d, pruned %d\n", (int)input_arrays.size (),
                          (int)xlist.size (), npruned);
        }

        if (dolmc) {
                selectUniqueArrays (xlist, output_arrays, verbose);

        } else {
                sort (xlist.begin (), xlist.end ()); 

                for (size_t i = 0; i < xlist.size (); i++) {
                        output_arrays.push_back (xlist[i]);
                }
        }
        if (verbose)
                myprintf ("  selecting %d/%d arrays\n", (int)output_arrays.size (), (int)xlist.size ());
}
