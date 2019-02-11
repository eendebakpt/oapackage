#include <sstream>
#include <string>
#include <algorithm>
#include <sstream>

#include "arraytools.h"

#include "arrayproperties.h"
#include "mathtools.h"
#include "strength.h"
#include "tools.h"
#include <errno.h>


#ifdef FULLPACKAGE
#include "bitarray/bit_array.h"
#include "lmc.h"
#include <iostream>
#endif

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

int sizeof_array_t () { return sizeof (array_t); }

#ifdef WIN32
#else
#include <fcntl.h>
int get_file_status (FILE *f) {
        int fd = fileno (f);
        return fcntl (fd, F_GETFL);
}

#endif

array_transformation_t array_transformation_t::operator* (const array_transformation_t b) const {

	array_transformation_t c(this->ad);

	const array_transformation_t &a = *this;

	const int nc = this->ad->ncols;

	// perform the rows permutations
	perform_inv_perm(b.rperm, c.rperm, this->ad->N, a.rperm);

	// perform the column permutations
	perform_inv_perm(b.cperm, c.cperm, nc, a.cperm);

	/* level permutations */
	for (colindex_t ci = 0; ci < ad->ncols; ci++) {
		levelperm_t l1 = b.lperms[a.cperm[ci]];
		levelperm_t l2 = a.lperms[ci];

		composition_perm(l1, l2, this->ad->s[ci], c.lperms[ci]);
	}

	return c;
}
int array_transformation_t::operator== (const array_transformation_t &t2) const {
        if (this->ad->N != t2.ad->N) {
                return 0;
        }

        int nrows = this->ad->N;
        int ncols = this->ad->ncols;

        if (this->ad->ncols != t2.ad->ncols) {
                return 0;
        }

        if (!std::equal (this->rperm, this->rperm + nrows, t2.rperm)) {
                return 0;
        }
        if (!std::equal (this->cperm, this->cperm + ncols, t2.cperm)) {
                return 0;
        }

        for (int c = 0; c < ncols; c++) {
                int nlevels = this->ad->s[c];
                if (nlevels != t2.ad->s[c])
                        return 0;
                if (!std::equal (this->lperms[c], this->lperms[c] + nlevels, t2.lperms[c])) {
                        return 0;
                }
        }

        return 1;
}

std::vector< int > array_transformation_t::rowperm () const {
        std::vector< int > ww (this->rperm, this->rperm + this->ad->N);
        return ww;
}

std::vector< int > array_transformation_t::colperm () const {
        std::vector< int > ww (this->cperm, this->cperm + this->ad->ncols);
        return ww;
}
std::vector< int > array_transformation_t::lvlperm (int c) const {
        if (c < 0 || c >= this->ad->ncols) {
                std::vector< int > ww;
                return ww;
        }
        std::vector< int > ww (this->lperms[c], this->lperms[c] + this->ad->s[c]);
        return ww;
}

void array_transformation_t::setrowperm (std::vector< int > rp) {
        if ((int)rp.size () != this->ad->N) {
                myprintf ("array_transformation_t::setrowperm: argument has wrong dimensions\n");
                return;
        }
        std::copy (rp.begin (), rp.end (), this->rperm);
}
void array_transformation_t::setcolperm (std::vector< int > colperm) {
        if ((int)colperm.size () != this->ad->ncols) {
                myprintf ("array_transformation_t::setrowperm: argument has wrong dimensions\n");
                return;
        }
        std::copy (colperm.begin (), colperm.end (), this->cperm);
}
void array_transformation_t::setlevelperm (int colindex, std::vector< int > lvlperm) {
        if (colindex < 0 || colindex >= this->ad->ncols) {
                myprintf ("array_transformation_t::setrowperm: argument has wrong dimensions\n");
                return;
        }
        if ((int)lvlperm.size () != this->ad->s[colindex]) {
                myprintf ("array_transformation_t::setrowperm: argument has wrong dimensions\n");
                return;
        }
        std::copy (lvlperm.begin (), lvlperm.end (), this->lperms[colindex]);
}

/**
 * @brief Print an array transformation to an output stream
 * @param out
 */
void array_transformation_t::show (std::ostream &out) const {
     
        if(this->ad==0) {
           out << "array transformation: no class defined" << std::endl;
           return;    
        }
        out << "array transformation: N " << ad->N;
        if (this->isIdentity ()) {
                out << ": identity transformation" << std::endl;
        } else {
                out << std::endl;
                out << "column permutation: ";
                print_perm (out, cperm, ad->ncols);

                out << "level perms:" << endl;
                for (colindex_t c = 0; c < ad->ncols; c++) {
                        print_perm (out, lperms[c], ad->s[c]);
                }

                out << "row permutation: ";
                print_perm (out, rperm, ad->N);
        }
}

void array_transformation_t::show () const {
#ifdef FULLPACKAGE
        std::stringstream ss;
        this->show (ss);
        myprintf ("%s", ss.str ().c_str ());
#endif
}

array_transformation_t::array_transformation_t () {
        ad = 0;
        rperm = 0;
        lperms = 0;
        cperm = 0;
}

/**
 * @brief Create an OA transformation
 * @param adp
 */
array_transformation_t::array_transformation_t (const arraydata_t *adp) {
        ad = new arraydata_t (*adp);
        allocate_data_structures ();
}
array_transformation_t::array_transformation_t (const arraydata_t &adp) {
        ad = new arraydata_t (adp);
        allocate_data_structures ();
}

array_transformation_t::array_transformation_t (const array_transformation_t &tt) {
        ad = new arraydata_t (*(tt.ad));

        allocate_data_structures ();

        // copy data
        std::copy (tt.rperm, tt.rperm + ad->N, rperm);
        std::copy (tt.cperm, tt.cperm + ad->ncols, cperm);
        for (colindex_t c = 0; c < ad->ncols; c++) {
                std::copy (tt.lperms[c], tt.lperms[c] + ad->s[c], lperms[c]);
        }
}

void array_transformation_t::reset () {
        init_perm (this->rperm, this->ad->N);
        init_perm< colindex_t > (this->cperm, ad->ncols);

        for (colindex_t c = 0; c < ad->ncols; c++) {
                init_perm (lperms[c], ad->s[c]);
        }
}

void array_transformation_t::allocate_data_structures () {
        if (ad==0)
               throw_runtime_exception("need valid arraydata_t structure");
        rperm = new_perm_init< rowindex_t > (ad->N);
        cperm = new_perm_init< colindex_t > (ad->ncols);

        lperms = new levelperm_t[ad->ncols];
        for (colindex_t c = 0; c < ad->ncols; c++) {
                lperms[c] = new_perm_init< array_t > (ad->s[c]);
        }
}

void array_transformation_t::free_data_structures () {
        if (ad == 0) {
                return;
        }

        delete_perm (rperm);
        delete_perm (cperm);
        for (colindex_t c = 0; c < ad->ncols; c++) {
                delete_perm (lperms[c]);
        }
        delete[] lperms;

        if (ad != 0) {
                delete ad;
        }
        ad=0;
        rperm=0;
        cperm=0;
        lperms=0;
}

/// Assignment operator
array_transformation_t &array_transformation_t::operator= (const array_transformation_t &tt) {
       if (tt.ad==0) {
        free_data_structures ();
        
        return *this;
       }
        
        free_data_structures ();

        ad = new arraydata_t (*(tt.ad));
        allocate_data_structures ();

        std::copy (tt.rperm, tt.rperm + ad->N, rperm);
        std::copy (tt.cperm, tt.cperm + ad->ncols, cperm);
        for (colindex_t c = 0; c < ad->ncols; c++) {
                std::copy (tt.lperms[c], tt.lperms[c] + ad->s[c], lperms[c]);
        }

        return *this;
}

/**
 * @brief array_transformation_t destructor
 */
array_transformation_t::~array_transformation_t () {
        delete_perm (rperm);
        delete_perm (cperm);
        if (lperms != 0) {
                for (colindex_t c = 0; c < ad->ncols; c++) {
                        delete_perm (lperms[c]);
                }
                delete[] lperms;
                lperms = 0;
        }

        delete ad;
}

bool array_transformation_t::isIdentity () const {
        for (int i = 0; i < ad->ncols; ++i) {
                if (cperm[i] != i) {
                        return 0;
                }
        }
        for (int i = 0; i < ad->N; ++i) {
                if (rperm[i] != i) {
                        return 0;
                }
        }
        for (int c = 0; c < ad->ncols; ++c) {
                for (int i = 0; i < ad->s[c]; ++i) {
                        if (lperms[c][i] != i) {
                                return 0;
                        }
                }
        }
        return 1;
}

/// apply transformation and show resulting array
void array_transformation_t::print_transformed (carray_t *source) const {
        array_t *a = clone_array (source, this->ad->N, this->ad->ncols);
        print_array (a, this->ad->N, this->ad->ncols);
        destroy_array (a);
}

int array_link::operator!= (const array_link &rhs_array) const {
	if (!this->equal_size(rhs_array)) {
		myprintf("array_link::operator== comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
			this->index, rhs_array.index, this->n_rows, this->n_columns, rhs_array.n_rows, rhs_array.n_columns);
		throw_runtime_exception("comparing arrays with different sizes");
	}

	return (!std::equal(array, array + n_rows * n_columns, rhs_array.array));
}

array_link array_transformation_t::apply(const array_link &al) const {
	array_link trx(al);
	this->apply(al.array, trx.array);
	return trx;
}

int array_link::operator< (const array_link &rhs) const {
	if (!this->equal_size(rhs)) {
		myprintf("array_link::operator< comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
			this->index, rhs.index, this->n_rows, this->n_columns, rhs.n_rows, rhs.n_columns);
		return 0;
	}
	return std::lexicographical_compare(array, array + n_rows * n_columns, rhs.array,
		rhs.array + n_rows * n_columns);
}

int array_link::equalsize(const array_link &rhs) const {
	return (this->n_rows == rhs.n_rows && this->n_columns == rhs.n_columns);
}

int array_link::operator== (const array_link &rhs_array) const {
	if (!this->equal_size(rhs_array)) {
		return 0;
	}
	return std::equal(array, array + n_rows * n_columns, rhs_array.array);
}

int array_link::operator> (const array_link &rhs) const {
	if (!this->equal_size(rhs)) {
		myprintf("array_link::operator> comparing arrays (%d %d) with different sizes: (%d,%d) (%d, %d)!\n",
			this->index, rhs.index, this->n_rows, this->n_columns, rhs.n_rows, rhs.n_columns);
		return 0;
	}

	return std::lexicographical_compare(rhs.array, rhs.array + n_rows * n_columns, array,
		array + n_rows * n_columns);
}

array_link &array_link::operator= (const array_link &rhs) { return deepcopy(rhs); }

array_link &array_link::shallowcopy(const array_link &rhs) {
	this->n_rows = rhs.n_rows;
	this->n_columns = rhs.n_columns;
	this->index = rhs.index;
	this->array = rhs.array;
	return *this;
}

array_link &array_link::deepcopy(const array_link &rhs) {
	this->n_rows = rhs.n_rows;
	this->n_columns = rhs.n_columns;
	this->index = rhs.index;
	// perform deep copy
	if (this->array) {
		destroy_array(this->array);
	}
	if (rhs.array == 0) {
		this->array = create_array(this->n_rows, this->n_columns);
	}
	else {
		this->array = clone_array(rhs.array, this->n_rows, this->n_columns);
	}
	return *this;
}

void array_transformation_t::apply (array_t *sourcetarget) const {
        array_t *tmp = create_array (ad->N, ad->ncols);
        copy_array (sourcetarget, tmp, ad->N, ad->ncols);
        this->apply (tmp, sourcetarget);
        destroy_array (tmp);
}

/** @brief Perform a row permutation
*
* @param source Source array
* @param target Target array
* @param perm Permutation to perform
* @param nrows Number of rows
* @param ncols Numer of columns
*/
inline void perform_inv_row_permutation (const array_t *source, array_t *target, rowperm_t perm, int nrows,
                                         int ncols) {
        for (int i = 0; i < ncols; i++)
                for (int j = 0; j < nrows; j++) {
                        target[nrows * i + j] = source[nrows * i + perm[j]];
                }
}

/** @brief Perform inverse column permutation on an array
 *
 * @param source
 * @param target
 * @param perm
 * @param nrows
 * @param ncols
 */
inline void perform_inv_column_permutation (const array_t *source, array_t *target, colperm_t perm, int nrows,
                                            int ncols) {
        for (int i = 0; i < ncols; i++) {
                memcpy (&target[i * nrows], &source[perm[i] * nrows], nrows * sizeof (array_t));
        }
}

void array_transformation_t::apply (const array_t *source, array_t *target) const {

        array_t *tmp = create_array (ad);

        /* column permutations */
        perform_inv_column_permutation (source, tmp, cperm, ad->N, ad->ncols);

        /* level permutations */
        for (colindex_t c = 0; c < ad->ncols; c++) {
#ifdef SAFELPERM
                safe_perform_level_perm (tmp + c * ad->N, ad->N, lperms[c], ad->s[c]);
#else
                perform_level_perm (tmp + c * ad->N, ad->N, lperms[c]);
#endif
        }

        /* row permutations */
        perform_inv_row_permutation (tmp, target, rperm, ad->N, ad->ncols);

        destroy_array (tmp);
}

/// initialize to a random transformation
void array_transformation_t::randomize () {
        /* row permutation */
        random_perm (rperm, ad->N);

        /* column permutation */
        for (int x = 0; x < ad->ncolgroups; x++) {
                random_perm (cperm + ad->colgroupindex[x], +ad->colgroupsize[x]);
        }

        /* level permutation */
        for (colindex_t c = 0; c < ad->ncols; c++) {
                random_perm (lperms[c], ad->s[c]);
        }
}

/// initialize to a random transformation
void array_transformation_t::randomizecolperm () {
        /* column permutation */
        for (int x = 0; x < ad->ncolgroups; x++) {
                random_perm (cperm + ad->colgroupindex[x], +ad->colgroupsize[x]);
        }
}

/// initialize to a random transformation
void array_transformation_t::randomizerowperm () { random_perm (this->rperm, ad->N); }
array_transformation_t array_transformation_t::inverse () const {
        array_transformation_t A (*this);

        invert_permutation (this->rperm, this->ad->N, A.rperm);
        invert_permutation (this->cperm, this->ad->ncols, A.cperm);

        /* level permutations */
        for (colindex_t ci = 0; ci < ad->ncols; ci++) {

                colindex_t cir = this->cperm[ci];

                invert_permutation (this->lperms[ci], this->ad->s[ci], A.lperms[cir]);

                if (ci < 0) {
                        myprintf ("ci %d: this->lperms[ci] ", ci);
                        print_perm (this->lperms[ci], this->ad->s[cir]);
                        myprintf ("cir %d: A->lperms[cir] ", cir);
                        print_perm (A.lperms[cir], this->ad->s[ci]);
                }
        }
        return A;
}

/**
 * Return factors of an array. These can be read from the array of the array is in LMC
 * form? Better is to use an oaconfig file...
 *
 */
void get_factors (array_t *s, carray_t *array, const rowindex_t N, const colindex_t ncols) {
        int count = 0;

        array_t max;
        for (colindex_t i = 0; i < ncols; i++) {
                max = 0;
                for (rowindex_t j = 0; j < N; j++)
                        if (array[count++] > max) {
                                max = array[count - 1];
                        }
                s[i] = max + 1;
        }
}

arraydata_t arraylink2arraydata (const array_link &al, int extracols, int strength) {
        int verbose = 0;

        if(al.min() < 0) {
         throw_runtime_exception("array should have positive integer values to convert to arraydata_t structure");    
        }
        // create arraydatya
        int ncols0 = al.n_columns;
        int N = al.n_rows;
        int ncols = ncols0 + extracols;
        std::vector< int > s (ncols);
        std::vector< int > smin (ncols);
        if (N > 0) {
                int ss = -1;
                for (int ik = 0; ik < ncols0; ik++) {
                        array_t *xx = std::max_element (al.array + N * ik, al.array + N * (ik + 1));
                        ss = *xx + 1;
                        s[ik] = ss;
                        array_t minval = *(std::min_element (al.array + N * ik, al.array + N * (ik + 1)));
                        smin[ik] = 0;

                        if (s[ik] < 0) {
                                myprintf ("arraylink2arraydata: column %d: input array should have elements ranging "
                                          "from 0 to s-1\n",
                                          ik);
                                s[ik] = 1;
                        }
                        if (smin[ik] < 0) {
                                myprintf ("arraylink2arraydata: column %d: input array should have elements ranging "
                                          "from 0 to s-1\n",
                                          ik);
                        }
                }
                for (int ik = ncols0; ik < ncols; ik++) {
                        s[ik] = ss; // repeat last factor value
                }
        }

        if (strength < 0) {
                strength = al.strength ();
        }
        arraydata_t ad (s, al.n_rows, strength, ncols);
        if (verbose) {
                myprintf ("arraylink2arraydata: ");
                ad.show ();
        }
        return ad;
}

/// helper function
void foldtest (jstruct_t &js, const array_link &al, int jj, int verbose) {
        const rowindex_t N = al.n_rows;
        const colindex_t k = al.n_columns;
        int *pp = new_comb_init< int > (jj);

        array_t **tmpcol = malloc2d< array_t > (jj + 1, N);
        std::fill (tmpcol[0], tmpcol[0] + N, 0);

        int fold;
        int cb = 0;
        int nc = ncombs (al.n_columns, jj);
        for (int x = 0; x < nc; x++) {
#ifdef FULLPACKAGE
                if (verbose >= 4) {
                        myprintf ("x %d: ", x);
                        std::cout << printfstring ("comb: ");
                        print_perm (pp, jj);
                }
#endif
                // update invalid columns
                for (int ii = cb; ii < jj; ii++) {
                        for (int r = 0; r < N; r++) {
                                tmpcol[ii + 1][r] = tmpcol[ii][r] + al.array[r + pp[ii] * N];
                        }
                }

                // calculate j-value from last column
                int jval = 0;
                int tmp = 0;
                for (int r = 0; r < N; r++) {
                        tmp = tmpcol[jj][r];
                        tmp %= 2;
                        jval += tmp;
                }
                jval = -(2 * jval - N);

                js.values[x] = jval;
                fold = next_combination_fold (pp, jj, k);
                cb = fold;

                if (verbose >= 2) {
                        myprintf ("comb %d: jval %d\n", x, jval);
                }
        }

        free2d (tmpcol, jj + 1);
        delete_comb (pp);
}

arraylist_t addConstant (const arraylist_t &lst, int value) {
        arraylist_t output_arrays (lst.size ());

        for (size_t i = 0; i < lst.size (); i++) {
                output_arrays[i] = lst[i] + value;
        }

        return output_arrays;
}

/** Return number of arrays with j_{2n+1}=0 for n<m */
std::vector< int > getJcounts (arraylist_t *arraylist, int N, int k, int verbose) {
        std::vector< int > aa (k + 1);

        for (int i = 0; i < (int)arraylist->size (); i++) {
#ifdef FULLPACKAGE
                if (verbose) {
                        if (i < 100 || i % 1000 == 0) {
                                std::cout << "## analyzing array " << i << "/" << arraylist->size () << std::endl;
                        }
                }
#endif
                array_link al = arraylist->at (i);

                int jl = 5;
                while (jl <= al.n_columns) {
                        jstruct_t js (al.n_rows, al.n_columns, jl);
                        if (0) {
                                arraylist_t *all = new arraylist_t ();
                                all->push_back (al);
                                std::vector< jstruct_t > xx = analyseArrays (*all, verbose, jl);
                                js = jstruct_t (xx[0]);
                                if (verbose >= 2) {
                                        myprintf (" old: ");
                                        xx[0].show ();
                                }
                                delete all;
                        } else {
                                foldtest (js, al, jl, verbose);
                        }
                        if (!js.allzero ()) {
                                if (verbose >= 3) {
                                        myprintf (" new: ");
                                        js.show ();
                                }

                                aa[jl]++;
                                if (verbose >= 3) {
                                        js.show ();
                                }
                                break;
                        }
                        jl += 2;
                        if (verbose >= 2) {
                                myprintf ("  jl %d\n", jl);
                        }
                }
                if (jl > al.n_columns) {
                        aa[0]++;
                }
        }
        return aa;
}

/**
 * @brief Default copy constructor
 * @param A
 */
array_link::array_link (const array_link &rhs) {
        array = 0;
        deepcopy (rhs);
}

#ifdef SWIGCODE
array_link::array_link (long *pymatinput, int nrows, int ncols) {
        this->index = INDEX_DEFAULT;
        this->n_columns = ncols;
        this->n_rows = nrows;
        this->array = create_array (this->n_rows, this->n_columns);
        int i = 0;
        for (int row = 0; row < this->n_rows; row++) {
                for (int col = 0; col < this->n_columns; col++) {
                        this->array[row + col * this->n_rows] = pymatinput[i];
                        i++;
                }
        }
}

#endif

array_link::array_link (Eigen::MatrixXd &m) {
        this->index = INDEX_DEFAULT;
        this->n_columns = m.cols ();
        this->n_rows = m.rows ();
        this->array = create_array (this->n_rows, this->n_columns);
        for (int i = 0; i < this->n_columns * this->n_rows; i++) {
                this->array[i] = m (i);
        }
}

/// create an array by permuting columns
array_link::array_link (const array_link &rhs, const std::vector< int > &colperm) : array (0) {
        init (rhs.n_rows, rhs.n_columns);
        for (int c = 0; c < rhs.n_columns; c++) {
                int cp = colperm[c];
                std::copy (rhs.array + n_rows * cp, rhs.array + n_rows * (cp + 1), this->array + c * n_rows);
        }
}


array_t array_link::at (const int index) const {
        if ( ! this->_valid_index(index) ) {
                throw std::out_of_range (printfstring("array_link error: index %dout of bounds\n", index));
        }
        return this->_at(index);
}

array_t array_link::_at (const int index) const { return this->array[index]; }

array_link array_link::clone () const {

        array_link al (*this);
        return al;
}

void array_link::setvalue (int r, int c, int val) {
        if ((r < 0) || (r >= this->n_rows) || (c < 0) || (c >= this->n_columns)) {
                myprintf ("array_link::setvalue: index out of bounds %d %d shape (%d %d)!!\n", r, c, this->n_rows,
                          this->n_columns);
                return;
        }

        this->array[r + this->n_rows * c] = val;
}

void array_link::setvalue (int r, int c, double val) {
        if ((r < 0) || (r >= this->n_rows) || (c < 0) || (c >= this->n_columns)) {
                myprintf ("array_link::setvalue: index out of bounds %d %d shape (%d %d)!!\n", r, c, this->n_rows,
                          this->n_columns);
                return;
        }

        this->array[r + this->n_rows * c] = val;
}

void array_link::_setvalue (int r, int c, int val) { this->array[r + this->n_rows * c] = val; }

void array_link::negateRow (rowindex_t row) {
          for (int c = 0; c < this->n_columns; c++) {
               this->atfast (row, c) *= -1;
          }
}
void array_link::show () const { myprintf ("array_link: index: %d, shape (%d, %d), data pointer %p\n", index, n_rows, n_columns, (void *)array); }

std::string array_link::showstr () const {
          std::stringstream s;
          s << "array_link: " << n_rows << ", " << n_columns << "";
          std::string rs = s.str ();
          return rs;
}
        
bool array_link::_valid_index (int index) const {
        if ((index < 0) || (index >= this->n_rows*this->n_columns)) {
                return false;
        }
        return true;
}

bool array_link::_valid_index (const rowindex_t r, const colindex_t c) const {
        if ((r < 0) || (r >= this->n_rows) || (c < 0) || (c >= this->n_columns)) {
                return false;
        }
        return true;
}

array_t array_link::_at (const rowindex_t r, const colindex_t c) const { return this->array[r + this->n_rows * c]; }

/**
 * @brief Element access
 */
array_t array_link::at (const rowindex_t r, const colindex_t c) const {
        if ( ! this->_valid_index(r, c) ) {
                throw std::out_of_range (printfstring("array_link::at: index out of bounds %d %d shape (%d %d)!!\n", r, c, this->n_rows,
                          this->n_columns));
        }

        return this->array[r + this->n_rows * c];
}

/**
 * @brief Element access
 */
array_t &array_link::at (const rowindex_t r, const colindex_t c) {
        if ( ! this->_valid_index(r, c) ) {
                throw std::out_of_range (printfstring("array_link::at: index out of bounds %d %d shape (%d %d)!!\n", r, c, this->n_rows,
                          this->n_columns));
        }

        return this->array[r + this->n_rows * c];
}

/**
 * @brief Element access
 */
void array_link::setconstant (array_t value) { std::fill (this->array, this->array + n_rows * n_columns, value); }

//! Create array link with clone of an array
array_link::array_link (const array_t *array, rowindex_t nrows, colindex_t ncols, int index_)
    : n_rows (nrows), n_columns (ncols), index (index_) {
        this->array = clone_array (array, nrows, ncols);

}

array_link::array_link (const std::vector< int > &values, rowindex_t nrows, colindex_t ncols, int index_)
    : n_rows (nrows), n_columns (ncols), index (index_) {
        if (values.size () != (size_t)nrows * ncols) {
                throw_runtime_exception ("array_link: size of vector does not match number of rows and columns\n");
        }
        this->array = create_array (nrows, ncols);
        std::copy (values.begin (), values.begin () + nrows * ncols, this->array);
}

//! Default constructor
array_link::array_link () : n_rows (-1), n_columns (0), index (INDEX_NONE), array (0) {
}

void array_link::init (rowindex_t r, colindex_t c) {
        if (array != 0) {
                destroy_array (array);
                array = 0;
        }
        n_rows = r;
        n_columns = c;

        this->array = create_array (n_rows, n_columns);
}

//! Default destructor
array_link::~array_link () {
        if (array != 0) {
                destroy_array (array);
                array = 0;
        }
}

array_link::array_link (rowindex_t nrows, colindex_t ncols, int index_)
    : n_rows (nrows), n_columns (ncols), index (index_) {
        this->array = create_array (nrows, ncols);
        this->setconstant(0);
}

array_link::array_link (rowindex_t nrows, colindex_t ncols, int index_, carray_t *data)
    : n_rows (nrows), n_columns (ncols), index (index_) {
        this->array = create_array (nrows, ncols);
        this->setarraydata (data, nrows * ncols);
}

bool array_link::columnEqual (int column_index, const array_link &rhs, int column_index_rhs) const {
        if (this->n_rows != rhs.n_rows)
                return false;
		if (column_index < 0 || column_index >= rhs.n_columns)
			throw std::out_of_range("index out of bound");
		if (column_index_rhs < 0 || column_index_rhs>= rhs.n_columns)
			throw std::out_of_range("index out of bound");
        const array_t *pl = this->array + this->n_rows * column_index;
        const array_t *pr = rhs.array + this->n_rows * column_index_rhs;

        for (int r = 0; r < this->n_rows; r++) {
                if (pl[r] != pr[r]) {
                        return false;
                }
        }
        return true;
}

int array_link::firstColumnDifference (const array_link &A) const {
        int r = -1, c = -1;
        this->firstDiff (A, r, c, 0);
        return c;
}

bool array_link::firstDiff(const array_link &A, int &r, int &c, int verbose) const {
	r = 0;
	c = 0;
	for (c = 0; c < this->n_columns; c++) {
		for (r = 0; r < this->n_rows; r++) {
			if (this->at(r, c) != A.at(r, c)) {
				if (verbose) {
					myprintf("first difference of array at %d, %d\n", r, c);
				}
				return true;
			}
		}
	}
	return false;
}

/// return md5 sum of array representation (as represented with 32bit int datatype in memory)
std::string array_link::md5 () const {
        if (sizeof (array_t) == 4) {
                // we have int as the representation type
                return ::md5 (this->array, this->n_columns * this->n_rows * sizeof (array_t));
        } else {
                short int nn = this->n_columns * this->n_rows;
                short int *x = new short int[nn];
                std::copy (array, array + nn, x);
                std::string m = ::md5 (x, nn * sizeof (short int));
                delete[] x;
                return m;
        }
}

void showArrayList (const arraylist_t &lst) {
        for (size_t i = 0; i < lst.size (); i++) {
                myprintf ("array %d:\n", (int)i);
                lst[i].showarray ();
        }
}

template < class atype >
/// write array to output stream
void write_array_format(std::ostream &ss, const atype *array, const int nrows, const int ncols, int width = 3) {
	assert(array != 0 || ncols == 0);

	int count;
	for (int j = 0; j < nrows; j++) {
		count = j;
		for (int k = 0; k < ncols; k++) {
			const char *s = (k < ncols - 1) ? " " : "\n";
			ss << std::setw(width) << array[count] << s;
			count += nrows;
		}
	}
}

#ifdef FULLPACKAGE

/**
 * @brief Print an array to a stream (formatted)
 * @param stream
 * @param A
 * @return
 */
std::ostream &operator<< (std::ostream &stream, const array_link &A) {
        write_array_format (stream, A.array, A.n_rows, A.n_columns);
        return stream;
}

#endif // FULLPACKAGE

array_link array_link::selectFirstRows (int nrows) const {
        mycheck (nrows >= 0, "array_link::selectFirstColumn: nrows<0\n");
        mycheck (nrows <= this->n_rows, "array_link::selectFirstRows: nrows %d too large\n", nrows);
        array_link d (nrows, this->n_columns, -1);
        for (int i = 0; i < this->n_columns; i++) {
                std::copy (this->array + i * this->n_rows, this->array + (i) * this->n_rows + nrows,
                           d.array + i * nrows);
        }
        return d;
}

array_link array_link::selectFirstColumns (int ncols) const {
        mycheck (ncols >= 0, "array_link::selectFirstColumn: ncols<0\n");
        mycheck (ncols <= this->n_columns,
                 "array_link::selectFirstColumn: ncols %d too large for array with %d columns\n", ncols,
                 this->n_columns);
        ncols = std::min (ncols, this->n_columns);
        array_link d (this->n_rows, ncols, -1);
        for (int i = 0; i < ncols; i++) {
                std::copy (this->array + i * this->n_rows, this->array + (i + 1) * this->n_rows,
                           d.array + i * this->n_rows);
        }
        return d;
}

array_link array_link::selectLastColumns (int ii) const {
        mycheck (ii >= 0, "array_link::selectFirstColumn: ii<0\n");
        array_link d (this->n_rows, ii, INDEX_DEFAULT);
        for (int j = 0; j < ii; j++) {
                int i = this->n_columns - ii + j;
                std::copy (this->array + i * this->n_rows, this->array + (i + 1) * this->n_rows,
                           d.array + j * this->n_rows);
        }
        return d;
}

array_link array_link::selectColumns (int c) const {
        array_link d (this->n_rows, 1, INDEX_DEFAULT);
        mycheck (c >= 0, "array_link::selectColumns: c<0\n");
        mycheck (c < this->n_columns, "array_link::selectColumns: c>=ncols\n");
        std::copy (this->array + c * this->n_rows, this->array + (c + 1) * this->n_rows, d.array);
        return d;
}

array_link array_link::selectColumns (std::vector< int > c) const {
        array_link d (this->n_rows, c.size (), INDEX_DEFAULT);
        for (int j = 0; j < (int)c.size (); j++) {
                int i = c[j];
                mycheck (i >= 0, "array_link::selectColumns: i<0\n");
                mycheck (i < this->n_columns, "array_link::selectColumns: i>ncols\n");
                std::copy (this->array + i * this->n_rows, this->array + (i + 1) * this->n_rows,
                           d.array + j * this->n_rows);
        }
        return d;
}

bool array_link::isSymmetric () const {
        int n = std::min ((int)this->n_columns, (int)this->n_rows);
        for (int c = 0; c < n; c++) {
                for (int r = 0; r < c; r++) {
                        if (this->at (c, r) != this->at (r, c))
                                return false;
                }
        }
        return true;
}
void array_link::makeSymmetric () {
        int n = std::min ((int)this->n_columns, (int)this->n_rows);
        for (int c = 0; c < n; c++) {
                for (int r = 0; r < c; r++) {
                        this->at (c, r) = this->at (r, c);
                }
        }
}
array_link array_link::deleteColumn (int index) const {
        mycheck (this->n_columns >= 2, "array_link::deleteColumn: array has <2 columns\n");
        array_link d (this->n_rows, this->n_columns - 1, -1);
        if (index < 0) {
                index = this->n_columns - 1;
        }
        for (int i = 0; i < index; i++) {
                std::copy (this->array + i * this->n_rows, this->array + (i + 1) * this->n_rows,
                           d.array + i * this->n_rows);
        }
        for (int i = index + 1; i < this->n_columns; i++) {
                std::copy (this->array + i * this->n_rows, this->array + (i + 1) * this->n_rows,
                           d.array + (i - 1) * this->n_rows);
        }
        return d;
}

array_link array_link::transposed () const {
        array_link d (this->n_columns, this->n_rows, -1);
        for (int i = 0; i < this->n_columns; i++)
                for (int j = 0; j < this->n_rows; j++) {
                        d.array[i + j * d.n_rows] = this->array[j + i * this->n_rows];
                }

        return d;
}

array_link array_link::randomperm () const {
        arraydata_t arrayclass = arraylink2arraydata (*this);
        array_transformation_t trans (&arrayclass);
        trans.randomize ();
        return trans.apply (*this);
}

array_link array_link::randomcolperm () const {
        arraydata_t arrayclass = arraylink2arraydata (*this);
        array_transformation_t trans (&arrayclass);
        trans.randomizecolperm ();
        return trans.apply (*this);
}

array_link array_link::randomrowperm () const {
        arraydata_t arrayclass = arraylink2arraydata (*this);
        array_transformation_t trans (&arrayclass);
        trans.randomizerowperm ();
        return trans.apply (*this);
}

/** Return example array
 *
 * \param idx Index of example array
 * \param verbose If postive then print description of the array to stdout
 *
 */
array_link exampleArray (int idx, int verbose) {
        if (idx == -1) {
                for (int i = 0; i < 1000; i++) {
                        array_link al = exampleArray (i, verbose);
                        if (al.n_columns == 0)
                                break;
                }
                return exampleArray (0);
        }

        std::string dstr = "";

        switch (idx) {
        default:
                myprintf ("exampleArray: no example array with index %d exists\n", idx);
                return array_link ();
                break;
        case 39: {
                dstr = "first LMC0 conference design in C(8,6)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (8, 6, 0);
                int tmp[] = {0, 1,  1,  1, 1, 1, 1,  1, 1, 0,  1, 1,  1, -1, -1, -1, 1, -1, 0,  1,  -1, 1, 1, -1,
                             1, -1, -1, 0, 1, 1, -1, 1, 1, -1, 1, -1, 0, -1, 1,  1,  1, 1,  -1, -1, 1,  0, 1, -1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 38: {
                dstr = "LMC0 conference design in C(30,3)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (30, 3, 0);
                int tmp[] = {0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1, 1,
                             1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1, 1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1,  1,  0,  1,  1,  1,  1, 1, 1,
                             -1, -1, -1, -1, -1, -1, -1, 1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 35: {
                dstr = "first double conference design in DC(20,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (20, 4, 0);
                int tmp[] = {0, 0,  1, 1, 1, 1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             1, -1, 0, 1, 1, 1,  1,  -1, -1, -1, -1, 0,  1,  1,  1,  1,  -1, -1, -1, -1,
                             1, -1, 1, 0, 1, -1, -1, 1,  1,  -1, -1, -1, 1,  1,  -1, -1, 0,  1,  1,  -1,
                             1, -1, 1, 1, 0, -1, -1, -1, -1, 1,  1,  -1, -1, -1, 1,  1,  -1, 1,  1,  0};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 36: {
                dstr = "second double conference design in DC(20,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (20, 4, 0);
                int tmp[] = {0, 0,  1, 1,  1, 1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             1, -1, 0, 1,  1, 1,  1,  -1, -1, -1, -1, 0,  1,  1,  1,  1,  -1, -1, -1, -1,
                             1, -1, 1, 0,  1, -1, -1, 1,  1,  -1, -1, -1, 1,  1,  -1, -1, 0,  1,  1,  -1,
                             1, -1, 1, -1, 1, 0,  -1, -1, -1, 1,  1,  -1, -1, -1, 1,  1,  1,  0,  1,  -1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 37: {
                dstr = "third double conference design in DC(20,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (20, 4, 0);
                int tmp[] = {0, 0,  1, 1,  1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             1, -1, 0, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  1,  1,  1,  1,  -1, -1, -1, -1,
                             1, -1, 1, 0,  1,  -1, -1, 1,  1,  -1, -1, -1, 1,  1,  -1, -1, 0,  1,  1,  -1,
                             1, -1, 1, -1, -1, 0,  1,  1,  -1, 1,  -1, -1, 1,  -1, 1,  -1, 1,  0,  -1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 32: {
                dstr = "first double conference design in DC(18,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (18, 4, 0);
                int tmp[] = {0, 0, 1, 1, 1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, -1,
                             0, 0, 1, 1, 1,  1,  -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1,
                             0, 0, 1, 1, -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1,
                             0, 0, 1, 1, -1, -1, -1, -1, 1,  1,  -1, -1, 1,  1,  1,  1,  -1, -1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 33: {
                dstr = "second double conference design in DC(18,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (18, 4, 0);
                int tmp[] = {0, 0, 1, 1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, -1,
                             0, 0, 1, 1,  1,  1,  -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1,
                             0, 0, 1, 1,  -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1,
                             0, 0, 1, -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 34: {
                dstr = "third double conference design in DC(18,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (18, 4, 0);
                int tmp[] = {0, 0,  1, 1, 1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, -1,
                             0, 0,  1, 1, 1,  1,  -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1,
                             0, 0,  1, 1, -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1,
                             1, -1, 0, 0, 1,  -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1, 1,  -1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 31: {
                dstr = "conference design in C(8,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (4, 8, 0);
                int tmp[] = {0, 1, 1,  1, 1, 0,  -1, -1, 1, 1,  0, -1, 1, 1,  1,  0,
                             1, 1, -1, 1, 1, -1, 1,  1,  1, -1, 1, -1, 1, -1, -1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                al = al.transposed ();
                return al;
                break;
        }
        case 30: {
                dstr = "conference design in C(8,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (8, 4, 0);
                int tmp[] = {0, 1,  1, 1, 1,  1, 1, 1,  1, 0,  1, 1,  1, -1, -1, -1,
                             1, -1, 0, 1, -1, 1, 1, -1, 1, -1, 1, -1, 0, 1,  -1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 0: {
                dstr = "array in OA(8,2, 2^2)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (8, 2, 0);
                std::vector< int > s;
                s.push_back (2);
                s.push_back (2);
                arraydata_t ad (s, 8, 2, 2);
                al.create_root (ad);
                return al;
                break;
        }
        case 1: {
                dstr = "array 3 in OA(16, 2, 2^5)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (16, 5, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
                             0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0,
                             1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 2: {
                dstr = "array 6 in OA(16, 2, 2^6)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (16, 6, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                             0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                             0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1,
                             0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 3: {
                dstr = "array ? in OA(32, 3, 2^7)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (32, 7, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                             1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
                             1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                             0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
                             0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
                             0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0,
                             1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                             0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
                             0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 4: {
                dstr = "array 4 in OA(16, 2, 2^7)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                // array 4 in OA(16, 2, 2^7)
                array_link al (16, 7, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                             1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
                             0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,
                             1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 5: {
                dstr = "array 0 in OA(24, 2, 4 3 2^a)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (24, 5, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 2, 2,
                             0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
                             1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1,
                             1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 6: {

                dstr = "array in OA(4, 2, 2^a)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (4, 3, 0);
                int tmp[] = {0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 7: {
                dstr = "array 0 in OA(4, 2, 2^a)?";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (4, 3, 0);
                int tmp[] = {0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 8: {
                dstr = "array in OA(40, 3, 2^7)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (40, 7, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                             0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
                             0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                             0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0,
                             0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0,
                             1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
                             0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1,
                             1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 9: {
                dstr = "array in A(40, 2^7), D-optimal";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (40, 7, 0);
                int tmp[] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                             0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1,
                             1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                             0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1,
                             0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
                             0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0,
                             0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0,
                             1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0,
                             1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1,
                             1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 10: {
                if (verbose) {
                        myprintf ("exampleArray %d: array in OA(9, 3^2)\n", idx);
                }
                array_link al (9, 3, 0);
                int tmp[] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 1, 2, 1, 1, 2, 0, 0, 2, 2, 0, 2, 0, 2, 1, 0, 1, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 11: {
                if (verbose) {
                        myprintf ("exampleArray %d: D-optimal array in OA(44, 2^8)\n", idx);
                }
                array_link al (44, 8, 0);
                int tmp[] = {1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1,
                             0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
                             1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
                             1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                             1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0,
                             1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0,
                             1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
                             1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1,
                             0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                             1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0,
                             0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1,
                             0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 12: {
                if (verbose) {
                        myprintf ("exampleArray %d: even-odd array OA(64, 2^13)\n", idx);
                }
                array_link al (64, 13, 0);
                int tmp[] = {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                    0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0,
                    1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
                    0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0,
                    1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1,
                    0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1,
                    0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1,
                    0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1,
                    0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1,
                    0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0,
                    1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1,
                    0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0,
                    1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
                    0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0,
                    1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 13: {
                dstr = "array in OA(25, 2^5)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (24, 5, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                             1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1,
                             0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1,
                             0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 14: {
                dstr = "design in D(28, 2^5), D-efficiency is low";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (28, 5, 0);
                int tmp[] = {1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1,
                             0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0,
                             1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1,
                             1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0,
                             1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 15: {
                dstr = "design in D(56, 2^10), D-efficiency is low";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }

                array_link al (56, 10, 0);
                int tmp[] = {
                    0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
                    1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
                    0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0,
                    0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
                    1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0,
                    0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0,
                    0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,
                    1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0,
                    1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0,
                    0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1,
                    0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0,
                    1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1,
                    0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0,
                    0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 16: {
                dstr = "array in OA(32, 2, 2^5)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (32, 5, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                             0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0,
                             0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                             0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1,
                             0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 17: {
                dstr = "unique array in OA(64, 4, 2^7)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (64, 7, 0);
                int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                             0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                             0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0,
                             0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                             0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
                             1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
                             0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                             0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                             0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0,
                             0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1,
                             0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 18: {
                dstr = "conference matrix of size 16, 7";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (7, 16, 0);
                int tmp[] = {0,  1,  1,  1, 1,  1,  1,  1,  0,  -1, -1, -1, -1, -1, 1,  1,  0,  -1, -1, -1, 1, 1,  1,
                             1,  0,  -1, 1, -1, 1,  1,  1,  1,  0,  -1, -1, 1,  1,  1,  -1, 1,  0,  1,  1,  1, -1, 1,
                             1,  -1, 0,  1, 1,  -1, 1,  -1, 1,  1,  1,  1,  -1, -1, 1,  1,  -1, 1,  -1, 1,  1, 1,  1,
                             -1, 1,  -1, 1, 1,  -1, -1, -1, 1,  -1, 1,  -1, 1,  -1, 1,  1,  -1, 1,  -1, -1, 1, 1,  1,
                             -1, -1, 1,  1, -1, 1,  1,  -1, -1, 1,  -1, 1,  1,  1,  -1, -1, -1, 1,  1,  -1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                al = al.transposed ();
                return al;
                break;
        }

        case 19: {
                dstr = "conference matrix of size 4, 3";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (3, 4, 0);
                int tmp[] = {0, 1, 1, 1, 0, -1, 1, 1, 0, 1, -1, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                al = al.transposed ();
                return al;
                break;
        }

        case 20: {
                dstr = "first LMC-0 double conference matrix in DC(24,3)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (24, 3, 0);
                int tmp[] = {0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, 1,  -1, 0,  1,  1,  1,  1,  1,  -1, -1, -1, -1,
                             -1, 0,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, 1,  -1, -1, 0,  1,  1,
                             -1, -1, 1,  1,  1,  -1, -1, 1,  1,  1,  -1, -1, -1, 0,  1,  1,  -1, -1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 21: {
                dstr = "second LMC-0 double conference matrix in DC(16,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (16, 4, 0);
                int tmp[] = {0,  0,  1,  1,  1, 1,  1, 1,  1,  -1, -1, -1, -1, -1, -1, -1, 1, -1, 0,  1, 1, 1,
                             -1, -1, -1, 0,  1, 1,  1, -1, -1, -1, 1,  -1, -1, 0,  1,  -1, 1, 1,  -1, 1, 1, -1,
                             -1, 0,  1,  -1, 1, -1, 1, -1, 1,  -1, 0,  -1, 1,  -1, -1, 0,  1, 1,  1,  -1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 22: {
                dstr = "LMC-0 double conference matrix in DC(32,4)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (32, 4, 0);
                int tmp[] = {0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1,  -1, 0,  1,  1,  1,  1,  1,  1,  1,  -1, -1,
                             -1, -1, -1, -1, -1, 0,  1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, 1,  -1,
                             -1, 0,  1,  1,  1,  -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1,
                             -1, 0,  1,  1,  1,  -1, -1, -1, 1,  -1, -1, -1, 0,  1,  -1, 1,  1,  -1, 1,  1,  -1, -1,
                             1,  1,  -1, 1,  1,  -1, -1, 1,  1,  -1, -1, 1,  1,  -1, -1, 0,  1,  -1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 23: {
                dstr = "LMC-0 double conference matrix in DC(32,6)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (32, 6, 0);
                int tmp[] = {0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1,  -1, 0,  1,  1,  1,  1,  1,  1,  1,  -1, -1,
                             -1, -1, -1, -1, -1, 0,  1,  1,  1,  1,  1,  1,  1,  -1, -1, -1, -1, -1, -1, -1, 1,  -1,
                             -1, 0,  1,  1,  1,  -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1,
                             -1, 0,  1,  1,  1,  -1, -1, -1, 1,  -1, -1, -1, 0,  1,  1,  1,  -1, -1, 1,  -1, -1, -1,
                             1,  1,  1,  1,  -1, -1, -1, 1,  1,  1,  -1, 1,  1,  1,  -1, 0,  -1, -1, 1,  -1, -1, -1,
                             -1, 0,  1,  -1, 1,  1,  -1, 1,  1,  -1, 1,  1,  -1, 1,  1,  -1, -1, 1,  -1, -1, 1,  1,
                             -1, -1, 1,  1,  0,  -1, 1,  -1, -1, 1,  -1, 1,  -1, 0,  1,  -1, -1, 1,  -1, 1,  1,  -1,
                             1,  1,  -1, 1,  -1, -1, 1,  -1, 1,  -1, 1,  -1, 0,  1,  -1, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 24: {
                dstr = "design in OA(64, 3, 2^16) (even-odd)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (64, 16, 0);
                int tmp[] = {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1,
                    0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                    1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1,
                    0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1,
                    0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0,
                    1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0,
                    1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0,
                    1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1,
                    0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1,
                    0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0,
                    1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1,
                    0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0,
                    1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0,
                    0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                    1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
                    0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                    1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1,
                    0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0,
                    1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1,
                    0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0,
                    1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 25: {
                dstr = "design in OA(64, 3, 2^16) (even-odd)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (64, 16, 0);
                int tmp[] = {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
                    0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0,
                    0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1,
                    0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1,
                    0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1,
                    0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0,
                    0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1,
                    0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0,
                    0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1,
                    0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1,
                    0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0,
                    0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0,
                    0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0,
                    0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1,
                    1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
                    1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
                    1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0,
                    1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1,
                    1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0,
                    1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0,
                    1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1,
                    1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0,
                    1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1,
                    1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0,
                    1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1,
                    1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1,
                    1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                    1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 26: {
                dstr = "design in OA(64, 3, 2^16) (even-odd)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (16, 64, 0);
                int tmp[] = {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1,
                    0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0,
                    0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0,
                    0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0,
                    0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1,
                    0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1,
                    0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1,
                    0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0,
                    0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1,
                    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0,
                    0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1,
                    0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0,
                    0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1,
                    0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                    1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0,
                    1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1,
                    1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0,
                    1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0,
                    1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0,
                    1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1,
                    1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0,
                    1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1,
                    1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
                    1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1,
                    1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1,
                    1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                    1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                al = al.transposed ();
                return al;
                break;
        }
        case 27: {
                dstr = "design in OA(64, 3, 2^16) (even-odd)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (16, 64, 0);
                int tmp[] = {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0,
                    0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1,
                    0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0,
                    0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1,
                    0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1,
                    0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0,
                    0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0,
                    0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1,
                    0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0,
                    0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1,
                    0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0,
                    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1,
                    0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
                    1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0,
                    1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1,
                    1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0,
                    1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0,
                    1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0,
                    1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                    1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0,
                    1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1,
                    1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1,
                    1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
                    1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                    1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0};

                al.setarraydata (tmp, al.n_rows * al.n_columns);
                al = al.transposed ();
                return al;
                break;
        }
        case 28: {
                dstr = "conference design in C(4, 3) in LMC0 form";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (4, 3, 0);
                int tmp[] = {0, 1, 1, 1, 1, 0, 1, -1, 1, -1, 0, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 29: {
                dstr = "conference design in C(4, 3)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (4, 3, 0);
                int tmp[] = {0, 1, 1, 1, 1, 0, 1, -1, 1, 1, -1, 0};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }

        case 40: {
                dstr = "first conference design in C(14, 5)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (14, 5, 0);
                int tmp[] = {0,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1, 1,  1, 0,  1,  1,  1,  1, 1,  1, -1, -1,
                             -1, -1, -1, -1, 1,  1, 0, 1,  1, -1, -1, -1, 1, 1,  1, -1, -1, -1, 1,  1, 1,  0, -1, 1,
                             -1, -1, 1,  -1, -1, 1, 1, -1, 1, 1,  1,  -1, 0, -1, 1, -1, -1, 1,  -1, 1, -1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 41: {
                dstr = "second conference design in C(14, 5)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (14, 5, 0);
                int tmp[] = {0,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1, 1, 1, 0,  1,  1,  1,  1, 1,  1, -1, -1,
                             -1, -1, -1, -1, 1,  1, 0, 1,  1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1,  1, 1,  0, -1, 1,
                             -1, -1, 1,  -1, -1, 1, 1, -1, 1, 1,  -1, -1, 1, 1, 0, -1, -1, 1,  -1, 1, -1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
        case 42: {
                dstr = "third conference design in C(14, 5)";
                if (verbose) {
                        myprintf ("exampleArray %d: %s\n", idx, dstr.c_str ());
                }
                array_link al (14, 5, 0);
                int tmp[] = {0,  1,  1,  1,  1,  1, 1, 1,  1, 1,  1,  1,  1, 1,  1, 0,  1,  1,  1,  1, 1,  1, -1, -1,
                             -1, -1, -1, -1, 1,  1, 0, 1,  1, -1, -1, -1, 1, 1,  1, -1, -1, -1, 1,  1, 1,  0, -1, 1,
                             -1, -1, 1,  -1, -1, 1, 1, -1, 1, 1,  -1, -1, 1, -1, 0, 1,  1,  -1, -1, 1, -1, 1};
                al.setarraydata (tmp, al.n_rows * al.n_columns);
                return al;
                break;
        }
		case 43: {
			dstr = "2x2 array with zeros and a singe value -1";
			if (verbose) {
				myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			}
			array_link al(2, 2, 0);
			int tmp[] = { 0,  0,  0,  -1 };
			al.setarraydata(tmp, al.n_rows * al.n_columns);
			return al;
			break;
		}
		 case 44: {
                        dstr = "D-optimal strength 3 ortogonal array in OA(40,3, 2^7)";
                        if (verbose) {
                                myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
                        }
                         array_link array (40, 7, 0); 
                         int array_data_tmp[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1,0,0,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1,1,0,0,1,0,1,0,1,0,1,1,0,1,0,1,1,1,1,0,0,0,0,1,1,1,0,0,1,1,0,0,0,1,1,0,0,0,0,1,1,1,0,0,1,0,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,0,1,1,1,0,1,0,0,0,1,0,1,0,1,0,0,1,1};
                         array.setarraydata (array_data_tmp, array.n_rows * array.n_columns);
                         return array;
                         break;
                }
		 case 45: {
			 dstr = "first conference design in C(20,8)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(20, 8, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,0,1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,0,1,-1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,0,1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,1,1,-1,-1,-1,1,-1,0,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,0,-1,1,1,1,1,-1,1,-1,-1,1,-1,1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 46: {
			 dstr = "second conference design in C(20,8)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(20, 8, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,0,1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,0,1,-1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,0,1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,1,1,-1,-1,-1,1,-1,0,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,1,0,1,-1,1,1,-1,1,-1,-1,-1,1,1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 47: {
			 dstr = "third conference design in C(20,8)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(20, 8, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,0,1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,0,1,-1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,0,1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,1,1,-1,-1,-1,1,-1,0,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,1,0,1,-1,1,1,-1,-1,1,1,-1,-1,1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 48: {
			 dstr = "last conference design in C(20,8)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(20, 8, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,0,1,1,-1,1,1,-1,-1,1,1,-1,-1,-1,1,1,1,-1,1,-1,1,-1,1,1,-1,0,-1,1,-1,-1,-1,1,1,-1,1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,-1,1,0,1,1,1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,1,-1,-1,-1,1,-1,-1,0,1,-1,1,1,1,-1,1,1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,1,0,1,-1,-1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 49: {
			 dstr = "array 4347 C(20,8)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(20, 8, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,0,1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,0,1,1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,0,1,-1,-1,-1,1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,0,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,-1,1,1,1,-1,-1,-1,-1,1,0,1,1,-1,1,-1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 50: {
			 dstr = "array 4506 C(20,8)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(20, 8, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,0,1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,0,1,1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,1,-1,-1,1,-1,0,1,-1,1,1,1,-1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,0,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,-1,1,1,-1,1,0,1,-1,-1,1,-1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 51: {
			 dstr = "first array in C(12,4)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(12, 4, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,-1,-1,-1,-1,-1,1,-1,0,1,1,-1,-1,1,1,1,-1,-1,1,-1,-1,0,-1,1,1,1,1,-1,1,-1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 52: {
			 dstr = "second array in C(12,4)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(12, 4, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,-1,-1,-1,-1,-1,1,-1,0,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,-1,1,-1,0,1,-1,1,-1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }
		 case 53: {
			 dstr = "third array in C(12,4)";
			 if (verbose) {
				 myprintf("exampleArray %d: %s\n", idx, dstr.c_str());
			 }
			 array_link array(12, 4, 0);
			 int array_data_tmp[] = { 0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,-1,-1,-1,-1,-1,1,-1,0,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,0,1,1,-1,-1 };
			 array.setarraydata(array_data_tmp, array.n_rows * array.n_columns);
			 return array;
		 }

        } // end of switch

        return array_link (1, 1, -1);
}

array_link array_link::reduceDOP () const {
        array_link d = reduceDOPform (*this, 0);
        return d;
}

MatrixFloat array_link::getEigenMatrix() const {
	int ncolumns = this->n_columns;
	int nrows = this->n_rows;
	MatrixFloat mymatrix = MatrixFloat::Zero(nrows, ncolumns);

	for (int column = 0; column < ncolumns; ++column) {
		int column_offset = column * nrows;
		for (int row = 0; row < nrows; ++row) {
			mymatrix(row, column) = this->array[row + column_offset];
		}
	}
	return mymatrix;
}

int array_link::columnGreater (int c1, const array_link &rhs, int rhs_column) const {

          if ((this->n_rows != rhs.n_rows) || c1 < 0 || rhs_column < 0 || (c1 > this->n_columns - 1)) {
               myprintf ("array_link::columnGreater: warning: comparing arrays with different sizes\n");
               return 0;
          }

          int n_rows = this->n_rows;
          return std::lexicographical_compare (rhs.array + rhs_column * n_rows, rhs.array + rhs_column * n_rows + n_rows,
                                             array + c1 * n_rows, array + c1 * n_rows + n_rows);
}
        
array_link array_link::reduceLMC () const {
        int strength = this->strength ();
        arraydata_t ad = arraylink2arraydata (*this, 0, strength);
        LMCreduction_t reduction (&ad);
        reduction.mode = OA_REDUCE;

        OAextend oaextend;
        oaextend.setAlgorithm (MODE_ORIGINAL, &ad);

        reduction.init_state = SETROOT;

        LMCcheck (*this, ad, oaextend, reduction);

        return reduction.getArray ();
}

symmetry_group array_link::row_symmetry_group () const {

        const int nc = this->n_columns;

        std::vector< mvalue_t< int > > row_elements (this->n_rows);
        for (int i = 0; i < this->n_rows; i++) {
                mvalue_t< int > &m = row_elements[i];
                m.values.resize (nc);

                for (int k = 0; k < nc; k++) {
                        m.values[k] = this->atfast (i, k);
                }
        }
        symmetry_group row_symmetry_group (row_elements, true, 0);
        return row_symmetry_group;
}

double array_link::nonzero_fraction () const {
        double nz = 0;
        long nn = this->n_columns * this->n_rows;
        for (int i = 0; i < nn; i++)
                if (this->array[i] != 0) {
                        nz++;
                }
        return nz / nn;
}

void array_link::clear() { std::fill(array, array + n_rows * n_columns, 0); }

void array_link::setcolumn (int target_column, const array_link &source_array, int source_column ) const {
          assert (target_column >= 0);
		  assert(target_column <= this->n_columns);
          assert (this->n_rows == source_array.n_rows);
          std::copy (source_array.array + source_column * this->n_rows, source_array.array + (source_column + 1) * this->n_rows,
                    this->array + this->n_rows * target_column);
}
        
bool array_link::is_conference (int number_of_zeros) const {
        if (!this->is_conference ()) {
                return false;
        };

        for (int c = 0; c < this->n_columns; c++) {
                int n = 0;
                for (int r = 0; r < this->n_rows; r++) {
                        if (this->atfast (r, c) == 0)
                                n++;
                }
                if (n != number_of_zeros)
                        return false;
        }
        return true;
}

bool array_link::is_orthogonal_array() const {
	if (this->min() < 0)
		return false;
	return true;
}

bool array_link::is_mixed_level() const {
	if (this->min() < 0)
		return false;
	array_t max_value = this->max();
	const int N = this->n_rows;
	int s = max_value + 1;
	
	int k;
	for (k = 0; k < this->n_columns; k++) {
		array_t *max_column_value = std::max_element(this->array + N * k, this->array + (N * (k + 1)));
		if (*max_column_value + 1 != s)
			break;
	}
	int is_mixed = (k < this->n_columns);
	return is_mixed;
}

bool array_link::is_conference () const {
        int m = this->min ();
        if (m < -1) {
                return false;
        }
        if (this->max () > 1) {
                return false;
        }
        return true;
}

bool array_link::is2level () const {
        int N = this->n_rows;
        for (int r = 0; r < this->n_rows; r++) {
                for (int c = 0; c < this->n_columns; c++) {
                        if (this->array[r + c * N] < 0) {
                                return false;
                        }
                        if (this->array[r + c * N] > 1) {
                                return false;
                        }
                }
        }
        return true;
}

void array_link::showproperties () const {
        myprintf ("array: %d rows, %d cols\n", this->n_rows, this->n_columns);
        if (this->min () >= 0) {
                myprintf ("  strength %d, rank %d\n", this->strength (), this->rank ());
        }
        myprintf ("  D-efficiency %.3f\n", this->Defficiency ());
        std::vector< double > gwlp = this->GWLP ();
#ifdef FULLPACKAGE
        myprintf ("  GWLP ");
        printf_vector (gwlp, "%d", " ");
        myprintf ("\n");
#endif
        return;
}

void array_link::debug () const { myprintf ("debug: %ld %p %p", (long)this->array, (void *)array, (void *)array); }
#ifdef SWIGCODE
void *array_link::data () {
        return ((void *)(this->array));
}
#else
#endif

bool array_link::equal_size(const array_link &rhs_array) const {
     if( (this->n_rows != rhs_array.n_rows) || (this->n_columns != rhs_array.n_columns) ) return false;
     else return true;
}

void array_link::showarraycompact () const {
        for (int r = 0; r < this->n_rows; r++) {
                for (int c = 0; c < this->n_columns; c++) {
                        myprintf ("%d", this->_at (r, c));
                }
                myprintf ("\n");
        }
}

void array_link::showarray () const {
        std::string out = this->showarrayString();
        myprintf("%s", out.c_str());
}

std::string array_link::showarrayString () const {
        std::stringstream ss;
        ss << "array:\n";
        write_array_format (ss, array, this->n_rows, this->n_columns);
        return ss.str ();
}

void perform_column_permutation (const array_link source, array_link &target, const std::vector< int > perm) {
        int ncols = source.n_columns;
        int nrows = source.n_rows;
        for (int i = 0; i < ncols; i++) {
                memcpy (target.array + (perm[i] * nrows), source.array + i * nrows, nrows * sizeof (array_t));
        }
}

void perform_row_permutation (const array_link source, array_link &target, const std::vector< int > perm) {
        int ncols = source.n_columns;
        int nrows = source.n_rows;
        for (int i = 0; i < ncols; i++)
                for (int j = 0; j < nrows; j++) {
                        target.array[nrows * i + perm[j]] = source.array[nrows * i + j];
                }
}

void array_link::create_root (const arraydata_t &ad, int fill_value) {
        if (!(ad.N <= this->n_rows)) {
                myprintf ("array_link::create_root: number of columns too small for root of size %d\n", ad.N);
                return;
        }

        if (!(ad.strength <= this->n_columns)) {
                myprintf ("array_link::create_root: number of columns (%d) too small for specified strength %d\n",
                          this->n_columns, ad.strength);
                return;
        }
        myassert (ad.strength <= this->n_columns, "strength invalid");
        std::fill (this->array, this->array + this->n_rows * this->n_columns, fill_value);
        ::create_root ((this->array), &ad);
}

/*!
  Creates the root of an OA according to the strength, the number of rows and the levels per column.
  \brief Create Root
  \param array Pointer to OA, where the root is placed in
  \param ad Parameter with data of the array
  */
void create_root (array_t *array, const arraydata_t *ad) {
        int steps = ad->N;
        for (colindex_t i = 0; i < ad->strength; i++) { /* loop over all root columns */
                steps /= ad->s[i]; // adjust nr of steps per collumn
                if (steps == 0) {
                        steps = 1;
                }
                array_t l = 0;
                int coloffset = i * ad->N;
                for (int j = 0; j < ad->N; j += steps) {  // big steps
                        for (int k = 0; k < steps; k++) { // small steps
                                array[coloffset + j + k] = l;
                        }
                        l++;
                        if (l == ad->s[i]) {
                                l = 0; // reset value, it has looped
                        }
                }
        }
}

/*!
  Creates the root of an OA. The root is appended to the current list of arrays
  \brief Create Root
  \param ad Pointer to arraydata_t structure
  \param solutions List of current solutions
  \sa create_root
  */
void create_root (const arraydata_t *ad, arraylist_t &solutions) {
        array_link cur_solution = array_link (ad->N, ad->strength, 1);
        create_root (cur_solution.array, ad);
        solutions.push_back (cur_solution);
}

std::vector< int > numberModelParametersConference(const array_link &conference_design)
{
	int n_columns = conference_design.n_columns;
	std::vector< int > number_parameters(4);
	number_parameters[0] = 1;
	number_parameters[1] = n_columns;
	number_parameters[2] = n_columns * (n_columns - 1) / 2;
	number_parameters[3] = n_columns;
	return number_parameters;
}
std::vector< int > numberModelParams(const array_link &array, int order)
{
	if (order > 0) {
		myprintf("numberModelParams: warning: order argument is not used any more\n");
                exit(0);
	}

	if (array.is_conference()) {
		return numberModelParametersConference(array);
	}
	if (array.is_orthogonal_array()) {
		int k = array.n_columns;
		std::vector< int > number_parameters(4);
		number_parameters[0] = 1;
		number_parameters[1] = array.n_columns;

		arraydata_t arrayclass = arraylink2arraydata(array, 0, 0);
		std::vector< int > s = arrayclass.factor_levels();
		std::vector< int > df = s;
		std::transform(df.begin(), df.end(), df.begin(), std::bind2nd(std::minus< int >(), 1.0));

		/* main effect contrasts */
		int mesize = std::accumulate(df.begin(), df.end(), 0);
		number_parameters[1] = mesize;

		/* 2fi*/
		int number_of_2factor_interactions = 0;
		for (int ii = 0; ii < k - 1; ii++) {
			for (int jj = ii + 1; jj < k; jj++) {
				number_of_2factor_interactions += df[ii] * df[jj];
			}
		}
		number_parameters[2] = number_of_2factor_interactions;

		/* quadratic effects*/
		int number_of_quadratic_interactions = 0;
		for (int ii = 0; ii < k - 1; ii++) {
			number_of_quadratic_interactions += df[ii] * df[ii];
		}
		number_parameters[3] = number_of_quadratic_interactions;

		return number_parameters;
	}
	throw_runtime_exception("array type invalid");

        return std::vector<int>();
}

MatrixFloat array_link::getModelMatrix (int order, int intercept, int verbose) const {
        int N = this->n_rows;
        MatrixFloat intcpt = MatrixFloat::Zero (N, 1);
        intcpt.setConstant (1);

        std::vector< int > np = numberModelParams (*this);
        if (verbose >= 1) {
                myprintf ("array_link::getModelMatrix numberModelParams: ");
                printf_vector (np, "%d", " ");
                myprintf ("\n");
        }

        std::pair< MatrixFloat, MatrixFloat > mmx = array2eigenModelMatrixMixed (*this, 0);

        if (verbose >= 2) {
                myprintf ("array_link::model matrices\n");
                eigenInfo (mmx.first, "mmx.first");
                eigenInfo (mmx.second, "mmx.second");
        }
        if (order == 0) {
                return intcpt;
        }
        if (order == 1) {
                MatrixFloat mm (N, np[0] + np[1]);
                mm << intcpt, mmx.first;
                return mm;
        }
        if (order == 2) {
                MatrixFloat mm (N, np[0] + np[1] + np[2]);
                mm << intcpt, mmx.first, mmx.second;
                return mm;
        }

        throw_runtime_exception ("array_link::getModelMatrix: order > 2 not supported!\n");
        return intcpt;
}

double array_link::CL2discrepancy () const { return ::CL2discrepancy (*this); }

array_t array_link::max () const {
        if (this->n_rows * this->n_columns == 0) {
                return 0;
        }
        return *std::max_element (this->array, this->array + this->n_rows * this->n_columns);
}
array_t array_link::min () const {
        if (this->n_rows * this->n_columns == 0) {
                return 0;
        }
        return *std::min_element (this->array, this->array + this->n_rows * this->n_columns);
}

bool array_link::foldover () const {
        std::vector< double > g = this->GWLP ();
        for (size_t i = 1; i < g.size (); i += 2) {
                if (g[i] != 0) {
                        return false;
                }
        }
        return true;
}

std::vector< double > array_link::GWLP (int truncate, int verbose) const { return ::GWLP (*this, verbose, truncate); }

/// convert array to second order interaction matrix in Eigen format
MatrixFloat array2eigenX2 (const array_link &al) {
        int k = al.n_columns;
        int n = al.n_rows;

        MatrixFloat mymatrix = MatrixFloat::Zero (n, k * (k - 1) / 2);

        // init interactions
        int ww = 0;
        for (int c = 0; c < k; ++c) {
                for (int c2 = 0; c2 < c; ++c2) {
                        int ci = c * n;
                        int ci2 = c2 * n;

                        for (int r = 0; r < n; ++r) {
                                mymatrix (r, ww) = (al.array[r + ci] + al.array[r + ci2]) % 2;
                        }
                        ww++;
                }
        }

        mymatrix.array () -= .5;
        mymatrix.array () *= 2;

        return mymatrix;
}

// helper function for Python interface
void eigen2numpyHelper (double *pymat1, int n, const MatrixFloat &m) {
        myprintf ("pymat1 %p\n", (void *)pymat1);
        std::copy (m.data (), m.data () + m.size (), pymat1);
}

void eigenInfo (const MatrixFloat m, const char *str, int verbose) {
        if (verbose == 1) {
                myprintf ("%s: %dx%d\n", str, (int)m.rows (), (int)m.cols ());
        }
        if (verbose == 2) {
                myprintf ("%s: %dx%d\n", str, (int)m.rows (), (int)m.cols ());
        }
}
void eigenInfoF (const Eigen::MatrixXf m, const char *str, int verbose) {
        if (verbose == 1) {
                myprintf ("%s: %dx%d\n", str, (int)m.rows (), (int)m.cols ());
        }
        if (verbose == 2) {
                myprintf ("%s: %dx%d\n", str, (int)m.rows (), (int)m.cols ());
        }
}

void print_eigen_matrix(const MatrixFloat matrix) {
     std::stringstream buffer;
     buffer << matrix <<  std::endl;
     myprintf("%s", buffer.str().c_str()); 
}

MatrixFloat array2eigenMainEffects (const array_link &al, int verbose) {
        std::pair< MatrixFloat, MatrixFloat > mm = array2eigenModelMatrixMixed (al, verbose);
        return mm.first;
}

// code from Eric Schoen, adapted to work for arrays of strength < 1
std::pair< MatrixFloat, MatrixFloat > array2eigenModelMatrixMixed (const array_link &al, int verbose) {
        if ( !(al.min () >= 0) )
             throw_runtime_exception("array cannot have negative elements");

        if (verbose >= 2) {
                myprintf ("array2eigenModelMatrixMixed: start");
        }

        int N = al.n_rows;
        int k = al.n_columns;
        arraydata_t arrayclass = arraylink2arraydata (al, 0, 0);
        std::vector< int > s = arrayclass.factor_levels ();

        std::vector< int > df = s;
        std::transform (df.begin (), df.end (), df.begin (), std::bind2nd (std::minus< int > (), 1.0));

        if (verbose >= 2) {
                arrayclass.show ();
                myprintf ("array2eigenModelMatrixMixed: N %d, k %d\n", N, k);
                myprintf ("df ");
                printf_vector (df, "%d ");
                myprintf ("\n");
        }
        MatrixFloat AA = al.getEigenMatrix ();

        /* main effect contrasts */
        int mesize = std::accumulate (df.begin (), df.end (), 0);
        if (verbose >= 2) {
                printfd ("main effects: size %dx%d\n", N, mesize);
        }

        std::vector< int > np = numberModelParams (al);
        MatrixFloat main_effects = MatrixFloat::Zero (N, mesize);

        int meoffset = 0;
        for (int column = 0; column < k; column++) {
                int md = df[column];
                MatrixFloat Z = MatrixFloat::Zero (N, md + 1); // large tmp buffer

                for (int ii = 0; ii < md + 1; ii++) {
                        for (int r = 0; r < N; r++) {
                                Z (r, ii) = AA (r, column) > (ii - 1);
                        }
                }

                // make Helmert contrasts (these are automatically orthogonal)
                for (int r = 0; r < N; r++) {
                        int array_value = AA (r, column);
                        Z (r, 0) = 1;
                        if (array_value > 0) {
                                Z (r, array_value) = array_value;
                        }
                        for (int q = 1; q < array_value; q++) {
                                Z (r, q) = 0;
                        }
                        for (int q = array_value + 1; q < md + 1; q++) {
                                Z (r, q) = -1;
                        }
                }

               if (verbose>=3) {
                    eigenInfo (Z, "Z (before normalization)\n");
                    print_eigen_matrix(Z);
               }
               for (int ii = 0; ii < md; ii++) {
                       
                        MatrixFloat tmp_norm = Z.col (ii + 1).transpose () * Z.col (ii + 1);

                        if (verbose >= 3) {
                             myprintf(" normalize: factor %.2f (%.2f, tmp_norm %.1f tmp_norm sqrt %.1f)\n", sqrt (double(N))/sqrt (double(tmp_norm (0, 0))), sqrt (double(N)), double(tmp_norm(0,0)), sqrt (double(tmp_norm (0, 0)) ) );
                        }
                        main_effects.col (meoffset + ii) =
                            sqrt (double(N)) * Z.col (ii + 1) / sqrt (double(tmp_norm (0, 0)));
                }

#ifdef FULLPACKAGE
                if (verbose >= 2) {
                        eigenInfo (Z, "Z");
                        print_eigen_matrix(Z);
                }
#endif
                meoffset += md;
        }

        /* 2fi */
        if (verbose >= 2) {
                myprintf ("2fi\n");
        }

        MatrixFloat tfi = MatrixFloat::Zero (N, np[2]);

        if (verbose >= 2) {
                printfd ("create 2fi\n");
        }
        int tel = 0;
        int n = al.n_columns;
        int po = 0, qo = 0; // offsets
        for (int column1 = 0; column1 < n - 1; column1++) {
                int n1 = df[column1];
                po = std::accumulate (df.begin (), df.begin () + column1, 0);

                for (int column2 = column1 + 1; column2 < n; column2++) {
                        int n2 = df[column2];

                        qo = std::accumulate (df.begin (), df.begin () + column2, 0);

                        if (verbose >= 3) {
                                printfd ("create 2fi: p0=o %d, qo %d\n", po, qo);
                        }

                        for (int pp = 0; pp < n1; pp++) {
                                for (int qq = 0; qq < n2; qq++) {
                                        tfi.col (tel) =
                                            main_effects.col (pp + po).cwiseProduct (main_effects.col (qq + qo));
                                        tel++;
                                }
                        }
                        qo += n2;
                }
                po += n1;
        }

        if (verbose >= 2) {
                myprintf ("done\n");
        }

        return std::pair< MatrixFloat, MatrixFloat > (main_effects, tfi);
}

MatrixFloat array2eigenX1 (const array_link &al, int intercept) {

        int k = al.n_columns;
        int n = al.n_rows;

        intercept = intercept > 0;
        MatrixFloat mymatrix = MatrixFloat::Zero (n, k + intercept);

        int ww = 0;
        if (intercept) {

                // init first column
                for (int r = 0; r < n; ++r) {
                        mymatrix (r, ww) = 1;
                }
                ww += 1;
        }
        // init array
        for (int c = 0; c < k; ++c) {
                int ci = c * n;
                for (int r = 0; r < n; ++r) {
                        mymatrix (r, ww + c) = al.array[r + ci];
                }
        }

        mymatrix.array () -= .5;
        mymatrix.array () *= 2;

        return mymatrix;
}

/// specialized version array2eigenModelMatrix returning integer valued matrix
Eigen::MatrixXi array2eigenModelMatrixInt (const array_link &al) {
        int k = al.n_columns;
        int n = al.n_rows;
        int m = 1 + k + k * (k - 1) / 2;

        if (n * k > 0) {
			myassert(al.max() < 2, "array should be 2-level array");
        }

        // create data in integer type (we are working with 2-level arrays, convert them later */
        Eigen::MatrixXi mymatrix = Eigen::MatrixXi::Zero (n, m);
        int *modelmatrix_raw_data = mymatrix.data ();

        // init first column
        int ww = 0;
        for (int r = 0; r < n; ++r) {
                mymatrix (r, ww) = 1;
        }

        // init array
        ww = 1;
        for (int c = 0; c < k; ++c) {
                int ci = c * n;
                std::copy (al.array + ci, al.array + ci + n, modelmatrix_raw_data + (ww + c) * n);
        }

        // init interactions
        ww = k + 1;

        for (int c = 0; c < k; ++c) {
                for (int c2 = 0; c2 < c; ++c2) {
                        int ci = c * n;
                        int ci2 = c2 * n;

                        for (int r = 0; r < n; ++r) {
                                modelmatrix_raw_data[r + ww * n] = (al.array[r + ci] + al.array[r + ci2]) % 2;
                        }
                        ww++;
                }
        }

        mymatrix.array () *= 2;
        mymatrix.array () -= 1;

        return mymatrix;
}

MatrixFloat array2eigenModelMatrix (const array_link &al) {
        return array2eigenModelMatrixInt (al).cast< eigenFloat > ();
}

double array_link::DsEfficiency (int verbose) const {
        if (!this->is2level ()) {
                return this->Defficiencies ()[1];
        }

        const array_link &al = *this;
        int k = al.n_columns;
        int n = al.n_rows;
        int n2fi = k*(k-1)/2;

        MatrixFloat X = array2eigenModelMatrix (al);
        MatrixFloat X2 = X.block(0,1+k, n, n2fi);
        MatrixFloat X02 (n, 1 + n2fi);
        X02 << X.block (0, 0, n, 1), X2;
        
        MatrixFloat matXtX = (X.transpose () * X / n);
        double f1 = matXtX.determinant ();
        double f2i = (X02.transpose () * X02 / n).determinant ();
        double Ds = 0;
        if (fabs (f1) < 1e-15) {
                if (verbose) {
                        myprintf ("DsEfficiency: f1 < 1e-15, setting Ds to zero\n");
                }
        } else {
                Ds = pow ((f1 / f2i), 1. / k);
        }
        if (verbose) {
                myprintf ("f1 %e, f2i %e, Ds %f\n", f1, f2i, Ds);
        }

        return Ds;
}

const double NaN = std::numeric_limits< double >::quiet_NaN ();

std::vector< double > array_link::Defficiencies (int verbose, int addDs0) const {
        const array_link &al = *this;

        if (this->min () < 0) {
                myprintf ("Defficiencies: error: input array should have elements ranging from 0 to s-1\n");
                std::vector< double > x (3 + addDs0);
                x[0] = NaN;
                x[1] = NaN;
                return x;
        }
        if (this->n_rows < 1) {
                myprintf ("Defficiencies: error: input array should have more than zero rows\n");
                std::vector< double > x (3 + addDs0);
                x[0] = NaN;
                x[1] = NaN;
                return x;
        }
        arraydata_t arrayclass = arraylink2arraydata (al);

        return ::Defficiencies (al, arrayclass, verbose, addDs0);
}

double array_link::Defficiency () const {
        if (!this->is2level ()) {
                return this->Defficiencies ()[0];
        }

        return ::Defficiency (*this);
}
double array_link::VIFefficiency () const { return ::VIFefficiency (*this); }
double array_link::Aefficiency () const {
        return ::Aefficiency (*this);
}

double array_link::Eefficiency () const { return ::Eefficiency (*this); }

std::vector< int > array_link::Fvalues (int jj) const {
        /// NOTE: this assumes the strength is at least 3
        jstruct_t js (*this, jj);
        std::vector< int > FF = js.calculateF ();
        return FF;
}

std::vector< int > array_link::FvaluesConference (int jj) const {
        if (! this->is_conference ())
             throw_runtime_exception("array is not a conference design");

        const int N = this->n_rows;
        jstructconference_t js (*this, jj);
        std::vector< int > FF = js.calculateF ();
        return FF;
}

std::vector< double > array_link::PECsequence (int verbose) const { return ::PECsequence (*this, verbose); }

std::vector< double > array_link::PICsequence(int verbose) const { return ::PICsequence(*this, verbose); }

/** Calculate J-characteristics of matrix (the values are signed)
 *
 * The actual calculation depends on the type of array (2-level or conference)
 * 
 * \param jj Number of columns
 * \returns Vector with calculated Jk-characteristics
 */
std::vector< int > array_link::Jcharacteristics (int jj) const {
        if (this->is2level ()) {
                return ::Jcharacteristics (*this, jj, 0);
        } else {
                if (this->min () == -1 && this->max () == 1) {
                        /// assume design is conference matrix
                        return ::Jcharacteristics_conference (*this, jj, 0);
                } else {
                        myprintf ("Jcharacteristics not implemented for this type of array\n");
                        return std::vector< int > ();
                }
        }
}

int jvalue_conference (const array_link &ar, const int J, const int *column_indices) {
        int jval = 0;

        for (rowindex_t r = 0; r < ar.n_rows; r++) {
                int tmp = 1;
                for (int i = 0; i < J; i++) {
                        tmp *= ar.atfast (r, column_indices[i]);
                }
                jval += tmp;
        }
        return (jval);
}

/// calculate J-characteristics for a conference design
std::vector< int > Jcharacteristics_conference (const array_link &al, int jj, int verbose) {
        myassert (al.max () == 1 && al.min () == -1, "array is not conference design");

        const int k = al.n_columns;
        const int nc = ncombs (k, jj);
        std::vector< int > vals (nc);

        int *column_indices = new_perm_init< int > (jj);

        for (int x = 0; x < nc; x++) {
                int jv = jvalue_conference (al, jj, column_indices); 
                vals[x] = jv;
                next_comb_s (column_indices, jj, k);
        }

        delete_perm (column_indices);

        return vals;
}

int array_link::rank () const { return arrayrank (*this); }

int array_link::strength () const {
        int t = -1;
        for (int i = 0; i <= this->n_columns; i++) {
                int has_strength = strength_check (*this, i);
                if (has_strength) {
                        t = i;
                } else {
                        break;
                }
        }
        return t;
}

int array_cmp (carray_p A, carray_p B, const rowindex_t r, const colindex_t c) {
        int count = 0;
        for (int cpos = 0; cpos < c; cpos++) {
                for (int rpos = 0; rpos < r; rpos++) {
                        if (A[count] < B[count]) {
                                return -1;
                        }
                        if (A[count] > B[count]) {
                                return 1;
                        }
                        count++;
                }
        }
        return 0;
}

int array_diff (carray_p A, carray_p B, const rowindex_t r, const colindex_t c, rowindex_t &rpos, colindex_t &cpos) {
        int count = 0;
        for (cpos = 0; cpos < c; cpos++) {
                for (rpos = 0; rpos < r; rpos++) {
                        if (A[count] < B[count]) {
                                return -1;
                        }
                        if (A[count] > B[count]) {
                                return 1;
                        }
                        count++;
                }
        }
        return 0;
}

/// create new arraydata_t object
arraydata_t::arraydata_t (const array_t *s_, rowindex_t N_, colindex_t t, colindex_t nc)
    : N (N_), ncols (nc), strength (t), order (ORDER_LEX), colgroupindex (0), colgroupsize (0) {
        s = new array_t[nc];
        memcpy ((void *)s, (const void *)s_, sizeof (array_t) * nc);
        complete_arraydata ();
}
arraydata_t::arraydata_t (const std::vector< int > factor_levels, rowindex_t N_, colindex_t t, colindex_t nc)
    : N (N_), ncols (nc), strength (t), order (ORDER_LEX), colgroupindex (0), colgroupsize (0) {
        s = new array_t[ncols];
        myassert (factor_levels.size () > 0, "array class should have at least 1 factor");
        if ((int)factor_levels.size () < nc) {
                myprintf ("arraydata_t: warning: in constructor: size of factor levels %d < number of columns %d, padding with factor %d\n", (int)factor_levels.size(), nc, factor_levels[factor_levels.size () - 1]);
                nc = factor_levels.size ();
                std::fill (s, s + ncols, factor_levels[factor_levels.size () - 1]);
        }
        std::copy (factor_levels.begin (), factor_levels.end (), s);
        complete_arraydata ();
}

/// instantiate function
template void array_link::setarraydata (const short int *tmp, int n);
template void array_link::setarraydata (const int *tmp, int n);
template void array_link::setarraydata (const long *tmp, int n);
template void array_link::setarraydata (const std::vector< short int > tmp, int n);
template void array_link::setarraydata (const std::vector< long > tmp, int n);

arraydata_t::arraydata_t (array_t s_, rowindex_t N_, colindex_t t, colindex_t nc)
    : N (N_), ncols (nc), strength (t), order (ORDER_LEX), colgroupindex (0), colgroupsize (0) {
        if (s_ < 1 || s_ > 100) {
                myprintf ("arraydata_t: level factors should be > 0 and < 100\n");
        }
        s = new array_t[nc];
        for (int i = 0; i < nc; i++) {
                s[i] = s_;
        }
        complete_arraydata ();
}

arraydata_t::arraydata_t ()
    : N (0), ncols (0), strength (0), s (0), order (ORDER_LEX), colgroupindex (0), colgroupsize (0) {
}

arraydata_t::arraydata_t (const arraydata_t *adp, colindex_t newncols) {
        if (adp == 0) {
                myprintf ("arraydata_t: error: pointer to arraydata_t is invalid\n");
                return;
        }

        N = adp->N;
        strength = adp->strength;
        ncols = newncols;
        order = adp->order;
        colgroupindex = 0;
        colgroupsize = 0;

        if (ncols > adp->ncols) {
                myprintf ("arraydata_t: warning: number of columns %d is too large (object %d cols)\n", ncols,
                          adp->ncols);
                s = new array_t[ncols];
                for (size_t i = 0; i < (size_t)ncols; i++) {
                        s[i] = 2;
                }
                memcpy (s, adp->s, sizeof (array_t) * adp->ncols);
        } else {
                s = new array_t[ncols];
                memcpy (s, adp->s, sizeof (array_t) * ncols);
        }
        complete_arraydata ();
}

arraydata_t::arraydata_t (const arraydata_t &adp)
    : N (adp.N), ncols (adp.ncols), strength (adp.strength), order (adp.order), colgroupindex (0), colgroupsize (0) {
        s = new array_t[ncols];
        memcpy (s, adp.s, sizeof (array_t) * ncols);
        complete_arraydata ();
}
arraydata_t::~arraydata_t () {
        delete[] s;
        delete[] colgroupindex;
        delete[] colgroupsize;
}

arraydata_t& arraydata_t::operator= (const arraydata_t &ad2) {
	this->N = ad2.N;
	this->strength = ad2.strength;
	this->ncols = ad2.ncols;
	this->order = ad2.order;
	if (s != 0) {
		delete[] s;
	}
	this->s = new array_t[this->ncols];

	if (ad2.s == 0) {
		myprintf("error: invalid arraydata_t structure\n");
	}
	std::copy(ad2.s, ad2.s + this->ncols, s);
	return *this;
}

int arraydata_t::operator== (const arraydata_t &ad2) {
	if (this->N != ad2.N) {
		return 0;
	}

	if (this->ncols != ad2.ncols) {
		return 0;
	}
	if (!std::equal(this->s, this->s + this->ncols, ad2.s)) {
		return 0;
	}
	if (this->strength != ad2.strength) {
		return 0;
	}
	if (this->order != ad2.order) {
		return 0;
	}

	return 1;
}

void arraydata_t::writeConfigFile (const char *file) const {
        arraydata_t ad = *this;
        colindex_t ncols = ad.ncols;
        std::ofstream outFile;

        outFile.open (file);
        if (!outFile) {
                myprintf ("writeConfigFile: unable to open file %s\n", file);
                throw_runtime_exception ("writeConfigFile: unable to open file");
        }

        /* read design specifications: runs, strength, number of factors */
        std::string str;
        outFile << "runs " << ad.N << std::endl;
        outFile << "strength " << ad.strength << std::endl;
        outFile << "nfactors " << ad.ncols << std::endl;

        for (int j = 0; j < ncols; j++) {
                outFile << ad.s[j];
                if (j < ncols - 1) {
                        outFile << " ";
                } else {
                        outFile << std::endl;
                }
        }
        outFile.close ();
}

void arraydata_t::show (int verbose) const {
        myprintf ("%s\n", showstr ().c_str ());
        if (verbose >= 2) {
                for (int i = 0; i < this->ncolgroups; i++) {
                        myprintf (" colgroupindex[%d] %d\n", i, this->colgroupindex[i]);
                        myprintf (" colgroup size[%d] %d\n", i, this->colgroupsize[i]);
                }
        }
}

std::string arraydata_t::showstr () const {
        std::stringstream ss;
        ss << printfstring ("arrayclass: N %d, k %d, strength %d, s ", this->N, this->ncols, this->strength);
        print_perm (ss, this->s, this->ncols);
        std::string s = ss.str ();
        s = s.substr (0, s.size () - 1);
        s += printfstring (", order %d", this->order);
        return s;
}

std::string arraydata_t::latexstr (int cmd, int series) const {
        std::stringstream ss;
        if (cmd) {
                ss << printfstring ("\\oadesign{%d}{%d}{", this->N, this->strength);
        } else {
                ss << printfstring ("\\mathrm{OA}(%d, %d, ", this->N, this->strength);
        }

        for (int i = 0; i < this->ncolgroups; ++i) {
                int cgi = this->colgroupindex[i];
                int s = this->s[cgi];
                if (series > 0 && i == this->ncolgroups - 1 && this->colgroupsize[i] > 1) {
                        ss << printfstring ("%d^{a}", s);
                } else {
                        ss << printfstring ("%d^{%d}", s, this->colgroupsize[i]);
                }
        }
        if (cmd) {
                ss << printfstring ("}");
        } else {
                ss << printfstring (")");
        }
        return ss.str ();
}

/**
 * @brief Return identifier string
 */
std::string arraydata_t::idstrseriesfull () const {
        std::string fname = "";
        fname += itos (N);

        for (int i = 0; i < this->ncolgroups; ++i) {
                int cgi = this->colgroupindex[i];
                int s = this->s[cgi];
                if (i == 0) {
                        fname += '.';
                } else {
                        fname += '-';
                }
                if (i == this->ncolgroups - 1 && this->colgroupsize[i] > 1) {
                        fname += printfstring ("%da", s);
                } else {
                        fname += printfstring ("%d^%d", s, this->colgroupsize[i]);
                }
        }
        fname += "-t" + printfstring ("%d", this->strength);
        return fname;
}

/// return full identifier string
std::string arraydata_t::fullidstr (int series) const {
        if (series) {
                return this->idstrseriesfull (); 
        } else {
                return this->idstr () + "-t" + printfstring ("%d", this->strength);
        }
}

/**
 * @brief Return identifier string
 */
std::string arraydata_t::idstr () const {
        std::string fname = "";
        fname += itos (N);

        for (int i = 0; i < ncols; i++) {
                if (i == 0) {
                        fname += '.';
                } else {
                        fname += '-';
                }
                fname += itos (s[i]);
        }
        return fname;
}

/** Calculate derived data such as the index and column groups from a design
 */
void arraydata_t::complete_arraydata () {
        const int verbose = 0;

		if (!this->is_factor_levels_sorted() ) {			
			myprintf("arraydata_t: warning: the factor levels of the structure are not sorted, this can lead to undefined behaviour\n");
			this->show();
		}
		
		if (!check_divisibility(this)) {
			myprintf("arraydata_t: warning: no orthogonal arrays exist with the specified strength %d and specified factor levels\n", this->strength);
		}

        if (verbose) {
                myprintf ("complete_arraydata: strength %d\n", this->strength);
                for (int i = 0; i < this->ncols; i++) {
                        myprintf ("complete_arraydata: k %d, s %d\n", i, s[i]);
                }
        }
        if (this->strength > this->ncols) {
                myprintf ("arraydata_t: warning strength %d > ncols %d, reducing strength\n", this->strength,
                          this->ncols);
                this->strength = this->ncols;
        }
        if (this->strength < 1) {
                if (verbose >= 2) {
                        myprintf ("arraydata_t: warning strength < 1\n");
                }
        }
        arraydata_t *ad = this;
        this->calculate_oa_index (ad->strength);

        /* calculate column group structure */
        symmetry_group sg (ad->factor_levels(), 0);

        ad->ncolgroups = sg.ngroups;
        ad->colgroupindex = new colindex_t[ad->ncolgroups + 1];
        std::copy (sg.gstart.begin (), sg.gstart.end (), ad->colgroupindex);
        ad->colgroupsize = new colindex_t[ad->ncolgroups + 1];
        std::copy (sg.gsize.begin (), sg.gsize.end (), ad->colgroupsize);
}

void arraydata_t::lmc_overflow_check () const {
        const arraydata_t *ad = this;
#ifdef FULLPACKAGE
        int nbits = 8 * sizeof (rowsort_value_t);
        rowsort_value_t val = 1;
        for (int i = 0; i < ad->ncols; i++) {
                if (ad->s[i] == 0) {
                        continue;
                }
                if (val != 0 && (std::numeric_limits< rowsort_value_t >::max () / (rowsort_value_t)ad->s[i]) < val) {
                        // multiplication would exceed range of unsigned
                        printfd ("error: LMC checks for %d columns would lead to integer overflow\n", i);
                        std::cout << "  column: " << i << ": max rowsort value " << val
                                  << printfstring (" (ncols %d, nbits %d)", ad->ncols, nbits) << std::endl;
                        std::cout << printfstring ("      (ncols %d, nbits %d, ad->s[i] %d)", ad->ncols, nbits,
                                                   ad->s[i])
                                  << std::endl;
                }
                val *= ad->s[i];
        }
#endif
}

array_link arraydata_t::randomarray (int strength, int ncols) const {
        if (ncols == -1) {
                ncols = this->ncols;
        }
        array_link al (this->N, this->ncols, -1);
        al.setconstant (0);

        for (int i = 0; i < this->ncols; i++) {
                int coloffset = this->N * i;
                array_t s = this->getfactorlevel (i);

                int step = floor (double(N) / s);
                if (strength == 1) {
                        for (int j = 0; j < s; j++) {
                                std::fill (al.array + coloffset + step * j, al.array + coloffset + step * (j + 1), j);
                        }
                        random_perm (al.array + coloffset, this->N);

                } else {
                        for (int r = 0; r < N; r++) {
                                al.array[r + coloffset] = fastrand () % s;
                        }
                }
        }
        return al;
}

std::vector< int > arraydata_t::factor_levels () const {
          std::vector< int > s (this->ncols);
          for (int i = 0; i < this->ncols; i++) {
               s[i] = this->s[i];
          }
          return s;
}
        
array_link arraydata_t::create_root (int n_columns, int fill_value) const {
        if (n_columns == -1)
                n_columns = this->strength;
        array_link al (this->N, n_columns, -1);
        al.create_root (*this);

        for (int i = this->strength; i < al.n_columns; i++) {
                for (int r = 0; r < this->N; r++) {
                        al.at (r, i) = fill_value;
                }
        }
        return al;
}

/// return the factor level for the specified column return -1 if the column index is invalid
int arraydata_t::getfactorlevel(int idx) const {
	if (idx < 0) {
		return -1;
	}
	if (idx >= this->ncols) {
		return -1;
	}
	return this->s[idx];
}

void arraydata_t::reset_strength(colindex_t t) {
	strength = t;
	delete[] colgroupindex;
	delete[] colgroupsize;
	complete_arraydata();
}

colindex_t arraydata_t::get_col_group(const colindex_t col) const {
	colindex_t j = 0;
	for (colindex_t i = 0; i < ncolgroups; i++) {
		if (colgroupindex[i] <= col) {
			j = i;
		}
		else {
			break;
		}
	}
	return j;
}

/// return True if the vector is sorting in descending order
bool is_sorted_descending(const std::vector<int> values) {
	for (size_t i = 0; i < values.size() - 1; i++) {
		if (values[i] < values[i + 1])
			return false;
	}
	return true;
}
bool arraydata_t::is_factor_levels_sorted() const
{
	std::vector<int> factor_levels = this->factor_levels();

	return is_sorted_descending(factor_levels);
}

bool arraydata_t::is2level () const {
        for (int i = 0; i < this->ncols; i++) {
                if (this->s[i] != 2) {
                        return false;
                }
        }
        return true;
}

bool arraydata_t::ismixed () const {
        colindex_t s = this->s[0];
        for (int i = 0; i < this->ncols; i++) {
                if (s != this->s[i]) {
                        return true;
                }
        }
        return false;
}

void arraydata_t::set_colgroups (const std::vector< int > splits) {
        delete[] colgroupsize;
        delete[] colgroupindex;

        this->ncolgroups = splits.size ();
        this->colgroupindex = new colindex_t[this->ncolgroups + 1];
        std::copy (splits.begin (), splits.end (), this->colgroupindex);
        this->colgroupsize = new colindex_t[this->ncolgroups + 1];
        for (int j = 0; j < (ncolgroups - 1); j++) {
                this->colgroupsize[j] = this->colgroupindex[j + 1] - this->colgroupindex[j];
        }
        this->colgroupsize[ncolgroups - 1] = this->ncols - this->colgroupindex[ncolgroups - 1];
}

void arraydata_t::set_colgroups (const symmetry_group &sg) {
        delete[] colgroupsize;
        delete[] colgroupindex;

        mycheck (sg.n == this->ncols, "arraydata_t::set_colgroups: invalid size sg.n %d, ncols %d", sg.n, this->ncols);
        this->ncolgroups = sg.ngroups;
        this->colgroupindex = new colindex_t[this->ncolgroups + 1];
        std::copy (sg.gstart.begin (), sg.gstart.end (), this->colgroupindex);
        this->colgroupsize = new colindex_t[this->ncolgroups + 1];
        std::copy (sg.gsize.begin (), sg.gsize.end (), this->colgroupsize);
}

void arraydata_t::show_colgroups () const {
          myprintf ("arraydata_t: colgroups: ");
          print_perm (this->colgroupindex, this->ncolgroups);
          myprintf ("                  size: ");
          print_perm (this->colgroupsize, this->ncolgroups);
}

void arraydata_t::calculate_oa_index (colindex_t strength) {
          int combs = 1;
          for (int i = 0; i < this->strength; i++) {
               combs *= this->s[i];
          }

		  if ( ! check_divisibility(this) ) {
			  this->oaindex = 0;
			  return;
		  }

          if (combs == 0 || (this->N % combs) != 0) {
               this->oaindex = 0;
          } else {
               this->oaindex = this->N / combs;
          }
}
        
/**
 * @brief Complete arraydata but treat last column as single column group
 */
void arraydata_t::complete_arraydata_splitn (int ns) {
        delete[] colgroupsize;
        delete[] colgroupindex;

        for (int i = 0; i < (this->ncols - 1); i++) {
                if (this->s[i] != this->s[i + 1]) {
                        myprintf ("complete_arraydata_splitn: s[%d]=%d, s[%d]=%d, NOT IMPLEMENTED\n", i, this->s[i],
                                  i + 1, this->s[i + 1]);
                        myassert(0, "not implemented");
                }
        }

        arraydata_t *ad = this;
        int combs = 1;
        for (int i = 0; i < ad->strength; i++) {
                combs *= ad->s[i];
        }
        ad->oaindex = ad->N / combs;
        ad->ncolgroups = 2;

        colgroupindex = new colindex_t[ad->ncolgroups + 1];
        colgroupsize = new colindex_t[ad->ncolgroups + 1];
        colgroupindex[0] = 0;
        colgroupindex[1] = ns;
        colgroupsize[0] = ns;
        colgroupsize[1] = ad->ncols - ns;
}

/**
 * @brief Complete arraydata but treat last column as single column group
 */
void arraydata_t::complete_arraydata_fixlast () {
        delete[] colgroupsize;
        delete[] colgroupindex;

        complete_arraydata ();

        if (colgroupsize[ncolgroups - 1] == 1) {
                /* last column already a single column group */
                return;
        } else {
                colgroupsize[ncolgroups - 1]--;
                colgroupindex[ncolgroups] = ncols - 1;
                colgroupsize[ncolgroups] = 1;
                ncolgroups = ncolgroups + 1;
        }
}

int jstructbase_t::maxJ () const {
        int vmax = 0;
        for (size_t i = 0; i < this->values.size (); i++) {
                int v = abs (this->values[i]);
                if (v > vmax) {
                        vmax = v;
                }
        }
        return vmax;
}

std::vector< int > jstructbase_t::calculateF () const {
        int nn = this->jvalues.size (); // floor ( ( double ) N/x ) +1;
        std::vector< int > F (nn);

        for (size_t i = 0; i < this->values.size (); i++) {
                int fi = abs (values[i]);
                int idx = jvalue2index.find (fi)->second;
                F[idx]++;
        }
        return F;
}

/* analyse arrays */

jstruct_t::jstruct_t () {
        this->nc = 0;
        this->abberration = -1;
}

int jstruct_t::maxJ () const {
        int vmax = 0;
        for (size_t i = 0; i < this->values.size (); i++) {
                int v = abs (this->values[i]);
                if (v > vmax) {
                        vmax = v;
                }
        }
        return vmax;
}

std::vector< int > jstruct_t::Fval (int strength) const {
        int x = pow ((double)2, strength + 1); 
        int nn = floor ((double)N / x) + 1;
        std::vector< int > Fv (nn);
        for (int i = 0; i < nn; i++) {
                Fv[i] = N - x * i;
        }
        return Fv;
}

std::vector< int > jstruct_t::calculateF (int strength) const {
        int Nmax = N;

        int x = pow (double(2), strength + 1); 
        int nn = floor ((double)N / x) + 1;
        std::vector< int > F (nn);

        for (int i = 0; i < nc; i++) {
                int fi = (N - abs (values[i])) / x;
                F[fi]++;
        }
        return F;
}

void jstruct_t::calculateAberration() {
	jstruct_t *js = this;
	js->abberration = 0;
	for (int i = 0; i < js->nc; i++) {
		js->abberration += js->values[i] * js->values[i];
	}
	js->abberration /= N * N;
}

void jstruct_t::calc (const array_link &al) {
        int *column_indices = new_perm_init< int > (jj);

        myassert (al.is2level (), "array is not 2-level array");
        for (int x = 0; x < this->nc; x++) {
                int jv = jvaluefast (al.array, al.n_rows, jj, column_indices);
                this->values[x] = jv;
                next_comb_s (column_indices, jj, k);
        }

        delete_perm (column_indices);
}

void jstructconference_t::calcJvalues(int N, int jj) {
	if (jj % 2 != 0) {
		throw_runtime_exception("calculation of J-characteristics for conference matrices only supported for even J");
	}
	if (jj != 4) {
		throw_runtime_exception("calculation of J-characteristics for conference matrices only supported for J=4");
	}
	int nn = floor(double(int((N - jj + 1) / 4))) + 1;
	this->jvalues = std::vector< int >(nn);
	this->jvalue2index.clear();
	for (size_t i = 0; i < jvalues.size(); i++) {
		int jval = (N - jj) - i * 4;
		jvalues[i] = jval;
		jvalue2index[jval] = i;
	}
}
void jstructconference_t::calc(const array_link &al) { values = Jcharacteristics_conference(al, this->jj); }

/// create J2 table as intermediate result for J-characteristic calculations
array_link createJdtable (const array_link &al) {
        const int nr = al.n_rows;
        const int nc = al.n_columns * al.n_columns;

        // fill double column table
        array_link dtable (nr, nc, -1);

        // loop over all column pairs
        int idx = 0;
        for (int i = 0; i < al.n_columns; i++) {
                for (int j = 0; j < al.n_columns; j++) {
                        // loop over all rows of original array
                        const array_t *p1 = al.array + al.n_rows * i;
                        const array_t *p2 = al.array + al.n_rows * j;
                        array_t *pout = dtable.array + idx * dtable.n_rows;
                        for (int x = 0; x < nr; x++) {
                                pout[x] = p1[x] ^ p2[x];
                        }
                        idx++;
                }
        }

        return dtable;
}

void jstruct_t::calcj4 (const array_link &al) {
        myassert (jj == 4, "jj should be 4");
        int *pp = new_perm_init< int > (jj);
        int ncolcombs = ncombs (k, jj);
        const int nr = al.n_rows;
        const int N = nr;

        array_link dtable = createJdtable (al);

        for (int x = 0; x < this->nc; x++) {
                const int idx1 = pp[0] + pp[1] * al.n_columns;
                const int idx2 = pp[2] + pp[3] * al.n_columns;
                int jv = 0;
                {
                        const array_t *o1 = dtable.array + dtable.n_rows * idx1;
                        const array_t *o2 = dtable.array + dtable.n_rows * idx2;
                        for (int xr = 0; xr < nr; xr++) {
                                int tmp = (o1[xr]) ^ (o2[xr]);
                                jv += tmp;
                        }
                }
                jv = 2 * jv - N;
                this->values[x] = jv;

                next_comb_s (pp, jj, k);
        }

        delete_perm (pp);
}

void jstruct_t::calcj5 (const array_link &al) {
        myassert (jj == 5, "jj should be 5");
        int *pp = new_perm_init< int > (jj);
        int ncolcombs = ncombs (k, jj);
        const int nr = al.n_rows;
        const int N = nr;

        array_link dtable = createJdtable (al);


        for (int x = 0; x < this->nc; x++) {
                const int idx1 = pp[0] + pp[1] * al.n_columns;
                const int idx2 = pp[2] + pp[3] * al.n_columns;
                const int idx3 = pp[4];
                int jv = 0;
                {
                        const array_t *o1 = dtable.array + dtable.n_rows * idx1;
                        const array_t *o2 = dtable.array + dtable.n_rows * idx2;
                        const array_t *o3 = al.array + N * idx3;
                        for (int xr = 0; xr < nr; xr++) {
                                int tmp = (o1[xr]) ^ (o2[xr]) ^ o3[xr];
                                jv += tmp;
                        }
                }
                jv = 2 * jv - N;
                this->values[x] = jv;
                next_comb_s (pp, jj, k);
        }

        delete_perm (pp);
}

/** Create an object to calculate J-characteristics
 *
 * @param al Array to use
 * @param jj Type of J-characteristics to calculate
 *
 **/
jstruct_t::jstruct_t (const array_link &al, int jj) {
        const int k = al.n_columns;
        const int N = al.n_rows;

        this->init (N, k, jj);
        if (jj == 4 && 1) {
                this->calcj4 (al);
        } else {
                if (jj == 5 && 1) {
                        this->calcj5 (al);
                } else {
                        this->calc (al);
                }
        }
        // calculate A value
        this->calculateAberration ();
}

jstruct_t::jstruct_t (const int N_, const int k_, const int jj_) { init (N_, k_, jj_); }

void jstruct_t::init (int N_, int k_, int jj_) {
        this->N = N_;
        this->k = k_;
        this->jj = jj_;

        if (this->jj < 0) {
                myprintf ("jstruct_j: J-characteristics for J<0 are undefined, setting J to 0\n");
                jj = 0;
        }
        if (this->jj > 20) {
                myprintf ("jstruct_j: J-characteristics for J>20 are not supported, setting J to 0\n");
                jj = 0;
        }

        this->nc = ncombs< long > (k_, jj_);
        values = std::vector< int > (nc);
        this->abberration = -1;
}

jstruct_t::jstruct_t (const jstruct_t &js) {
        N = js.N;
        k = js.k;
        jj = js.jj;
        nc = js.nc;
        abberration = js.abberration;
        values = std::vector< int > (nc);
        std::copy (js.values.begin (), js.values.begin () + nc, values.begin ());
}

jstruct_t &jstruct_t::operator= (const jstruct_t &rhs) {
        this->N = rhs.N;
        this->k = rhs.k;
        this->jj = rhs.jj;
        this->nc = rhs.nc;

        this->abberration = rhs.abberration;
        values = std::vector< int > (nc);
        std::copy (rhs.values.begin (), rhs.values.begin () + nc, values.begin ());

        return *this;
}

jstruct_t::~jstruct_t () {
}

std::string jstruct_t::showstr () {
        std::stringstream ss;
        ss << "jstruct_t: " << printfstring ("N %d, jj %d", N, jj);
        return ss.str ();
}

int jstruct_t::allzero() const {
	for (int i = 0; i < this->nc; ++i) {
		if (this->values[i] != 0) {
			return 0;
		}
	}
	return 1;
}

void jstruct_t::show () const {
#ifdef FULLPACKAGE
        std::cout << "jstruct_t: " << printfstring ("N %d, jj %d, values ", N, jj);
        for (int x = 0; x < this->nc; x++) {
                cout << printfstring (" %d", values[x]);
        }
        std::cout << std::endl;
#endif
}

void jstruct_t::showdata () {
        for (size_t x = 0; x < this->values.size (); x++) {
                myprintf("%s", printfstring (" %d", values[x]).c_str());
        }
		myprintf("\n");
}

std::string jstructbase_t::showstr () {
        std::string s = "jstruct_t: " + printfstring ("jj %d, values ", jj);
        return s;
}
void jstructbase_t::show () {
#ifdef FULLPACKAGE
        std::cout << "jstruct_t: " << printfstring ("jj %d, values ", jj);
        for (size_t x = 0; x < this->values.size (); x++) {
                std::cout << printfstring (" %d", values[x]);
        }
        std::cout << std::endl;
#endif
}

void jstructbase_t::showdata (int verbose) {
#ifdef FULLPACKAGE

        for (size_t x = 0; x < this->values.size (); x++) {
                std::cout << printfstring (" %d", values[x]);
        }
        std::cout << std::endl;

        if (verbose >= 2) {
                typedef std::map< int, int >::const_iterator it_type;
                for (it_type iterator = this->jvalue2index.begin (); iterator != this->jvalue2index.end ();
                     iterator++) {
                        myprintf ("this->jvalue2index[%d]=%d\n", iterator->first, iterator->second);
                }
        }
#endif
}

/** Calculate J-value for an array
 *
 * The array should be a 2-factor array
 *
 * @param array_link Array
 * @param J Number of columns
 * @param pp Indices of columns to use
 * @return j-value
 */
int jvalue (const array_link &ar, const int J, const int *pp) {
        int jval = 0;

        for (rowindex_t r = 0; r < ar.n_rows; r++) {
                int tmp = 0;
                for (int i = 0; i < J; i++) {
                        tmp += ar.at (r, pp[i]);
                }
                tmp %= 2;
                tmp *= 2;
                tmp--;
                jval -= tmp;
        }
        return (jval);
}

symmdata::symmdata (const array_link &al, int minlen) {
        orig = al;
        rowvalue = al;
        const int N = al.n_rows;

        for (int c = 0; c < al.n_columns; c++) {
                array_t *rvc = rowvalue.array + c * N;
                rvc[0] = 0;
                for (int r = 1; r < al.n_rows; r++) {
                        if (al.atfast (r, c) != al.atfast (r - 1, c)) {
                                rvc[r] = rvc[r - 1] + 1;
                                ;
                        } else {
                                if (c > 0) {
                                        if (rowvalue.atfast (r, c - 1) != rowvalue.atfast (r - 1, c - 1)) {
                                                rvc[r] = rvc[r - 1] + 1;
                                        } else {
                                                rvc[r] = rvc[r - 1];
                                        }
                                } else {
                                        rvc[r] = rvc[r - 1];
                                }
                        }
                }
        }

        ft = array_link (2 * al.n_rows + 2, al.n_columns, -1);
        ft.setconstant (0);
        size_t nfrow = ft.n_rows - 1;
        for (int c = 0; c < al.n_columns; c++) {
                array_t v = rowvalue.at (0, c);
                int nf = 0;
                int prevr = 0;
                carray_t *rvc = rowvalue.array + c * N;

                for (int r = 1; r < al.n_rows; r++) {
                        if (rvc[r] != v) {
                                if ((r - prevr) >= minlen) {
                                        ft.atfast (2 * nf, c) = prevr;
                                        ft.atfast (2 * nf + 1, c) = r;
                                        nf++;
                                }
                                prevr = r;
                                v = rvc[r];
                        }
                }
                ft.atfast (2 * nf, c) = prevr;
                ft.atfast (2 * nf + 1, c) = al.n_rows;
                nf++;
                ft.atfast (nfrow, c) = nf;
        }
}

/** @brief Calculate J-characteristic of a 2-level array for a column combination
*
* We assume the array has values 0 and 1
*/
int fastjX (const array_t *array, rowindex_t N, const int J, const colindex_t *pp) {
        int jval = 0;

        const array_t *cp[MAXCOLS];

        for (int i = 0; i < J; i++) {
                cp[i] = array + N * pp[i];
        }

        for (rowindex_t r = 0; r < N; r++) {
                array_t tmp = 0;
                for (int i = 0; i < J; i++) {
                        tmp += cp[i][r]; 
                }
                tmp %= 2;
                jval += tmp;
        }
        jval = 2 * jval - N;
        return (jval);
}


int jvaluefast (const array_t *array, rowindex_t N, const int J, const colindex_t *column_indices) {
        array_t tmpval[MAXROWS];

        std::fill_n (tmpval, N, 0);
        fastJupdate (array, N, J, column_indices, tmpval);
        int jval = 0;
        for (rowindex_t r = 0; r < N; r++) {
                jval += tmpval[r] % 2;
        }
        jval = 2 * jval - N;
        return (jval);
}

/** Analyse a list of arrays
 *
 * For each array the j-values are calculated
 *
 */
vector< jstruct_t > analyseArrays (const arraylist_t &arraylist, const int verbose, const int jj) {
        if (verbose) {
                myprintf ("analyseArrays (j-values): %ld arrays, jj %d\n", (long)arraylist.size (), jj);
        }

        vector< jstruct_t > results;
        results.reserve (arraylist.size ());
        jstruct_t *js;

        for (unsigned int ii = 0; ii < arraylist.size (); ii++) {

                const array_link ll = arraylist.at (ii);

                int *pp = new_perm_init< int > (jj);

                js = new jstruct_t (ll, jj);

                if (verbose >= 2) {
                        std::vector< int > FF = js->calculateF ();
                        myprintf ("F%d (high to low): ", jj);
                        display_vector (FF);
						myprintf("\n");
                }
                delete_perm (pp);

                results.push_back (*js);
                delete js;
        }

        return results;
}

/* reading and writing of arrays */

void write_array_latex (FILE *fid, carray_t *array, const int nrows, const int ncols);

/**
 * @brief Write an array to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void write_array (FILE *fid, carray_t *array, const int nrows, const int ncols) {
        int count;
        stringstream ss;

        for (int j = 0; j < nrows; j++) {
                count = j;
                for (int k = 0; k < ncols; k++) {
                        const char *s = (k < ncols - 1) ? " " : "\n";
                        ss << array[count] << s;
                        count += nrows;
                }
        }
        fprintf (fid, "%s", ss.str ().c_str ());
}

/**
 * @brief Write array file in LaTeX format
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void write_array_latex (FILE *fid, carray_t *array, const int nrows, const int ncols) {
        int count;
        stringstream ss;

        ss << "\\begin{tabular}{";
        for (int x = 0; x < ncols; x++) {
                ss << 'c';
        }
        ss << "}" << endl;

        for (int j = 0; j < nrows; j++) {
                count = j;
                for (int k = 0; k < ncols; k++) {
                        const char *s = (k < ncols - 1) ? " & " : " \\\\ \n";
                        ss << array[count] << s;
                        count += nrows;
                }
        }
        ss << "\\end{tabular}" << endl;

        fprintf (fid, "%s", ss.str ().c_str ());
}

int append_arrays (FILE *fid, arraylist_t &arrays, int startidx = 0) {
        arraylist_t::iterator it;

        for (it = arrays.begin (); it != arrays.end (); it++) {
                fprintf (fid, "%i\n", startidx++);
                write_array (fid, it->array, it->n_rows, it->n_columns);
        }

        return 0;
}

template < class TypeIn, class TypeOut >
/// Write array to binary blob of selected datatype
void writeblob (const TypeIn *src, int n, FILE *fid) {
        TypeOut *dst = new TypeOut[n];

        for (int i = 0; i < n; i++) {
                dst[i] = src[i];
        }
        fwrite ((const void *)dst, sizeof (TypeOut), n, fid);

        delete[] dst;
}

template < class TypeIn, class TypeOut >
/// Convert binary blob to datatype
void readblob (TypeOut *dst, int n, FILE *fid) {
        TypeIn *src = new TypeIn[n];
        int r = fread ((void *)src, sizeof (TypeIn), n, fid);

        for (int i = 0; i < n; i++) {
                dst[i] = src[i];
        }

        delete[] src;
}
#ifdef USEZLIB
template < class TypeIn, class TypeOut >
/// Convert binary blob to datatype
void readblob (TypeOut *dst, int n, gzFile gzfid) {
        TypeIn *src = new TypeIn[n];
        gzread (gzfid, src, sizeof (TypeIn) * n);

        for (int i = 0; i < n; i++) {
                dst[i] = src[i];
        }

        delete[] src;
}
#endif

/**
 * @brief Read an array from a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void read_array (FILE *fid, array_t *array, const int nrows, const int ncols) {
        const int maxbuf = 1024;
        char buf[maxbuf];
        int count;
        for (int j = 0; j < nrows; j++) {
                count = j;

                char *r = fgets (buf, maxbuf, fid);
                stringstream ss (stringstream::in | stringstream::out);
                ss << buf;

                for (int k = 0; k < ncols; k++) {
                        ss >> array[count];
                        count += nrows;
                }
        }
}

bool file_exists (const std::string filename) { return file_exists (filename.c_str ()); }

bool oa_file_exists (const char *filename) {
        std::string s = filename;
        return oa_file_exists (s);
}

bool oa_file_exists (const std::string filename) {
        if (file_exists (filename.c_str ())) {
                return true;
        }
        std::string gzfile = filename + ".gz";
        if (file_exists (gzfile.c_str ())) {
                return true;
        }
        return false;
}

#ifdef HAVE_BOOST
#include "boost/filesystem.hpp"
#endif

/// return true if the specified file exists
bool file_exists (const char *filename) {
#ifdef HAVE_BOOST
        return boost::filesystem::exists (filename);
#else
        FILE *fid = fopen (filename, "r");
        if (fid == 0) {
                return false;
        } else {
                fclose (fid);
                return true;
        }
#endif
}

#ifdef FULLPACKAGE

bool arrayfile_t::isbinary () const {
        return (this->mode == ABINARY || this->mode == ABINARY_DIFF || this->mode == ABINARY_DIFFZERO);
}

int arrayfile_t::swigcheck () const {
#ifdef SWIGCODE
          if (sizeof (int) != 4) {
               fprintf (stderr, "arrayfile_t: error int is not 32-bit?!\n");
          }
          return 1;
#else
          return 0;
#endif
}

/// return string describing the object
std::string arrayfile_t::showstr () const {
          if (this->isopen ()) {
               std::string modestr;
               switch (mode) {
               case ALATEX:
                         modestr = "latex";
                         break;
               case ATEXT:
                         modestr = "text";
                         break;
               case ABINARY:
                         modestr = "binary";
                         break;
               case ABINARY_DIFF:
                         modestr = "binary_diff";
                         break;
               case ABINARY_DIFFZERO:
                         modestr = "binary_diffzero";
                         break;
               case AERROR:
                         modestr = "invalid";
                         break;
               default:
                         modestr = "error";
                         myprintf ("arrayfile_t: showstr(): no such mode\n");
                         break;
               }

               int na = narrays;
               if (this->rwmode == WRITE) {
                         na = narraycounter;
               }

               std::string s = printfstring ("file %s: %d rows, %d columns, %d arrays", filename.c_str (),
                                             nrows, ncols, na);
               s += printfstring (", mode %s, nbits %d", modestr.c_str (), nbits);
               return s;
          } else {
               std::string s = "file " + filename + ": invalid file";
               return s;
          }
}
        
void arrayfile_t::append_array (const array_link &a, int specialindex) {
        int r;
        if (!this->isopen ()) {
                return;
        }

        arrayfile_t *afile = this;
        int index;
        if (specialindex == -1) {
                index = a.index;
        } else {
                index = specialindex;
        }

        switch (afile->mode) {
        case arrayfile::ATEXT:
                fprintf (afile->nfid, "%i\n", index);
                write_array (afile->nfid, a.array, a.n_rows, a.n_columns);
                break;
        case arrayfile::ALATEX:
                fprintf (afile->nfid, "%i\n", index);
                write_array_latex (afile->nfid, a.array, a.n_rows, a.n_columns);
                break;
        case arrayfile::ABINARY: {
                r = fwrite (&index, 1, sizeof (int32_t), afile->nfid);
                if (r != sizeof (int32_t) * 1) {
                        printfd ("error during write to file\n");
                }
                afile->write_array_binary (a.array, a.n_rows, a.n_columns);
                break;
        }
        case arrayfile::ABINARY_DIFF: {
                fwrite (&index, 1, sizeof (int32_t), afile->nfid);
                afile->write_array_binary_diff (a);
                break;
        }
        case arrayfile::ABINARY_DIFFZERO: {
                afile->write_array_binary_diffzero (a);
                break;
        }
        default:
                std::cout << "warning: arrayfile_t::append_array: no such mode " << afile->mode << std::endl;
                break;
        }
        narraycounter++;
}

/// return true if file is open
int arrayfile_t::isopen () const {
#ifdef USEZLIB
        return (this->nfid != 0 || this->gzfid != 0);
#else
        return (this->nfid != 0);
#endif
}

int arrayfile_t::append_arrays (const arraylist_t &arrays, int startidx) {
         myassert (this->nfid!=0, "file pointer is zero");

        for (arraylist_t::const_iterator it = arrays.begin (); it != arrays.end (); it++) {
                this->append_array (*it, startidx);
                if (startidx >= 0)
                        startidx++;
        }
        return 0;
}

arrayfile_t::arrayfile_t (const std::string fname, int nrows, int ncols, int narrays_, arrayfilemode_t m, int nb) {
        this->verbose = 0;
        this->nfid = 0;
        this->narrays = -1;
        this->narraycounter = -1;
        this->rwmode = READ;
        this->mode = ATEXT;

#ifdef USEZLIB
        this->gzfid = 0;
#endif

        createfile (fname, nrows, ncols, narrays_, m, nb);
}

int arrayfile_t::headersize () const { return 8 * sizeof (int32_t); }

int arrayfile_t::barraysize () const {
          int num = sizeof (int32_t);

          switch (this->nbits) {
          case 8:
               num += nrows * ncols;
               break;
          case 32:
               num += nrows * ncols * 4;
               break;
          case 1: {
               word_addr_t num_of_words = nwords (nrows * ncols);
               num += sizeof (word_t) * num_of_words;
          } break;
          default:
               myprintf ("error: number of bits undefined\n");
               break;
          }
          return num;
}

size_t arrayfile_t::afwrite (void *ptr, size_t t, size_t n) {
          if (this->nfid == 0) {
               myprintf ("afwrite: not implemented, we cannot write compressed files\n");
               return 0;
          }
          return fwrite (ptr, t, n, nfid);
}
        
size_t arrayfile_t::afread (void *ptr, size_t sz, size_t cnt) {
#ifdef USEZLIB
        size_t r;
        if (iscompressed) {
                assert (this->gzfid);
                r = gzread (this->gzfid, ptr, sz * cnt);
                r /= sz;
        } else {
                r = fread (ptr, sz, cnt, this->nfid);
        }
        if (r == 0) {
                // myprintf("could not read from file");
        }
#else
        size_t r;
        r = fread (ptr, sz, cnt, this->nfid);
#endif
        return r;
}

int arrayfile_t::seek (int pos) {
        if (mode != ABINARY) {
                myprintf ("error: seek only possible for binary files\n");
                return -1;
        }

        if (pos < 0 || pos > 1000000000) {
                myprintf ("error: position specified is invalid\n");
                return -1;
        }
        if (nfid) {
                fseek (nfid, this->headersize () + this->barraysize () * pos, SEEK_SET);
        }
#ifdef USEZLIB
        if (gzfid) {
                gzseek (gzfid, this->headersize () + this->barraysize () * pos, SEEK_SET);
        }
#endif

        narraycounter = pos;
        return pos;
}

int arrayfile_t::read_array_binary_zero (array_link &a) {

        // no index is written or read (to save disk space
        a.index = array_link::INDEX_NONE;
        int index = array_link::INDEX_NONE;

        if (a.n_columns != diffarray.n_columns && diffarray.n_columns != -1 && diffarray.n_columns != 0) {
                myprintf ("arrayfile_t::read_array_binary_zero: error different number of columns %d %d\n",
                          diffarray.n_columns, a.n_columns);
                return array_link::INDEX_ERROR;
        }

        int16_t nrest;
        int result = afread (&nrest, sizeof (int16_t), 1);
        if (result != 1) {
                // error reading array, we could have reached the end of the file
                if (this->narrays == -1) {
                } else {
                        myprintf ("arrayfile_t::read_array_binary_zero: error reading array: index %d, result %d\n",
                                  index, result);
                }
                index = array_link::INDEX_ERROR;
                return index;
        }

        int nrows = a.n_rows;
        int ncols = a.n_columns;
        int ngood = ncols - nrest;
        if (diffarray.n_columns > 0) {
                copy_array (diffarray.array, a.array, nrows, ngood);
        } else {
                // diffarray not initialized yet...
        }
        this->read_array_binary (a.array + nrows * ngood, nrows, nrest);

        if (a.n_columns == diffarray.n_columns) {
                // update array

                int N = a.n_rows;
                for (int i = 0; i < nrest; i++) {
                        for (int j = 0; j < N; j++) {
                                int idx = j + (i + ngood) * N;
                                a.array[idx] += diffarray.array[idx];
                                a.array[idx] = a.array[idx] % 2;
                        }
                }
        }
        diffarray = a;
        return index;
}

arraylist_t arrayfile_t::readarrays (int nmax, int verbose) {
        arraylist_t ll;

        if (!this->isopen ()) {
                if (verbose <= 1) {
                        myprintf ("arrayfile_t::readarrays: array file is not open (file %s)\n",
                                  this->filename.c_str ());
                }

                return ll;
        }

        for (int i = 0; i < nmax; i++) {

                array_link x = this->readnext ();
                if (x.index == array_link::INDEX_ERROR) {
                        break;
                }
                ll.push_back (x);
        }
        return ll;
}

array_link arrayfile_t::readnext () {
        array_link al (this->nrows, this->ncols, array_link::INDEX_DEFAULT);
        this->read_array (al);
        return al;
}

void arrayfile_t::flush () {
        if (this->nfid != 0) {
                fflush (nfid);
        }
#ifdef USEZLIB
        if (this->gzfid != 0) {
                // gzflush(gzfid, Z_SYNC_FLUSH);
        }
#endif
}
int arrayfile_t::read_array (array_link &a) {
        int32_t index;
        switch (this->mode) {
        case arrayfile::ATEXT:
                index = this->read_array (a.array, a.n_rows, a.n_columns);
                a.index = index;
                break;
        case arrayfile::ABINARY:
                index = this->read_array (a.array, a.n_rows, a.n_columns);
                a.index = index;
                break;
        case arrayfile::ABINARY_DIFFZERO:
                index = this->read_array_binary_zero (a);
                a.index = index;
                break;
        case arrayfile::ABINARY_DIFF: {
                int result = afread (&index, sizeof (int32_t), 1);
                if (result != 1) {
                        // error reading array, we could have reached the end of the file
                        if (this->narrays == -1) {
                        } else {
                                myprintf ("error reading array: index %d, result %d\n", index, result);
                        }
                        index = array_link::INDEX_ERROR;
                        a.index = index;
                        break;
                }

                int32_t nrest;
                result = afread (&nrest, 1, sizeof (int32_t));

                int nrows = a.n_rows;
                int ncols = a.n_columns;
                int ngood = ncols - nrest;
                copy_array (diffarray.array, a.array, nrows, ngood);

                this->read_array_binary (a.array + nrows * ngood, nrows, nrest);
                a.index = index;
                diffarray = a;
        } break;
        default:
                index = -1;
                break;
        }

        return index;
}

int arrayfile_t::read_array (array_t *array, const int nrows, const int ncols) {
        int index = -10;
        if (nrows != this->nrows) {
                myprintf ("arrayfile_t::read_array: nrows unequal %d %d\n", nrows, this->nrows);
        }
        if (ncols != this->ncols) {
                myprintf ("arrayfile_t::read_array: ncols unequal\n");
        }

        switch (this->mode) {
        case arrayfile::ATEXT: {
                int r = fscanf (nfid, "%d\n", &index);
                ::read_array (nfid, array, nrows, ncols);
                break;
        }
        case arrayfile::ABINARY: {
                int result = afread (&index, sizeof (int32_t), 1);
                if (result != 1) {
                        index = -1;
                        myprintf ("arrayfile_t::read_array: error: could not read array index: result %d, index %d\n",
                                  result, index);
                        break;
                }
                if (index == -1) {
                        myprintf ("updating index (-1 is invalid)\n");
                        index = 0;
                }

                this->read_array_binary (array, nrows, ncols);
        } break;
        default:
                myprintf ("arrayfile_t::read_array: error: no such mode %d\n", this->mode);
                break;
        }

        return index;
}

void arrayfile_t::finisharrayfile() {
	if (this->mode == ATEXT) {
		fprintf(this->nfid, "-1\n");
	}
	this->closefile();
}

void arrayfile_t::setVerbose(int v) { this->verbose = v; }

void arrayfile_t::writeheader () {
        assert (this->nfid);

        if (this->mode == arrayfile::ABINARY || this->mode == arrayfile::ABINARY_DIFF ||
            this->mode == arrayfile::ABINARY_DIFFZERO) {
                int magic = 65;
                int32_t reserved = 0;

                if (this->mode == arrayfile::ABINARY_DIFFZERO) {
                        // NOTE: needed because in ABINARY_DIFFZERO we diff modulo 2
                        if (this->nbits != 1)
                                myprintf ("not implemented...!\n");
                        assert (this->nbits == 1);
                }
                fwrite ((const void *)&magic, sizeof (int), 1, this->nfid);
                fwrite ((const void *)&this->nbits, sizeof (int), 1, this->nfid);
                fwrite ((const void *)&this->nrows, sizeof (int), 1, this->nfid);
                fwrite ((const void *)&this->ncols, sizeof (int), 1, this->nfid);
                fwrite ((const void *)&this->narrays, sizeof (int), 1, this->nfid);

                if (this->mode == arrayfile::ABINARY) {
                        reserved = 1001;
                } else {
                        if (this->mode == arrayfile::ABINARY_DIFF) {
                                reserved = 1002;
                        } else {
                                reserved = 1003;
                        }
                }
                fwrite ((const void *)&reserved, sizeof (int), 1, this->nfid);
                reserved = 9999;
                fwrite ((const void *)&reserved, sizeof (int), 1, this->nfid);
                fwrite ((const void *)&reserved, sizeof (int), 1, this->nfid);
        } else {
                fprintf (this->nfid, "%i %i %i\n", this->ncols, this->nrows, this->narrays);
        }
}

//! Create new array link object, clone an array
array_link::array_link (const array_t *array, rowindex_t nrows, colindex_t ncolsorig, colindex_t ncols,
                        int index_ = -1)
    : n_rows (nrows), n_columns (ncols), index (index_) {
        this->array = create_array (nrows, ncols);
        memcpy (this->array, array, nrows * ncolsorig * sizeof (array_t)); // FIX: replace by copy_array
}

std::string printfstring(const char *message, ...) {
	char buf[32 * 1024];

	va_list va;
	va_start(va, message);
	vsprintf(buf, message, va);
	va_end(va);

	std::string str(buf);
	return str;
}

/**
 * @brief Read file with design of OA
 * @param file
 * @return
 */
arraydata_t *readConfigFile (const char *file) {
        //***open config file***
        colindex_t N, strength, ncols;
        array_t *s;
        ifstream inFile;

        inFile.open (file);
        if (!inFile) {
                myprintf ("readConfigFile: unable to open file %s\n", file);
                return 0;
        }

        /* read design specifications: runs, strength, number of factors */
        string str;
        inFile >> str >> N;
        assert (str.compare ("runs "));
        inFile >> str >> strength;
        assert (str.compare ("strength "));
        inFile >> str >> ncols;
        assert (strcmp (str.c_str (), "nfactors "));
        if (N > 10000 || N < 1) {
                myprintf ("readConfigFile: file %s: invalid number of runs %d\n", file, N);
                return 0;
        }
        if (strength > 1000 || strength < 1) {
                myprintf ("readConfigFile: file %s: invalid strength %d\n", file, strength);
                return 0;
        }
        if (ncols > 1000 || ncols < 1) {
                myprintf ("readConfigFile: file %s: invalid ncols %d\n", file, ncols);
                return 0;
        }
        s = (array_t *)malloc (ncols * sizeof (array_t));
        for (int j = 0; j < ncols; j++) {
                inFile >> s[j];
                if ((s[j] < 1) || (s[j] > 25)) {
                        myprintf ("warning: number of levels specified is %d\n", s[j]);
						throw_runtime_exception("incorrect specification of factor levels");
                }
        }
        inFile.close ();

        arraydata_t *ad = new arraydata_t (s, N, strength, ncols);
        ad->lmc_overflow_check ();
        free (s);
        return ad;
}

long nArrays (const char *fname) {
        arrayfile_t af (fname, 0);
        long n = af.narrays;
        af.closefile ();
        return n;
}

void arrayfileinfo(const char *filename, int &number_of_arrays, int &number_of_rows, int &number_of_columns) {
	arrayfile_t af(filename, 0);
	number_of_arrays = af.narrays;
	number_of_rows = af.nrows;
	number_of_columns = af.ncols;
	af.closefile();
}

void appendArrays(const arraylist_t &arrays_to_append, arraylist_t &dst) {
	for (arraylist_t::const_iterator it = arrays_to_append.begin(); it < arrays_to_append.end(); ++it) {
		dst.push_back(*it);
	}
}

void convert_array_file(std::string input_filename, std::string output_filename, arrayfile::arrayfilemode_t output_format, int verbose)
{
	arrayfile_t input_array_file(input_filename.c_str());
	int nc = input_array_file.ncols;
	int nr = input_array_file.nrows;
	array_link al(nr, nc, -1);

	if (!input_array_file.isopen()) {
		fprintf(stderr, "oaconvert: could not open input file, aborting...\n");
		throw_runtime_exception(printfstring("convert_array_file: could not open input file %s", input_filename.c_str()));
	}
	if (verbose)
		printf("oaconvert: output mode %d, nr %d nc %d\n", output_format, nr, nc);

	{
		// streaming mode
		int narrays = input_array_file.narrays;
		int nb = 8; // we have not read any arrays so far, so nb is hard to predict
		if (output_format == ABINARY_DIFFZERO)
			nb = 1;

		arrayfile_t afout(output_filename.c_str(), nr, nc, narrays, output_format, nb);

		if (input_array_file.narrays < 0) {
			narrays = arrayfile_t::NARRAYS_MAX;
			printf("warning: untested code! (number of arrays undefined)\n");
		}

		int index;
		for (long i = 0; i < narrays; i++) {
			if ((i % 10000 == 0 && verbose) || (verbose >= 3)) {
				log_print(QUIET, "oaconvert: loading arrays: %d/%d\n", i, input_array_file.narrays);
			}

			int g = input_array_file.read_array(al);
			if (g < 0) {
				if (verbose)
					printf("   oaconvert: read_array returned index %d, end of file\n", g);
				break;
			}

			afout.append_array(al);
		}

		afout.finisharrayfile();
	}

}

arraylist_t readarrayfile (const char *fname, int verbose, int *setcols) {
        arraylist_t v;
        readarrayfile (fname, &v, verbose, setcols);
        return v;
}


int readarrayfile (const char *fname, arraylist_t *arraylist, int verbose, int *setcols, int *setrows,
                   int *setbits) {
        if (arraylist == 0) {
                arraylist = new arraylist_t;
        }

        arrayfile_t *afile = new arrayfile_t (fname, verbose);
        if (setcols != 0) {
                (*setcols) = afile->ncols;
        }
        if (setrows != 0) {
                (*setrows) = afile->nrows;
        }
        if (setbits != 0) {
                (*setbits) = afile->nbits;
        }

        if (!afile->isopen ()) {
                if (verbose) {
                        myprintf ("readarrayfile: problem with file %s\n", fname);
                        myprintf (" readarrayfile %d\n", (int)arraylist->size ());
                }

                return 0;
        }

        int narrays = afile->narrays;
        if (afile->narrays < 0) {
                narrays = arrayfile_t::NARRAYS_MAX;
        }

        array_link *alink;
        long i;
        int index;
        for (i = 0; i < narrays; i++) {

                if (i % 10000 == 0 && verbose) {
                        log_print (QUIET, "readarrayfile: loading arrays: %d/%d\n", i, afile->narrays);
                }

                if (verbose >= 3) {
                        myprintf ("readarrayfile: i %ld\n", i);
                }
                alink = new array_link (afile->nrows, afile->ncols, i + 1);
                index = afile->read_array (*alink);

                // NOTE: for array files we have used -1 to read arrays to the end of file,
                if (index < 0) {
                        myprintf ("readarrayfile: index %d, problem?\n", index);
                        delete alink;
                        break;
                }
                if (verbose >= 4) {
                        alink->showarray ();
                }
                arraylist->push_back (*alink);
                delete alink;
        }

        delete afile;
        return i;
}

int writearrayfile (const char *fname, const arraylist_t *arraylist, arrayfile::arrayfilemode_t mode, int nrows,
                    int ncols) {
        int nb = 8; // default: char

        if (arraylist->size () == 0) {
                if (mode == arrayfile::ABINARY_DIFFZERO) {
                        // special case: diffzero should always have 1 bit files
                        nb = 1;
                }
                if (nrows <= 0) {
                        myprintf ("writearrayfile: warning: empty list, using nrows %d, ncols %d\n", nrows, ncols);
                }
        } else {
                nrows = arraylist->at (0).n_rows;
                ncols = arraylist->at (0).n_columns;

                nb = arrayfile_t::arrayNbits (arraylist->at (0));
        }

        arrayfile_t *afile = new arrayfile_t (fname, nrows, ncols, arraylist->size (), mode, nb);

        if (!afile->isopen ()) {
                myprintf ("writearrayfile: problem with file %s\n", fname);
                return 0;
        }

        int i = afile->append_arrays (*arraylist, 1); 
        afile->finisharrayfile ();
        delete afile;

        return i;
}

int writearrayfile(const char *filename, const arraylist_t &arraylist, arrayfile::arrayfilemode_t mode, int nrows,
	int ncols) {
	return writearrayfile(filename, &arraylist, mode, nrows, ncols);
}

arrayfile_t::arrayfile_t () {
        this->verbose = 0;

        this->nfid = 0;
        this->filename = "";
#ifdef USEZLIB
        this->gzfid = 0;
#endif

        this->mode = ATEXT;

        this->rwmode = READ;

        this->nrows = 0;
        this->ncols = 0;
        this->narrays = 0;
        this->narraycounter = 0;
}

void arrayfile_t::createfile (const std::string fname, int nrows, int ncols, int narrays, arrayfilemode_t m, int nb) {
        this->closefile ();

        this->verbose = 0;
        this->filename = fname;

        this->iscompressed = 0;
#ifdef USEZLIB
        this->gzfid = 0;
#endif
        this->nrows = nrows;
        this->ncols = ncols;
        this->narrays = narrays;
        this->narraycounter = 0;
        this->mode = m;
        this->rwmode = WRITE;
        this->nbits = nb; // only used in binary mode

        if (strcmp (fname.c_str (), "<standard output>") == 0) {
                this->nfid = stdout;
        } else {
                this->nfid = fopen (fname.c_str (), "w+b");
                if (this->nfid == 0) {
                        myprintf ("arrayfile_t: ERROR: problem opening %s\n", fname.c_str ());
                        return;
                }
        }
        writeheader ();
}

arrayfile_t::arrayfile_t (const std::string fnamein, int verbose) {
        int warngz = verbose >= 2;
#ifdef SWIG
        swigcheck ();
#endif

        narraycounter = -1; // make invalid

        this->verbose = verbose;

        std::string fname = fnamein;
        this->rwmode = arrayfile::READ;
        this->filename = fname;
        this->narrays = -1; // initialize to -1

        if (verbose >= 3) {
                printfd ("start\n");
        }

        std::string gzname = fname + ".gz";
        if (!file_exists (fname.c_str ()) && file_exists (gzname.c_str ())) {
                if (verbose && warngz) {
                        myprintf ("  file %s does not exist, but gzipped file does\n", fname.c_str ());
                }
                this->filename = gzname;
                fname = gzname;
        }

#ifdef USEZLIB
        // simple check for compressed files
        if (verbose >= 2) {
                printfd ("zlib file %s\n", fname.c_str ());
        }
        if (fname.substr (fname.find_last_of (".") + 1) == "gz") {
                this->iscompressed = 1;
                this->gzfid = gzopen (fname.c_str (), "rb");
                this->nfid = 0;
                if (verbose >= 2) {
                        myprintf ("   opened file |%s| in compressed mode: nfid %ld, gzfid %ld\n", fname.c_str (),
                                  (long)this->nfid, (long)this->gzfid);
                }
        } else {
                this->iscompressed = 0;
                this->nfid = fopen (fname.c_str (), "r+b");

                this->gzfid = 0;
                if (verbose >= 2) {
                        myprintf ("   opened file |%s|: nfid %ld, gzfid %ld\n", fname.c_str (), (long)this->nfid,
                                  (long)this->gzfid);
                }
        }
#else
        this->iscompressed = 0;
        this->nfid = fopen (fname.c_str (), "r+b");
        this->gzfid = 0;
#endif

        if (verbose >= 2) {
                myprintf (" nfid %ld, gzfid %ld, isopen %d \n", (long)nfid, (long)gzfid, this->isopen ());
        }

        if (!this->isopen ()) {
                if (verbose) {
                        myprintf ("problem opening file %s (iscompressed %d)\n", fname.c_str (), this->iscompressed);
                }
                this->closefile ();
                return;
        }
        int magic;
        int result = afread (&magic, 1, sizeof (int32_t));
        if (result == 0) {
                this->closefile ();
                return;
        }
        if (this->nfid) {
                fclose (this->nfid);
                this->nfid = 0;
        }
#ifdef USEZLIB
        if (this->gzfid) {
                gzclose (this->gzfid);
                this->gzfid = 0;
        }
#endif

        if (verbose >= 2) {
                myprintf ("arrayfile_t: reading array file %s, magic %d\n", fname.c_str (), magic);
        }

        // read the header
        if (magic == 65) {
// we have a file in binary format
#ifdef USEZLIB
                if (iscompressed) {
                        this->gzfid = gzopen (fname.c_str (), "rb");
                } else {
                        this->nfid = fopen (fname.c_str (), "r+b");
                }
#else
                this->nfid = fopen (fname.c_str (), "r+b");
#endif

                int result = afread (&magic, 1, sizeof (int32_t));
                result = afread (&this->nbits, sizeof (int32_t), 1);
                assert (result == 1);
                result = afread (&this->nrows, sizeof (int32_t), 1);
                assert (result == 1);
                result = afread (&this->ncols, sizeof (int32_t), 1);
                assert (result == 1);
                result = afread (&this->narrays, sizeof (int32_t), 1);

                int binary_mode;

                int reserved;
                result = afread (&binary_mode, sizeof (int32_t), 1);
                assert (result == 1);
                result = afread (&reserved, sizeof (int32_t), 1);
                result = afread (&reserved, sizeof (int32_t), 1);

                switch (binary_mode) {
                case 1001:
                        this->mode = arrayfile::ABINARY;
                        break;
                case 1002:
                        this->mode = arrayfile::ABINARY_DIFF;
                        break;
                case 1003:
                        this->mode = arrayfile::ABINARY_DIFFZERO;
                        break;
                case 0:
                        if (verbose) {
                                myprintf ("  arrayfile_t::arrayfile_t: legacy file format, file %s!\n",
                                          this->filename.c_str ());
                        }
                        this->mode = arrayfile::ABINARY;
                        break;
                default:
                        this->mode = arrayfile::AERROR;
                        myprintf ("arrayfile_t::arrayfile_t : error opening binary file (binary_mode %d)\n",
                                  binary_mode);
                        break;
                }

                if (this->mode == arrayfile::ABINARY_DIFFZERO && this->nbits != 1) {
                        myprintf ("arrayfile_t: error for mode ABINARY_DIFFZERO we need 1 bit: bits %d\n",
                                  this->nbits);
                        this->mode = AERROR;
                        this->closefile ();
                        return;
                }

                if (result != 1) {
                        myprintf ("open binary file: wrong count in afread! %d\n", result);
                }
        } else {
                if (iscompressed) {
                        if (verbose) {
                                myprintf ("   compressed file: file or corrupt or gzipped text file, cannot read "
                                          "compressed files in TEXT mode..\n");
                        }
                        this->mode = AERROR;
                        this->closefile ();
                        return;
                } else {
                        this->nfid = fopen (fname.c_str (), "rb");
                        this->gzfid = 0;
                }

                this->mode = arrayfile::ATEXT;

                char buf[1];
                buf[0] = -1;
                int r = fread (buf, sizeof (char), 1, this->nfid);
                if (buf[0] < 48 || buf[0] > 57 || r < 0) {
                        if (verbose >= 1) {
                                myprintf ("   problem opening file %s (iscompressed %d)\n", fname.c_str (),
                                          this->iscompressed);
                        }
                        this->closefile ();
                        return;
                }
                r = fseek (this->nfid, 0, SEEK_SET);

                r = fscanf (this->nfid, "%i %i %i\n", &this->ncols, &this->nrows, &this->narrays);
                this->nbits = 0;
                if (verbose >= 2) {
                        myprintf ("arrayfile_t: text mode: header %d %d %d\n", this->ncols, this->nrows,
                                  this->narrays);
                }
                if (this->nrows < 0 || this->nrows > 20000 || this->ncols < 0 || this->ncols > 10000 ||
                    this->narrays < 0) {
                        if (verbose >= 1) {
                                myprintf ("   problem opening file %s (iscompressed %d)\n", fname.c_str (),
                                          this->iscompressed);
                        }
                        this->nfid = 0;
                        this->gzfid = 0;
                }
        }

        narraycounter = 0;

        if (verbose >= 2) {
                myprintf ("   arrayfile_t::arrayfile_t: mode %d, nrows %d, ncols %d, narrays %d\n", mode, this->nrows,
                          this->ncols, this->narrays);
        }
}

void arrayfile_t::closefile () {
        if (verbose >= 2) {
                myprintf ("arrayfile_t::closefile(): nfid %ld\n", (long)nfid);
        }

        if (!this->isopen ()) {
                if (verbose >= 2) {
                        myprintf ("arrayfile_t::closefile(): file already closed\n");
                }
                return;
        }

        if (verbose >= 3) {
                myprintf ("arrayfile_t::closefile: narrays: %d\n", narrays);
                myprintf ("arrayfile_t::closefile: narraycounter: %d\n", narraycounter);
                myprintf ("arrayfile_t::closefile: rwmode: %d\n", rwmode);
        }

        // for a binary file update the number of arrays stored
        updatenumbers ();

        // close file handles
        if (this->nfid != 0) {
                fclose (this->nfid);
                this->nfid = 0;
        }
#ifdef USEZLIB
        if (this->gzfid != 0) {
                gzclose (this->gzfid);
                this->gzfid = 0;
        }
#endif
        // delete afile;
}

void arrayfile_t::updatenumbers () {
        if (narraycounter >= 0 && narrays == -1 && (this->rwmode == WRITE || this->rwmode == READWRITE) &&
            this->isbinary ()) {
                if (verbose >= 2) {
                        myprintf ("arrayfile_t: updating numbers %d->%d\n", narrays, narraycounter);
                }
                if (verbose >= 3) {
                        myprintf ("arrayfile_t: nfid %ld\n", long(nfid));
                }
                if (nfid != 0) {
                        long pos = ftell (nfid);
                        int r = fseek (nfid, 4 * sizeof (int32_t), SEEK_SET);
                        r = this->afwrite (&narraycounter, sizeof (int32_t), 1);
                        if (verbose >= 2) {
                                myprintf ("   arrayfile_t::updatenumbers: result of write: %d\n", (int)r);
                        }
                        fseek (nfid, pos, SEEK_SET); // place back pointer
                }
        }
}

arrayfile_t::~arrayfile_t () {
#ifdef SWIG
        swigcheck ();
#endif

        if (verbose >= 2) {
                myprintf ("arrayfile_t: destructor: filename %s, nfid %ld, narraycounter %d, this->rwmode %d\n",
                          filename.c_str (), (long)nfid, narraycounter, this->rwmode);
        }

        closefile ();
}


void arrayfile_t::read_array_binary (array_t *array, const int nrows, const int ncols) {
        switch (this->nbits) {
        case 1: {

                // construct bit array
                BIT_ARRAY *bitarr = bit_array_create (nrows * ncols);
                word_addr_t num_of_words = nwords (nrows * ncols);
                int result = afread (bitarr->words, num_of_words, sizeof (word_t));
                assert (result);

                // fill bit array
                for (int i = 0; i < nrows * ncols; i++) {
                        array[i] = bit_array_get_bit_nocheck (bitarr, i);
                }
                bit_array_free (bitarr);
                break;
        }
        case 8:
#ifdef USEZLIB
                if (this->iscompressed) {
                        readblob< char, array_t > (array, nrows * ncols, this->gzfid);
                } else
#endif
                        readblob< char, array_t > (array, nrows * ncols, this->nfid);
                break;
        case 32:
#ifdef USEZLIB
                if (this->iscompressed) {
                        readblob< int32_t, array_t > (array, nrows * ncols, this->gzfid);
                } else
#endif
                        readblob< int32_t, array_t > (array, nrows * ncols, this->nfid);
                break;
        default:
                myprintf (" no such number of bits supported\n");
                break;
        }
}


void arrayfile_t::write_array_binary_diff (const array_link &A) {
        myassert (this->rwmode == WRITE || this->rwmode == READWRITE, "error: arrayfile_t not in write mode");
        myassert (A.max () <= 1, "arrayfile_t::write_array_binary_diff: array is not binary");

        int ngood = 0;

        rowindex_t N = this->nrows;
        const int num = N * sizeof (array_t);

        for (int i = 0; i < diffarray.n_columns; i++) {
                if (!memcmp (this->diffarray.array + N * i, A.array + N * i, num)) {
                        ngood++;

                } else {
                        break;
                }
        }

        int32_t nwrite = A.n_columns - ngood;
        array_link rest = A.selectLastColumns (nwrite);

        int n = fwrite ((const void *)&nwrite, sizeof (int32_t), 1, nfid);

        this->write_array_binary (rest);

        // update with previous array
        this->diffarray = A;
}

int arrayfile_t::getnbits () { return nbits; }

arrayfile::arrayfilemode_t arrayfile_t::parseModeString (const std::string format) {
          arrayfile::arrayfilemode_t mode = arrayfile::ATEXT;
          if (format == "AUTO" || format == "A") {
               mode = arrayfile::A_AUTOMATIC;

          } else {
               if (format == "BINARY" || format == "B") {
                         mode = arrayfile::ABINARY;
               } else {
                         if (format == "D" || format == "DIFF") {
                              mode = arrayfile::ABINARY_DIFF;
                         } else {
                              if (format == "Z" || format == "DIFFZERO") {
                                        mode = arrayfile::ABINARY_DIFFZERO;
                              } else {
                                        if (format == "AB" || format == "AUTOBINARY") {
                                                  mode = arrayfile::A_AUTOMATIC_BINARY;
                                        } else {
                                                  mode = arrayfile::ATEXT;
                                        }
                              }
                         }
               }
          }
          return mode;
}

/**
 * @brief Write an array in binary mode to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void arrayfile_t::write_array_binary_diffzero (const array_link &A) {
        myassert (this->rwmode == WRITE || this->rwmode == READWRITE, "error: arrayfile_t not in write mode");
        int ngood = 0;

        rowindex_t N = this->nrows;
        const int num = N * sizeof (array_t);

        for (int i = 0; i < diffarray.n_columns; i++) {
                if (!memcmp (this->diffarray.array + N * i, A.array + N * i, num)) {
                        ngood++;

                } else {
                        break;
                }
        }

        int16_t nwrite = A.n_columns - ngood;
        array_link rest = A.selectLastColumns (nwrite);

        int n = fwrite ((const void *)&nwrite, sizeof (int16_t), 1, nfid);

        if (diffarray.n_columns > 0) {

                array_link diffrest = this->diffarray.selectLastColumns (nwrite);

                array_link z = rest - diffrest;
                for (int i = 0; i < z.n_columns * z.n_rows; i++) {
                        z.array[i] = (2 + z.array[i]) % 2;
                }

                this->write_array_binary (z);
        } else {
                array_link z = rest;
                this->write_array_binary (z);
        }

        // update with previous array
        this->diffarray = A;
}

/// Write an array in binary mode to a file
void arrayfile_t::write_array_binary (const array_link &A) { 
     write_array_binary (A.array, A.n_rows, A.n_columns);     
}

/**
 * @brief Write an array in binary mode to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void arrayfile_t::write_array_binary (carray_t *array, const int nrows, const int ncols) {
        int m = 0;
        for (int i = 0; i < nrows * ncols; ++i)
                if (array[i] > m) {
                        m = array[i];
                }
        if (this->nbits == 1) {
                if (m > 1) {
                        myprintf ("arrayfile_t::write_array_binary: ERROR!\n");
                }
        }

        if (sizeof (array_t) == sizeof (int32_t) && this->nbits == 32) {
                int n = fwrite ((const void *)array, sizeof (int32_t), nrows * ncols, nfid);
        } else {
                switch (this->nbits) {
                case 8:
                        writeblob< array_t, char > (array, nrows * ncols, nfid);
                        break;
                case 32:
                        writeblob< array_t, int32_t > (array, nrows * ncols, nfid);
                        break;
                case 1: {
                        // construct bit array
                        BIT_ARRAY *bitarr = bit_array_create (nrows * ncols);
                        // fill bit array
                        for (int i = 0; i < nrows * ncols; i++) {
                                if (array[i]) {
                                        bit_array_set_bit (bitarr, i);
                                } else {
                                        bit_array_clear_bit (bitarr, i);
                                }
                        }
                        word_addr_t num_of_words =
                            nwords (nrows * ncols); // printf("num_of_words: %d\n", (int)num_of_words);

                        fwrite (bitarr->words, num_of_words, sizeof (word_t), this->nfid);

                        // printf ( "1-bit write: %ld bytes for array with %d elements\n", num_of_words*sizeof ( word_t
                        // ), nrows*ncols );
                        bit_array_free (bitarr);
                }

                break;
                default:
                        myprintf ("error: number of bits not supported\n");
                        break;
                }
        }
}

/** Read header for binary data file. Return true if valid header file
*
* The header consists of 4 integers: 2 magic numbers, then the number of rows and columns
*/
bool readbinheader(FILE *fid, int &nr, int &nc) {
	if (fid == 0) {
		return false;
	}

	double h[4];
	int nn = fread(h, sizeof(double), 4, fid);
	nr = (int)h[2];
	nc = (int)h[3];

	// myprintf("readbinheader: nn %d magic %f %f %f %f check %d %d\number_of_arrays", nn, h[0], h[1], h[2], h[3],
	// h[0]==30397995, h[1]==12224883);
	bool valid = false;

	// check 2 numbers of the magic header
	if (nn == 4 && h[0] == 30397995 && h[1] == 12224883) {
		return true;
	}

	return valid;
}

/// Write header for binary data file
void writebinheader(FILE *fid, int nr, int nc) {
	double h[4];
	// write 2 numbers of the magic header
	h[0] = 30397995;
	h[1] = 12224883;
	h[2] = nr;
	h[3] = nc;
	fwrite(h, sizeof(double), 4, fid);
}

/// Write a vector of vector elements to binary file
void vectorvector2binfile(const std::string fname, const std::vector< std::vector< double > > vals,
	int writeheader, int na) {
	FILE *fid = fopen(fname.c_str(), "wb");

	if (fid == 0) {
		fprintf(stderr, "vectorvector2binfile: error with file %s\n", fname.c_str());

		throw_runtime_exception("vectorvector2binfile: error with file");
	}

	if (na == -1) {
		if (vals.size() > 0) {
			na = vals[0].size();
		}
	}
	if (writeheader) {
		writebinheader(fid, vals.size(), na);
	}
	else {
		myprintf("warning: legacy file format\n");
	}
	for (unsigned int i = 0; i < vals.size(); i++) {
		const std::vector< double > x = vals[i];
		if ((int)x.size() != na) {
			myprintf("error: writing incorrect number of elements to binary file\n");
		}
		for (unsigned int j = 0; j < x.size(); j++) {
			fwrite(&(x[j]), sizeof(double), 1, fid);
		}
	}
	fclose(fid);
}

int writearrayfile (const char *filename, const array_link &array, arrayfile::arrayfilemode_t mode) {
        arraylist_t s;
        s.push_back (array);
        return writearrayfile (filename, s, mode);
}

/// append a single array to an array file. creates a new file if no file exists
int append_arrayfile (const char *fname, const array_link al) {
        int dverbose = 0;
        int nb = 8; // default: char
        int nrows = -1, ncols = -1;
        arrayfilemode_t mode = arrayfile::ABINARY;

        arraylist_t arraylist;
        arraylist.push_back (al);

        arrayfile_t *afile = new arrayfile_t (fname);

        if (dverbose) {
                myprintf ("\n### appendarrayfile: opened array file: %s\n", afile->showstr ().c_str ());
        }

        if (!afile->isopen ()) {
                myprintf ("appendarrayfile: creating new array file %s\n", fname);
                {
                        nrows = al.n_rows;
                        ncols = al.n_columns;
                        nb = arrayfile_t::arrayNbits (al);
                }

                afile = new arrayfile_t (fname, nrows, ncols, -1, mode, nb);
        } else {
                if (dverbose) {
                        printfd ("file is at position %d\n", ftell (afile->nfid));
                }

        }

        afile->seek (afile->narrays);

        if (!afile->isopen ()) {
                myprintf ("writearrayfile: problem with file %s\n", fname);
                return 0;
        }

        int i = afile->append_arrays (arraylist, 1); 
        afile->narrays = -1;
        afile->rwmode = READWRITE;

        afile->finisharrayfile ();
        delete afile;

        return i;
}

void selectArrays (const std::string filename, std::vector< int > &idx, arraylist_t &rl, int verbose) {
        arrayfile_t af (filename, 0);
        array_link al (af.nrows, af.ncols, -1);
        if (af.mode == ABINARY) {
                for (std::vector< int >::iterator it = idx.begin (); it < idx.end (); ++it) {
                        if (verbose) {
                                myprintf ("selectArrays: idx %d\n", *it);
                        }
                        af.seek (*it);
                        af.read_array (al);
                        rl.push_back (al);
                }
        } else {
                // check whether list is sorted
                indexsort vv (idx);

                if (vv.issorted ()) {
                        int cpos = 0;
                        for (std::vector< int >::iterator it = idx.begin (); it < idx.end (); ++it) {
                                int pos = *it;
                                if (verbose) {
                                        myprintf ("selectArrays: idx %d, cpos %d\n", pos, cpos);
                                }
                                if (pos < cpos) {
                                        myprintf ("selectArrays: arrayfile in text mode and negative seek, "
                                                  "aborting!!! %d %d\n",
                                                  pos, cpos);
                                        return;
                                }
                                int nsk = pos - cpos;
                                if (verbose) {
                                        myprintf ("selectArrays: skipping %d arrays\n", nsk);
                                }
                                for (int j = 0; j < (nsk); j++) {
                                        af.read_array (al);
                                        if (verbose >= 3) {
                                                std::vector< double > tmp = al.GWLP ();
                                                myprintf ("  gwlp: ");
                                                display_vector (tmp);
                                                myprintf ("\n");
                                        }
                                        cpos++;
                                }
                                af.read_array (al);
                                cpos++;
                                rl.push_back (al);
                        }
                } else {
                        if (verbose >= 2) {
                                myprintf ("selectArrays: text file with unsorted indices, UNTESTED CODE!\n");
                        }
                        if (verbose) {
                                myprintf ("selectArrays: no sorted indices!\n");
                        }
                        if (verbose >= 2) {
                                cout << "idx: ";
                                display_vector< int > (idx);
                                std::cout << std::endl;
                        }
                        // not sorted!
                        int nn = vv.indices.size ();
                        std::vector< int > sidx = vv.sorted (idx);
                        if (verbose >= 2) {
                                std::cout << "sidx: ";
                                display_vector< int > (sidx);
                                std::cout << std::endl;
                        }
                        arraylist_t tmp;
                        if (verbose >= 2) {
                                std::cout << "vv.indices: ";
                                display_vector< int > (vv.indices);
                                std::cout << std::endl;
                        }
                        selectArrays (filename, sidx, tmp, 0);
                        if ((int)tmp.size () != nn) {
                                myprintf ("ERROR!\n");
                                return;
                        }
                        for (int j = 0; j < nn; j++) {
                                int x = vv.indices[j];
                                rl.push_back (tmp[x]);
                        }
                }
        }
}

array_link selectArrays (const std::string filename, int ii) {
        arrayfile_t af (filename);
        array_link al (af.nrows, af.ncols, -1);
        if (ii < 0) {
                myprintf ("selectArrays: error: negative index\n");
                return al;
        }
        if (af.mode == ABINARY) {
                af.seek (ii);
                af.read_array (al);
        } else {
                for (int i = 0; i < ii; i++) {
                        af.read_array (al);
                }
                af.read_array (al);
        }
        return al;
}

#endif // FULLPACKAGE, related to arrayfile_t

void selectArrays (const arraylist_t &input_list, std::vector< int > &idx, arraylist_t &output_list) {
        for (std::vector< int >::iterator it = idx.begin (); it < idx.end (); it++) {
                output_list.push_back (input_list.at (*it));
        }
}
void selectArrays (const arraylist_t &input_list, std::vector< long > &idx, arraylist_t &output_list) {
        for (std::vector< long >::iterator it = idx.begin (); it < idx.end (); it++) {
                if ((*it) >= 0 && ((*it) < (int)input_list.size ())) {
                        output_list.push_back (input_list.at (*it));
                } else {
                        myprintf ("selectArrays: index %ld out of bounds!\n", *it);
                }
        }
}

arraylist_t selectArrays (const arraylist_t &input_list, std::vector< int > &idx) {
        arraylist_t rl;
        for (std::vector< int >::iterator it = idx.begin (); it < idx.end (); ++it) {
                int val = *it;
                if (val >= 0 && val <= (int)input_list.size ()) {
                        rl.push_back (input_list.at (val));
                } else {
                        myprintf ("selectArrays: error: index out of bounds: index %d, size %d\n", val,
                                  (int)input_list.size ());
                }
        }
        return rl;
}

arraylist_t selectArrays (const arraylist_t &input_list, std::vector< long > &idx) {
        arraylist_t rl;
        for (std::vector< long >::iterator it = idx.begin (); it < idx.end (); ++it) {
                int val = *it;
                if (val >= 0 && val <= (int)input_list.size ()) {
                        rl.push_back (input_list.at (val));
                } else {
                        myprintf ("selectArrays: error: index out of bounds: index %d, size %ld\n", val,
                                  (long)input_list.size ());
                }
        }
        return rl;
}

array_link array_link::operator+ (array_t value) const {
        array_link tmp = (*this);
        for (int i = 0; i < tmp.n_columns * tmp.n_rows; i++) {
                tmp.array[i] += value;
        }
        return tmp;
}
array_link array_link::operator- (array_t value) const {
        array_link tmp = (*this);
        for (int i = 0; i < tmp.n_columns * tmp.n_rows; i++) {
                tmp.array[i] -= value;
        }
        return tmp;
}

array_link array_link::operator+ (const array_link &array) const {
        assert (this->equalsize (array));
        array_link tmp = (*this);
        for (int i = 0; i < tmp.n_columns * tmp.n_rows; i++) {
                tmp.array[i] += array.array[i];
        }
        return tmp;
}

array_link array_link::operator- (const array_link &array) const {
        assert (this->equalsize (array));
        array_link tmp = (*this);
        for (int i = 0; i < tmp.n_columns * tmp.n_rows; i++) {
                tmp.array[i] -= array.array[i];
        }
        return tmp;
}

array_link array_link::operator* (const array_link &rhs) const {
        assert (this->equalsize (rhs));
        array_link tmp = (*this);
        for (int i = 0; i < tmp.n_columns * tmp.n_rows; i++) {
                tmp.array[i] *= rhs.array[i];
        }
        return tmp;
}

array_link array_link::operator* (array_t val) const {
          array_link al (*this);
          int NN = this->n_rows * this->n_columns;
          for (int i = 0; i < NN; i++) {
               al.array[i] *= val;
          }
          return al;
}

array_link array_link::operator*= (array_t val) {
          int NN = this->n_rows * this->n_columns;
          for (int i = 0; i < NN; i++) {
               this->array[i] *= val;
          }
          return *this;
}

array_link array_link::operator+= (array_t val) {
          int NN = this->n_rows * this->n_columns;
          for (int i = 0; i < NN; i++) {
               this->array[i] += val;
          }
          return *this;
}
array_link array_link::operator-= (array_t val) {
          int NN = this->n_rows * this->n_columns;
          for (int i = 0; i < NN; i++) {
               this->array[i] -= val;
          }
          return *this;
}
        
/// stack two arrays together
array_link vstack (const array_link &al, const array_link &b) {
        assert (al.n_columns == b.n_columns);
        array_link v (al.n_rows + b.n_rows, al.n_columns, array_link::INDEX_NONE);
        int N1 = al.n_rows;
        int N2 = b.n_rows;
        int N = N1 + N2;
        for (int c = 0; c < al.n_columns; c++) {
                std::copy (al.array + c * N1, al.array + (c + 1) * N1, v.array + c * N);
                std::copy (b.array + c * N2, al.array + (c + 1) * N2, v.array + c * N + N1);
        }
        return v;
}

array_link hstack (const array_link &array, const conference_column &column) {
        assert (array.n_rows == (int)column.size ());
        array_link v (array.n_rows, array.n_columns + 1, array_link::INDEX_NONE);
        std::copy (array.array, array.array + array.n_columns * array.n_rows, v.array);
        std::copy (column.begin (), column.end (), v.array + v.n_rows * array.n_columns);
        return v;
}

array_link hstack (const array_link &array_left, const array_link &array_right) {
        assert (array_left.n_rows == array_right.n_rows);
        array_link output_array (array_left.n_rows, array_left.n_columns + array_right.n_columns, array_link::INDEX_NONE);
        std::copy (array_left.array, array_left.array + array_left.n_columns * array_left.n_rows, output_array.array);
        std::copy (array_right.array, array_right.array + array_right.n_columns * array_left.n_rows, output_array.array + output_array.n_rows * array_left.n_columns);
        return output_array;
}

array_link hstacklastcol (const array_link &al, const array_link &b) {
        array_link v (al.n_rows, al.n_columns + 1, array_link::INDEX_NONE);
        std::copy (al.array, al.array + al.n_columns * al.n_rows, v.array);
        size_t offset = al.n_rows * (b.n_columns - 1);
        std::copy (b.array + offset, b.array + offset + al.n_rows, v.array + v.n_rows * al.n_columns);
        return v;
}

/// concatenate two columns
conference_column vstack(const conference_column &column_top, const conference_column &column_bottom) {
	conference_column output_column(column_top.size() + column_bottom.size());

	std::copy(column_top.begin(), column_top.end(), output_column.begin());
	std::copy(column_bottom.begin(), column_bottom.end(), output_column.begin() + column_top.size());
	return output_column;
}

void conference_transformation_t::show (int verbose) const {
        myprintf ("row permutation: ");
        print_perm_int (rperm);
        myprintf ("  row flips: ");
        print_perm_int (rswitch);
        myprintf ("column permutation: ");
        print_perm_int (cperm);
        myprintf ("  col flips: ");
        print_perm_int (cswitch);
}

/// helper function to invert a permutation and sign switch
void invert_perm_switch (const std::vector< int > perm, const std::vector< int > switchin, std::vector< int > &permout,
                         std::vector< int > &switchout) {
        invert_permutation (perm, permout);

        for (size_t i = 0; i < perm.size (); i++) {
                switchout[permout[i]] = switchin[i];
        }
}

conference_transformation_t conference_transformation_t::inverse () const {
        conference_transformation_t I (nrows, ncols);

        invert_perm_switch (rperm, rswitch, I.rperm, I.rswitch);
        invert_perm_switch (cperm, cswitch, I.cperm, I.cswitch);
        return I;
}

array_link conference_transformation_t::apply (const array_link &al) const {
        array_link alx (al.n_rows, al.n_columns, al.index);
        array_link tmp (al.n_rows, al.n_columns, al.index);

        /* column permutations */
        perform_column_permutation (al, tmp, cperm);
        perform_row_permutation (tmp, alx, rperm);

        /* sign flips */
        for (size_t r = 0; r < (size_t)nrows; r++) {
                for (size_t c = 0; c < (size_t)ncols; c++) {
                        alx.array[r + nrows * c] *= rswitch[r] * cswitch[c];
                }
        }

        return alx;
}

int conference_transformation_t::operator== (const conference_transformation_t &rhs) const {
        if (this->ncols != rhs.ncols)
                return 0;
        if (this->nrows != rhs.nrows)
                return 0;

        if (this->rperm != rhs.rperm)
                return 0;
        if (this->cperm != rhs.cperm)
                return 0;
        if (this->rswitch != rhs.rswitch)
                return 0;
        if (this->cswitch != rhs.cswitch)
                return 0;

        return 1;
}

conference_transformation_t conference_transformation_t::operator* (const conference_transformation_t &rhs) const {
        const int N = this->nrows;
        const int ncols = this->ncols;

        conference_transformation_t c (N, ncols);

        const conference_transformation_t &lhs = *this;

        // perform the rows permutations
        composition_perm (rhs.rperm, lhs.rperm, c.rperm);

        // perform the column permutations
        composition_perm (rhs.cperm, lhs.cperm, c.cperm);

        /* rowsign switches */
        for (rowindex_t ri = 0; ri < N; ri++) {
                int riz = rhs.rperm[ri];
                int rix = c.rperm[ri];
                c.rswitch[rix] = lhs.rswitch[rix] * rhs.rswitch[riz];
        }

        /* column sign switches */
        for (colindex_t ci = 0; ci < ncols; ci++) {
                int ciz = rhs.cperm[ci];
                int cix = c.cperm[ci];
                c.cswitch[cix] = lhs.cswitch[cix] * rhs.cswitch[ciz];
        }

        return c;
}

void conference_transformation_t::init (int nr, int nc) {
        this->nrows = nr;
        this->ncols = nc;

        this->rperm.resize (nr);
        this->cperm.resize (nc);
        rswitch.resize (nr);
        cswitch.resize (nc);

        reset ();
}
void conference_transformation_t::reset () {
        for (size_t i = 0; i < rperm.size (); i++) {
                rperm[i] = i;
        }
        for (size_t i = 0; i < cperm.size (); i++) {
                cperm[i] = i;
        }
        std::fill (rswitch.begin (), rswitch.end (), 1);
        std::fill (cswitch.begin (), cswitch.end (), 1);
}

conference_transformation_t::conference_transformation_t (const conference_transformation_t &T) {
        this->nrows = T.nrows;
        this->ncols = T.ncols;

        this->rperm = T.rperm;
        this->rswitch = T.rswitch;

        this->cperm = T.cperm;
        this->cswitch = T.cswitch;
}

conference_transformation_t::conference_transformation_t () { init (1, 1); }
conference_transformation_t::conference_transformation_t (int nr, int nc) { init (nr, nc); }

conference_transformation_t::conference_transformation_t (const array_link &al) { init (al.n_rows, al.n_columns); }

bool conference_transformation_t::isIdentity () const {
        for (int i = 0; i < ncols; ++i) {
                if (cperm[i] != i) {
                        return 0;
                }
        }
        for (int i = 0; i < nrows; ++i) {
                if (rperm[i] != i) {
                        return 0;
                }
        }
        for (int c = 0; c < ncols; ++c) {
                if (cswitch[c] != 1) {
                        return 0;
                }
        }
        for (int c = 0; c < nrows; ++c) {
                if (rswitch[c] != 1) {
                        return 0;
                }
        }
        //	myprintf("isIdentity: yes\n");
        return 1;
}

/// initialize to a random column permutation
void conference_transformation_t::randomizecolperm () { random_perm (cperm); }

/// initialize to a random row permutation
void conference_transformation_t::randomizerowperm () { random_perm (rperm); }

void conference_transformation_t::randomizecolflips () {
        for (size_t i = 0; i < this->cswitch.size (); i++) {
                this->cswitch[i] = 2 * (rand () % 2) - 1;
        }
}

void conference_transformation_t::randomizerowflips () {
        for (size_t i = 0; i < this->rswitch.size (); i++) {
                this->rswitch[i] = 2 * (rand () % 2) - 1;
        }
}

/// initialize to a random transformation
void conference_transformation_t::randomize () {
        /* row and columns permutations */
        random_perm (rperm);
        random_perm (cperm);

        /* sign switches */
        for (size_t x = 0; x < rswitch.size (); x++) {
                rswitch[x] = 2 * (rand () % 2) - 1;
        }
        for (size_t x = 0; x < cswitch.size (); x++) {
                cswitch[x] = 2 * (rand () % 2) - 1;
        }
}

int arrayInList (const array_link &al, const arraylist_t &ll, int verbose) {
        for (size_t jj = 0; jj < ll.size (); jj++) {
                const array_link &alx = ll[jj];
                if (alx == al.selectFirstColumns (alx.n_columns)) {
                        if (verbose) {
                                myprintf ("arrayInList: found array at position %ld\n", jj);
                        }
                        return jj;
                }
        }
        return -1;
}

int arrayInFile (const array_link &al, const char *afile, int verbose) {
        bool found = false;
        arrayfile_t af (afile, 1);
        long aidx = -1;
        for (long jj = 0; jj < af.narrays; jj++) {

                array_link alx = af.readnext ();
                if (alx == al.selectFirstColumns (alx.n_columns)) {
                        if (verbose) {
                                myprintf ("arrayInFile: found array at position %ld\n", jj);
                        }
                        found = true;
                        aidx = jj;
                        break;
                }
        }
        if (found) {
                return aidx;
        } else {
                if (verbose) {
                        myprintf ("arrayInFile: could not find array\n");
                }
                return -1;
        }
        return -1;
}

// kate: indent-mode cstyle; indent-width 5; replace-tabs on;
