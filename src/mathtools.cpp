#include <numeric>

#include "mathtools.h"

void set_srand (unsigned int s) { srand (s); }

template < class Type > void symmetry_group::init (const std::vector< Type > vals, bool ascendingx, int verbose) {

        if (verbose >= 2) {
                myprintf ("symmetry_group::init: %ld elements: ", vals.size ());
                for (size_t i = 0; i < vals.size (); i++)
                        std::cout << vals[i] << ", ";
                myprintf ("\n");
        }

        n = vals.size ();
        ascending = ascendingx;

        if (verbose >= 2)
                myprintf ("symmetry_group::init: check sorting\n");

        indexsort is (vals.size ());

        if (ascending)
                is.sort (vals);
        else
                is.sortdescending (vals);

        if (verbose >= 2) {
                if (ascending) {
                        if (!is.issorted ()) {
                                myprintf ("symmetry_group: input group was not sorted!\n");
                                is.show ();
                                myprintf ("\n");
                        }
                } else {
                        if (!is.issorted ()) {
                                myprintf ("symmetry_group: input group was not sorted!\n");
                                is.show ();
                                myprintf ("\n");
                        }
                }
        }
        // calc group
        int nsg = 0;
        Type prev; 

        prev = std::numeric_limits< Type >::quiet_NaN ();
        /* count number of symmetry groups */
        for (int i = 0; i < n; i++) {
                if (vals[i] != prev || i == 0) {
                        nsg++;
                }
                prev = vals[i];
        }

        if (verbose) {
                myprintf ("symmetry_group: %d elements, %d groups\n", n, nsg);
        }

        this->ngroups = nsg;

        // initialize structures
        gidx.resize (n);
        gstart.resize (nsg + 1);
        gsize.resize (nsg + 1);

        /* find starting positions of symmetry group */
        if (verbose >= 2)
                myprintf ("symmetry_group::init:  find starting positions\n");

        prev = std::numeric_limits< Type >::quiet_NaN ();
        nsg = 0;
        int i;
        for (i = 0; i < n; i++) {
                if (vals[i] != prev || i == 0) {
                        gstart[nsg] = i;
                        nsg++;
                }
                gidx[i] = nsg - 1;
                prev = vals[i];
        }
        gstart[nsg] = n; /* add dummy element */

        for (int i = 0; i < nsg; i++) {
                gsize[i] = gstart[i + 1] - gstart[i];
        }
        gsize[nsg] = 0;
}

symmetry_group::symmetry_group (const std::vector< float > &vals, bool ascending, int verbose) {
        this->init< float > (vals, ascending, verbose);
}

symmetry_group::symmetry_group (const std::vector< double > &vals, bool ascending, int verbose) {
        this->init< double > (vals, ascending, verbose);
}

symmetry_group::symmetry_group (const std::vector< short int > &vals, bool ascending, int verbose) {
        this->init (vals, ascending, verbose);
}

symmetry_group::symmetry_group (const std::vector< unsigned int > &vals, bool ascending, int verbose) {
        this->init (vals, ascending, verbose);
}

symmetry_group::symmetry_group (const std::vector< int > &vals, bool ascending, int verbose) {
        this->init< int > (vals, ascending, verbose);
}

symmetry_group::symmetry_group (const std::vector< mvalue_t< double > > &vals, bool ascending, int verbose) {
        this->init (vals, ascending, verbose);
}

symmetry_group::symmetry_group (const std::vector< mvalue_t< int > > &vals, bool ascending, int verbose) {
        this->init (vals, ascending, verbose);
}

symmetry_group &symmetry_group::operator= (const symmetry_group &sgx) {
        gidx = sgx.gidx;
        gstart = sgx.gstart;
        gsize = sgx.gsize;
        ngroups = sgx.ngroups;
        n = sgx.n;
        ascending = sgx.ascending;

        return *this;
}

symmetry_group::symmetry_group (const symmetry_group &sgx) {
        gidx = sgx.gidx;
        gstart = sgx.gstart;
        gsize = sgx.gsize;
        ngroups = sgx.ngroups;
        n = sgx.n;
        ascending = sgx.ascending;
}
symmetry_group::symmetry_group () {
        ngroups = 0;
        n = 0;
        ascending = 0;
}

/// representation function (for python interface)
std::string symmetry_group::__repr__() const {
	std::stringstream ss;
	ss << "symmetry group: " << n << " elements, " << ngroups << " subgroups: ";
	for (int i = 0; i < ngroups; i++)
		ss << gsize[i] << " ";

	std::string s = ss.str();
	return s;
}

/// show the symmetry group
void symmetry_group::show(int verbose) const {
	myprintf("symmetry group: %d elements, %d subgroups: ", n, ngroups);
	for (int i = 0; i < ngroups; i++)
		myprintf("%d ", gsize[i]);
	myprintf("\n");

	if (verbose >= 2) {
		myprintf("gidx: ");
		for (int i = 0; i < n; i++)
			myprintf("%d, ", gidx[i]);
		myprintf("\n");
		myprintf("gstart: ");
		for (int i = 0; i < ngroups; i++)
			myprintf("%d, ", gstart[i]);
		myprintf("\n");
		myprintf("gsize: ");
		for (int i = 0; i < ngroups; i++)
			myprintf("%d, ", gsize[i]);
		myprintf("\n");
	}
}

void symmetry_group_walker::show (int verbose) const {
        myprintf ("symmetry_group_walker: ");
        if (verbose >= 2)
                myprintf ("\n");
        for (size_t i = 0; i < (size_t)perms.size (); i++) {
                if (verbose >= 2) {
                        myprintf ("  block %ld: ", i);
                        print_perm (perms[i]);
                } else {
                        print_perm (perms[i], 100, false);
                        myprintf (" ");
                }
        }
        if (verbose == 1)
                myprintf ("\n");
}

bool symmetry_group_walker::nextsub (int g) {
        bool of = next_perm (perms[g]);
        if (of && g > 0)
                of = nextsub (g - 1);

        return of;
}

std::vector< int > symmetry_group_walker::fullperm () const {
        std::vector< int > w (sg.n);
        for (int i = 0; i < sg.n; i++)
                w[i] = i;

        std::vector< int > ww (sg.n);
        for (size_t j = 0; j < perms.size (); j++) {
                for (size_t k = 0; k < perms[j].size (); k++) {
                        int offset = sg.gstart[j];
                        ww[offset + k] = w[offset + perms[j][k]];
                }
        }
        return ww;
}

/* Random number generators */

int g_seed = 123;
void seedfastrand (int s) { g_seed = s; }
int fastrand () {
        g_seed = (214013 * g_seed + 2531011);
        return (g_seed >> 16) & 0x7FFF;
}

/// calculate a random integer modulo K
int fastrandK (int K) { return fastrand () % K; }

int Combinations::ncombscachemax = 0;
long **Combinations::ncombsdata = 0;

void Combinations::initialize_number_combinations (int N) {
        if (N <= number_combinations_max_n ())
                return;

#pragma omp critical
        {
                const int rowsize = N + 1;
                const int nrows = N + 1;
                if (ncombsdata != 0) {
                        delete[] ncombsdata[0];
                        delete[] ncombsdata;
                }

                ncombsdata = new long *[nrows];

                ncombsdata[0] = new long[nrows * rowsize];

                int offset = 0;
                for (int i = 0; i < nrows; i++) {
                        ncombsdata[i] = ncombsdata[0] + offset;
                        offset += rowsize;
                }

                for (int i = 0; i < nrows; i++) {
                        for (int j = 0; j < rowsize; j++) {
                                ncombsdata[i][j] = ncombs (i, j);
                        }
                        ncombscachemax = N;
                }
        }
}

Combinations::~Combinations() {
	// we define the empty destructor because swig otherwise generates an error
}

int Combinations::number_combinations_max_n () { return Combinations::ncombscachemax; }

long Combinations::number_combinations (int n, int k) {
#ifdef SWIGPYTHON
        myassert (Combinations::ncombsdata != 0);
#endif
        return Combinations::ncombsdata[n][k];
}

