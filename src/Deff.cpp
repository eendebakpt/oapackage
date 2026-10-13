/*! \file Deff.cpp
 *  \brief Contains functions to optimize designs
 *
 */

#include <algorithm>
#include <cmath>
#include <memory>

#include <Eigen/Dense>

#include "printfheader.h"
#include "arrayproperties.h"
#include "arraytools.h"

#ifdef DOOPENMP
#include "omp.h"
#endif

#include "Deff.h"

#ifndef myprintf
#define myprintf printf
#endif

#ifdef MAINMEX
#define MATLABUPDATE
#else
#ifdef MATLAB_MEX
#define MATLABUPDATE mexEvalString ("drawnow");
#else
#define MATLABUPDATE
#endif
#endif

const double NaN = std::numeric_limits< double >::quiet_NaN();

DoptimReturn Doptimize (const arraydata_t &arrayclass, int nrestartsmax, std::vector< double > alpha, int verbose,
                        coordinate_exchange_method_t method, int niter, double maxtime, int nabort) {
        if (method == DOPTIM_AUTOMATIC)
                method = DOPTIM_UPDATE;

        double t0 = get_time_ms ();
        std::vector< std::vector< double > > dds;
        arraylist_t AA;

        bool abort = false;
        int nrestarts = 0;

#ifdef DOOPENMP
#pragma omp parallel for num_threads(4) schedule(dynamic, 1)
#endif
        for (int i = 0; i < nrestartsmax; i++) {
                if (abort)
                        continue;

#ifdef DOOPENMP
#pragma omp critical
#endif
                {
                        if (verbose && (i % 500 == 0 || i == nrestartsmax - 1)) {
                                myprintf ("Doptimize: iteration %d/%d\n", i, nrestartsmax);
                                MATLABUPDATE // helper macro from matlab interface
                        }
                }

                array_link al = arrayclass.randomarray (1);

                array_link A = optimDeff (al, arrayclass, alpha, verbose >= 2, method, niter, nabort);
                std::vector< double > dd = A.Defficiencies ();
                if (verbose >= 2) {
#ifdef DOOPENMP
#pragma omp critical
#endif
                        {
                                double score = scoreD (dd, alpha);
                                myprintf ("Doptimize: iteration %d/%d: %f %f %f: score %.3f\n", i, nrestartsmax, dd[0],
                                          dd[1], dd[2], score);
                        }
                }

#ifdef DOOPENMP
#pragma omp critical
#endif
                {
                        AA.push_back (A);
                        dds.push_back (dd);
                        nrestarts++;
                }

                if ((get_time_ms () - t0) > maxtime) {
                        if (verbose)
                                myprintf ("max running time exceeded, aborting\n");
#pragma omp critical
                        { abort = true; }
                }
        }

        // loop is complete
        DoptimReturn a = {dds, AA, nrestarts, nrestarts};
        return a;
}

double scoreD (const std::vector< double > efficiencies, const std::vector< double > weights) {
        double v = 0;
        for (size_t i = 0; i < efficiencies.size (); i++)
                v += efficiencies[i] * weights[i];
        return v;
}

DoptimReturn DoptimizeMixed (const arraylist_t &array_list, const arraydata_t &arrayclass, const std::vector< double > alpha,
                             int verbose, int nabort) {
        const size_t nn = array_list.size ();
        double t0 = get_time_ms ();
        std::vector< std::vector< double > > dds (nn);
        arraylist_t AA (nn);

        bool abort = false;
        int nimproved = 0;

        coordinate_exchange_method_t method1 = DOPTIM_SWAP;
        coordinate_exchange_method_t method2 = DOPTIM_UPDATE;
        const int niter = 600000;
        if (nabort < 0)
                nabort = arrayclass.N * arrayclass.ncols * 20 + 500;

        if (verbose >= 3)
                myprintf ("DoptimizeMixed: nabort %d\n", nabort);
#ifdef DOOPENMP
#pragma omp parallel for num_threads(4) schedule(dynamic, 1)
#endif
        for (size_t i = 0; i < nn; i++) {
                if (abort)
                        continue;

#ifdef DOOPENMP
#pragma omp critical
#endif
                {
                        if (verbose && (i % 100 == 0 || i == nn - 1)) {
                                myprintf ("DoptimizeMixed: iteration %d/%d\n", (int)i, (int)nn);
                                MATLABUPDATE // helper macro from matlab interface
                        }
                }

                const array_link &array = array_list[i];
                double score0 = scoreD (array.Defficiencies (), alpha);

                array_link alu = optimDeff (array, arrayclass, alpha, verbose >= 3, method1, niter, nabort);
                double score1 = scoreD (alu.Defficiencies (), alpha);

                array_link alu2 = optimDeff (alu, arrayclass, alpha, verbose >= 3, method2, niter, 0);
                double score2 = scoreD (alu2.Defficiencies (), alpha);

                dds[i] = alu2.Defficiencies ();
                AA[i] = alu2;

                if (score2 > score0) {
                        if (verbose >= 2)
                                myprintf ("DoptimizeMixed: array %d/%d: improve %.6f -> %.6f -> %.6f\n", (int)i,
                                          (int)nn, score0, score1, score2);

                        nimproved = nimproved + 1;
                }
        }

        double dt = get_time_ms () - t0;
        if (verbose) {
                myprintf ("DoptimizeMixed: improved %d/%ld arrays, %.2f [s]\n", nimproved, (long)nn, dt);
        }

        // loop is complete
        DoptimReturn a = {dds, AA, (int)nn, (int)nimproved};
        return a;
}

namespace {

/** Fill one row of the second order model matrix of a 2-level design
 *
 * The column convention is identical to the one used in array2eigenModelMatrixInt: intercept, main effects,
 * followed by the interactions (c, c2) with c2 < c in lexicographic order.
 */
inline void modelrow_2level (const array_link &A, int row, int k, double *x) {
        const int N = A.n_rows;
        const array_t *a = A.array + row;
        x[0] = 1;
        for (int c = 0; c < k; c++) {
                x[1 + c] = 2 * (double)a[c * N] - 1;
        }
        int ww = 1 + k;
        for (int c = 0; c < k; c++) {
                for (int c2 = 0; c2 < c; c2++) {
                        x[ww++] = 2 * (double)((a[c * N] + a[c2 * N]) % 2) - 1;
                }
        }
}

/** Cholesky factor and log-determinant of a principal submatrix of X^T X with rank-1 updates
 *
 * If p rows of the model matrix X change, the matrix M = X^T X changes by a matrix of rank at most 2p. The determinant
 * of the updated matrix follows from the matrix determinant lemma, using triangular solves with the Cholesky factor.
 * When the change is accepted, the factor is updated with rank-1 updates and downdates. Both operations cost O(n^2),
 * instead of the O(N n^2 + n^3) of a direct calculation. The Cholesky updates are numerically stable, so rounding
 * errors do not accumulate as fast as with updates of the inverse matrix.
 */
struct DeffSubsystem {
        /// indices of the columns of the full model matrix in this subsystem
        std::vector< int > idx;
        /// number of rows of the design
        int N = 0;
        /// lower triangular Cholesky factor of the current matrix M = L L^T
        Eigen::MatrixXd L;
        /// natural logarithm of the determinant of the current matrix
        double logdet = 0;
        /// true if the current matrix is singular (or could not be factored reliably)
        bool singular = true;
        /// number of rank updates applied since the last direct calculation
        int nupdates = 0;

        // data of the pending candidate update
        Eigen::MatrixXd U, Z;
        std::vector< double > work;
        int q = 0;
        double pending_logdet = 0;
        bool pending_singular = true;
        /// true if the pending candidate is too ill-conditioned for the incremental calculation
        bool pending_unreliable = false;

        int size () const { return (int)idx.size (); }

        /// normalized determinant: det(M/N)
        static double normalized_determinant (double logdet, int n, int N) { return exp (logdet - n * log ((double)N)); }

        /// calculate the Cholesky factor and the determinant directly from the full model matrix X
        void init (const Eigen::MatrixXd &X) {
                N = (int)X.rows ();
                const int n = size ();
                nupdates = 0;
                singular = true;

                Eigen::MatrixXd Xs (N, n);
                for (int j = 0; j < n; j++) {
                        Xs.col (j) = X.col (idx[j]);
                }
                Eigen::MatrixXd M = Xs.transpose () * Xs;
                Eigen::LLT< Eigen::MatrixXd > llt (M);
                if (llt.info () != Eigen::Success) {
                        return;
                }
                L = llt.matrixL ();
                logdet = 2 * L.diagonal ().array ().log ().sum ();
                if (normalized_determinant (logdet, n, N) < 1e-15) {
                        // same criterion as in Defficiencies: check the rank with a full pivoting LU
                        Eigen::FullPivLU< Eigen::MatrixXd > lu_decomp (Xs);
                        if (lu_decomp.rank () < n) {
                                return;
                        }
                }
                U.resize (n, 4);
                Z.resize (n, 4);
                work.resize (n);
                singular = false;
        }

        /** Evaluate the determinant of the matrix after a change of p rows of the model matrix
         *
         * The new matrix is M + U W^T with U = [xnew_1, xold_1, ...] and W = [xnew_1, -xold_1, ...]. The determinant
         * lemma gives det(M + U W^T) = det(M) det(I + W^T M^{-1} U) with W^T M^{-1} U = (L^{-1} W)^T (L^{-1} U).
         *
         * \param p Number of changed rows (1 or 2)
         * \param xold Full model matrix rows before the change
         * \param xnew Full model matrix rows after the change
         */
        void candidate (int p, const double *const *xold, const double *const *xnew) {
                const int n = size ();
                q = 2 * p;
                for (int i = 0; i < p; i++) {
                        for (int j = 0; j < n; j++) {
                                U (j, 2 * i) = xnew[i][idx[j]];
                                U (j, 2 * i + 1) = xold[i][idx[j]];
                        }
                }
                // forward substitution Z = L^{-1} U (a plain loop is much faster than the blocked Eigen solver for small n)
                const double *Ldata = L.data ();
                for (int col = 0; col < q; col++) {
                        double *z = Z.data () + col * n;
                        std::copy (U.data () + col * n, U.data () + (col + 1) * n, z);
                        for (int k = 0; k < n; k++) {
                                const double *Lk = Ldata + k * n;
                                const double zk = z[k] / Lk[k];
                                z[k] = zk;
                                for (int i = k + 1; i < n; i++) {
                                        z[i] -= Lk[i] * zk;
                                }
                        }
                }

                double detS = 0;
                if (q == 2) {
                        const double a = Z.col (0).dot (Z.col (0));
                        const double b = Z.col (0).dot (Z.col (1));
                        const double c = Z.col (1).dot (Z.col (1));
                        detS = (1 + a) * (1 - c) + b * b;
                } else {
                        Eigen::Matrix4d S;
                        for (int i = 0; i < q; i++) {
                                const double sign = (i % 2 == 0) ? 1 : -1;
                                for (int j = 0; j < q; j++) {
                                        S (i, j) = sign * Z.col (i).dot (Z.col (j)) + (i == j);
                                }
                        }
                        detS = S.determinant ();
                }

                pending_unreliable = false;
                pending_singular = false;
                if (detS < 1e-8) {
                        // the new matrix is singular or nearly singular: let the direct calculation decide
                        pending_unreliable = true;
                        pending_singular = true;
                        return;
                }
                pending_logdet = logdet + log (detS);
        }

        /// rank-1 update of the Cholesky factor: L L^T + x x^T
        void rank1_update (double *x) {
                const int n = size ();
                double *Ldata = L.data ();
                for (int k = 0; k < n; k++) {
                        double *Lk = Ldata + k * n;
                        const double lkk = Lk[k];
                        const double r = sqrt (lkk * lkk + x[k] * x[k]);
                        const double c = r / lkk;
                        const double s = x[k] / lkk;
                        Lk[k] = r;
                        for (int i = k + 1; i < n; i++) {
                                Lk[i] = (Lk[i] + s * x[i]) / c;
                                x[i] = c * x[i] - s * Lk[i];
                        }
                }
        }

        /// rank-1 downdate of the Cholesky factor: L L^T - x x^T. Returns false if the result is not positive definite
        bool rank1_downdate (double *x) {
                const int n = size ();
                double *Ldata = L.data ();
                for (int k = 0; k < n; k++) {
                        double *Lk = Ldata + k * n;
                        const double lkk = Lk[k];
                        const double r2 = lkk * lkk - x[k] * x[k];
                        if (r2 <= 0) {
                                return false;
                        }
                        const double r = sqrt (r2);
                        const double c = r / lkk;
                        const double s = x[k] / lkk;
                        Lk[k] = r;
                        for (int i = k + 1; i < n; i++) {
                                Lk[i] = (Lk[i] - s * x[i]) / c;
                                x[i] = c * x[i] - s * Lk[i];
                        }
                }
                return true;
        }

        /// apply the pending candidate update
        void accept () {
                if (pending_singular) {
                        singular = true;
                        return;
                }
                const int n = size ();
                // first the updates, then the downdates, so that the intermediate matrices stay positive definite
                for (int i = 0; i < q; i += 2) {
                        std::copy (U.col (i).data (), U.col (i).data () + n, work.data ());
                        rank1_update (work.data ());
                }
                for (int i = 1; i < q; i += 2) {
                        std::copy (U.col (i).data (), U.col (i).data () + n, work.data ());
                        if (!rank1_downdate (work.data ())) {
                                singular = true;
                                return;
                        }
                }
                logdet = 2 * L.diagonal ().array ().log ().sum ();
                nupdates++;
        }
};

/** Incremental calculation of the D-efficiencies of a 2-level design during coordinate exchange
 *
 * The efficiencies D, Ds and D1 (see Defficiencies) are determined by the determinants of three principal
 * submatrices of X^T X, with X the second order model matrix. The class keeps the inverses of these matrices and
 * updates them when a few rows of the design change. Whenever a matrix is singular, or an update is numerically
 * unreliable, the class falls back to the direct calculation of Defficiencies. After a number of updates the
 * inverses are recalculated from scratch to limit the accumulation of rounding errors.
 */
class Deff2levelUpdater {
      public:
        /** Create the incremental calculation for a design
         *
         * \param arrayclass Class of the design
         * \param A Initial design
         * \param alpha Weights of the efficiencies D, Ds and D1 in the score. Efficiencies with zero weight are not calculated
         * \param refresh_interval Number of accepted updates after which the state is recalculated directly
         */
        Deff2levelUpdater (const arraydata_t &arrayclass, const array_link &A, const std::vector< double > &alpha,
                           int refresh_interval = 64)
            : arrayclass (arrayclass), N (A.n_rows), k (A.n_columns), n2fi (k * (k - 1) / 2), m (1 + k + n2fi),
              refresh_interval (refresh_interval) {
                if (alpha.size () >= 3) {
                        need_full = alpha[0] != 0 || alpha[1] != 0;
                        need_02 = alpha[1] != 0;
                        need_1 = alpha[2] != 0;
                }
                for (int i = 0; i < m; i++) {
                        full.idx.push_back (i);
                }
                sub02.idx.push_back (0);
                for (int i = 1 + k; i < m; i++) {
                        sub02.idx.push_back (i);
                }
                for (int i = 0; i < 1 + k; i++) {
                        sub1.idx.push_back (i);
                }
                // a matrix with fewer rows than columns is always singular
                full_possible = need_full && N >= m;
                sub1_possible = need_1 && N >= 1 + k;

                rowbuffer.resize (4 * m);
                for (int i = 0; i < 2; i++) {
                        xold[i] = &rowbuffer[i * m];
                        xnew[i] = &rowbuffer[(2 + i) * m];
                }
                reset (A);
        }

        /// calculate the state directly from the design
        void reset (const array_link &A) {
                Eigen::MatrixXd X = array2eigenModelMatrix (A);
                fallback = false;
                if (full_possible) {
                        full.init (X);
                        fallback = fallback || full.singular;
                        if (need_02) {
                                sub02.init (X);
                                fallback = fallback || sub02.singular;
                        }
                }
                if (sub1_possible) {
                        sub1.init (X);
                        fallback = fallback || sub1.singular;
                }
                nresets++;
        }

        /// record the rows of the design that are about to change
        void begin (const array_link &A, const int *rows, int p) {
                this->p = p;
                for (int i = 0; i < p; i++) {
                        this->rows[i] = rows[i];
                        modelrow_2level (A, rows[i], k, xold[i]);
                }
        }

        /// calculate the efficiencies of the changed design
        std::vector< double > evaluate (const array_link &A) {
                pending_fallback = fallback;
                if (!pending_fallback) {
                        for (int i = 0; i < p; i++) {
                                modelrow_2level (A, rows[i], k, xnew[i]);
                        }
                        if (full_possible) {
                                full.candidate (p, xold, xnew);
                                pending_fallback = full.pending_unreliable;
                                if (!full.pending_singular && need_02) {
                                        sub02.candidate (p, xold, xnew);
                                        pending_fallback = pending_fallback || sub02.pending_unreliable;
                                }
                        }
                        if (sub1_possible) {
                                sub1.candidate (p, xold, xnew);
                                pending_fallback = pending_fallback || sub1.pending_unreliable;
                        }
                }
                if (pending_fallback) {
                        nfallback++;
                        return ::Defficiencies (A, arrayclass);
                }

                // efficiencies that are not needed for the score are reported as zero
                std::vector< double > efficiencies (3, 0.);
                if (full_possible && !full.pending_singular) {
                        efficiencies[0] = exp (full.pending_logdet / m) / N;
                        if (need_02) {
                                efficiencies[1] = exp ((full.pending_logdet - sub02.pending_logdet) / k) / N;
                        }
                }
                if (sub1_possible && !sub1.pending_singular) {
                        efficiencies[2] = exp (sub1.pending_logdet / (k + 1)) / N;
                }
                return efficiencies;
        }

        /// accept the pending change of the design
        void accept (const array_link &A) {
                if (pending_fallback) {
                        reset (A);
                        return;
                }
                bool refresh = false;
                if (full_possible) {
                        full.accept ();
                        refresh = refresh || full.singular || full.nupdates >= refresh_interval;
                        if (need_02) {
                                sub02.accept ();
                                refresh = refresh || sub02.singular;
                        }
                }
                if (sub1_possible) {
                        sub1.accept ();
                        refresh = refresh || sub1.singular || sub1.nupdates >= refresh_interval;
                }
                if (refresh) {
                        reset (A);
                }
        }

        /// number of direct calculations performed
        int nresets = 0;
        /// number of candidates evaluated with the direct calculation
        int nfallback = 0;

      private:
        const arraydata_t &arrayclass;
        const int N, k, n2fi, m;
        const int refresh_interval;
        bool need_full = true, need_02 = true, need_1 = true;
        bool full_possible = false;
        bool sub1_possible = false;
        /// true if the incremental state is not valid and the direct calculation is used
        bool fallback = true;
        bool pending_fallback = true;

        DeffSubsystem full, sub02, sub1;

        int p = 0;
        int rows[2] = {0, 0};
        std::vector< double > rowbuffer;
        double *xold[2];
        double *xnew[2];
};

/// restore the array to the state before the coordinate exchange step
inline void undo_exchange (array_link &A, coordinate_exchange_method_t optimmethod, int r, int c, int r2, int c2,
                           array_t o, array_t o2) {
        switch (optimmethod) {
        case DOPTIM_SWAP:
                A._setvalue (r, c, o);
                A._setvalue (r2, c2, o2);
                break;
        case DOPTIM_FLIP:
        case DOPTIM_UPDATE:
                A._setvalue (r, c, o);
                break;
        case DOPTIM_NONE:
                break;
        case DOPTIM_AUTOMATIC:
                throw_runtime_exception ("coordinate_exchange_method_t needs to be set");
        }
}

} // namespace

double check_incremental_Defficiencies (const arraydata_t &arrayclass, int nsteps, int verbose) {
        if (!arrayclass.is2level ()) {
                throw_runtime_exception ("check_incremental_Defficiencies: only implemented for 2-level designs");
        }
        const int N = arrayclass.N;
        const int k = arrayclass.ncols;
        array_link A = arrayclass.randomarray (0);
        std::vector< double > alpha (3, 1.);
        Deff2levelUpdater updater (arrayclass, A, alpha, 1000000); // large refresh interval to test accumulation of errors

        double maxdiff = 0;
        for (int ii = 0; ii < nsteps; ii++) {
                int r = fastrandK (N);
                int c = fastrandK (k);
                int r2 = fastrandK (N);
                int c2 = fastrandK (k);
                array_t o = A._at (r, c);
                array_t o2 = A._at (r2, c2);
                bool swap = (ii % 2 == 1) && (o != o2);
                int rows[2] = {r, r2};
                updater.begin (A, rows, (swap && r != r2) ? 2 : 1);
                if (swap) {
                        A._setvalue (r, c, o2);
                        A._setvalue (r2, c2, o);
                } else {
                        A._setvalue (r, c, 1 - o);
                }
                std::vector< double > dd = updater.evaluate (A);
                std::vector< double > dd_direct = Defficiencies (A, arrayclass);
                for (int j = 0; j < 3; j++) {
                        maxdiff = std::max (maxdiff, fabs (dd[j] - dd_direct[j]));
                }
                if (verbose >= 2) {
                        myprintf ("check_incremental_Defficiencies: step %d: incremental %.8f %.8f %.8f direct %.8f %.8f %.8f\n", ii, dd[0],
                                  dd[1], dd[2], dd_direct[0], dd_direct[1], dd_direct[2]);
                }
                // accept roughly 2 out of 3 changes
                if (ii % 3 != 0) {
                        updater.accept (A);
                } else {
                        if (swap) {
                                A._setvalue (r, c, o);
                                A._setvalue (r2, c2, o2);
                        } else {
                                A._setvalue (r, c, o);
                        }
                }
        }
        if (verbose) {
                myprintf ("check_incremental_Defficiencies: N %d, k %d, %d steps: max deviation %e (%d direct calculations)\n", N, k, nsteps,
                          maxdiff, updater.nresets + updater.nfallback);
        }
        return maxdiff;
}

array_link optimDeff (const array_link &A0, const arraydata_t &arrayclass, const std::vector< double > alpha,
                      int verbose, coordinate_exchange_method_t optimmethod, int niter, int nabort) {
        const int N = arrayclass.N;
        const int k = arrayclass.ncols;

        if (nabort <= 0) {
			/// select a sane default
                if (arrayclass.is2level ()) {
                        nabort = (int)((N * k) * 2.5 + 1); // factor 2 to compensate for equal error switches
                } else {
                        nabort = (N * k) * (arrayclass.s[0]) + 2;
                }
        }

        std::vector< int > factor_levels = arrayclass.factor_levels ();
        if (optimmethod == DOPTIM_UPDATE) {
                if (arrayclass.is2level ())
                        optimmethod = DOPTIM_FLIP;
        }
        array_link A = A0;
        symmetry_group sg = symmetry_group (factor_levels);

        std::vector< int > gidx = sg.gidx;

        /// for 2-level designs the efficiencies are updated incrementally
        const bool incremental = arrayclass.is2level () && A.n_rows == N && A.n_columns == k && A.max () <= 1 && A.min () >= 0;
        std::unique_ptr< Deff2levelUpdater > updater;
        if (incremental) {
                updater.reset (new Deff2levelUpdater (arrayclass, A, alpha));
        }

        std::vector< double > dd0 = A0.Defficiencies ();

        if (verbose) {
                myprintf ("optimDeff: initial D-efficiency %.4f\n", dd0[0]);
        }

        // initialize score
        double current_score = scoreD (dd0, alpha);

        /// differences in score below this tolerance are considered to be rounding errors
        const double score_tolerance = 1e-10;

        int last_update = 0; // index of last change to array

        // initialize arary with random permutation
        const int nn = N * k;
        std::vector< int > updatepos = permutation< int > (nn);
        my_random_shuffle (updatepos.begin (), updatepos.end ());
        int updateidx = 0;

        for (int ii = 0; ii < niter; ii++) {
                // select random row and column
                int r = updatepos[updateidx] % N;
                int c = updatepos[updateidx] / N;
                updateidx = (updateidx + 1) % (nn);

                int r2 = fastrandK (N);

                // make sure column is proper column group
                int c2 = sg.gstart[sg.gidx[c]] + fastrandK (sg.gsize[gidx[c]]);

                // get values
                array_t o = A._at (r, c);
                array_t o2 = A._at (r2, c2);

                // update

                if (optimmethod == DOPTIM_SWAP && o == o2)
                        continue;

                if (incremental) {
                        int rows[2] = {r, r2};
                        updater->begin (A, rows, (optimmethod == DOPTIM_SWAP && r2 != r) ? 2 : 1);
                }

                switch (optimmethod) {
                case DOPTIM_SWAP: // swap
                        A._setvalue (r, c, o2);
                        A._setvalue (r2, c2, o);
                        break;
                case DOPTIM_UPDATE: { // random update
                        int val = fastrandK (factor_levels[c]);
                        A._setvalue (r, c, val);
                        break;
                }
                case DOPTIM_FLIP: // flip
                        A._setvalue (r, c, 1 - o);
                        break;
                case DOPTIM_NONE:
                        break;
                default:
                        myprintf ("optimDeff: no such optimization method\n");
                        break;
                }
                // evaluate
                std::vector< double > dd = incremental ? updater->evaluate (A) : A.Defficiencies ();
                double next_score = scoreD (dd, alpha);

                if (verbose >= 4) {
                        myprintf ("  optimDeff: switch %d, %d: %d with %d, %d: %d: d %.4f -> dn %.4f \n", r, c, r2, c2,
                                  o, o2, current_score, next_score);
                }

                // switch back if necessary

                if (next_score >= current_score - score_tolerance) {
                        if (incremental) {
                                updater->accept (A);
                        }
                        if (next_score > current_score + score_tolerance) {
                                last_update = ii;

                                if (verbose >= 3)
                                        myprintf ("optimDeff: ii %d: %.6f -> %.6f\n", ii, current_score, next_score);
                        } else {
                                if (verbose >= 2) {
                                        myprintf ("optimDeff: equal ii %d: %.6f -> %.6f\n", ii, current_score, next_score);
                                }
                        }
                        current_score = next_score;
                } else {
                        // restore to original
                        undo_exchange (A, optimmethod, r, c, r2, c2, o, o2);
                }

                // check progress
                if ((ii - last_update) > nabort) {
                        if (verbose >= 2)
                                myprintf ("optimDeff: early abort ii %d, lc %d: %d\n", ii, last_update, (ii - last_update));
                        break;
                }
        }

        std::vector< double > dd = A.Defficiencies ();
        double dn = scoreD (dd, alpha);

        if (verbose) {
                myprintf ("optimDeff: final score %.4f, final D-efficiency %.4f\n", dn, dd[0]);
        }

        return A;
}

array_link optimDeff2level (const array_link &A0, const arraydata_t &arrayclass, std::vector< double > alpha,
                            int verbose, int optimmethod, int niter, int nabort) {
        int N = arrayclass.N;
        int k = arrayclass.ncols;
        std::vector< int > s = arrayclass.factor_levels ();

        if (!arrayclass.is2level ()) {
                myprintf ("optimDeff2level: error arrayclass is not 2-level\n");
                return A0;
        }

        if (nabort <= 0) {
                nabort = N * k + 2;
        }

        if (optimmethod == DOPTIM_UPDATE) {
                if (arrayclass.is2level ())
                        optimmethod = DOPTIM_FLIP;
        }
        array_link A = A0;

        std::vector< double > dd0 = A0.Defficiencies ();

        if (verbose) {
                myprintf ("optimDeff: initial D-efficiency %.4f\n", dd0[0]);
        }

        // initialize arary with random permutation
        std::vector< int > updatepos = permutation< int > (N * k);
        my_random_shuffle (updatepos.begin (), updatepos.end ());
        int updateidx = 0;

        // initialize score
        double current_score = scoreD (dd0, alpha);

        int last_update = 0; // index of last change to array

        int nn = updatepos.size ();

        for (int ii = 0; ii < niter; ii++) {
                // select random row and column
                int r = fastrandK (N);
                int c = fastrandK (k);
                updateidx = (updateidx + 1) % (nn);

                int r2 = fastrandK (N);
                int c2 = fastrandK (k);

                // get values
                array_t o = A._at (r, c);
                array_t o2 = A._at (r2, c2); // no extra error checking

                // update
                if (optimmethod == DOPTIM_SWAP && o == o2)
                        continue;

                switch (optimmethod) {
                case DOPTIM_SWAP: // swap
                        A._setvalue (r, c, o2);
                        A._setvalue (r2, c2, o);
                        break;
                case DOPTIM_UPDATE: { // random update
                        int val = fastrandK (s[c]);
                        A._setvalue (r, c, val);
                        break;
                }
                case DOPTIM_FLIP: // flip
                        A._setvalue (r, c, 1 - o);
                        break;
                case DOPTIM_NONE:
                        break;
                default:
                        myprintf ("optimDeff: no such optimization method\n");
                        break;
                }
                // evaluate
                std::vector< double > dd = A.Defficiencies ();
                double next_score = scoreD (dd, alpha);

                if (verbose >= 4) {
                        myprintf ("  optimDeff: switch %d, %d: %d with %d, %d: %d: d %.4f -> dn %.4f \n", r, c, r2, c2,
                                  o, o2, current_score, next_score);
                }

                // switch back if necessary
                if (next_score >= current_score) {
                        if (next_score > current_score) {
                                last_update = ii;
                                if (verbose >= 3)
                                        myprintf ("optimDeff: ii %d: %.6f -> %.6f\n", ii, current_score, next_score);
                                current_score = next_score;
                        } else {
                                myprintf ("optimDeff: equal ii %d: %.6f -> %.6f\n", ii, current_score, next_score);
                        }
                } else {
                        // restore to original
                        switch (optimmethod) {
                        case DOPTIM_SWAP:
                                A._setvalue (r, c, o);
                                A._setvalue (r2, c2, o2);
                                break;
                        case DOPTIM_FLIP:
                        case DOPTIM_UPDATE:
                                A._setvalue (r, c, o);
                                break;
                        case DOPTIM_NONE:
                                break;
                        }
                }

                // check progress
                if ((ii - last_update) > nabort) {
                        if (verbose >= 2)
                                myprintf ("optimDeff: early abort ii %d, lc %d: %d\n", ii, last_update, (ii - last_update));
                        break;
                }
        }

        std::vector< double > dd = A.Defficiencies ();
        double dn = scoreD (dd, alpha);

        if (verbose) {
                myprintf ("optimDeff: final score %.4f, final D-efficiency %.4f\n", dn, dd[0]);
        }

        return A;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
