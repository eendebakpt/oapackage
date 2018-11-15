#include <algorithm>
#include <math.h>

#include "arraytools.h"
#include "mathtools.h"
#include "tools.h"
#include "lmc.h"

#include "conference.h"

/// check transformation inverse. return 0 if test is good
int checkConferenceInverse(const array_link &al) {
	conference_transformation_t T1(al);
	T1.randomize();

	conference_transformation_t T1i = T1.inverse();
	conference_transformation_t II = T1i * T1;

	myassert(II.isIdentity(), "unittest error: inverse of conference matrix transformation\n");
	myassert(al == II.apply(al), "unittest error: inverse of conference matrix transformation\n");

	return 0;
}

/// check composition operator. returns 0 if test id good
int checkConferenceComposition(const array_link &al, int verbose) {
	conference_transformation_t T1(al);
	T1.randomize();

	conference_transformation_t T2(al);
	T2.randomize();

	conference_transformation_t T3 = T2 * T1;

	array_link al1 = T1.apply(al);
	array_link al1t2 = T2.apply(al1);
	array_link al3 = T3.apply(al);

	if (verbose) {
		printfd("checkTransformationComposition: transforms\n");
		T1.show();
		T2.show();
		T3.show();
		printfd("checkTransformationComposition: arrays\n");
		al.showarray();
		al1.showarray();
		al1t2.showarray();
		al3.showarray();
	}

	myassert(al3 == al1t2, "unittest error: composition of conference transformations\n");

	return 0;
}

void test_conference_candidate_generators(int verbose = 1) {
	int N = 20;
	conference_t ct(N, 4, 1);
	ct.ctype = conference_t::DCONFERENCE;
	ct.j3zero = 1;

	array_link al = exampleArray(35, 1);

	CandidateGeneratorDouble cgenerator(array_link(), ct);
	cgenerator.verbose = 0;
	for (int i = 0; i < 2; i++) {
		{
			const std::vector< conference_column > &cl = cgenerator.generateCandidates(al);
			myassert(cl.size() == 3, "unittest error: inverse of array transformation\n");

			if (verbose >= 2) {
				printfd("generated %d\n", cl.size());
				cgenerator.showCandidates(2);
			}
		}
	}

	// generator for conference matrices
	{
		const int example_idx = 39;
		array_link al = exampleArray(example_idx, 1);
		conference_t ct(al.n_rows, al.n_columns + 4, 0);

		if (verbose >= 2) {
			myprintf("checking generator on array:\n");
			al.showarray();
		}
		myassert(al.is_conference());
		myassert(al.min() == -1);

		int filterj2 = 1;
		int filtersymminline = 1;
		int averbose = verbose;
		std::vector< conference_column > ccX = generateSingleConferenceExtensions(al, ct, -1, averbose, 1, filterj2,
			ct.j3zero, filtersymminline);
		if (verbose >= 2) {

			showCandidates(ccX);
			printf("\n-----------\n");
		}
		myassert(ccX.size() == 2, "number of candidnates generated");
		{
			CandidateGenerator cgenerator(array_link(), ct);
			int kz = maxz(al) + 1;
			cgenerator.verbose = verbose;
			std::vector< conference_column > ee = cgenerator.generateCandidatesZero(al, kz);
			printf("ee.size() %d\n", (int)ee.size());
			myassert(ee.size() == 1, "number of candidnates generated");
			if (verbose >= 2) {
				cgenerator.showCandidates(2);
				printf("generateCandidatesZero: %d\n-------------\n", (int)ee.size());
			}
		}
	}
}

/// check transformation inverse. return 0 if test is good
int checkTransformationInverse(const array_link &al) {
	arraydata_t adataX = arraylink2arraydata(al);
	array_transformation_t T1(&adataX);
	T1.randomize();

	array_transformation_t T1i = T1.inverse();
	array_transformation_t II = T1i * T1;

	myassert(II.isIdentity(), "unittest error: inverse of array transformation\n");

	return 0;
}

bool testLMC0checkDC(const array_link &al, int verbose = 1) {
	const int niter = 20; // number of iterations used in check

	myassert(al.is_conference(2));

	// perform LMC0 test
	lmc_t r = LMC0checkDC(al, verbose >= 2);

	if (verbose >= 2) {
		printf("testLMC0checkDC: result %d (LMC_LESS %d, LMC_MORE %d)\n", (int)r, LMC_LESS, LMC_MORE);
	}
	if (r == LMC_MORE)
		return true;

	// array is in LMC0 format, perform random transformations
	for (int jj = 0; jj < niter; jj++) {
		conference_transformation_t tr(al);
		tr.randomizecolperm();
		tr.randomizerowperm();
		tr.randomizecolflips();

		array_link alx = tr.apply(al);

		if (al != alx) {
			lmc_t r = LMC0checkDC(alx, verbose >= 2);

			if (verbose >= 2) {
				printf("testLMC0checkDC: randomized array: result %d (LMC_LESS %d, LMC_MORE %d)\n",
					(int)r, LMC_LESS, LMC_MORE);
			}

			if (r != LMC_LESS) {
				printf("testLMC0checkDC: error?: LMC0checkDC on randomized array did not return "
					"LMC_LESS\n");
				return false;
			}
		}
		else {
			if (verbose >= 2)
				printf("testLMC0checkDC: randomized array resulted in same array\n");
		}
	}
	return true;
}
