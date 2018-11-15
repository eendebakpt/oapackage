#include <algorithm>
#include <math.h>

#include "arraytools.h"
#include "mathtools.h"
#include "tools.h"
#include "lmc.h"

#include "conference.h"



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
