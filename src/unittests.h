/** \file unittests.h

 \brief Contains unit tests
 */

#pragma once

int unittest_reduceConferenceTransformation();

int unittest_nautynormalform(const array_link &al, int verbose);

int checkTransformationComposition(const array_link &al, int verbose = 0);

void test_array_manipulation(int verbose = 1);

int checkConferenceComposition(const array_link &al, int verbose = 0);

void test_conference_candidate_generators(int verbose = 1);

int checkTransformationInverse(const array_link &al);

int checkConferenceInverse(const array_link &array);

bool testLMC0checkDC(const array_link &al, int verbose = 1);
