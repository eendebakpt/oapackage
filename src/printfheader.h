#pragma once

#ifdef RPACKAGE
// do not include stdio or iostream for an R package
#include "R_ext/Print.h"
#define myprintf Rprintf
#else
#define FULLPACKAGE 1
#include <stdio.h>
#define myprintf printf
#endif

