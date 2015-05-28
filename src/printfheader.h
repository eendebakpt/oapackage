#pragma once

#ifdef RPACKAGE
//#include <stdio.h>
#include "R_ext/Print.h"
#define myprintf Rprintf
//#define FULLPACKAGEX 0
#else
#define FULLPACKAGE 1
//#define FULLPACKAGEX 1
#include <stdio.h>
#define myprintf printf
#endif

