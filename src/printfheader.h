#pragma once

#ifdef SWIGCODE

#include <Python.h>
inline void pyprintf (const char *message, ...) {
        char buf[32 * 1024];

        va_list va;
        va_start (va, message);
        vsprintf (buf, message, va);
        va_end (va);

        PyObject *f = PySys_GetObject ((char *)"stdout");
        if (f == 0) {
                printf ("error: could not get Python stdout object\n");
                return;
        }
        PyFile_WriteString (buf, f);

        return;
}
#define myprintf pyprintf

#else
// the full package is build (no restrictions from R or Python)
#define FULLPACKAGE 1
#include <stdio.h>
#define myprintf printf
#endif // SWIGCODE

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
