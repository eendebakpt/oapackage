#pragma once

#ifdef RPACKAGE
  // do not include stdio or iostream for an R package
  #include "R_ext/Print.h"
  #define myprintf Rprintf
#else
  #ifdef SWIGCODE

#include <Python.h>
inline void pyprintf ( const char *message, ... )
{
    char buf[32 * 1024];

    va_list va;
    va_start ( va, message );
    vsprintf ( buf, message, va );
    //vsnprintf;
    
    va_end ( va );

  PyObject *f = PySys_GetObject((char *)"stdout");
  if (f==0) {
    printf("error: could not get Python stdout object\n");
   return;
  }
  PyFile_WriteString(buf, f);
  
//  printf("pyprintf called with |%s|\n", buf);
    return ;
}

/*
inline void pyprintf(const std::string msg) {
  PyObject *f = PySys_GetObject((char *)"stdout");
  PyFile_WriteString(msg.c_str(), f);
}

inline void pyprintf(const char *msg) {
  PyObject *f = PySys_GetObject((char *)"stdout");
  PyFile_WriteString(msg, f);
}
*/
    #define myprintf pyprintf

#else
    #define FULLPACKAGE 1
    #include <stdio.h>
    #define myprintf printf
  #endif
#endif

