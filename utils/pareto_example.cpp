/** \file pareto.cpp

Class for calculation the Pareto optimal elements from a set of multi-valued objects


Copyright (c) 2013, Pieter Eendebak <pieter.eendebak@gmail.com>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <stdio.h>
#include <stdlib.h>
#include "pareto.h"


int main ( int argc, char **argv )
{
	printf("Example program for calculating Pareto optimnal elements\n");
	printf("    author: <pieter.eendebak@gmail.com>\n");
	
   /* initialize random seed: */
#ifdef WIN32
   srand ( 0 );
#else
	srand(time(NULL));
#endif

   Pareto<long, int> pareto;

   // select 12 random objects with each 4 values
   const int nelem=12;
   const int nvals=3;

   const int nmax=4;
   std::vector<long> x ( nvals );
   printf ( "adding 12 random elements\n" );
   for ( int i=0; i<nelem; i++ ) {
      for ( int j=0; j<nvals; j++ )
         x[j]=rand() % nmax;

      pareto.addvalue ( x, i );
      printf ( "  element %d: ", i );
      detail::display_vector ( x );
      printf ( "\n" );
   }

   printf ( "The pareto optimal values are:\n" );
   pareto.show ( 3 );
}
