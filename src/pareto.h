/** \file pareto.h

\brief Class for calculation the Pareto optimal elements from a set of multi-valued objects


Copyright (c) 2013, Pieter Eendebak <pieter.eendebak@gmail.com>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following
disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#pragma once



#include <stdio.h>
#ifdef SWIGCODE
#include "printfheader.h"
#else
#define myprintf printf
#endif

#include <deque>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdlib.h>
#include <vector>

namespace detail {
template < class atype >
/// generic function to print a std::vector
void display_vector (const std::vector< atype > &v, const char *sep = " ") {
        std::stringstream buffer;
        std::copy (v.begin (), v.end (), std::ostream_iterator< atype > (buffer, sep));
        myprintf ("%s", buffer.str ().c_str ());
}
template <>
/// specialized function to print a std::vector
inline void display_vector (const std::vector< int > &v, const char *sep) {
        for (size_t i = 0; i < v.size (); i++) {
                myprintf ("%d", v[i]);
                if (i < v.size () - 1) {
                        myprintf ("%s", sep);
                }
        }
}
template <>
/// specialized function to print a std::vector
inline void display_vector (const std::vector< long > &v, const char *sep) {
        // myprintf("long case");
        for (size_t i = 0; i < v.size (); i++) {
                myprintf ("%ld", v[i]);
                if (i < v.size () - 1) {
                        myprintf ("%s", sep);
                }
        }
}

template <>
/// specialized function to print a std::vector
inline void display_vector (const std::vector< double > &v, const char *sep) {
        // myprintf("int case");
        for (size_t i = 0; i < v.size (); i++) {
                myprintf ("%f", v[i]);
                if (i < v.size () - 1) {
                        myprintf ("%s", sep);
                }
        }
}
} // end of namespace

template < class ValueType, class IndexType >
/// helper class for the Pareto class to hold elements
struct pareto_element {

        typedef std::vector< ValueType > pValue;

        pValue value;
        std::vector< IndexType > indices;

        /// return true of the argument element dominates this value
        bool dominates (pValue v) {
                for (size_t i = 0; i < v.size (); i++) {
                        if (value[i] < v[i]) {
                                return false;
                        }
                }
                return true;
        }
		/// return true of the argument element is dominated by this value
		bool isdominated (pValue v) {
                for (size_t i = 0; i < v.size (); i++) {
                        if (value[i] > v[i]) {
                                return false;
                        }
                }
                return true;
        }
        /// return true of the argument element is equal to this element
        bool equal (pValue v) {
                for (size_t i = 0; i < v.size (); i++) {
                        if (value[i] != v[i]) {
                                return false;
                        }
                }
                return true;
        }
};

template < class ValueType, class IndexType >
/** @brief Class to the calculate Pareto optimal elements.
 *
 * The class is templated by the type of values to be compared and an index type. The index type is used to index the
 * elements.
 *
 * For elements added to the Pareto structure larger is better.
 */
class Pareto {
      public:
		/// type for values of Pareto elements
        typedef std::vector< ValueType > pValue;

		/// a pareto element consists of a pair (value, index)
        typedef pareto_element< ValueType, IndexType > pElement;

		/// Verbosity level
        int verbose;
        /// contains a list of all Pareto optimal elements
        std::deque< pareto_element< ValueType, IndexType > > elements;

        /// Create an empty Pareto class
        Pareto () : verbose (1){};
        ~Pareto (){};

        /// return the total number of Pareto optimal values
        int number () const { return elements.size (); }
        /// return the total number Pareto optimal objects
        int numberindices () const {
                int t = 0;
                for (size_t i = 0; i < elements.size (); i++) {
                        t += elements[i].indices.size ();
                }
                return t;
        }

        std::string __repr__ () const {
                std::string ss =
                    printfstring ("Pareto: %zu optimal values, %zu elements\n", elements.size (), numberindices ());
                return ss;
        }

		/// show a Pareto element 
		static void showvalue (const pValue p) { detail::display_vector (p, "; "); }
        /// show the current set of Pareto optimal elements
        void show (int verbose = 1) {
                if (verbose == 0) {
                        return;
                }
                myprintf ("Pareto: %ld optimal values, %d objects\n", (long)elements.size (), numberindices ());
                if (verbose >= 2) {
                        for (size_t i = 0; i < elements.size (); i++) {
                                myprintf ("value %d: ", (int)i);
                                fflush (stdout);
                                detail::display_vector (elements[i].value, "; ");
                                myprintf ("\n");
                                if (verbose >= 3) {
                                        myprintf ("  indices: ");
                                        detail::display_vector (elements[i].indices, ", ");
                                        myprintf ("\n");
                                }
                        }
                }
        }

        /// return all indices of the Pareto optimal elements as a std::deque
        std::deque< IndexType > allindicesdeque () const {
                std::deque< IndexType > lst;
                for (size_t i = 0; i < elements.size (); i++) {
                        lst.insert (lst.end (), elements[i].indices.begin (), elements[i].indices.end ());
                }
                return lst;
        }

        /// return all indices of the Pareto optimal elements
        std::vector< IndexType > allindices () const {
                std::vector< IndexType > lst;
                for (size_t i = 0; i < elements.size (); i++) {
                        lst.insert (lst.end (), elements[i].indices.begin (), elements[i].indices.end ());
                }
                return lst;
        }

        /// return the values of all Pareto optimal elements
        std::vector< pValue > allvalues () const {
                std::vector< pValue > lst;
                for (size_t i = 0; i < this->elements.size (); i++) {
                        lst.push_back (this->elements[i].value);
                }
                return lst;
        }

        /// add a new element
        bool addvalue (const pValue value, const IndexType idx) {
                size_t ii = 0;
                while (ii < elements.size ()) {
                        if (verbose >= 4) {
                                myprintf ("Pareto::addvalue: compare new element to element %d\n", (int)ii);
                        }
                        if (elements[ii].dominates (value)) {
                                if (elements[ii].equal (value)) {
                                        elements[ii].indices.push_back (idx);
                                        if (verbose >= 3) {
                                                myprintf ("Pareto::addvalue: new pareto item (same value)\n");
                                        }
                                        return true;
                                } else {
                                        // not a pareto element, so continue
                                        if (verbose >= 3) {
                                                myprintf ("Pareto::addvalue: not pareto\n");
                                                myprintf (" new elememnt : ");
                                                this->showvalue (value);
                                                myprintf ("\n");
                                                myprintf ("  dominated by: ");
                                                this->showvalue (elements[ii].value);
                                                myprintf ("\n");
                                        }
                                        return false;
                                }
                        }
                        if (elements[ii].isdominated (value)) {
                                // element ii is dominated by the new element, we remove element ii
                                if (verbose >= 2) {
                                        printf ("Pareto::addvalue: removing element\n");
                                        if (verbose >= 3) {
                                                myprintf ("  new element : ");
                                                this->showvalue (value);
                                                myprintf ("\n");
                                                myprintf ("removing %3d : ", (int)ii);
                                                this->showvalue (elements[ii].value);
                                                myprintf ("\n");
                                        }
                                }
                                elements.erase (elements.begin () + ii);
                        } else {
                                if (verbose >= 4) {
                                        myprintf ("Pareto:addvalue: ii++\n");
                                }
                                ii++;
                        }
                }

                // we have a new pareto element: add it to the current set
                pElement pareto_element;
                pareto_element.value = value;
                pareto_element.indices.push_back (idx);
                this->elements.push_back (pareto_element);
                if (verbose >= 2) {
                        myprintf ("Pareto: addvalue: new pareto item (new value), total is %ld\n",
                                  (long)this->elements.size ());
                }
                return true;
        }
};

