/** \file pareto.h

\brief Class for calculation the Pareto optimal elements from a set of multi-valued objects


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

#pragma once

#ifndef PARETO_H
#define PARETO_H

#ifdef RPACKAGE
// not implemented...
//#define myPrintf Rprintf
#else
#include <stdio.h>
//#define myPrintf printf
#endif

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <deque>
#include <iterator>

namespace detail
{
template <class atype>
/// generic function to print a std::vector
void display_vector ( const std::vector<atype> &v, const char *sep = " " )
{
    std::copy ( v.begin(), v.end(), std::ostream_iterator<atype> ( std::cout, sep ) );
}
}



template <class ValueType, class IndexType>
/// helper class for the Pareto class
struct pareto_element {

    typedef std::vector<ValueType> pValue;

    pValue value;
    std::vector<IndexType> indices;

    /// return true of the argument element dominates this value
    bool dominates ( pValue v )
    {
        for ( size_t i=0; i<v.size(); i++ ) {
            if ( value[i] < v[i] )
                return false;
        }
        return true;
    }
    bool isdominated ( pValue v )
    {
        for ( size_t i=0; i<v.size(); i++ ) {
            if ( value[i] > v[i] )
                return false;
        }
        return true;
    }
    /// return true of the argument element is equal to this element
    bool equal ( pValue v )
    {
        for ( size_t i=0; i<v.size(); i++ ) {
            if ( value[i] != v[i] )
                return false;
        }
        return true;
    }

};

template <class ValueType, class IndexType>
/** @brief Class to the calculate Pareto optimal elements.
 *
 * The class is templated by the type of values to be compared and an index type. The index type is used to index the elements.
 *
 */
class Pareto
{
public:
    typedef std::vector<ValueType> pValue;
    typedef pareto_element<ValueType, IndexType> pElement;

    int verbose;
    /// contains a list of all Pareto optimal elements
    std::deque<pareto_element<ValueType, IndexType> > elements;

    /// constructor
    Pareto() : verbose ( 1 ) {};
    ~Pareto() {};

    /// return the total number of Pareto optimal values
    int number() const
    {
        return elements.size();
    }
    /// return the toal number Pareto optimal objects
    int numberindices() const
    {
        int t=0;
        for ( size_t i=0; i<elements.size(); i++ ) {
            t+= elements[i].indices.size() ;
        }
        return t;
    }

    std::string __repr__() const
    {
        std::string ss = printfstring ( "Pareto: %zu optimal values, %zu elements\n", elements.size(), numberindices() );
        return ss;
    }

    static void showvalue(const pValue p)
    {
        detail::display_vector ( p, "; " );

    }
    /// show the current set of Pareto optimal elements
    void show ( int verbose=1 )
    {
        if ( verbose==0 )
            return;
        printf ( "Pareto: %ld optimal values, %d objects\n", (long)elements.size(),  numberindices()  );
        if ( verbose>=2 ) {
            for ( size_t i=0; i<elements.size(); i++ ) {
                printf ( "value %d: ", (int)i );
                detail::display_vector ( elements[i].value, "; " );
                printf ( "\n" );
                if ( verbose>=3 ) {
                    printf ( "  indices: " );
                    detail::display_vector ( elements[i].indices, ", " );
                    printf ( "\n" );
                }
            }
        }
    }

    /// return all indices of the Pareto optimal elements as a std::deque
    std::deque<IndexType> allindicesdeque() const
    {
        std::deque<IndexType> lst;
        for ( size_t i=0; i<elements.size(); i++ ) {
            lst.insert ( lst.end(), elements[i].indices.begin(), elements[i].indices.end() );
        }
        return lst;
    }

    /// return all indices of the Pareto optimal elements
    std::vector<IndexType> allindices() const
    {
        std::vector<IndexType> lst;
        for ( size_t i=0; i<elements.size(); i++ ) {
            lst.insert ( lst.end(), elements[i].indices.begin(), elements[i].indices.end() );
        }
        return lst;
    }

    /// return all Paretop optimal elements
    std::vector<pValue> allvalues() const
    {
        std::vector<pValue> lst;
        for ( size_t i=0; i<this->elements.size(); i++ ) {
            lst.push_back(this->elements[i].value);
        }
        return lst;
    }

    /// add a new element
    bool addvalue ( const pValue val, const IndexType idx )
    {
        size_t ii=0;
        while ( ii<elements.size() ) {
            if ( verbose>=4 )
                printf ( "Pareto::addvalue: compare new element to element %d\n", (int) ii );
            if ( elements[ii].dominates ( val ) ) {
                if ( elements[ii].equal ( val ) ) {
                    elements[ii].indices.push_back ( idx );
                    if ( verbose>=3 )
                        printf ( "Pareto::addvalue: new pareto item (same value)\n" );
                    return true;
                } else {
                    // not a pareto element, so continue
                    if ( verbose>=3 ) {
                        printf ( "Pareto::addvalue: not pareto\n" );
			printf(" new elememnt : "); this->showvalue(val); printf("\n");
			printf("  dominated by: "); this->showvalue(elements[ii].value); printf("\n");
		    }
                    return false;
                }
            }
            if ( elements[ii].isdominated ( val ) ) {
                // element ii is dominated by the new element, we remove element ii
                if ( verbose>=2 ) {
                    printf ( "Pareto::addvalue: removing element\n" );
	      if ( verbose>=3 ) {
			printf("  new element : "); this->showvalue(val); printf("\n");
			printf("removing %3d : ", (int)ii); this->showvalue(elements[ii].value); printf("\n");
		    }

	    }
                elements.erase ( elements.begin() +ii );
            } else {
                if ( verbose>=4 )
                    printf ( "Pareto:addvalue: ii++\n" );
                ii++;
            }
        }

        // we have a new pareto element: add it to the current set
        pElement p;
        p.value = val;
        p.indices.push_back ( idx );
        this->elements.push_back ( p );
        if ( verbose>=2 )
            printf ( "Pareto: addvalue: new pareto item (new value), total is %ld\n", (long)this->elements.size() );
        return true;
    }
};


#endif
