/** \file mathtools.h

 \brief Contains declarations of mathematical tools (mostly inline and template)

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include "printfheader.h"
#ifdef FULLPACKAGE
#include <iostream>
#endif
#include <limits>
#include <algorithm> /* defines max and min template functions */
#include <assert.h>
#include <cstdlib>
#include <vector>


#include <string.h>
#include <deque>
#include <cmath>
#include <sstream>      // std::stringstream

#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <numeric>
#include <iterator>


#if defined(_MSC_VER)
#pragma warning(disable: 4800)
#pragma warning(disable: 4996)
#pragma warning (disable : 4068 )
#pragma warning (disable : 4267 )

// add round function, since visual studio does not define it...
inline double round ( double x )
{
    return floor ( x + 0.5 );
}
#endif

#ifdef SWIG
//%ignore    mvalue_t::operator<<;
//%rename(Complex_add_dc) operator+(double, const Complex &);
//%ignore    mvalue_t::operator<<( std::ostream& stream, const mvalue_t<W>& mval );
#endif

#include <queue>

// from: http://codereview.stackexchange.com/questions/13979/simple-object-pool-template-container-in-c
template <class TYPE>
/// Class to make a pool of objects that can be re-used
class object_pool
{
    std::vector< TYPE *> pool;
    //std::queue<size_t> avail;
    int maxpoolsize; /// maximum size of the pool
    int nn;		/// number newly created objects
    int rr;		/// re-used object count

public:
    int verbose;	// for debugging

public:
    typedef typename std::vector<TYPE>::iterator       iterator;
    typedef typename std::vector<TYPE>::const_iterator const_iterator;

    /// constructor
    object_pool()
    {
        nn=0;
        rr=0;
        pool.reserve ( 1000 );
        maxpoolsize=100000;
        verbose=0;
    }

    /// assume the pool is filled with pointers: remove them all and call the destructor
    void reset()
    {
        while ( pool.size() > 0 ) {
            TYPE *t = pool.back();
            pool.pop_back();
            delete t;
        }
    }

    iterator       begin()
    {
        return pool.begin();
    }
    iterator       end()
    {
        return pool.end();
    }
    const_iterator begin() const
    {
        return pool.begin();
    }
    const_iterator end()   const
    {
        return pool.end();
    }

    TYPE&       operator[] ( std::size_t index )
    {
        return pool[index];
    }
    TYPE&       at ( std::size_t index )
    {
        return pool.at ( index );
    }

    TYPE const& operator[] ( std::size_t index ) const
    {
        return pool[index];
    }
    TYPE const& at ( std::size_t index ) const
    {
        return pool.at ( index );
    }

    size_t size() const
    {
        return pool.size();
    }

    TYPE* New()
    {
        if ( ( pool.size() ) ==0 ) { // no reusable object
            nn++;
            if ( nn%10000==0 ) {
                //   myprintf("object_pool::New(): allocate object # %d\n", nn );
            }
            if ( verbose ) {
                myprintf ( "  object_pool::New(): allocate object # %d\n", nn );
            }
            TYPE *t  = new TYPE();
            return t;
        } else {
            rr++;
            if ( rr%5000==0 ) {
                //myprintf("object_pool::New(): re-use object (nn %d, rr %d)\n", nn, rr);
            }
            if ( verbose ) {
                myprintf ( "  object_pool::New(): re-use object (nn %d, rr %d)\n", nn, rr );
            }
            TYPE *t = pool.back();
            pool.pop_back();
            return t;
        }

    }
    void Delete ( TYPE *t )
    {
        pool.push_back ( t );
        //if (pool.size() % 100000==0) myprintf("object_pool::Delete() stored object %zu\n", pool.size() );
        if ( verbose || 0 ) {
            myprintf ( "  object_pool::Delete() stored object %ld\n", ( long ) pool.size() );
        }

    }


};

template <class numtype>
/// lightweight array class
class larray
{
public:
    numtype *d;
    int n;

    larray()
    {
        d=0;
        n=-1;
    }

    larray ( const numtype *data,int nn ) : d(0)
    {
#ifdef OADEBUG
        if ( nn<=0 )
            myprintf ( "larray: constructor from pointer: nn %d\n", nn );
#endif
        alloc ( nn );
        std::copy ( data, data+nn, this->d );
    }

    larray ( int nn ) : d(0)
    {
#ifdef OADEBUG
        if ( nn<=0 )
            myprintf ( "larray: constructor nn %d\n", nn );
#endif
        alloc ( nn );
    }

    numtype * begin()
    {
        return d;
    }

    larray ( const larray &rhs )
    {
        if ( rhs.n<0 ) {
            //myprintf("larray: copy constructor: rhs.n<0! %d\n", rhs.n);
            //this->d[-40]=0;
            this->n=-1;
            this->d=0;
            return;
        }
        alloc ( rhs.n );
        std::copy ( rhs.d, rhs.d+n, this->d );
        //for(int i=0; i<n; i++) this->d[i]=rhs.d[i];
    }

    ~larray()
    {
        //  myprintf("larray::~larray n %d d %ld\n", this->n, (long)this->d );
        if ( d!=0 ) {
            delete [] d;
            n=-3;
        }
        this->d=0;
        //  myprintf("  --> larray::~larray n %d d %ld\n", this->n, (long)this->d );

    }

    void resize ( size_t n )
    {
        if ( this->n!= ( int ) n )
            alloc ( n );
    }
    size_t size() const
    {
        return n;
    }
    //Copy assignment operator
    larray &operator= ( const larray &rhs )
    {
        //  myprintf("larray::operator= n %d %d\n", this->n, rhs.n );
        if ( this->n != rhs.n ) {
            release();
            alloc ( rhs.n );
        } else {
            // myprintf("larray::operator= n %d %d (no re-allocation)\n", this->n, rhs.n );
        }
        //myprintf("larray::operator= after alloc: n %d rhs.n %d\n", this->n, rhs.n );
        std::copy ( rhs.d, rhs.d+n, this->d );

        //for(int i=0; i<n; i++) this->d[i]=rhs.d[i];

        return *this;
    }

    //Copy assignemnt operator
    larray &operator= ( const std::vector<numtype> &rhs )
    {
        //myprintf("larray::operator= (from vector) n %d %zu\n", this->n, rhs.size() );
        int nn =rhs.size();
        release();
        // myprintf("  nn %d\n", nn);
        this->alloc ( nn );
        // myprintf("  alloc done nn %d\n", nn);
        for ( int i=0; i<nn; i++ ) {
            this->d[i]=rhs[i];
        }
        return *this;
    }

    numtype operator[] ( int i ) const
    {
        //  if(i<0) myprintf("larray::operator[] i %d\n", i);
        //   if(i>=n) myprintf("larray::operator[] i %d (n %d)\n", i, n);
        return d[i];
    }

    numtype& at ( size_t i )
    {
        return d[i];
    }

    numtype& operator[] ( int i )
    {
        return d[i];
    }
    bool operator== ( const larray &rhs ) const
    {
        //  myprintf("larray::operator==\n");
        if ( this->n != rhs.n )
            return false;
        for ( int i=0; i<n; i++ ) {
            if ( this->d[i]!=rhs.d[i] )
                return false;
        }
        return true;
    }

    bool operator!= ( const larray &rhs ) const
    {
        return ! ( *this == rhs );
    }
    larray addelement ( numtype v ) const
    {
        larray l ( this->n+1 );
        for ( int i=0; i<n; i++ ) {
            l.d[i]=this->d[i];
        }
        // myprintf("addelement: l.n %d\n", l.n);
        l.d[l.n-1]=v;
        return l;
    }
private:
    void alloc ( int nn )
    {
        if ( nn<=0 )
            myprintf ( "larray: alloc %d\n", nn );
        this->n=nn;
        this->d=new numtype[nn];
    }
    void release()
    {
        // myprintf("larray::release %d, d %ld\n", n, (long)d);
        if ( d!=0 ) {
            delete [] d;
            n=-2;
        }
        this->d=0;
    }
};

template <class NumericType>
/** @brief Multi-value type
 *
 * This object represents a multi-valued object. The objects are ordered using lexicographic ordering.
 */
struct mvalue_t {
public:
    /// vector containing the values
    std::vector<NumericType> v;
    //int direction;
    enum direction_t {HIGH, LOW};
    /// value representing the ordering used
    direction_t d;

    mvalue_t() : d ( HIGH ) {};
    ~mvalue_t() {};

    mvalue_t ( NumericType m, direction_t dd = HIGH )
    {
        v.push_back ( m );
        d= dd;
    }
    mvalue_t ( std::vector<NumericType> vv, direction_t dd = HIGH )
    {
        d= dd;
        this->v = vv;
    }

    template <class T>
    mvalue_t ( std::vector<T> vv , direction_t dd = HIGH )
    {
        d= dd;
        v.clear();
        v.resize ( vv.size() );
        for ( size_t ii =0; ii<vv.size(); ii++ ) {
            this->v[ii] = vv[ii];
        }
    }

    size_t size() const
    {
        return this->v.size();
    }

//Copy assignemnt operator
    mvalue_t &operator= ( const mvalue_t &rhs )
    {
        this->d = rhs.d;
        this->v = rhs.v;
        return *this;
    }
    /// comparison operator
    bool operator== ( const mvalue_t &rhs ) const
    {
        if ( this->v.size() !=rhs.size() )
            return 0;

        for ( size_t i=0; i<this->v.size(); i++ ) {
            if ( v[i]!=rhs.v[i] ) {
                return 0;
            }
        }
        return true;
    }
    bool operator!= ( const mvalue_t &rhs ) const
    {
        if ( this->v.size() !=rhs.size() )
            return 1;

        for ( size_t i=0; i<this->v.size(); i++ ) {
            if ( v[i]!=rhs.v[i] )
                return 1;
        }
        return 0;
    }
    bool operator<= ( const mvalue_t &rhs ) const
    {
        return ! rhs.operator> ( *this );
    }

    bool operator< ( const mvalue_t &rhs ) const
    {
        bool val=0;
        if ( d==HIGH )
            val = (bool)worse ( rhs );
        else
            val= (bool)better ( rhs );
        return val;
    }
    bool operator> ( const mvalue_t &rhs ) const
    {
        bool val=0;
        if ( d==HIGH )
            val= (bool)better ( rhs );
        else
            val = (bool)worse ( rhs );
        // if (dverbose) myprintf("mvalue_t: operator>: %d\n", val);
        return val;
    }
    bool operator>= ( const mvalue_t &rhs ) const
    {
        //        if (dverbose) myprintf("mvalue_t: operator<=");

        return ! rhs.operator< ( *this );
    }

    void show_integer() const {
		for(size_t i=0; i<this->v.size(); i++) {
			myprintf("%ld", (long)this->v[i] );
			if (i<this->v.size()-1 )
				myprintf(",");
		}
	}
	
    template <class W>
    friend std::ostream& operator<< ( std::ostream& stream, const mvalue_t<W>& mval );

#ifdef SWIGCODE
	std::string __repr__x()
	{
		std::stringstream s;
		for(size_t i=0; i<this->v.size(); i++) {
			s << v[i];
			if (i<this->v.size()-1 )
				s <<(",");
		}
		
		return s.str();
	}
#endif	
private:
    int equal ( const mvalue_t &rhs ) const
    {
        for ( size_t i=0; i<this->size(); i++ ) {
            if ( this->v[i]!=rhs.v[i] )
                return 0;
        }
        return 1;
    }

    /** return true if the argument element is larger than this value
	 * 
	 * The comparision is from left to right.
	 * 
	 */    
    int better ( const mvalue_t &rhs ) const
    {
        for ( size_t i=0; i<this->v.size(); i++ ) {
            if ( v[i]>rhs.v[i] )
                return 1;
            if ( v[i]<rhs.v[i] )
                return 0;
        }
        return 0;
    }
    int worse ( const mvalue_t &rhs ) const
    {
        for ( size_t i=0; i<this->v.size(); i++ ) {
            if ( v[i]<rhs.v[i] )
                return 1;
            if ( v[i]>rhs.v[i] )
                return 0;
        }
        return 0;
    }

};

template <class NumericType>
std::ostream& operator << ( std::ostream& stream, const mvalue_t<NumericType> & mval )
{
    std::copy ( mval.v.begin(), mval.v.end(), std::ostream_iterator<NumericType> ( stream, " " ) );
    return stream;
}


template <class Type>
/// Return maximum element of a std::vector
Type vectormax ( const std::vector<Type> &v, Type defaultvalue )
{
    if ( v.size() ==0 )
        return defaultvalue;
    else {
        typename std::vector<Type>::const_iterator p = std::max_element ( v.begin(), v.end() );
        return *p;
    }
}

template <class Type>
/// Return minimum element of a std::vector
Type vectormin ( const std::vector<Type> &v, Type defaultvalue )
{
    if ( v.size() ==0 )
        return defaultvalue;
    else {
        typename std::vector<Type>::const_iterator p = std::min_element ( v.begin(), v.end() );
        return *p;
    }
}



template<class NumType>
/// calculate cumulative sum of a vector
std::vector<NumType> cumsum ( const std::vector<NumType> x )
{
    // initialize the result vector
    std::vector<NumType> res ( x.size() );
    std::partial_sum ( x.begin(), x.end(), res.begin() );
    return res;
}

template<class NumType>
/// calculate cumulative sum of a vector with added zero
std::vector<NumType> cumsum0 ( const std::vector<NumType> x )
{
    // initialize the result vector
    std::vector<NumType> res ( x.size() +1 );
    res[0]=0;
    std::partial_sum ( x.begin(), x.end(), res.begin() +1 );
    return res;
}

template<class Type, class InputType>
std::vector<Type> cumsum0(std::vector<InputType> s)
{
    std::vector<Type> c(s.size()+1);
    c[0]=0;
    for(int i=0; i<s.size(); i++) c[i+1]=c[i]+s[i];
}


template<class NumType>
/// create permutation of specified length
std::vector<NumType> permutation ( int n )
{
    std::vector<NumType> p ( n );
    for ( int i=0; i<n; i++ )
        p[i]=i;
    return p;
}

template<class NumType, class NumTypeIn>
/// convert array given by pointer to std::vector
std::vector<NumType> array2vector ( const NumTypeIn *x, int len )
{
    std::vector<NumType> w;
    w.assign ( x, x + len );
    return w;
}

template<class NumType, class NumTypeIn>
larray<NumType> array2larray ( const NumTypeIn *x, int len )
{
    larray<NumType> w ( len );
    std::copy ( x, x+len, w.d );
    return w;
}

/*! Prints a permutation to output stream
  \param out Output stream
  \param s Pointer to start of array
  \param len Length of array to be printed
  \param maxlen (optional) Maximum length to print
  \brief Print permutation
  */
template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( std::ostream &out, const permutationType *s, const int len, const int maxlen = 256 )
{
    out << "{";

    int plen = std::min ( len, maxlen );
    for ( int i = 0; i < plen-1 ; i++ )
        out << s[i] << ",";

    if ( len==0 ) {
        // corner case
        out << "}\n";

    } else {
        if ( plen<len )
            out  << s[plen-1] << ",...}\n";
        else
            out << s[plen-1] << "}\n";
    }
}
template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( std::ostream &out, const larray<permutationType> s, const int maxlen = 256, const bool ret = true )
{
    out << "{";
    const int len = s.size();
    int plen = std::min ( len, maxlen );
    for ( int i = 0; i < plen-1 ; i++ )
        out << s[i] << ",";

    if ( len==0 ) {
        // corner case
        out << "}\n";

    } else {
        if ( plen<len )
            out  << s[plen-1] << ",...}";
        else
            out << s[plen-1] << "}";
    }
    if ( ret ) {
        out << "\n";
    }
}
template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( std::ostream &out, const std::vector<permutationType> s, const int maxlen = 256, const bool ret = true )
{
    int len = s.size();

    out << "{";

    int plen = std::min ( len, maxlen );
    for ( int i = 0; i < plen-1 ; i++ )
        out << s[i] << ",";

    if ( len==0 ) {
        // corner case
        out << "}";

    } else {
        if ( plen<len )
            out  << s[plen-1] << ",...}";
        else
            out << s[plen-1] << "}";
    }
    if ( ret ) {
        out << "\n";
    }
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm_int ( const std::vector<permutationType> s, const int maxlen = 256, const bool ret = true )
{
    int len = s.size();
    int plen = std::min ( len, maxlen );

    myprintf("{");

    for ( int i = 0; i < plen-1 ; i++ )
        myprintf("%d,", s[i]);

    if ( len==0 ) {
        // corner case
        myprintf("}");

    } else {
        if ( plen<len )
            myprintf("%d,...",  s[plen-1] );
        else
            myprintf("%d}",  s[plen-1] );
    }
    if ( ret ) {
		myprintf("\n");
    }
}

#ifdef FULLPACKAGE

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
/// print permutation with string in front 
static void print_perm ( const char *msg,  const std::vector<permutationType> s, const int maxlen = 256, const bool ret = true )
{
	myprintf("%s: ", msg);
    print_perm ( std::cout, s, maxlen, ret );
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( const larray<permutationType> s, const int maxlen = 256, const bool ret = true )
{
    print_perm ( std::cout, s, maxlen, ret );
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( const std::vector<permutationType> s, const int maxlen = 256, const bool ret = true )
{
    print_perm ( std::cout, s, maxlen, ret );
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( const permutationType *s, const int len, const int maxlen = 256 )
{
    print_perm ( std::cout, s, len, maxlen );
}

#else
// dummy values
template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( const larray<permutationType> s, const int maxlen = 256, const bool ret = true )
{
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( const std::vector<permutationType> s, const int maxlen = 256, const bool ret = true )
{
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
static void print_perm ( const permutationType *s, const int len, const int maxlen = 256 )
{
}
#endif

#define print_comb print_perm


template <class numtype>
/**
 * Compare two arrays and return whether equal or not.
 * @param A Pointer to array
 * @param B Pointer to array
 * @param r Number of rows
 * @param c Number of columns
 * @return
 */
int compare_matrix ( const numtype *A, const numtype *B, int r, int c )
{
    for ( int x=0; x<r; x++ )
        for ( int y=0; y<c; y++ ) {
            if ( A[x+y*r]!=B[x+y*r] ) {
                myprintf ( "arrays unequal: %d, %d\n", x, y );
                return 0;
            }
        }
    myprintf ( "arrays equal\n" );
    return 1;
}


/*!
 *	A small function that calculates the factorial of a number. Can be inlined if the compiler decides
 *	it is faster. Returns one if the number is smaller or equal to 1
 	\brief Calculates factorial
 *	\param f	Number to calculate the factorial of
 */
template <class Type>
inline Type fact ( const Type f )
{
    //Factorial
    if ( f <= 1 )
        return 1;
    Type  sol = 1;
    for ( int i = f; i > 1; i-- ) {
#ifdef OAOVERFLOW
        if ( sol>std::numeric_limits<Type>::max() /100 ) {
            myprintf ( "fact: f %ld, i %d:  %ld, %ld\n", ( long ) f, i, ( long ) sol, ( long ) std::numeric_limits<Type>::max() );
        }
#endif
        sol *= i;
    }
    return sol;
}

/*!
  This inline factorial function is the same as the standard factorial calculation, except that the
  return type is generic
  \brief Calculates factorial of type numtype
  \param f number to calculate factorial of
  */
template <class numtype, class argtype>
static inline numtype factorial ( const argtype f )
{
    numtype	sol = 1;

    if ( f <= 1 )
        return 1;
    for ( argtype i = 2; i <= f; i++ ) {
#ifdef OAOVERFLOW
        if ( sol > std::numeric_limits<numtype>::max() / ( 4*i ) ) {
            std::cout << "factorial: possible numberic overflow: f " << f << ", i " << i << ":  " << sol << ", " <<  std::numeric_limits<int>::max() << std::endl;
        }
#endif
        sol *= i;
    }
    return sol;
}


/*!
  The number of combinations is calculated using the an addapted formula
  \brief Calculates number of combinations
  \param n Total number of entries to choose from
  \param k Number of entries in a certain combination
  */
template <class Type>
inline Type ncombs ( const Type n, const Type k )
{
    register int i;
    Type sol = 1;	///n!/(k! * (n-k)!) = (n - k + 1) * ..... * n/k!
    for ( i = n - k + 1; i <= n; i++ ) // since n-k > k usually
        sol *= i;
    return sol/fact ( k );
}

int ncombscacheNumber();

/// initialize datastructure, this function is not thread safe
void initncombscache(int N);

/** return number of combinations from previously calculated results
 *
 *  The results should be calculated with initncombscache
 */
long ncombscache(int n, int k);

template <class Type>
/// calculate using multiplicative formula, see http://en.wikipedia.org/wiki/Binomial_coefficient#Computing_the_value_of_binomial_coefficients
inline Type ncombsm ( const Type &n, const Type &k )
{
    //if (k > n) throw std::logic_error("k can not be larger than n");
    Type result = 1;
    for ( Type i = 1; i <= k; i++ ) {
        result *= n- ( k-i ); //separating * and / allows us to use integer * and /, do not combine
        result /= i;
    }
    return result;
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
/*
* Note: this function does check for the last permutation, but does not act on this
*/
bool next_perm ( std::vector<permutationType> &s )
{
    int len=s.size();
    int	i = len - 1, j = len;
    //permutationType	tmp;

    while ( i>0 && ( s[i-1] >= s[i] ) )
        i--;

    if ( i==0 ) {
        // last permutation reached!
        // myprintf("last permutation reached\n");
        for ( int k=0; k< ( len/2 ); k++ ) {
            std::swap ( s[k], s[len-k-1] );
        }
        return true;
    }

    while ( s[j-1] <= s[i-1] )
        j--;

    std::swap ( s[j - 1], s[i - 1] );	// swap values at positions (i-1) and (j-1)

    j = len-1;
    while ( i < j ) {
        std::swap ( s[j], s[i] ); // swap values at positions i and j
        i++;
        j--;
    }

    return false;
}

/* Random number generators */

// return random integer
int fastrand();
void seedfastrand ( int s );

// return random integer in range 0 to k-1
int fastrandK ( int k );

#ifdef RPACKAGE
// R packages are not allowed to use rand
#define myrand fastrand
#else
/// set the random number seed using srand
void set_srand ( unsigned int s );
#define myrand rand
#endif

#ifdef RPACKAGE
template<typename myRandomAccessIterator>
inline void
my_random_shuffle ( myRandomAccessIterator myfirst, myRandomAccessIterator mylast )
{
    // concept requirements
    //__glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<
    //  _RandomAccessIterator>)
    //__glibcxx_requires_valid_range(__first, __last);

    if ( myfirst != mylast )
        for ( myRandomAccessIterator __i = myfirst + 1; __i != mylast; ++__i ) {
            myRandomAccessIterator __j = myfirst + fastrand() % ( ( __i - myfirst ) + 1 );
            if ( __i != __j )
                std::iter_swap ( __i, __j );
        }
}
#else
#define my_random_shuffle std::random_shuffle
#endif

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
/*
* Specialized function for permutations of length two
*/
inline void next_perm_twoperm ( permutationType *s, const int len )
{
    std::swap ( s[0], s[1] );
}

template <class permutationType>	/* permtype should be a numeric type, i.e. int or long */
/*
* Note: this function does check for the last permutation, but does not act on this
*/
void next_perm ( permutationType *s, const int len )
{
    int	i = len - 1, j = len;
    permutationType	tmp;

    while ( i>0 && ( s[i-1] >= s[i] ) )
        i--;

    if ( i==0 ) {
        // last permutation reached!
        return;
    }

    while ( s[j-1] <= s[i-1] )
        j--;

    tmp = s[j - 1];		// swap values at positions (i-1) and (j-1)
    s[j - 1] = s[i - 1];
    s[i - 1] = tmp;

    j = len-1;
    while ( i < j ) {
        tmp=s[j];
        s[j]=s[i];
        s[i]=tmp; // swap values at positions i and j
        i++;
        j--;
    }
}

template <class numtype, class objecttype>	/* permtype should be a numeric type, i.e. int or long */
/**
 * See also http://en.wikipedia.org/wiki/Permutation#Numbering_permutations
 * @param
 * @param s
 * @param len
 * @return
 */
numtype* permutationLex ( numtype k, objecttype *s, numtype n )
{
    numtype fact= factorial<numtype,int> ( n-1 );
    objecttype tempj, temps;

    for ( int j=0; j<n-1; j++ ) {
#ifdef OADEBUG
        if ( fact==0 )
            myprintf ( "division by zero: j %d, fact %d, k %d\n", j, fact, k );
#endif
        tempj = ( k/ fact ) % ( n - j );
        //printf("j %d, fact %d, tempj %d\n", j, fact, tempj);
        temps = s[j+ tempj];
        for ( int i=j+ tempj; i>=j+1; i-- ) {
            s[i]= s[i- 1];      // shift the chain right
        }
        s[j]= temps;
        fact= fact/ ( n- ( j+1 ) );
    }
    return s;
}



template <class objecttype, class numtype>
/** Create random permutation using Fisher-Yates shuffle, or Knuth shuffle
 */
void random_perm ( objecttype *s, numtype len )
{
    for ( numtype i=0; i<len-1; i++ ) {
        numtype j = i+myrand() % ( len-i );
        std::swap ( s[i], s[j] );
    }
}

template <class objecttype>
/** Create random permutation using Fisher-Yates shuffle, or Knuth shuffle
 */
void random_perm ( std::vector<objecttype> &s)
{
    int len=s.size();
    for ( int i=0; i<len-1; i++ ) {
        int j = i+myrand() % ( len-i );
        std::swap ( s[i], s[j] );
    }
}

template <class numtype>
//! @brief Create a new combination and initialize
inline numtype *new_comb_init ( int len )
{
    numtype *comb = new numtype [len];
    for ( int i=0; i<len; i++ )
        comb[i] = i;
    return comb;
}

template <class numtype>
//! @brief Delete combination
inline void delete_comb ( numtype *comb )
{
    delete [] comb;
}


template <class numtype>
/**
 * Initialize a combination
 *
 * @param comb Pointer to combination array
 * @param k Number of elements
 * @param n Numbers to choose from
 * @return Number of combinations possible
 */
inline int init_comb ( numtype *comb, int k, int n )
{
    for ( int i=0; i<k; i++ )
        comb[i]=i;
    return ncombs ( n,k );
}


/*!
  Gives combination number k based on an algorithm from wikipedia.
  \brief Gives combination nr k
  \param comb Pointer to combination
  \param k Number of the current combination
  \param n Number of elements in combination
  */
template <class numtype>
numtype next_combination ( numtype *comb, int k, int n )
{
    //all possible combinations of n out of k numbers
    int             i;// = k - 1;
    const int       offset = n - k + 1;
    //comb[k - 1] = n - 1;
    i = k - 1;
    comb[i]++;
    while ( ( comb[i] >= offset + i ) && ( i > 0 ) ) {
        //myprintf("next_combination: while 1: i %d, comb[i] %d\n", i, comb[i]);
        i--;
        comb[i]++;
    }

    if ( comb[0] > n - k )
        return 0; /* No more combinations can be generated */

    /* comb now looks like (…, x, n, n, n, …, n).
       Turn it into (…, x, x + 1, x + 2, …) */
    for ( i++; i < k; i++ )
        comb[i] = comb[i - 1] + 1;

    return 1;
}

template <class numtype>
numtype next_combination_fold ( numtype *comb, int k, int n )
{
    //all possible combinations of n out of k numbers
    int             i;// = k - 1;
    const int       offset = n - k + 1;
    //comb[k - 1] = n - 1;
    i = k - 1;
    comb[i]++;
    while ( ( comb[i] >= offset + i ) && ( i > 0 ) ) {
        //printf("next_combination: while 1: i %d, comb[i] %d\n", i, comb[i]);
        i--;
        comb[i]++;

    }

    int fold = i;

    if ( comb[0] > n - k )
        return 0; /* No more combinations can be generated */

    /* comb now looks like (…, x, n, n, n, …, n).
       Turn it into (…, x, x + 1, x + 2, …) */
    for ( i++; i < k; i++ )
        comb[i] = comb[i - 1] + 1;

    return fold;
}

inline void print_combinations ( int n, int k )
{
    int *comb = new_comb_init<int> ( k );
    int nc = ncombs ( n, k );
    for ( int i=0; i<nc; i++ ) {
        print_comb ( comb, k );
        next_combination<int> ( comb, k, n );
    }
    delete_comb ( comb );
}

/* code related to permutations */
template <class numtype>
/**
 * @brief Check whether a permutation is ordered or not
 * @param perm
 * @param len
 * @return Return 1 of an ordered permutation, 0 otherwise
 */
int perm_is_ordered ( numtype *perm, int len )
{
    for ( int i=0; i<len-1; i++ ) {
        if ( perm[i]>perm[i+1] )
            return 0;
    }
    return 1;
}

template <class numtype>
//! @brief Create a new permutation
numtype *new_perm ( int len )
{
    return ( numtype * ) malloc ( sizeof ( numtype ) *len );
}

template <class numtype>
//! @brief Create a new permutation
numtype *clone_perm ( numtype *source, int len )
{
    numtype *perm = new_perm<numtype> ( len );
    copy_perm ( source, perm, len );
    return perm;
}

template <class numtype>
//! @brief Delete a permutation
inline void delete_perm ( numtype *perm )
{
    free ( perm );
}

template <class numtype>
/**
 * @brief Invert a permutation
 * @param perm Permutation as integer type std::vector
 * @return New permutation that is the inverse of the argument
 */
std::vector<numtype> invert_permutation ( const std::vector<numtype> perm )
{
    std::vector<numtype> iperm(perm.size());
    for ( size_t x=0; x<perm.size(); x++ )
        iperm[perm[x]]=x;
    return iperm;
}

template <class numtype>
/**
 * @brief Invert a permutation
 * @param perm Permutation as integer type std::vector
 * @param permout Output permutation
 */
void invert_permutation ( const std::vector<numtype> perm, std::vector<numtype> &iperm )
{
    iperm.resize(perm.size());
    for ( size_t x=0; x<perm.size(); x++ )
        iperm[perm[x]]=x;
}

template <class numtype>
/**
 * @brief Invert a permutation
 * @param perm Pointer to permutation
 * @param len
 * @return Pointer to new permutation that is the inverse of the argument
 */
void invert_permutation ( numtype *perm, int len, numtype *iperm )
{
    for ( int x=0; x<len; x++ )
        iperm[perm[x]]=x;
}

template <class numtype>
/**
 * @brief Invert a permutation
 * @param perm Pointer to permutation
 * @param len
 * @return Pointer to new permutation that is the inverse of the argument
 */
numtype *invert_permutation ( numtype *perm, int len )
{
    numtype *iperm = new_perm<numtype> ( len );

    for ( int x=0; x<len; x++ )
        iperm[perm[x]]=x;
    return iperm;
}

#ifdef SAFELPERM
template <class numtype>
//! Perform level permutation with bounds check
inline numtype safe_lperm ( numtype val, const numtype *lperm, int n )
{
    if ( val<0 || val>=n ) {
        return val;
    }
    return lperm[val];
}


template <class numtype>
//! Perform level permutation with bounds check
inline void safe_perform_level_perm ( numtype * src, int n, const numtype * perm, const int pmax )
{
    for ( int i=0; i<n; i++ ) {
        src[i]=safe_lperm<numtype> ( src[i], perm, pmax );
    }
}

template <class numtype>
//! Perform level permutation with bounds check
inline void safe_perform_level_perm ( const numtype * src, numtype * dst, int n, const numtype * perm, const int pmax )
{
    for ( int i=0; i<n; i++ ) {
        dst[i]=safe_lperm<numtype> ( src[i], perm, pmax );
    }
}
#endif

template <class numtype>
/**
 * Perform level permutation on an array
 * @param src Pointer to array
 * @param n Length of array
 * @param perm Permutation to perform
 */
inline void perform_level_perm ( numtype * src, int n, const numtype * perm )
{
    for ( int i=0; i<n; i++ ) {
        src[i]=perm[src[i]];
    }
}

template <class numtype>
/**
 * @brief Calculate composition of 2 permutations
 *
 * Calculates C = B \circ A
 * @param A
 * @param B
 * @param n
 * @param C
 */
inline void composition_perm ( const numtype * A, const  numtype* B, int n, numtype * C )
{
    for ( int i=0; i<n; i++ ) {
        C[i] = B[A[i]];
    }
}

template <class numtype>
/**
 * @brief Calculate composition of 2 permutations
 *
 * Calculates C = B \circ A
 * @param A
 * @param B
 * @param n
 * @param C
 */
inline void composition_perm ( const std::vector<numtype> &A, const std::vector<numtype> &B,  std::vector<numtype> &C )
{
    for ( size_t i=0; i<A.size(); i++ ) {
        C[i] = B[A[i]];
    }
}

template <class object, class numtype>
/**
 * @brief Perform a permutation on a set of objects
 * @param src
 * @param target
 * @param n
 * @param perm
 */
inline void perform_perm ( const object *const src, object *const target, const int n, const numtype * perm )
{
    for ( int i=0; i<n; i++ ) {
        target[perm[i]]=src[i];
    }
}

/// Perform a permutation on a set of objects
template <class object, class numtype>
inline std::vector<object> perform_perm ( const std::vector<object> src, const std::vector<numtype> perm )
{
    //myprintf("perform_perm: src %d, perm %d\n", src.size(), perm.size() );
    std::vector<object> dst ( src.size() );
    for ( size_t i=0; i<perm.size(); i++ ) {
        dst[perm[i]]=src[i];
    }
    return dst;
}

/// Perform inverse permutation
template <class object, class numtype>
inline void perform_inv_perm ( const std::vector<object> src, std::vector<object> &target, const int n, const std::vector<numtype> perm )
{
    for ( int i=0; i<n; i++ ) {
        target[i]=src[perm[i]];
    }
}

/// Perform inverse permutation
template <class object, class numtype>
inline void perform_inv_perm ( const object *const src, object *const target, const int n, const numtype * perm )
{
    for ( int i=0; i<n; i++ ) {
        target[i]=src[perm[i]];
    }
}

/// Perform inverse permutation
template <class object, class numtype>
inline void perform_inv_perm ( const object *const src, object *const target, const int n, const std::vector<numtype> &perm )
{
    for ( int i=0; i<n; i++ ) {
        target[i]=src[perm[i]];
    }
}

template <class numtype>
/**
 * @brief Perform a permutation on a set of data elements
 * @param src
 * @param target
 * @param n
 * @param perm
 */
inline void perform_level_perm ( const numtype *const src, numtype *const target, const int n, const numtype * perm )
{
    for ( int i=0; i<n; i++ ) {
        target[i]=perm[src[i]];
    }
}

template <class numtype>
/**
 * Initialiaze a permutation
 * @param perm
 * @param len
 */
void init_perm ( std::vector<numtype> &perm)
{
    for ( size_t i=0; i<perm.size(); i++ )
        perm[i]=i;
}

template <class numtype>
/**
 * Initialiaze a permutation
 * @param perm
 * @param len
 */
void init_perm ( numtype *perm, int len )
{
    for ( int i=0; i<len; i++ )
        perm[i]=i;
}

// Initialize sign permutation
template <class numtype>
/**
 * Initialiaze a permutation
 * @param perm
 * @param len
 */

void init_signperm ( std::vector<numtype> &signperm)
{
    for ( size_t i=0; i<signperm.size(); i++ )
        signperm[i]=1;
}


template <class numtype>
bool compare_perm ( const numtype *permA, const numtype *permB, int len )
{
    return std::equal ( permA, permA+len, permB );
}


template <class numtype>
/// copy a permuntation
inline void copy_perm ( const numtype *source, numtype *target, int len )
{
    memcpy ( target, source, sizeof ( numtype ) *len );
}

template <class numtype, class outtype>
/// initialize a permutation and return the number of permutations
inline outtype init_perm_n ( numtype *perm, int len )
{
    for ( int i=0; i<len; i++ )
        perm[i]=i;
    return factorial<outtype> ( len );
}


template <class numtype>
numtype *new_perm_init ( int len )
{
    numtype *perm = new_perm<numtype> ( len );
    init_perm ( perm, len );

    return perm;
}

//#include <iostream>
//#include <iterator>

template <typename _ForwardIterator>
inline bool issorted ( _ForwardIterator first, const _ForwardIterator last )
{
    if ( first == last )
        return true;
//    _ForwardIterator first = __first;
    _ForwardIterator __next = first;
    for ( ++__next; __next != last; first = __next, ++__next )
        if ( *__next < *first )
            return false;
    return true;

}

/* templates for valueindex */

template <class returntype, class basetype, class numtype>
returntype* new_valueindex ( const basetype *bases, const numtype n )
{
    //myprintf("init_valueindex: n %d\n", n);
    returntype *valueindex = new returntype [n * sizeof ( returntype )];
    assert ( valueindex!=0 );
#ifdef OADEBUG
    if ( n==0 ) {
        myprintf ( "valueindex of size 0\n" );
        exit ( 0 );
    }
#endif
    valueindex[n - 1] = 1;

    for ( int i = n - 2; i >= 0; i-- )
        valueindex[i] = valueindex[i + 1] * bases[i + 1];

    //return valueindex[0]*bases[0];	// return max value
    return valueindex;
}

template <class numtype>
numtype* init_valueindex_forward ( numtype *valueindex, const numtype *bases, const numtype n )
{
    valueindex[0] = 1;

    for ( int i = 0; i< ( n-1 ); i++ )
        valueindex[i+1] = valueindex[i] * bases[i];

    return valueindex;
}

template <class numtype>
numtype* init_valueindex ( numtype *valueindex, const numtype *bases, const numtype n )
{
    valueindex[n - 1] = 1;

    for ( int i = n - 2; i >= 0; i-- )
        valueindex[i] = valueindex[i + 1] * bases[i + 1];

    return valueindex;
}

template <class Type >
/// Helper class
class sort_indices
{
private:
    Type *mparr;
public:
    sort_indices ( Type *parr ) : mparr ( parr ) {}
    bool operator() ( int i, int j )
    {
        return mparr[i]<mparr[j];
    }
};

template <class ContainerType >
/// Helper class
class sort_indices_container
{
private:
    const ContainerType *mparr;
    bool ascending;
public:
    sort_indices_container ( const ContainerType *parr, bool ascendingx = true ) : mparr ( parr ), ascending ( ascendingx ) {}
    bool operator() ( int i, int j )
    {
        if ( ascending )
            return mparr->at ( i ) < mparr->at ( j );
        else
            return mparr->at ( i ) > mparr->at ( j );
    }
};

template <class Type >
/// Helper class
class sort_indices_deque
{
private:
    const std::deque<Type> *mparr;
    bool ascending;
public:
    sort_indices_deque ( const std::deque<Type> *parr, bool ascendingx = true ) : mparr ( parr ), ascending ( ascendingx ) {}
    bool operator() ( int i, int j )
    {
        if ( ascending )
            return mparr->at ( i ) <mparr->at ( j );
        else
            return mparr->at ( i ) >mparr->at ( j );
    }
};

template <class Type >
/// Helper class
class sort_indices_vector
{
private:
    const std::vector<Type> *mparr;
    bool ascending;
public:
    sort_indices_vector ( const std::vector<Type> *parr, bool ascendingx = true ) : mparr ( parr ), ascending ( ascendingx ) {}
    bool operator() ( int i, int j )
    {
        if ( ascending )
            return mparr->at ( i ) <mparr->at ( j );
        else
            return mparr->at ( i ) >mparr->at ( j );
    }
};

/** @brief Class to sort data without moving the data in memory
 *
 * The data is sorted by using a list of indices. A stable sort is being used.
 *
 */
class indexsort
{
public:
    std::vector<int> indices;
private:
    int n;
public:
    indexsort ( int nn ) : n ( nn )
    {
        indices.resize ( n );
        for ( int i=0; i<n; i++ )
            indices[i]=i;
    }

    template<class Type>
    /// Constructor for deque class
    indexsort ( const std::deque<Type> &vals )
    {
        init ( vals );
    }

    template<class Type>
    /// Constructor for vector class
    indexsort ( const std::vector<Type> &vals )
    {
        init ( vals );
    }

    template<class Type>
    /// initialize sorting structure with specified values
    void init ( const std::deque<Type> &vals )
    {
        n = vals.size();
        indices.resize ( n );
        for ( int i=0; i<n; i++ )
            indices[i]=i;
        this->sort ( vals );
    }
    template<class Type>
    /// initialize sorting structure with specified values
    void init ( const std::vector<Type> &vals )
    {
        n = vals.size();
        indices.resize ( n );
        for ( int i=0; i<n; i++ )
            indices[i]=i;
        this->sort ( vals );
    }

    template<class Type>
    /// sort values and store the indices
    void sort ( const Type *vals )
    {
        std::stable_sort ( indices.begin(), indices.end(), sort_indices<Type> ( vals ) );
    }
    template<class Type>
    /// sort values and store the indices
    void sort ( const std::vector<Type> &vals )
    {
        std::stable_sort ( indices.begin(),indices.end(), sort_indices_vector<Type> ( &vals ) );
    }
    template<class Type>
    /// sort values and store the indices
    void sort ( const std::deque<Type> &vals )
    {
        std::stable_sort ( indices.begin(),indices.end(), sort_indices_deque<Type > ( &vals ) );
    }

    template<class Type>
    /// sort values and store the indices
    void sortdescending ( const std::vector<Type> &vals )
    {
        std::stable_sort ( indices.begin(),indices.end(), sort_indices_vector<Type> ( &vals, false ) );
    }
    void show() const
    {
        for ( int i=0; i<n; i++ )
            myprintf ( "%d ", indices[i] );
    }

    template<class Type>
    /// return array sorted using the order from the indexsort structure
    std::vector<Type> sorted ( const std::vector<Type> &vals ) const
    {
        std::vector<Type> s ( n );
        //myprintf("here\n");
        for ( int i=0; i<n; i++ )
            s[i] = vals[indices[i]];
        //myprintf("here 2\n");
        return s;
    }
    template<class ContainerType>
    /// return array sorted using the order from the indexsort structure
    ContainerType sorted ( const ContainerType &vals ) const
    {
        ContainerType s ( n );
        for ( int i=0; i<n; i++ )
            s[i] = vals[indices[i]];
        return s;
    }

    /// Returns true of the data is sorted ascending
    bool issorted() const
    {
        for ( int i=0; i<n-1; i++ ) {
            if ( indices[i]>indices[i+1] )
                return false;
        }
        return true;
    }
    /// Returns true of the data is sorted descending
    bool issorteddescending() const
    {
        for ( int i=0; i<n-1; i++ ) {
            if ( indices[i]<indices[i+1] )
                return false;
        }
        return true;
    }

};

template <class numtype>
std::vector<int> argsort(const std::vector<numtype> vv)
{
    indexsort idx(vv);
    return idx.indices;
}

#ifdef FULLPACKAGE
#include "InfInt.h"
#endif

#include <limits>

/** @brief Class to describe the symmetry group of a list of elements
 *
 * The class assumes the list is sorted. The symmetry group is then a direct product of full permutation groups.
 *
 * We do not implement this using templates because we want to export to Python.
 */
class symmetry_group
{
public:
    std::vector<int> gidx;
    std::vector<int> gstart;	/// start of the subgroups
    std::vector<int> gsize; 	/// size of the subgroups
    int ngroups; /// number of subgroups
    int n; /// number of elements
    bool ascending; /// ordering of elements

public:
    symmetry_group ( const std::vector<int> &vals, bool ascending = true, int verbose=0 );
    symmetry_group ( const std::vector<double> &vals, bool ascending = true,int verbose=0 );
    symmetry_group ( const std::vector<float> &vals, bool ascending = true,int verbose=0 );
    symmetry_group ( const std::vector<short int> &vals, bool ascending = true,int verbose=0 );
    symmetry_group ( const std::vector<unsigned int> &vals, bool ascending = true,int verbose=0 );
    symmetry_group ( const std::vector<mvalue_t<double> > &vals, bool ascending = true,int verbose=0 );
    symmetry_group ( const std::vector<mvalue_t<int> > &vals, bool ascending = true,int verbose=0 );
    symmetry_group ( const symmetry_group &sgx ); /// copy constructor
    symmetry_group ( ); /// default constructor

    //	symmetry_group(const std::vector<int> vals, int verbose=0);
private:
    template<class Type>
    void init ( const std::vector<Type> vals, bool ascending = true, int verbose=0 );
public:

    symmetry_group& operator= ( const symmetry_group &sgx ); /// copy assignment

#ifdef FULLPACKAGE

    typedef long perm_return_type;
//#define perm_return_type long
    /// return size of the group of all permutations respecting the symmetry
    perm_return_type permsize() const
    {
        perm_return_type  s = 1;
        for ( int i=0; i<ngroups; i++ ) {
            perm_return_type f = factorial<long> ( gsize[i] );

            //myprintf("i %d: f %ld, s %ld s int %d, max %ld\n", i,f, s, (int)(s), std::numeric_limits< long long>::max() );
            if ( f != 0 && ( std::numeric_limits< perm_return_type>::max() / f ) < s ) {
                // multiplication would exceed range of unsigned
                myprintf ( "symmetry_group::init: group size outside limits, please use permsize_large\n" );
                throw;
                return -1;
            }

            s *= f;
        }
        return s;
    }
#endif

#ifdef FULLPACKAGE
    /// return size of the group of all permutations respecting the symmetry
    InfInt permsize_large() const
    {
        InfInt s = 1;
        for ( int i=0; i<ngroups; i++ ) {
            long f = factorial<long> ( gsize[i] );

            s *= f;
        }
        return s;
    }
#endif

    /// list with indices set to check for symmetry reductions
    std::vector<int> checkIndices() const
    {
        std::vector<int> check_indices(n);

        // never check first index
        for (int row=1; row<n; row++ ) {
            if (this->gidx[row]==this->gidx[row-1] )
                check_indices[row]=1;
        }
        return check_indices;
    }


    /// representation function (for python interface)
    std::string __repr__() const
    {
        std::stringstream ss;
        ss << "symmetry group: " << n << " elements, " << ngroups << " subgroups: ";
        for ( int i=0; i<ngroups; i++ )
            ss << gsize[i] << " ";
        // ss << std::endl;

        std::string s = ss.str();
        return s;
    }

    /// show the symmetry group
    void show ( int verbose=1 ) const
    {
        myprintf ( "symmetry group: %d elements, %d subgroups: ", n, ngroups );
        for ( int i=0; i<ngroups; i++ )
            myprintf ( "%d ", gsize[i] );
        myprintf ( "\n" );

        if ( verbose>=2 ) {
            myprintf ( "gidx: " );
            for ( int i=0; i<n; i++ )
                myprintf ( "%d, ", gidx[i] );
            myprintf ( "\n" );
            myprintf ( "gstart: " );
            for ( int i=0; i<ngroups; i++ )
                myprintf ( "%d, ", gstart[i] );
            myprintf ( "\n" );
            myprintf ( "gsize: " );
            for ( int i=0; i<ngroups; i++ )
                myprintf ( "%d, ", gsize[i] );
            myprintf ( "\n" );

        }
    }
    // void symmetry_group< int >(std::vector< int > vals, int verbose);

};


/// class to walk over all element of a symmetry group
class symmetry_group_walker
{
public:
    symmetry_group sg;
    std::vector< std::vector<int> > perms;

    symmetry_group_walker ( symmetry_group sg )
    {
        this->sg = sg;
        perms.resize ( sg.ngroups );

        init();
    }

    void show ( int verbose=1 ) const;


    void init()
    {
        perms.resize ( sg.ngroups );
        for ( size_t i=0; i< ( size_t ) sg.ngroups; i++ ) {
            perms[i].resize ( sg.gsize[i] );
            for ( size_t j=0; j< ( size_t ) sg.gsize[i]; j++ )
                perms[i][j]=j;
        }
    }

    bool next()
    {
        int g = perms.size()-1;
        return nextsub ( g );
    }

    std::vector<int> fullperm() const;


protected:
    bool nextsub ( int g );

};


template<class Type, class IndexType>
/// Permute a std::vector
std::vector<Type> permute ( const std::vector<Type> x, const std::vector<IndexType> indices )
{
    std::vector<Type> y ( x.size() );
    for ( size_t i=0; i<x.size(); i++ ) {
        //myprintf ( " permute %d: y[%d]=x[%d]=%f\n", i, i, indices[i], x[indices[i]] );
        y[i]=x[indices[i]];
    }
    return y;
}

template<class Type, class IndexType>
/// Permute a std::vector with inverse permutation
std::vector<Type> permuteback ( const std::vector<Type> x, const std::vector<IndexType> indices )
{
    std::vector<Type> y ( x.size() );
    for ( int i=0; i<x.size(); i++ )
        y[indices[i]]=x[i];
    return y;
}

template <class numtype, class itype>
/**
 * @brief Calculate symmetry groups of a list of integers under permutations
 * @param vec
 * @param n
 * @param idx
 * @param gstart
 * @param gsize
 * @return Number of groups found
 */
int symm_group_index_plain ( const numtype *vec, const int n, itype *& idx, itype *& gstart, itype *& gsize )
{
    int i;
    int nsg=0;
    numtype prev;

    prev = -1;
    /* count number of symmetry groups */
    for ( i=0; i<n; i++ ) {
        if ( vec[i]!=prev ) {
            nsg++;
        }
        prev=vec[i];
    }

    /* initialize data structure */
    idx = new itype[n];
    gstart = new itype[nsg+1];
    gsize = new itype[nsg+1];

    /* find starting positions of symmetry group */
    prev = -1;
    nsg = 0;
    for ( i=0; i<n; i++ ) {
        //myprintf("symm_group_index: nsg: %d i %d\n", nsg, i);
        if ( vec[i]!=prev ) {
            gstart[nsg]=i;
            nsg++;
        }
        idx[i]=nsg-1;
        prev=vec[i];
    }
    gstart[nsg]=i;	/* add dummy element */

    for ( i=0; i<nsg; i++ ) {
        gsize[i] = gstart[i+1]-gstart[i];
    }

    //printf("symm_grp:\n"); print_perm(vec,n); print_perm(idx, n); print_perm(gstart, nsg);
    return nsg;
}

/// Power of two integers
inline long ipow ( long x, long y )
{
    long r = 1;
    for ( int j=0; j<y; j++ )
        r*=x;

    return r;
}

/// Power of two unsigned integers
inline unsigned int ipow ( unsigned int x, unsigned int p )
{
    if ( p == 0 )
        return 1;
    if ( x == 0 )
        return 0;

    unsigned int r = 1;
    for ( ;; ) {
        if ( p & 1 )
            r *= x;
        if ( ( p >>= 1 ) == 0 )
            return r;
        x *= x;
    }
}

/// -1 to the power n (integer)
inline int powmo ( int n )
{
    return ( n%2==0 ) ? 1: -1;
}

template <class Type>
/// calculate value of Krawtchouk polynomial
inline Type krawtchouk ( Type j, Type x, Type n, Type s, int verbose=0 )
{
    Type val=0;

    for ( Type i=0; i<=j; i++ ) {
        //printf("       ncombs(%d, %d) %d\n", x,i, ncombs(x, i) );

        // TODO: implement dynamic programming for this function
        //val += powmo(i) * std::pow ((double)s-1, j-i) * ncombs(x, i) * ncombs(n-x, j-i);
        val += powmo ( i ) * ipow ( s-1, j-i ) * ncombs ( x, i ) * ncombs ( n-x, j-i );
        if ( verbose ) {
            Type tt= std::pow ( double ( -1 ), ( double ) i ) * std::pow ( ( double ) ( s-1 ), ( double ) ( j-i ) ) * ncombs ( x, i ) * ncombs ( n-x, j-i );
            myprintf ( "    krawtchouk(%d, %d, %d, %d) term %d: %d=%d*%d*%d*%d\n", ( int ) j, ( int ) x, ( int ) n, ( int ) s, ( int ) i, ( int ) tt, ( int ) std::pow ( ( double )-1, ( double ) i ), ( int ) std::pow ( ( double ) s-1, ( double ) ( j-i ) ) , ( int ) ncombs ( x, i ), ( int ) ncombs ( n-x, j-i ) );
            myprintf ( "   ncombs(%d, %d) = %d \n", ( int ) ( n-x ), ( int ) ( j-i ), ( int ) ncombs ( n-x, j-i ) );
        }
    }
    return val;
}

template <class Type>
/// calculate value of Krawtchouk polynomial
inline Type krawtchouksCache ( Type j, Type x, Type n )
{
    Type val=0;

    for ( Type i=0; i<=j; i++ ) {
        val += powmo ( i ) * ncombscache ( x, i ) * ncombscache ( n-x, j-i );
    }
    //printf("krawtchouks: %d, %d, %d, %ld\n", (int) j, (int)x, (int)n, val);
    return val;
}

template <class Type>
/// calculate value of Krawtchouk polynomial
inline Type krawtchouks ( Type j, Type x, Type n )
{
    Type val=0;

    for ( Type i=0; i<=j; i++ ) {
        val += powmo ( i ) * ncombs ( x, i ) * ncombs ( n-x, j-i );
    }
    //printf("krawtchouks: %d, %d, %d, %ld\n", (int) j, (int)x, (int)n, val);
    return val;
}

#include <Eigen/SVD>
#include <Eigen/Core>

template < class Type>
/// return the condition number of a matrix
double conditionNumber(const Eigen::Matrix<Type,-1,-1> A)
{
    Eigen::JacobiSVD<Eigen::Matrix<Type,-1,-1> > svd(A);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    return cond;
}

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
