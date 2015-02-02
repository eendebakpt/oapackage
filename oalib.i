/* File: example.i */
%module(docstring="Python Orthogonal Array interface") oalib

// see http://www.swig.org/Doc1.3/Python.html#Python_nn65
//%feature("autodoc", "docstring")

// basic features, see http://realmike.org/blog/2010/07/18/python-extensions-in-cpp-using-swig/
%include "std_string.i"
%include "std_vector.i"
%include "std_deque.i"

#define NEWINTERFACE

//%feature("shadow") array_link::getarray() %{
//def getarray(self, *args):
//  print('getting array: size %d %d' % (self.n_rows, self.n_columns))
//  x=self.getarraydata(self.n_rows*self.n_columns)
//  return x.reshape((self.n_columns, self.n_rows)).transpose()
//  #$action
//%}

//%template(_int_list) std::vector< int >;


//%feature("pythonprepend") array_link::clear() %{
//  print('gr')
//%}

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
#ifdef SWIGPYTHON
import_array();
#endif
%}

%apply (double* ARGOUT_ARRAY2, int DIM1, int DIM2) {(double* pymat2, int nrows, int ncols)}
%apply (int* ARGOUT_ARRAY2, int DIM1, int DIM2) {(int* pymat2, int nrows, int ncols)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* pymat1, int n)}
%apply (array_t* ARGOUT_ARRAY1, int DIM1) {(array_t* pymat1, int n)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* rangevec, int n)}
%apply (int* ARGOUT_ARRAY1, int DIM1, double xx) {(int* rangevec, int n, double xx)}
%apply (int* ARGOUT_ARRAY2, int DIM1, int DIM2) {(array_t* pymat2, int nrows, int ncols)}


%{
#ifdef WIN32
#else
#include <stdint.h>
#endif
#include "oaoptions.h"
#include "mathtools.h"
#include "arraytools.h"
#include "tools.h"
#include "md5.h"
#include "pareto.h"
#include "arrayproperties.h"
#include "extend.h"
#include "lmc.h"
#ifdef OADEV
#include "oadevelop.h"
#endif
%}

%extend array_link {
%insert("python") %{
def showarray(self):
  """ Show array"""
  # overridden to fix problems with ipython
  print(self.showarrayS() )
  
def getarray(self, verbose=0, *args):
  if verbose:
      print('getting array: size %d %d' % (self.n_rows, self.n_columns))
  x=self.getarraydata( int(self.n_rows*self.n_columns) )
  return x.reshape((self.n_columns, self.n_rows)).transpose()
  #$action
def setarray(self, X, verbose=0):
  self.init(X.shape[0], X.shape[1])
  self.index=-1
  iv = intVector(X.T.astype(int).flatten().tolist())
  #iv =_oalib.intVector(X.flatten().tolist())
  self.setarraydata(iv, X.size)
def __getitem__(self,index):
  if type(index)==int:
      if index<0 or index > self.n_rows*self.n_columns:
        raise IndexError('index out of bounds')
      return self.at(index)
  else:
      if len(index)==2:
        # FIXME: error checking
        a=index[0]
        b=index[1]
        if a<0 or a >= self.n_rows:
          raise IndexError('index out of bounds')
        if b<0 or b >= self.n_columns:
          raise IndexError('index out of bounds')
        return self.at(a, b)	  
      else:
        raise IndexError('invalid index')
def __setitem__(self,index, value):
  if type(index)==int:
      if index<0 or index > self.n_rows*self.n_columns:
        raise IndexError('index out of bounds')
      self.setvalue(index, 0, value) 
  else:
      if len(index)==2:
        a=index[0]
        b=index[1]
        if a<0 or a >= self.n_rows:
          raise IndexError('index out of bounds')
        if b<0 or b >= self.n_columns:
          raise IndexError('index out of bounds')
        self.setvalue(a, b, value)
      else:
        raise IndexError('invalid index')
%}
}


%feature("autodoc", "1");

// http://www.swig.org/Doc2.0/Python.html#Python_nn47
//%include "carrays.i"
//%array_class(int, intArray);

%include "cpointer.i"
//%pointer_functions(int, intp);


// ignore variable argument length functions
%ignore printfstring;   

// rename problem names
%rename(__lt__) ::operator<;
%rename(__gt__) ::operator>;

#pragma SWIG nowarn=454

namespace std {
   %template(arraylist_t) deque<array_link>; // arraylist_t
   %template(jstructArray) vector<jstruct_t>; // results list
   %template(intVector) vector<int>;
   %template(longVector) vector<long>;
   %template(longDeque) deque<long>;
   %template(doubleVector) vector<double>;
};


// prevent memory leaks
#ifdef OADEV

//%newobject extend_arraylist;
//arraylist_t & extend_arraylist(arraylist_t & alist, arraydata_t &fullad,   OAextend const &oaextend);
#endif


//%newobject readarrayfile;
//arraylist_t & readarrayfile(const char *fname, int verbose=0, int *setcols = 0); 

%newobject readarrayfile;
arraylist_t readarrayfile(const char *fname, int verbose=0, int *setcols = 0); 


// do this before the real arraylink is included...

#ifdef SWIGPYTHON
%pythoncode %{
import numpy
%}
#endif

%include "oaoptions.h"
%include "mathtools.h"
%include "arraytools.h"
%include "tools.h"
%include "arrayproperties.h"
%include "md5.h"
%include "pareto.h"
%include "extend.h"
%include "lmc.h"
#ifdef OADEV
%include "oadevelop.h"
#endif

%template(mvalue_t_long) mvalue_t<long>;
%template(mvalue_t_double) mvalue_t<double>;
%template(ParetoLong) Pareto<mvalue_t<long>,long>;
%template(ParetoMLong) Pareto<mvalue_t<long>,long>;
%template(ParetoDoubleLong) Pareto<double,long>;
%template(ParetoElementLong) pareto_element<mvalue_t<long>,long>;
%template(ParetoMElementLong) pareto_element<mvalue_t<long>,long>;
%template(mvalueVector) std::vector<mvalue_t<long> >;
%template(DequeParetoElementLong) std::deque<pareto_element<mvalue_t<long>,long> >;
//%template(DequeParetoLong) std::deque< Pareto< mvalue_t< long >,long >;
//%template(GWLPvalueVector2) std::vector<GWLPvalue>; 
%template(GWLPvalueVector) std::vector< mvalue_t<double> >;

/* representation functions */

%extend arraydata_t {
public:
    std::string __repr__() {
      return $self->showstr();
    }
} 

%extend array_link {
public:
    std::string __repr__() {
      return $self->showstr();
    }
} 

%extend arrayfile::arrayfile_t {
public:
    std::string __repr__() {
      return $self->showstr();
    }
} 

%extend std::deque<array_link> {
public:
    std::string __repr__() {
      return printfstring("arraylink array with %d elements", $self->size() ); // $self->showstr();
    }
}

%extend jstruct_t {
public:
    std::string __repr__() {
      return $self->showstr();
    }
} 


// Full Doxygen documentation
//%include "./swig_doc.i"


#ifdef SWIGPYTHON
// Add module docstring
%pythoncode 
%{
__doc__ = """
Python Orthogonal Array Interface 2
"""
%}
#endif



//%module darray
%inline %{

double iarray_get(int *a, int index) {
  return a[index];
}

// note: the size of the array is not easy to pass to the C function
// see: http://www.scipy.org/Cookbook/SWIG_Memory_Deallocation
// see also: http://stackoverflow.com/questions/2209395/in-python-how-to-access-a-uint163-array-wrapped-by-swig-i-e-unwrap-a-pyswigo

%}

