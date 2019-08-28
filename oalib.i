/* File: oalib.i
 * 
 * Defines the Python interface to the OApackage
 *
 */
%module(docstring="Python Orthogonal Array interface") oalib


// basic features, see http://realmike.org/blog/2010/07/18/python-extensions-in-cpp-using-swig/
%include "std_string.i"
%include "std_vector.i"
%include "std_deque.i"
%include "std_pair.i"
%include "exception.i"

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply ( long* IN_ARRAY2, int DIM1, int DIM2 ) { (long* pymatinput, int nrows, int ncols) }
%apply ( double* IN_ARRAY2, int DIM1, int DIM2 ) { (double* pymatdoubleinput, int nrows, int ncols) }

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* pymat1, int nrows)}
%apply (int* ARGOUT_ARRAY2, int DIM1, int DIM2) {(int* pymat2, int nrows, int ncols)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* pymat1, int n)}
%apply (array_t* ARGOUT_ARRAY1, int DIM1) {(array_t* pymat1, int n)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* rangevec, int n)}

%{
#include <utility>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <numpy/arrayobject.h>

#ifdef _WIN32
#else
#include <stdint.h>
#endif
#include "printfheader.h"
#include "oaoptions.h"
#include "mathtools.h"
#include "arraytools.h"
#include "tools.h"
#include "md5.h"
#include "pareto.h"
#include "arrayproperties.h"
#include "extend.h"
#include "lmc.h"
#include "Deff.h"
#include "graphtools.h"
#include "evenodd.h"
#include "conference.h"
#include "unittests.h"
#ifdef OADEV
#include "oadevelop.h"
#endif
%}

%typemap(in) Eigen::MatrixXd (Eigen::MatrixXd inputEigen)
{
  /* note that Eigen is column-major by default and numpy is row major by default */
  int rows = 0;
  int cols = 0;

  PyArrayObject *pp = (PyArrayObject *)($input);
  rows = PyArray_DIM(pp,0);
  cols = PyArray_DIM(pp,1);

  PyArrayObject* temp;
  PyArg_ParseTuple($input, "O", &temp);  

  inputEigen.resize(rows,cols);
  inputEigen.fill(0);

  double *  values = ((double *) PyArray_DATA( pp ));
  for (long int i = 0; i != rows; ++i){
      for(long int j = 0; j != cols; ++j){
          inputEigen(i,j) = values[i*rows+j];
      }
  }  

}

// see http://sourceforge.net/p/swig/mailman/message/32490448/
//http://mail.scipy.org/pipermail/numpy-discussion/2013-February/065637.html
%typemap(out) Eigen::VectorXd 
{
    npy_intp dims[1] = {$1.size()};
    PyObject* array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double* data = ((double *) (PyArray_DATA( (PyArrayObject*) array ) ) );
    for (int i = 0; i != dims[0]; ++i){
        *data++ = $1.data()[i];
    }
    $result = array;
}
%typemap(out) Eigen::MatrixXd 
{
  /* note that Eigen is column-major by default and numpy is row major by default */

  const int verbose=0;
  if (verbose) {
    printf("typemap out for Eigen::MatrixXd: \n");
    eigenInfo($1);
    }
    Eigen::MatrixXd mt = $1.transpose();
    
    npy_intp dims[2] = {$1.rows(), $1.cols() };
    PyObject* array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double* data = ((double *) (PyArray_DATA( (PyArrayObject*) array ) ) );
    for (int i = 0; i != dims[0]*dims[1]; ++i){
        *data++ = mt.data()[i];
    }
    $result = array;
}

%pythoncode %{
import sys
import numpy as np
import copy

def reduceGraphNauty(G, colors=None, verbose=1):
  """ Return vertex transformation reducing array to normal form
  
  The reduction is calculated using `Nauty <http://users.cecs.anu.edu.au/~bdm/nauty/>`_

  Args:
      G (numpy array or array_link) :	the graph in incidence matrix form
      colors (list or None): an optional vertex coloring
  Returns:
      v: relabelling of the vertices
  """
  if isinstance(G, np.ndarray):
      al=array_link()
      al.setarray(G)
  else:
      al = copy.copy(G)
  if colors is None:
    colors = [0] * G.shape[0]
  v = _oalib.reduceNauty ( al, colors, verbose )
  return v
  
def transformGraphMatrix(G, tr, verbose=1):
    """ Apply a vertex permutation to a graph 
    
    Arguments
    ---------
      G : Numpy array
	the graph in incidence matrix form
      tr : list 
	the vertex transformation as a list
	
    Returns
    -------
      The transformed graph
      
    """
    al=array_link()
    al.setarray(G)
    alt = _oalib.transformGraph(al, tr, verbose)  
    return np.array(alt)
  
%}


/* Convert from C --> Python */
%typemap(out) void * {
    $result = PyLong_FromVoidPtr($1);
}

%extend array_link {
%insert("python") %{

def __getattr__(self, attr):
    if attr=='__array_interface__':
      a = dict()
      a['version']=3
      a['shape']=(self.n_rows, self.n_columns)
      sizeofdata=_oalib.sizeof_array_t()
      a['typestr']='<i%d' % sizeofdata # sizeof(array_t)
      a['data']=(self.data(), True)
      # convert from the OAP column-major style to Numpy row-major style
      a['strides']=(sizeofdata, sizeofdata*self.n_rows)
      return a
    else:
      raise AttributeError("%r object has no attribute %r" %
                         (self.__class__, attr))

@property
def shape(self):
    return (self.n_rows, self.n_columns)

@property
def size(self):
    return self.n_rows*self.n_columns
                             
def showarray(self):
  """ Show array """
  # overridden to fix problems with ipython
  sys.stdout.write(self.showarrayString())

def getarray(self, verbose=0, *args):
  """ Return Numpy style array """
  if verbose:
      print('getting array: size %d %d' % (self.n_rows, self.n_columns))
  x=self.getarraydata( int(self.n_rows*self.n_columns) )
  return x.reshape((self.n_columns, self.n_rows)).transpose()

def setarray(self, X, verbose=0):
  """ Update the array link object with a Numpy array

  Args:
     X (numpy array): array to be copied to the object
  """
  self.init(X.shape[0], X.shape[1])
  self.index=-1
  iv = intVector(X.T.astype(int).flatten().tolist())
  self.setarraydata(iv, X.size)

def _slice2range(self, slice, max_value):
    """ Convert a python slice object to a range """
    if isinstance(slice, int):
        return [slice]
    if slice.start  is None:
        start = 0
    else:
        start = slice.start
    if slice.stop is None:
        stop = max_value
    else:
        stop = slice.stop
    if slice.step is None:
        step = 1
    else:
        step = slice.step
    return list(range(start, stop, step))

def _ranges2subarray(self, row_range, col_range):
      """ From a list of row element and a list of column element construct a submatrix """
      al=array_link(len(row_range), len(col_range), array_link.INDEX_DEFAULT )
      for ii, row in enumerate(row_range):
          for jj, col in enumerate(col_range):
              al[ii, jj]=self.at(row, col)
      return al

def __getitem__(self, index):
  """ Return element of array """
  if type(index)==int:
      if index<0 or index > self.n_rows*self.n_columns:
        raise IndexError('index out of bounds')
      return self.at(index)
  elif isinstance(index, slice):
      indices=self._slice2range(index, self.n_rows*self.n_columns)
      return np.array( [self.at(a) for a in indices])
  else:
      if len(index)==2:
          index0=index[0]
          index1=index[1]
          if isinstance(index0, int) and isinstance(index1, int):
            if index0<0 or index0 >= self.n_rows:
              raise IndexError('index out of bounds')
            if index1<0 or index1 >= self.n_columns:
              raise IndexError('index out of bounds')
            return self.at(index0, index1)	  
          elif isinstance(index0, int) and isinstance(index1, slice):
              row_range=[index0]
              col_range=self._slice2range(index1, self.n_columns)
              
              return self._ranges2subarray(row_range, col_range)
          elif isinstance(index0, slice) and isinstance(index1, int):
              row_range=self._slice2range(index0, self.n_rows)
              col_range=[index1]
              
              return self._ranges2subarray(row_range, col_range)
          elif isinstance(index0, slice) and isinstance(index1, slice):
              row_range=self._slice2range(index0, self.n_rows)
              col_range=self._slice2range(index1, self.n_columns)
              
              return self._ranges2subarray(row_range, col_range)
          else:
              raise NotImplementedError('slice indexing not supported')
      else:
        raise IndexError('invalid index')

def __setitem__(self, index, value):
  """ Set specified value at specified index in the array """
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

// see http://www.swig.org/Doc1.3/Python.html#Python_nn65
//%feature("autodoc", "docstring")


%feature("autodoc", "1");
// to generate the oadoxy.i:
// doxygen Doxyfile; python doxy2swig.py docs/xml/index.xml oadoxy.i
// see also: http://www.enricozini.org/2007/tips/swig-doxygen-docstring/

%include "oadoxy.i"


// see http://www.swig.org/Doc2.0/Python.html#Python_nn47

%include "cpointer.i"
%include "std_map.i"


// ignore variable argument length functions
%ignore printfstring;   

// rename problem names
%rename(__lt__) ::operator<;
%rename(__gt__) ::operator>;

#pragma SWIG nowarn=454

namespace std {
   %template(arraylist_t) deque<array_link>; 
   %template(jstructArray) vector<jstruct_t>; 
   %template(uint8Vector) std::vector<unsigned char>;
   %template(charVector) std::vector<signed char>;
   %template(intVector) std::vector<int>;
   %template(longVector) std::vector<long>;
   %template(longDeque) deque<long>;
   %template(doubleVector) std::vector<double>;
   %template(stringVector) std::vector<std::string>;
   
   %template(map_int_long) std::map<int, long>;
   
};

%exception array_link::selectFirstColumns {
  try {
    $action
  } catch (std::runtime_error& e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(e.what()));
  }
}
%exception mycheck_handler {
  try {
    $action
  } catch (std::runtime_error& e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(e.what()));
  }
}
%exception throw_runtime_exception {
  try {
    $action
  } catch (std::runtime_error& e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(e.what()));
  }
}

%exception {
  try {
    $action
  } 
  SWIG_CATCH_STDEXCEPT // catch std::exception
  catch (...) {
     SWIG_exception_fail(SWIG_UnknownError, "Unknown exception");
  }
}

// prevent memory leaks
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
%include "Deff.h"
%include "pareto.h"
%include "extend.h"
%include "lmc.h"
%include "Deff.h"
%include "graphtools.h"
%include "evenodd.h"
%include "conference.h"
%include "unittests.h"
#ifdef OADEV
%include "oadevelop.h"
#endif

%template(pairDoptimize) std::pair< std::vector< std::vector<double> > ,arraylist_t>;
%template(pairGraphColors) std::pair< array_link  , std::vector<int>  >;
%template(pairEigenMatrix) std::pair< MatrixFloat  , MatrixFloat >;


%extend mvalue_t<double> {
%insert("python") %{

def __getattr__(self, attr):
    if attr=='__array_interface__':
      a = dict()
      a['version']=3
      a['shape']=(self.size(), )
      sizeofdata=_oalib.sizeof_double()
      a['typestr']='<f%d' % sizeofdata
      a['data']=(np.array(self.values), True)
      return a
    else:
      raise AttributeError("%r object has no attribute %r" %
                         (self.__class__, attr))
                         
%}
}

%template(mvalue_t_long) mvalue_t<long>;
%template(mvalue_t_double) mvalue_t<double>;
%template(ParetoLongLong) Pareto<long,long>;
%template(ParetoMultiLongLong) Pareto<mvalue_t<long>,long>;
%template(ParetoMultiDoubleLong) Pareto<mvalue_t<double>,long>;
%template(ParetoDoubleLong) Pareto<double,long>;
%template(ParetoElementLong) pareto_element<mvalue_t<long>,long>;
%template(ParetoMElementLong) pareto_element<mvalue_t<long>,long>;
%template(vector_mvalue_t_double) std::vector<mvalue_t<double> >;
%template(vector_mvalue_t_int) std::vector<mvalue_t<int> >;
%template(vector_mvalue_t_long) std::vector<mvalue_t<long> >;
%template(DequeParetoElementLong) std::deque<pareto_element<mvalue_t<long>,long> >;


%pythoncode %{
# for legacy reasons and for name consistency
GWLPvalueVector = vector_mvalue_t_double
mvalueVector = vector_mvalue_t_long
#%template(mvalueVector) std::vector<mvalue_t<long> >;
%}


%template(conference_columnVector) std::vector< conference_column >;
%template(calculateArrayParetoJ5) calculateArrayParetoJ5<array_link>;
%template(calculateArrayParetoJ5int) calculateArrayParetoJ5<int>;
%template(calculateArrayParetoJ5long) calculateArrayParetoJ5<long>;
%template(vector_vector_double) std::vector< std::vector<double> >;

/* representation functions */

%extend arraydata_t {
public:
    std::string __repr__() {
      return $self->showstr();
    }
} 

%extend array_transformation_t {
public:
    std::string __repr__() {
      if($self->ad!=0)
	return printfstring("array_transformation_t: transformation for array of size %d x %d", $self->ad->N, $self->ad->ncols);
      else
	return printfstring("array_transformation_t: no class defined");      
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
      return printfstring("list of array_link objects with %d elements", $self->size() ); 
    }
}

%extend jstruct_t {
public:
    std::string __repr__() {
      return $self->showstr();
    }
} 

#ifdef SWIGPYTHON
// Add module docstring
%pythoncode  
%{
__doc__ = """
Python Orthogonal Array Interface 
"""
%}
#endif


