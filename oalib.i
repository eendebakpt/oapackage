/* File: example.i */
%module(docstring="Python Orthogonal Array interface") oalib


// basic features, see http://realmike.org/blog/2010/07/18/python-extensions-in-cpp-using-swig/
%include "std_string.i"
%include "std_vector.i"
%include "std_deque.i"

#define NEWINTERFACE

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

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* pymat1, int nrows)}
%apply (double* ARGOUT_ARRAY2, int DIM1, int DIM2) {(double* pymat2, int nrows, int ncols)}
%apply (int* ARGOUT_ARRAY2, int DIM1, int DIM2) {(int* pymat2, int nrows, int ncols)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* pymat1, int n)}
%apply (array_t* ARGOUT_ARRAY1, int DIM1) {(array_t* pymat1, int n)}
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* rangevec, int n)}
%apply (int* ARGOUT_ARRAY1, int DIM1, double xx) {(int* rangevec, int n, double xx)}
%apply (int* ARGOUT_ARRAY2, int DIM1, int DIM2) {(array_t* pymat2, int nrows, int ncols)}

// enable keyword arguments in interface
//%feature ("kwargs")

%include "std_pair.i"

%{
#include <utility>
#include <Eigen/Core>
#include <Eigen/Dense>
//#include <Eigen/Dense>

//#include <Python.h>
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
#ifdef OADEV
#include "oadevelop.h"
#endif
%}

/* Instantiate a few different versions of the template */
//%template(EigenMatrix) Eigen::MatrixXd;

/*
%pythoncode %{ 
def eigen2numpy(m):
  print('convert eigen to numpy matrix')
  
  m.rows()

%}
*/



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
          // std::cout << "data " << data[i] << std::endl;
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

def reduceGraphNauty(G, colors, verbose=1):
  """ Reduce vertex transformation reducing array to normal form
  
  The reduction is calculated using `Nauty <http://pallini.di.uniroma1.it/>`_

  Arguments
  ---------
      G : Numpy array
	the graph in incidence matrix form
      colors : list
	an optional vertex coloring
  """
  
  al=array_link()
  al.setarray(G)
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




%extend array_link {
%insert("python") %{
#__array_interface__ = None

def __getattr__(self, attr):
    if attr=='__array_interface__':
      a = dict()
      a['version']=3
      a['shape']=(self.n_rows, self.n_columns)
      a['typestr']='<i2'
      a['data']=(self.data(), False)
      # convert from the OAP column-major style to Numpy row-major style
      a['strides']=(2, 2*self.n_rows)
      return a
    else:
      raise AttributeError("%r object has no attribute %r" %
                         (self.__class__, attr))
                         
def showarray(self):
  """ Show array"""
  # overridden to fix problems with ipython
  #print(self.showarrayS(), end='',flush=True)	# not valid in python2
  sys.stdout.write(self.showarrayS())

def getarray(self, verbose=0, *args):
  """ Return Numpy style array """
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
  """ Return element of array """
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

// see http://www.swig.org/Doc1.3/Python.html#Python_nn65
//%feature("autodoc", "docstring")


%feature("autodoc", "1");
// to generate the oadoxy.i:
// doxygen Doxyfile
// python2 doxy2swig.py xml/index.xml oadoxy.i
// see also: http://www.enricozini.org/2007/tips/swig-doxygen-docstring/

%include "oadoxy.i"


// http://www.swig.org/Doc2.0/Python.html#Python_nn47

%include "cpointer.i"
%include "std_map.i"


// ignore variable argument length functions
%ignore printfstring;   

// rename problem names
%rename(__lt__) ::operator<;
%rename(__gt__) ::operator>;

#pragma SWIG nowarn=454

namespace std {
   %template(arraylist_t) deque<array_link>; // arraylist_t
   %template(jstructArray) vector<jstruct_t>; // results list
   %template(uint8Vector) std::vector<unsigned char>;
   %template(charVector) std::vector<signed char>;
   %template(intVector) std::vector<int>;
   %template(longVector) std::vector<long>;
   %template(longDeque) deque<long>;
   %template(doubleVector) std::vector<double>;
   %template(stringVector) std::vector<std::string>;
   
   %template(map_int_long) std::map<int, long>;
   
};


// prevent memory leaks

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

//%include "printfheader.h"
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
#ifdef OADEV
%include "oadevelop.h"
#endif

%template(pairDoptimize) std::pair< std::vector< std::vector<double> > ,arraylist_t>;
%template(pairGraphColors) std::pair< array_link  , std::vector<int>  >;

%template(mvalue_t_long) mvalue_t<long>;
%template(mvalue_t_double) mvalue_t<double>;
%template(ParetoLong) Pareto<mvalue_t<long>,long>;
%template(ParetoMLong) Pareto<mvalue_t<long>,long>;
%template(ParetoDoubleLong) Pareto<double,long>;
%template(ParetoElementLong) pareto_element<mvalue_t<long>,long>;
%template(ParetoMElementLong) pareto_element<mvalue_t<long>,long>;
%template(mvalueVector) std::vector<mvalue_t<long> >;
%template(DequeParetoElementLong) std::deque<pareto_element<mvalue_t<long>,long> >;
//%template(GWLPvalueVector2) std::vector<GWLPvalue>; 
%template(GWLPvalueVector) std::vector< mvalue_t<double> >;

%template(cpermVector) std::vector< cperm >;

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
      return printfstring("array_transformation_t: transformation for array of size %d x %d", $self->ad->N, $self->ad->ncols);
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

