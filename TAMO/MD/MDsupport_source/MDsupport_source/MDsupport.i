// Copyright (2005) Whitehead Institute for Biomedical Research
// All Rights Reserved
// Author: David Benjamin Gordon

// For typemaps related to argument input (in,ignore,default,arginit,check), replace
// $source by $input and $target by $1.   For typemaps related to return values (out,
// argout,ret,except), replace $source by $1 and $target by $result.  See the file
// Doc/Manual/Typemaps.html for complete details.

%include stl.i
/* instantiate the required template specializations */
namespace std {
    %template(IntVector)    vector<int>;
    %template(DoubleVector) vector<double>;
}

// A typemap for handling any int [][] array
%typemap(in) double [ANY][ANY] {
  int i,j;
  for (i = 0; i < $dim0; i++) {
    for (j = 0; j < $dim1; j++) {
      $1[i][j] = *($input+$dim1*i+j); */
    }
  }
}

%typemap(python,in) dlist {
  /* Convert a Python List into an array of values */
  int i,sz;
  sz = PyList_Size($input);
  double *arr = new double[sz+1];
  for (i = 0; i < sz; i++) {
    arr[i] = (double) PyFloat_AsDouble(PyList_GetItem($input,i));
  }
  arr[i] = 1e20;
  $1 = arr;
}

%typemap(out) flist {
  int i,j;
  int len;
  PyObject *theFloat;
  PyObject *theList;
 
  /* Major Kludge -- crazy values determine end of array */
  for (len=0; $1[len]>-999; len++);
  theList = PyList_New(len);
  for (i=0; i<len; i++) {
    theFloat = PyFloat_FromDouble($1[i]);
    PyList_SetItem(theList, i, theFloat);
  }
  $result = theList;
};

%typemap(python, in) PyObject * {
  $1 = $input;
}

%typemap(python, out) PyObject* {
  $result = $1;
}


%pythoncode %{
def Motif2c_PSSM(motif):
    c_PSSM = SeqMat(motif.width)
    LjT = zip(['A','C','G','T'],[0,1,2,3])
    for i in range(motif.width):
        for L,j in LjT:
            c_PSSM.set(i,j, motif.ll[i][L])
    c_PSSM.compute_ambig()
    return c_PSSM
%}


%module MDsupport
%{
/* Put header files here (optional) */
/* #include "MDsupport.h" */
#include "MDsupport.h"
%}

%include "MDsupport.h"
