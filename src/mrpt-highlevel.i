//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

// SWIG interface file to generate python module from
%module mrpt 

%{
#define SWIG_FILE_WITH_INIT
#include "mrpt-highlevel.h"
%}

%include "typemaps.i"
//get numpy typemaps
%include "numpy.i"

%init %{
import_array();
%}

//apply numpy typemaps to array accessor functions
//these typemaps don't seem to work with unsigned instead of int as DIM_TYPE
%apply (int* DIM1, double** ARGOUTVIEW_ARRAY1) {(int *outN_k, double** outArray1)};
%apply (int* DIM1, int** ARGOUTVIEW_ARRAY1) {(int *outN_k, int** outArray1)};
%apply (int* DIM1, double** ARGOUTVIEW_ARRAY1) {(int* outM, double** outArray1)};
%apply (int* DIM1, double** ARGOUTVIEW_ARRAY1) {(int* outK, double** outArray1)};
%apply (int* DIM1, int* DIM2, int** ARGOUTVIEW_ARRAY2 ) {(int* outK, int* outM, int** outArray2)};
%apply (int* DIM1, int** ARGOUTVIEW_ARRAY1) {(int *outM, int** outArray1)};
%apply (int* DIM1, int* DIM2, double** ARGOUTVIEW_ARRAY2 ) {(int* outK, int* outM, double** outArray2)};
%apply (int* DIM1, int* DIM2, double** ARGOUTVIEW_ARRAY2 ) {(int* outL, int* outM, double** outArray2)};
%apply (int* DIM1, int* DIM2, int** ARGOUTVIEW_ARRAY2 ) {(int* outK, int* outL, int** outArray2)};
%apply (int* DIM1, int* DIM2, int** ARGOUTVIEW_ARRAY2 ) {(int* outL, int* outM, int** outArray2)};

//this is for
// double entropyDifference(double,double, double&)
//to turn the double& into a second return value
%apply double* OUTPUT {double& outIntegrationError};


//SWIG-typemap to use a list of (utf-8) python strings instead of argc, argv
%typemap(in) (int argc, char **argv) {
  /* Check if is a list */
    if (PyList_Check($input)) {
        int i;
        $1 = PyList_Size($input);
        $2 = (char **) malloc(($1+1)*sizeof(char *));
        for (i = 0; i < $1; i++) {
            PyObject *o = PyList_GetItem($input,i);
            if (PyString_Check(o)) {
                $2[i] = PyString_AsString(PyList_GetItem($input,i));
            } else if (PyUnicode_Check(o)) {
                //utf-8 support (assumes char == 1 byte)
                $2[i] = (char *)PyUnicode_AsEncodedString(o, "utf-8", "Error ~");
                $2[i] = (char *)PyBytes_AS_STRING($2[i]);
            } else {
                PyErr_SetString(PyExc_TypeError,"list must contain (unicode) strings");
                free($2);
                return NULL;
            }
        }
        $2[i] = 0;
    } else {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
    }
}

%typemap(freearg) (int argc, char **argv) {
  free((char *) $2);
}


%include "mrpt-highlevel.h"
