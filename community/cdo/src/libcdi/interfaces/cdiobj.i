%module CdiObj
%{
#define SWIG_FILE_WITH_INIT
#include "cdi.hpp"
%}
%include "stl.i"
namespace std {
   %template(IntVector)    vector<int>;
   %template(DoubleVector) vector<double>;
   %template(DoubleDoubleVector) vector< vector<double> >;
   %template(StringVector) vector<string>;
   %template(VarsVector)   vector<CdiVariable>;
   %template(VarsMap)      map<string,CdiVariable>;
   %template(VarsByCode)   map<int,CdiVariable>;
   %template(TaxesMap)     map<int,CdiTaxis>;
   %template(ZaxesMap)     map<int,CdiZaxis>;
   %template(GridsMap)     map<int,CdiGrid>;
}
%include "cdi.hpp"
