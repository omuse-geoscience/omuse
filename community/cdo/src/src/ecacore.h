/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef ECA_H_
#define ECA_H_


typedef enum {
  NONE,
  MEAN,
  PERCENT_OF_TIME,
  PERCENT_OF_TOTAL_AMOUNT
}
ECA_EPILOG;

typedef void (*ECA_FUNC_1)(field_t *, double);
typedef void (*ECA_FUNC_2)(field_t *, field_t);
typedef void (*ECA_FUNC_3)(field_t *, field_t, double);

/**
 * Structure describing a processing request of the form
 * 
 * o = F3(F2(a * F1(i) + b))
 * 
 * where i and o denote the input and output fields, and
 * F1, F2 and F3 are field operators.
 * 
 * The structure contains the following elements:
 * 
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * f1        the 1st field operator
 * f1arg     the argument of the 1st field operator
 * f2        the 2nd field operator
 * f3        the 3rd field operator
 * mulc      the multiplier a
 * addc      the addend b
 * epilog    the final operation carried out after processing
 */ 
typedef struct {
  const char *name;      
  const char *longname;  
  const char *units;   
  ECA_FUNC_1  f1;
  double      f1arg;
  ECA_FUNC_2  f2;
  ECA_FUNC_2  f3;
  double      mulc;
  double      addc;
  ECA_EPILOG  epilog;
}
ECA_MAJOR_REQUEST_ELEMENT_1;

/**
 * Structure describing a processing request of the form
 * 
 * o = H3(H2(H1(i)))
 * 
 * where i and o denote the input and output fields, and
 * H1, H2 and H3 are field operators.
 * 
 * The structure contains the following elements:
 * 
 * name      the name of the output variable
 * longname  the longname of the output variable
 * h1        the 1st field operator
 * h1arg     the argument of the 1st field operator
 * h2        the 2nd field operator
 * h3        the 3rd field operator
 */ 
typedef struct {
  const char *name;
  const char *longname;
  const char *units;
  ECA_FUNC_1  h1;
  double      h1arg;
  ECA_FUNC_2  h2;
  ECA_FUNC_2  h3;
}
ECA_MINOR_REQUEST_ELEMENT_1;


typedef struct {
  ECA_MAJOR_REQUEST_ELEMENT_1 var1;
  ECA_MINOR_REQUEST_ELEMENT_1 var2;
}
ECA_REQUEST_1;


/**
 * Structure describing a processing request of the form
 * 
 * o = F5(F4(F3(F1(i1), F2(i2))))
 * 
 * where i1, i2 and o denote the input and output fields,
 * and F1, F2, F3, F4 and F3 are field operators.
 * 
 * The structure contains the following elements:
 * 
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * f1        the 1st field operator
 * f1arg     the argument of the 1st field operator
 * f2        the 2nd field operator
 * f2arg     the argument of the 2nd field operator
 * f3        the 3rd field operator
 * f4        the 4th field operator
 * f5        the 5th field operator
 * f5arg     the argument of the 5th field operator
 * epilog    the final operation carried out after processing
 */ 
typedef struct {
  const char *name;
  const char *longname;
  const char *units;
  ECA_FUNC_1  f1;
  double      f1arg;
  ECA_FUNC_1  f2;
  double      f2arg;
  ECA_FUNC_2  f3;
  ECA_FUNC_2  f4;
  ECA_FUNC_3  f5;
  double      f5arg;
  ECA_EPILOG  epilog;
}
ECA_MAJOR_REQUEST_ELEMENT_2;


/**
 * Structure describing a processing request of the form
 * 
 * o = H2(H1(i))
 * 
 * where i and o denote the input and output fields, and
 * H1, and H2 are field operators.
 * 
 * The structure contains the following elements:
 * 
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * h1        the 1st field operator
 * h1arg     the argument of the 1st field operator
 * h2        the 2nd field operator
 */ 
typedef struct {
  const char *name;
  const char *longname;
  const char *units;
  ECA_FUNC_1  h1;
  double      h1arg;
  ECA_FUNC_2  h2;
}
ECA_MINOR_REQUEST_ELEMENT_2;


typedef struct {
  ECA_MAJOR_REQUEST_ELEMENT_2 var1;
  ECA_MINOR_REQUEST_ELEMENT_2 var2;
}
ECA_REQUEST_2;


/**
 * Structure describing a processing request of the form
 * 
 * o = F3(F1(i1), F2(i2))
 * 
 * where i1, i2 and o denote the input and output fields,
 * and F1, F2 and F3 are field operators.
 * 
 * The structure contains the following elements:
 * 
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * f1        the 1st field operator
 * f2        the 2nd field operator
 * f3        the 3rd field operator
 */ 
typedef struct {
  const char *name;
  const char *longname;
  const char *units;
  ECA_FUNC_2  f1;
  ECA_FUNC_2  f2;
  ECA_FUNC_2  f3;
}
ECA_REQUEST_3;


/**
 * Structure describing a GSL-like processing request. The structure
 * contains the following elements:
 * 
 * name       the name of the 1st output variable
 * longname   the longname of the 1st output variable
 * units      the units of the 1st output variable
 * name2      the name of the 2nd output variable
 * longname2  the longname of the 2nd output variable
 * units2     the units of the 2nd output variable
 * name3      the name of the 3rd output variable
 * longname3  the longname of the 3rd output variable
 * units3     the units of the 3rd output variable
 * s1         the 1st field selector
 * s1arg      argument of the 1st field selector
 * s2         the 2nd field selector
 * s2arg      argument of the 2nd field selector
 * consecutiveDays  the number od concecutive days
 */ 
typedef struct {
  const char *name;
  const char *longname;
  const char *units;
  const char *name2;
  const char *longname2;
  const char *units2;
  ECA_FUNC_1  s1;
  double      s1arg;
  ECA_FUNC_1  s2;
  double      s2arg;
  ECA_FUNC_1  s3;
  double      s3arg;
  int         consecutiveDays;
}
ECA_REQUEST_4;


/**
 * Function processing a request of type 1.
 * 
 * @param request the processing request
 */ 
void eca1(const ECA_REQUEST_1 *request);

/**
 * Function processing a request of type 2.
 * 
 * @param request the processing request
 */ 
void eca2(const ECA_REQUEST_2 *request);

/**
 * Function processing a request of type 3.
 * 
 * @param request the processing request
 */ 
void eca3(const ECA_REQUEST_3 *request);

/**
 * Function processing a request of type 4.
 * 
 * @param request the processing request
 */ 
void eca4(const ECA_REQUEST_4 *request);

#endif /*ECA_H_*/
