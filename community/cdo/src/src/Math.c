/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Math       abs             Absolute value
      Math       sqr             Square
      Math       sqrt            Square root
      Math       exp             Exponential
      Math       ln              Natural logarithm
      Math       log10           Base 10 logarithm
      Math       sin             Sine
      Math       cos             Cosine
      Math       tan             Tangent
      Math       asin            Arc sine
      Math       acos            Arc cosine
      Math       atan            Arc tangent
      Math       pow             Power
      Math       reci            Reciprocal
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Math(void *argument)
{
  enum {ABS, FINT, FNINT, SQR, SQRT, EXP, LN, LOG10, SIN, COS, TAN, ASIN, ACOS, ATAN, POW, RECI};
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int gridsize;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss, nmiss2;
  int i;
  int number;
  double missval1;
  double *array1, *array2;
  double rc = 0;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("abs",   ABS,   0, NULL);
  cdoOperatorAdd("int",   FINT,  0, NULL);
  cdoOperatorAdd("nint",  FNINT, 0, NULL);
  cdoOperatorAdd("sqr",   SQR,   0, NULL);
  cdoOperatorAdd("sqrt",  SQRT,  0, NULL);
  cdoOperatorAdd("exp",   EXP,   0, NULL);
  cdoOperatorAdd("ln",    LN,    0, NULL);
  cdoOperatorAdd("log10", LOG10, 0, NULL);
  cdoOperatorAdd("sin",   SIN,   0, NULL);
  cdoOperatorAdd("cos",   COS,   0, NULL);
  cdoOperatorAdd("tan",   TAN,   0, NULL);
  cdoOperatorAdd("asin",  ASIN,  0, NULL);
  cdoOperatorAdd("acos",  ACOS,  0, NULL);
  cdoOperatorAdd("atan",  ATAN,  0, NULL);
  cdoOperatorAdd("pow",   POW,   0, NULL);
  cdoOperatorAdd("reci",  RECI,  0, NULL);
 
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  if ( operfunc == POW )
    {
      operatorInputArg("value");
      rc = parameter2double(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  missval1 = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          number   = vlistInqVarNumber(vlistID1, varID);

          if ( number == CDI_REAL )
            {
              switch ( operfunc )
                {
                case ABS:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : fabs(array1[i]);
                  break;
                case FINT:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : (int)(array1[i]);
                  break;
                case FNINT:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : round(array1[i]);
                  break;
                case SQR:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : array1[i]*array1[i];
                  break;
                case SQRT:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : SQRT(array1[i]);
                  break;
                case EXP:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : exp(array1[i]);
                  break;
                case LN:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 ? missval1 : log(array1[i]);
                  break;
                case LOG10:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 ? missval1 : log10(array1[i]);
                  break;
                case SIN:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : sin(array1[i]);
                  break;
                case COS:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : cos(array1[i]);
                  break;
                case TAN:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : tan(array1[i]);
                  break;
                case ASIN:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < -1
                              || array1[i] > 1 ? missval1 : asin(array1[i]);
                  break;
                case ACOS:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < -1
                              || array1[i] > 1 ? missval1 : acos(array1[i]);
                  break;
                case ATAN:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : atan(array1[i]);
                  break;
                case POW:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : pow(array1[i], rc);
                  break;
                case RECI:
                  for ( i = 0; i < gridsize; i++ )
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array1[i], 0.) ? missval1 : 1/array1[i];
                  break;
                default:
                  cdoAbort("operator not implemented!");
                  break;
                }
            }
          else
            {
              switch ( operfunc )
                {
                case SQR:
                  for ( i = 0; i < gridsize; i++ )
                    {
                      array2[i*2]   = array1[i*2]*array1[i*2] + array1[i*2+1]*array1[i*2+1];
                      array2[i*2+1] = 0;
                    }
                  break;
                default:
                  cdoAbort("operator not implemented for complex numbers!");
                  break;
                }
            }

          nmiss2 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval1) ) nmiss2++;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss2);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
