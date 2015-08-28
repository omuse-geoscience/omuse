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

/*
   This module contains the following operators:

      Wct     wct          Compute the windchill temperature (degree C)
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static const char WCT_NAME[]     = "wind_chill_temperature";
static const char WCT_LONGNAME[] = "Windchill temperature describes the fact that low temperatures are felt to be even lower in case of wind. It is based on the rate of heat loss from exposed skin caused by wind and cold. It is calculated according to the empirical formula: 33 + (T - 33) * (0.478 + 0.237 * (SQRT(ff*3.6) - 0.0124 * ff * 3.6)) with T  = air temperature in degree Celsius, ff = 10 m wind speed in m/s. Windchill temperature is only defined for temperatures at or below 33 degree Celsius and wind speeds above 1.39 m/s. It is mainly used for freezing temperatures.";
static const char WCT_UNITS[]    = "Celsius";

static const int FIRST_VAR = 0;


static double windchillTemperature(double t, double ff, double missval)
{
  static const double tmax = 33.0; 
  static const double vmin = 1.39; /* minimum wind speed (m/s) */
  
  return ff < vmin || t > tmax ? missval : tmax + (t - tmax) * (0.478 + 0.237 * (sqrt(ff * 3.6) - 0.0124 * ff * 3.6));
}


static void farexpr(field_t *field1, field_t field2, double (*expression)(double, double, double))
{
  int   i, len;
  const int     grid1    = field1->grid;
  const int     nmiss1   = field1->nmiss;
  const double  missval1 = field1->missval;
  double       *array1   = field1->ptr;
  const int     grid2    = field2.grid;
  const int     nmiss2   = field2.nmiss;
  const double  missval2 = field2.missval;
  const double *array2   = field2.ptr;

  len = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) )  
	  array1[i] = missval1;
	else
	  array1[i] = expression(array1[i], array2[i], missval1);
    }
  else
    {
      for ( i = 0; i < len; i++ )
        array1[i] = expression(array1[i], array2[i], missval1);  
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}

   
void *Wct(void *argument)
{
  int streamID1, streamID2, streamID3;
  int gridsize;
  int nrecs, nrecs2, recID;
  int tsID;
  int gridID, zaxisID;
  int varID1, varID2, varID3;
  int levelID1, levelID2;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  field_t field1, field2;

  cdoInitialize(argument);
  cdoOperatorAdd("wct", 0, 0, NULL);

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);

  vlistCompare(vlistID1, vlistID2, CMP_DIM);
  
  gridsize = vlistGridsizeMax(vlistID1);
  
  field_init(&field1);
  field_init(&field2);

  field1.ptr = (double*) malloc(gridsize*sizeof(double));
  field2.ptr = (double*) malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d", vlistNtsteps(vlistID1), vlistNtsteps(vlistID2));

  vlistID3 = vlistCreate();
  gridID   = vlistInqVarGrid(vlistID1, FIRST_VAR);
  zaxisID  = vlistInqVarZaxis(vlistID1, FIRST_VAR);
  varID3   = vlistDefVar(vlistID3, gridID, zaxisID, TSTEP_INSTANT);

  taxisID3 = taxisCreate(TAXIS_RELATIVE);
  taxisDefTunit(taxisID3, TUNIT_MINUTE);
  taxisDefCalendar(taxisID3, CALENDAR_STANDARD);
  taxisDefRdate(taxisID3, 19550101);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  vlistDefVarName(vlistID3, varID3, WCT_NAME);
  vlistDefVarLongname(vlistID3, varID3, WCT_LONGNAME);
  vlistDefVarUnits(vlistID3, varID3, WCT_UNITS);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      nrecs2 = streamInqTimestep(streamID2, tsID);
      if ( nrecs2 == 0 )
        cdoAbort("Input streams have different number of timesteps!");

      taxisCopyTimestep(taxisID3, taxisID1);
      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID1, &levelID1);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  streamInqRecord(streamID2, &varID2, &levelID2);
	  streamReadRecord(streamID2, field2.ptr, &field2.nmiss);
	  
	  if ( varID1 != varID2 || levelID1 != levelID2 )
	    cdoAbort("Input streams have different structure!");
	    
          if ( varID1 != FIRST_VAR )
            continue;
            
	  field1.grid    = vlistInqVarGrid(vlistID1, varID1);
	  field1.missval = vlistInqVarMissval(vlistID1, varID1);

	  field2.grid    = vlistInqVarGrid(vlistID2, varID2);
	  field2.missval = vlistInqVarMissval(vlistID2, varID2);

	  farexpr(&field1, field2, windchillTemperature);
	  
	  streamDefRecord(streamID3, varID3, levelID1);
	  streamWriteRecord(streamID3, field1.ptr, field1.nmiss);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
