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

      Mastrfu    mastrfu         Mass stream function
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


static
void mastrfu(int gridID, int zaxisID, double *array1, double *array2, int nmiss, double missval)
{
  int ilev, ilat, n;
  double fact =  4*atan(1.0) * 6371000 / 9.81;
  char units[CDI_MAX_NAME];

  int nlat = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);
  double *phi    = (double*) malloc(nlat*sizeof(double));
  double *dummy  = (double*) malloc(nlat*sizeof(double));
  double *cosphi = (double*) malloc(nlat*sizeof(double));
  double *plevel = (double*) malloc(nlev*sizeof(double));
  double **field1 = (double**) malloc(nlev*sizeof(double*));
  double **field2 = (double**) malloc(nlev*sizeof(double*));

  zaxisInqLevels(zaxisID, plevel);

  // gaussaw(phi, dummy, nlat);
  
  gridInqYvals(gridID, phi);
  gridInqYunits(gridID, units);

  if ( memcmp(units, "degree", 6) == 0 )
    for ( ilat = 0; ilat < nlat; ilat++ ) phi[ilat] *= DEG2RAD;

  for ( ilat = 0; ilat < nlat; ilat++ ) phi[ilat] = sin(phi[ilat]);

  for ( ilat = 0; ilat < nlat; ilat++ )
    cosphi[ilat] = sqrt(1.0 - phi[ilat]*phi[ilat]);

  for ( ilev = 0; ilev < nlev; ilev++ )
    {
      field1[ilev] = array1 + ilev*nlat;
      field2[ilev] = array2 + ilev*nlat;
    }

  for ( ilev = 0; ilev < nlev; ilev++ )
    for ( ilat = 0; ilat < nlat; ilat++ )
      field2[ilev][ilat] = 0.0;

  if ( nmiss == 0 )
    {
      for ( ilev = nlev-1; ilev >= 0; ilev-- )
	for ( n = ilev; n < nlev-1; n++ )
	  for ( ilat = 0; ilat < nlat; ilat++ )
	    {
	      field2[ilev][ilat] += fact*(field1[n][ilat]+field1[n+1][ilat])*cosphi[ilat]*(plevel[n]-plevel[n+1]);
	    }
    }
  else
    {
      for ( ilat = 0; ilat < nlat; ilat++ )
	for ( ilev = nlev-1; ilev >= 0; ilev-- )
	  for ( n = ilev; n < nlev-1; n++ )
	    {
	      if ( DBL_IS_EQUAL(field1[n][ilat], missval) )
		{
		  field2[ilev][ilat] = missval;
		  break;
		}
	      else
		field2[ilev][ilat] += fact*(field1[n][ilat]+field1[n+1][ilat])*cosphi[ilat]*(plevel[n]-plevel[n+1]);
	    }
    }

  free(field2);
  free(field1);
  free(plevel);
  free(cosphi);
  free(dummy);
  free(phi);
}


void *Mastrfu(void *argument)
{
  int nrecs;
  int recID, varID, levelID;
  int offset;
  int nmiss, nmiss1;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);

  int nvars = vlistNvars(vlistID1);
  if ( nvars != 1 ) cdoAbort("This operator works only with one variable!");

  int code = vlistInqVarCode(vlistID1, 0);
  if ( code > 0 && code != 132 ) cdoWarning("Unexpected code %d!", code);

  double missval = vlistInqVarMissval(vlistID1, 0);

  int zaxisID = vlistInqVarZaxis(vlistID1, 0);
  if ( zaxisInqType(zaxisID) != ZAXIS_PRESSURE &&
       zaxisInqType(zaxisID) != ZAXIS_GENERIC )
    {
      char longname[CDI_MAX_NAME];
      zaxisInqLongname(zaxisID, longname);
      cdoWarning("Unexpected vertical grid %s!", longname);
    }

  int gridID = vlistInqVarGrid(vlistID1, 0);
  if ( gridInqXsize(gridID) > 1 ) cdoAbort("Grid must be a zonal mean!");

  int gridsize = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);

  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  vlistDefVarCode(vlistID2, 0, 272);
  vlistDefVarName(vlistID2, 0, "mastrfu");
  vlistDefVarLongname(vlistID2, 0, "mass stream function");
  vlistDefVarUnits(vlistID2, 0, "kg/s");
  vlistDefVarDatatype(vlistID2, 0, DATATYPE_FLT32);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  double *array1 = (double*) malloc(gridsize*nlev*sizeof(double));
  double *array2 = (double*) malloc(gridsize*nlev*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      nmiss = 0;
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  offset  = gridsize*levelID;
	  streamReadRecord(streamID1, array1+offset, &nmiss1);
	  nmiss += nmiss1;
	}

      mastrfu(gridID, zaxisID, array1, array2, nmiss, missval);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  varID = 0;
	  levelID = recID;
	  streamDefRecord(streamID2, varID,  levelID);
	  offset  = gridsize*levelID;
	  streamWriteRecord(streamID2, array2+offset, nmiss);     
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  cdoFinish();

  return (0);
}
