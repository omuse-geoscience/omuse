/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Brockmann Consult
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

      Seascount   seascount         Seasonal counts
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Seascount(void *argument)
{
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs, nrecords;
  int varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  int i;
  int year, month, day, seas, seas0 = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nwpv; // number of words per value; real:1  complex:2
  int *recVarID, *recLevelID;
  int newseas, oldmon = 0, newmon;
  field_t **vars1 = NULL;
  field_t field;
  int season_start;

  cdoInitialize(argument);

  cdoOperatorAdd("seascount", 0, 0, NULL);

  season_start = get_season_start();

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;

  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  vars1 = field_malloc(vlistID1, FIELD_PTR);

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      newseas = FALSE;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
        {
          vdate = taxisInqVdate(taxisID1);
          vtime = taxisInqVtime(taxisID1);
	  cdiDecodeDate(vdate, &year, &month, &day);

	  newmon = month;
	  if ( season_start == START_DEC && newmon == 12 ) newmon = 0;

          seas = month_to_season(month);

          if ( nsets == 0 )
            {
              seas0 = seas;
              oldmon = newmon;
            }

          if ( newmon < oldmon ) newseas = TRUE;

          if ( (seas != seas0) || newseas ) break;

          oldmon = newmon;

          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(streamID1, &varID, &levelID);

              if ( tsID == 0 )
                {
                  recVarID[recID]   = varID;
                  recLevelID[recID] = levelID;
                }

	      nwpv     = vars1[varID][levelID].nwpv;
              gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

              if ( nsets == 0 )
                {
                  for ( i = 0; i < nwpv*gridsize; i++ )
                    vars1[varID][levelID].ptr[i] = vars1[varID][levelID].missval;
		  vars1[varID][levelID].nmiss = gridsize;
                }

              streamReadRecord(streamID1, field.ptr, &field.nmiss);
              field.grid    = vars1[varID][levelID].grid;
              field.missval = vars1[varID][levelID].missval;

              farcount(&vars1[varID][levelID], field);
            }

          vdate0 = vdate;
          vtime0 = vtime;
          nsets++;
          tsID++;
        }

      if ( nrecs == 0 && nsets == 0 ) break;

      taxisDefVdate(taxisID2, vdate0);
      taxisDefVtime(taxisID2, vtime0);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
        {
          varID   = recVarID[recID];
          levelID = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
        }

      if ( nrecs == 0 ) break;
      otsID++;
    }

  field_free(vars1, vlistID1);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
