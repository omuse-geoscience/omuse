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

      Sort sortcode  Sort by code number
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


typedef struct
{
  int        nmiss;
  int        levelID;
  double     level;
}
levinfo_t;

typedef struct
{
  int        varID;
  int        nlevs;
  int        code;
  char       name[CDI_MAX_NAME];
  levinfo_t *levInfo;
}
varinfo_t;

static
int cmpvarcode(const void *s1, const void *s2)
{
  int cmp = 0;
  const varinfo_t *x = (const varinfo_t *) s1;
  const varinfo_t *y = (const varinfo_t *) s2;
  /*
  printf("%d %d  %d %d\n", x->code, y->code, x, y);
  */
  if      ( x->code < y->code ) cmp = -1;
  else if ( x->code > y->code ) cmp =  1;

  return (cmp);
}

static
int cmpvarname(const void *s1, const void *s2)
{
  const varinfo_t *x = (const varinfo_t *) s1;
  const varinfo_t *y = (const varinfo_t *) s2;

  return (strcmp(x->name, y->name));
}

static
int cmpvarlevel(const void *s1, const void *s2)
{
  int cmp = 0;
  const levinfo_t *x = (const levinfo_t *) s1;
  const levinfo_t *y = (const levinfo_t *) s2;

  if      ( x->level < y->level ) cmp = -1;
  else if ( x->level > y->level ) cmp =  1;

  return (cmp);
}

static
int cmpvarlevelrev(const void *s1, const void *s2)
{
  int cmp = 0;
  const levinfo_t *x = (const levinfo_t *) s1;
  const levinfo_t *y = (const levinfo_t *) s2;

  if      ( x->level > y->level ) cmp = -1;
  else if ( x->level < y->level ) cmp =  1;

  return (cmp);
}

static
void setNmiss(int varID, int levelID, int nvars, varinfo_t *varInfo, int nmiss)
{
  int vindex, lindex;
  int nlevs;

  for ( vindex = 0; vindex < nvars; vindex++ )
    if ( varInfo[vindex].varID == varID ) break;

  if ( vindex == nvars ) cdoAbort("Internal problem; varID not found!");

  nlevs = varInfo[vindex].nlevs; 
  for ( lindex = 0; lindex < nlevs; lindex++ )
    if ( varInfo[vindex].levInfo[lindex].levelID == levelID ) break;

  if ( lindex == nlevs ) cdoAbort("Internal problem; levelID not found!");

  varInfo[vindex].levInfo[lindex].nmiss = nmiss;
}


void *Sort(void *argument)
{
  int SORTCODE, SORTNAME, SORTLEVEL;
  int operatorID;
  int streamID1, streamID2;
  int nrecs;
  int tsID, recID, varID, levelID, zaxisID;
  int vindex, lindex;
  int nvars, nlevs, offset;
  int vlistID1, vlistID2;
  int gridsize;
  int nmiss;
  double *single;
  double **vardata = NULL;
  varinfo_t *varInfo;
  int taxisID1, taxisID2;
  int (*cmpvarlev)(const void *, const void *) = cmpvarlevel;

  cdoInitialize(argument);

  SORTCODE  = cdoOperatorAdd("sortcode",  0, 0, NULL);
  SORTNAME  = cdoOperatorAdd("sortname",  0, 0, NULL);
  SORTLEVEL = cdoOperatorAdd("sortlevel", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  if ( operatorID == SORTLEVEL && operatorArgc() == 1 )
    {
      int iarg = parameter2int(operatorArgv()[0]);
      if ( iarg < 0 ) cmpvarlev = cmpvarlevelrev;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  /*
  if ( operatorID == SORTCODE )
      vlistSortCode(vlistID2);
   else if ( operatorID == SORTNAME )
      ;
   else if ( operatorID == SORTLEVEL )
      ;
  */

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars   = vlistNvars(vlistID1);

  varInfo = (varinfo_t*) malloc(nvars*sizeof(varinfo_t));
  for ( varID = 0; varID < nvars; ++varID )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      varInfo[varID].nlevs = nlevs;
      varInfo[varID].levInfo = (levinfo_t*) malloc(nlevs*sizeof(levinfo_t));
    }

  vardata = (double**) malloc(nvars*sizeof(double*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata[varID] = (double*) malloc(gridsize*nlevs*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      varInfo[varID].varID = varID;
	      varInfo[varID].code  = vlistInqVarCode(vlistID1, varID);
	      vlistInqVarName(vlistID1, varID, varInfo[varID].name);
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      varInfo[varID].levInfo[levelID].levelID = levelID;
	      varInfo[varID].levInfo[levelID].level   = zaxisInqLevel(zaxisID, levelID);
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata[varID] + offset;

	  streamReadRecord(streamID1, single, &nmiss);

	  setNmiss(varID, levelID, nvars, varInfo, nmiss);
	  // varInfo[varID].levInfo[levelID].nmiss = nmiss;
	}

      if ( tsID == 0 )
	{
	  if ( cdoVerbose )
	    for ( vindex = 0; vindex < nvars; vindex++ )
	      {
		nlevs = varInfo[vindex].nlevs;
		for ( lindex = 0; lindex < nlevs; ++lindex )
		  printf("sort in: %d %s %d %d %g\n",
			 vindex, varInfo[vindex].name, varInfo[vindex].code, varInfo[vindex].nlevs, varInfo[vindex].levInfo[lindex].level);
	      }

	  if      ( operatorID == SORTCODE )
	    qsort(varInfo, nvars, sizeof(varinfo_t), cmpvarcode);
	  else if ( operatorID == SORTNAME )
	    qsort(varInfo, nvars, sizeof(varinfo_t), cmpvarname);
	  else if ( operatorID == SORTLEVEL )
	    {
	      for ( vindex = 0; vindex < nvars; vindex++ )
		{
		  nlevs = varInfo[vindex].nlevs;
		  qsort(varInfo[vindex].levInfo, nlevs, sizeof(levinfo_t), cmpvarlev);
		}
	    }

	  if ( cdoVerbose )
	    for ( vindex = 0; vindex < nvars; vindex++ )
	      {
		nlevs = varInfo[vindex].nlevs;
		for ( lindex = 0; lindex < nlevs; ++lindex )
		  printf("sort out: %d %s %d %d %g\n",
			 vindex, varInfo[vindex].name, varInfo[vindex].code, varInfo[vindex].nlevs, varInfo[vindex].levInfo[lindex].level);
	      }
	}

      for ( vindex = 0; vindex < nvars; vindex++ )
	{
	  varID = varInfo[vindex].varID;
	  nlevs = varInfo[vindex].nlevs;
	  for ( lindex = 0; lindex < nlevs; ++lindex )
	    {
	      levelID = varInfo[vindex].levInfo[lindex].levelID;
	      nmiss   = varInfo[vindex].levInfo[lindex].nmiss;

	      if ( tsID == 0 || vlistInqVarTsteptype(vlistID1, varID) != TSTEP_CONSTANT )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
		  offset   = gridsize*levelID;
		  single   = vardata[varID] + offset;

		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, single, nmiss);
		}
	    }
	}

      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  for ( varID = 0; varID < nvars; varID++ ) free(vardata[varID]);
  free(vardata);

  for ( vindex = 0; vindex < nvars; vindex++ ) free(varInfo[vindex].levInfo);
  free(varInfo);

  cdoFinish();

  return (0);
}
