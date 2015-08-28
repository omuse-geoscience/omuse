#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_ext.h"
#include "stream_fcommon.h"
#include "varscan.h"
#include "datetime.h"
#include "extra.h"
#include "vlist.h"


#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID

#define SINGLE_PRECISION  4
#define DOUBLE_PRECISION  8

#if defined (HAVE_LIBEXTRA)


typedef struct {
  int param;
  int level;
} extcompvar_t;

static
int extInqDatatype(int prec, int number)
{
  int datatype;

  if ( number == 2 )
    {
      if ( prec == DOUBLE_PRECISION ) datatype = DATATYPE_CPX64;
      else                            datatype = DATATYPE_CPX32;
    }
  else
    {
      if ( prec == DOUBLE_PRECISION ) datatype = DATATYPE_FLT64;
      else                            datatype = DATATYPE_FLT32;
    }

  return (datatype);
}

static
void extDefDatatype(int datatype, int *prec, int *number)
{

  if ( datatype != DATATYPE_FLT32 && datatype != DATATYPE_FLT64 &&
       datatype != DATATYPE_CPX32 && datatype != DATATYPE_CPX64 )
    datatype = DATATYPE_FLT32;

  if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    *number = 2;
  else
    *number = 1;

  if ( datatype == DATATYPE_FLT64 || datatype == DATATYPE_CPX64 )
    *prec = DOUBLE_PRECISION;
  else 
    *prec = SINGLE_PRECISION;
}

/* not used
int extInqRecord(stream_t *streamptr, int *varID, int *levelID)
{
  int status;
  int fileID;
  int icode, ilevel;
  int zaxisID = -1;
  int header[4];
  int vlistID;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  *varID   = -1;
  *levelID = -1;

  status = extRead(fileID, extp);
  if ( status != 0 ) return (0);

  extInqHeader(extp, header);

  icode  = header[1];
  ilevel = header[2];

  *varID = vlistInqVarID(vlistID, icode);

  if ( *varID == UNDEFID ) Error("Code %d undefined", icode);

  zaxisID = vlistInqVarZaxis(vlistID, *varID);

  *levelID = zaxisInqLevelID(zaxisID, (double) ilevel);

  return (1);
}
*/

void extReadRecord(stream_t *streamptr, double *data, int *nmiss)
{
  int vlistID, fileID;
  int status;
  int recID, vrecID, tsID;
  off_t recpos;
  int header[4];
  int varID, gridID;
  int i, size;
  double missval;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;
  tsID    = streamptr->curTsID;
  vrecID  = streamptr->tsteps[tsID].curRecID;
  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  recpos  = streamptr->tsteps[tsID].records[recID].position;
  varID   = streamptr->tsteps[tsID].records[recID].varID;

  fileSetPos(fileID, recpos, SEEK_SET);

  status = extRead(fileID, extp);
  if ( status != 0 )
    Error("Failed to read EXTRA record");

  extInqHeader(extp, header);
  extInqDataDP(extp, data);

  missval = vlistInqVarMissval(vlistID, varID);
  gridID  = vlistInqVarGrid(vlistID, varID);
  size    = gridInqSize(gridID);

  streamptr->numvals += size;

  *nmiss = 0;
  if ( vlistInqVarNumber(vlistID, varID) == CDI_REAL )
    {
      for ( i = 0; i < size; i++ )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
  else
    {
      for ( i = 0; i < 2*size; i+=2 )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
}


void extCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  streamFCopyRecord(streamptr2, streamptr1, "EXTRA");
}


void extDefRecord(stream_t *streamptr)
{
  int gridID;
  int header[4];
  int pdis, pcat, pnum;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  gridID   = streamptr->record->gridID;

  cdiDecodeParam(streamptr->record->param, &pnum, &pcat, &pdis);
  header[0] = streamptr->record->date;
  header[1] = pnum;
  header[2] = streamptr->record->level;
  header[3] = gridInqSize(gridID);

  extDefDatatype(streamptr->record->prec, &extp->prec, &extp->number);

  extDefHeader(extp, header);
}


void extWriteRecord(stream_t *streamptr, const double *data)
{
  int fileID = streamptr->fileID;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  extDefDataDP(extp, data);
  extWrite(fileID, extp);
}

static
void extAddRecord(stream_t *streamptr, int param, int level, int xysize,
		  size_t recsize, off_t position, int prec, int number)
{
  int leveltype;
  int gridID = CDI_UNDEFID;
  int levelID = 0;
  int tsID, recID, varID;
  record_t *record;
  grid_t grid;
  int vlistID;

  vlistID = streamptr->vlistID;
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamptr, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  (*record).size     = recsize;
  (*record).position = position;
  (*record).param     = param;
  (*record).ilevel   = level;

  memset(&grid, 0, sizeof(grid_t));
  grid.type  = GRID_GENERIC;
  grid.size  = xysize;
  grid.xsize = xysize;
  grid.ysize = 0;
  grid.xvals = NULL;
  grid.yvals = NULL;
  gridID = varDefGrid(vlistID, &grid, 0);
  /*
  if ( level == 0 ) leveltype = ZAXIS_SURFACE;
  else              leveltype = ZAXIS_GENERIC;
  */
  leveltype = ZAXIS_GENERIC;

  varAddRecord(recID, param, gridID, leveltype, 0, level, 0, 0, 0,
	       extInqDatatype(prec, number), &varID, &levelID, TSTEP_INSTANT, 0, 0, -1,
               NULL, NULL, NULL, NULL, NULL, NULL);

  (*record).varID   = (short)varID;
  (*record).levelID = (short)levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d gridID = %d levelID = %d",
	    varID, gridID, levelID);
}

static
void extScanTimestep1(stream_t *streamptr)
{
  int header[4];
  int status;
  int fileID;
  int rxysize = 0;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID;
  long recsize;
  off_t recpos;
  int nrecords, nrecs, recID;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  extcompvar_t compVar, compVar0;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  streamptr->curTsID = 0;

  tsID  = tstepsNewEntry(streamptr);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamptr->fileID;

  nrecs = 0;
  while ( TRUE )
    {
      recpos = fileGetPos(fileID);
      status = extRead(fileID, extp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 1;
	  break;
	}
      recsize = fileGetPos(fileID) - recpos;

      extInqHeader(extp, header);

      vdate   = header[0];
      vtime   = 0;
      rcode   = header[1];
      rlevel  = header[2];
      rxysize = header[3];

      param = cdiEncodeParam(rcode, 255, 255);

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}
      else
	{
	  datetime.date = vdate;
	  datetime.time = vtime;
	  compVar.param = param;
          compVar.level = rlevel;
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      compVar0.param  = streamptr->tsteps[0].records[recID].param;
	      compVar0.level = streamptr->tsteps[0].records[recID].ilevel;

	      if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) == 0 ) break;
	    }
	  if ( recID < nrecs ) break;
	  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) )
	    Warning("Inconsistent verification time for code %d level %d", rcode, rlevel);
	}

      nrecs++;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", nrecs, (int)recpos, rcode, rlevel, vdate, vtime);

      extAddRecord(streamptr, param, rlevel, rxysize, (size_t)recsize, recpos, extp->prec, extp->number);
    }

  streamptr->rtsteps = 1;

  cdi_generate_vars(streamptr);

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  taxis->type  = TAXIS_ABSOLUTE;
  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  vlist_check_contents(vlistID);

  nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
        (record_t *)xrealloc(streamptr->tsteps[0].records, (size_t)nrecords * sizeof (record_t));
    }

  streamptr->tsteps[0].recIDs = (int *)xmalloc((size_t)nrecords * sizeof (int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for ( varID = 0; varID < streamptr->nvars; varID++ )
	    {
	      vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	    }
	}
    }
}

static
int extScanTimestep2(stream_t *streamptr)
{
  int header[4];
  int status;
  int fileID;
  // int rxysize = 0;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  int tsID;
  int varID;
  off_t recpos = 0;
  int nrecords, nrecs, recID, rindex;
  int nextstep;
  taxis_t *taxis;
  int vlistID;
  extcompvar_t compVar, compVar0;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  streamptr->curTsID = 1;

  fileID  = streamptr->fileID;
  vlistID = streamptr->vlistID;

  tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpected timestep %d", tsID+1);

  taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  nrecords = streamptr->tsteps[0].nallrecs;
  streamptr->tsteps[1].recIDs = xmalloc((size_t)nrecords * sizeof (int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position =
	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     =
	streamptr->tsteps[0].records[recID].size;
    }

  for ( rindex = 0; rindex <= nrecords; rindex++ )
    {
      recpos = fileGetPos(fileID);
      status = extRead(fileID, extp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

      extInqHeader(extp, header);

      vdate  = header[0];
      vtime  = 0;
      rcode  = header[1];
      rlevel = header[2];
      // rxysize = header[3];

      param = cdiEncodeParam(rcode, 255, 255);

      if ( rindex == 0 )
	{
	  taxis->type  = TAXIS_ABSOLUTE;
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;
	}

      compVar.param = param;
      compVar.level = rlevel;
      nextstep = FALSE;
      for ( recID = 0; recID < nrecords; recID++ )
	{
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) == 0 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  nextstep = TRUE;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = TRUE;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}
	      break;
	    }
	}
      if ( recID == nrecords )
	{
	  Warning("Code %d level %d not found at timestep %d", rcode, rlevel, tsID+1);
	  return (CDI_EUFSTRUCT);
	}

      if ( nextstep ) break;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", rindex+1, (int)recpos, rcode, rlevel, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
      compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

      if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
    }

  nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  varID = streamptr->tsteps[tsID].records[recID].varID;
          vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  return (0);
}


int extInqContents(stream_t *streamptr)
{
  int fileID;
  int status = 0;

  fileID = streamptr->fileID;

  streamptr->curTsID = 0;

  extScanTimestep1(streamptr);

  if ( streamptr->ntsteps == -1 ) status = extScanTimestep2(streamptr);

  fileSetPos(fileID, 0, SEEK_SET);

  return (status);
}

static
long extScanTimestep(stream_t *streamptr)
{
  int header[4];
  int status;
  int fileID;
  // int rxysize = 0;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  off_t recpos = 0;
  int recID;
  int rindex, nrecs = 0;
  extcompvar_t compVar, compVar0;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;
  /*
  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }
  */

  int tsID  = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *)xmalloc((size_t)nrecs * sizeof (int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      for ( rindex = 0; rindex <= nrecs; rindex++ )
	{
	  recpos = fileGetPos(fileID);
	  status = extRead(fileID, extp);
	  if ( status != 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

	  extInqHeader(extp, header);

	  vdate  = header[0];
	  vtime  = 0;
	  rcode  = header[1];
	  rlevel = header[2];
	  // rxysize = header[3];

	  param = cdiEncodeParam(rcode, 255, 255);

	  // if ( rindex == nrecs ) break; gcc-4.5 internal compiler error
	  if ( rindex == nrecs ) continue;
	  recID = streamptr->tsteps[tsID].recIDs[rindex];

	  if ( rindex == 0 )
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;
	    }

	  compVar.param  = param;
          compVar.level  = rlevel;
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	      Error("Invalid, unsupported or inconsistent record structure!");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d%8d%4d%8d%8d%6d", rindex, (int)recpos, rcode, rlevel, vdate, vtime);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamptr);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = 1;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return (streamptr->ntsteps);
}


int extInqTimestep(stream_t *streamptr, int tsID)
{
  int nrecs;
  long ntsteps;

  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);

  ntsteps = UNDEFID;
  while ( ( tsID + 1 ) > streamptr->rtsteps && ntsteps == UNDEFID )
    ntsteps = extScanTimestep(streamptr);

  if ( tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID )
    {
      nrecs = 0;
    }
  else
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return (nrecs);
}


void extReadVarDP(stream_t *streamptr, int varID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int levID, nlevs, gridID, gridsize;
  off_t recpos, currentfilepos;
  int header[4];
  int tsid;
  int recID;
  int i;
  double missval;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  /* NOTE: tiles are not supported here! */
  nlevs    = streamptr->vars[varID].recordTable[0].nlevs;
  missval  = vlistInqVarMissval(vlistID, varID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsid     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  for (levID = 0; levID < nlevs; levID++)
    {
      /* NOTE: tiles are not supported here! */
      recID = streamptr->vars[varID].recordTable[0].recordID[levID];
      recpos = streamptr->tsteps[tsid].records[recID].position;
      fileSetPos(fileID, recpos, SEEK_SET);
      extRead(fileID, extp);
      extInqHeader(extp, header);
      extInqDataDP(extp, &data[levID*gridsize]);
    }
  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  if ( vlistInqVarNumber(vlistID, varID) == CDI_REAL )
    {
      for ( i = 0; i < nlevs*gridsize; i++ )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
  else
    {
      for ( i = 0; i < 2*nlevs*gridsize; i+=2 )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
}


void extReadVarSliceDP(stream_t *streamptr, int varID, int levID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int nlevs, gridID, gridsize;
  off_t recpos, currentfilepos;
  int header[4];
  int tsid;
  int recID;
  int i;
  double missval;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  /* NOTE: tiles are not supported here! */
  nlevs    = streamptr->vars[varID].recordTable[0].nlevs;
  missval  = vlistInqVarMissval(vlistID, varID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsid     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d",
	     nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  /* NOTE: tiles are not supported here! */
  recID = streamptr->vars[varID].recordTable[0].recordID[levID];
  recpos = streamptr->tsteps[tsid].records[recID].position;
  fileSetPos(fileID, recpos, SEEK_SET);
  extRead(fileID, extp);
  extInqHeader(extp, header);
  extInqDataDP(extp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  if ( vlistInqVarNumber(vlistID, varID) == CDI_REAL )
    {
      for ( i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
  else
    {
      for ( i = 0; i < 2*gridsize; i+=2 )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
}


void extWriteVarDP(stream_t *streamptr, int varID, const double *data)
{
  int fileID;
  int levID, nlevs, gridID, gridsize;
  int zaxisID;
  double level;
  int header[4];
  int tsID;
  int vlistID;
  int pdis, pcat, pnum;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  tsID     = streamptr->curTsID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  nlevs    = zaxisInqSize(zaxisID);

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  header[0] = streamptr->tsteps[tsID].taxis.vdate;
  header[1] = pnum;
  header[3] = gridInqSize(gridID);

  extDefDatatype(vlistInqVarDatatype(vlistID, varID), &extp->prec, &extp->number);

  for ( levID = 0;  levID < nlevs; levID++ )
    {
      level = zaxisInqLevel(zaxisID, levID);

      header[2] = (int) level;
      extDefHeader(extp, header);
      extDefDataDP(extp, &data[levID*gridsize]);
      extWrite(fileID, extp);
    }
}


void extWriteVarSliceDP(stream_t *streamptr, int varID, int levID, const double *data)
{
  int fileID;
  int gridID;
  int zaxisID;
  double level;
  int header[4];
  int tsID;
  int vlistID;
  int pdis, pcat, pnum;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  tsID     = streamptr->curTsID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  level    = zaxisInqLevel(zaxisID, levID);

  if ( CDI_Debug )
    Message("gridID = %d zaxisID = %d", gridID, zaxisID);

  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  header[0] = streamptr->tsteps[tsID].taxis.vdate;
  header[1] = pnum;
  header[2] = (int) level;
  header[3] = gridInqSize(gridID);

  extDefDatatype(vlistInqVarDatatype(vlistID, varID), &extp->prec, &extp->number);

  extDefHeader(extp, header);
  extDefDataDP(extp, data);
  extWrite(fileID, extp);
}

#endif /* HAVE_LIBEXTRA */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
