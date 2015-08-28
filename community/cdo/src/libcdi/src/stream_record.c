#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdf_int.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"


void recordInitEntry(record_t *record)
{
  record->position = CDI_UNDEFID;
  record->size     = 0;
  record->param    = 0;
  record->ilevel   = CDI_UNDEFID;
  record->used     = FALSE;
  record->varID    = CDI_UNDEFID;
  record->levelID  = CDI_UNDEFID;
  memset(record->varname, 0, sizeof(record->varname));
  memset(&record->tiles, 0, sizeof(record->tiles));
}


int recordNewEntry(stream_t *streamptr, int tsID)
{
  int recordID = 0;
  int recordSize = streamptr->tsteps[tsID].recordSize;
  record_t *records = streamptr->tsteps[tsID].records;
  /*
    Look for a free slot in record.
    (Create the table the first time through).
  */
  if ( ! recordSize )
    {
      int i;
      recordSize = 1;   /*  <<<<----  */
      records = (record_t *)xmalloc((size_t)recordSize * sizeof (record_t));
      if ( records == NULL )
	{
          Message("recordSize = %d", recordSize);
	  SysError("Allocation of records failed");
	}

      for ( i = 0; i < recordSize; i++ )
	records[i].used = CDI_UNDEFID;
    }
  else
    {
      while ( recordID < recordSize )
	{
	  if ( records[recordID].used == CDI_UNDEFID ) break;
	  recordID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( recordID == recordSize )
    {
      int i;

      recordSize = 2*recordSize;
      records    = (record_t *)xrealloc(records,
                                        (size_t)recordSize * sizeof (record_t));
      if ( records == NULL )
	{
          Message("recordSize = %d", recordSize);
	  SysError("Reallocation of records failed");
	}
      recordID = recordSize/2;

      for ( i = recordID; i < recordSize; i++ )
	records[i].used = CDI_UNDEFID;
    }

  recordInitEntry(&records[recordID]);

  records[recordID].used = 1;

  streamptr->tsteps[tsID].recordSize = recordSize;
  streamptr->tsteps[tsID].records    = records;

  return (recordID);
}

static
void cdiInitRecord(stream_t *streamptr)
{
  streamptr->record = (Record *) malloc(sizeof(Record));

  streamptr->record->param      = 0;
  streamptr->record->level      = 0;
  streamptr->record->date       = 0;
  streamptr->record->time       = 0;
  streamptr->record->gridID     = 0;
  streamptr->record->buffer     = NULL;
  streamptr->record->buffersize = 0;
  streamptr->record->position   = 0;
  streamptr->record->varID      = 0;
  streamptr->record->levelID    = CDI_UNDEFID;
}


void streamInqRecord(int streamID, int *varID, int *levelID)
{
  check_parg(varID);
  check_parg(levelID);

  stream_t *streamptr = stream_to_pointer(streamID);
  stream_check_ptr(__func__, streamptr);

  cdiDefAccesstype(streamID, TYPE_REC);

  if ( ! streamptr->record ) cdiInitRecord(streamptr);

  int tsID   = streamptr->curTsID;
  int rindex = streamptr->tsteps[tsID].curRecID + 1;

  if ( rindex >= streamptr->tsteps[tsID].nrecs )
    Error("record %d not available at timestep %d", rindex+1, tsID+1);

  int recID  = streamptr->tsteps[tsID].recIDs[rindex];

  if ( recID == -1 || recID >= streamptr->tsteps[tsID].nallrecs )
    Error("Internal problem! tsID = %d recID = %d", tsID, recID);

  *varID   = streamptr->tsteps[tsID].records[recID].varID;
  int lindex = streamptr->tsteps[tsID].records[recID].levelID;

  int isub = subtypeInqActiveIndex(streamptr->vars[*varID].subtypeID);
  *levelID = streamptr->vars[*varID].recordTable[isub].lindex[lindex];

  if ( CDI_Debug )
    Message("tsID = %d, recID = %d, varID = %d, levelID = %d\n", tsID, recID, *varID, *levelID);

  streamptr->curTsID = tsID;
  streamptr->tsteps[tsID].curRecID = rindex;
}

/*
@Function  streamDefRecord
@Title     Define the next record

@Prototype void streamDefRecord(int streamID, int varID, int levelID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.

@Description
The function streamDefRecord defines the meta-data of the next record.
@EndFunction
*/
void streamDefRecord(int streamID, int varID, int levelID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  stream_check_ptr(__func__, streamptr);

  int tsID = streamptr->curTsID;

  if ( tsID == CDI_UNDEFID )
    {
      tsID++;
      streamDefTimestep(streamID, tsID);
    }

  if ( ! streamptr->record ) cdiInitRecord(streamptr);

  int vlistID = streamptr->vlistID;
  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int param   = vlistInqVarParam(vlistID, varID);
  int level   = (int) zaxisInqLevel(zaxisID, levelID);

  streamptr->record->varID    = varID;
  streamptr->record->levelID  = levelID;
  streamptr->record->param    = param;
  streamptr->record->level    = level;
  streamptr->record->date     = streamptr->tsteps[tsID].taxis.vdate;
  streamptr->record->time     = streamptr->tsteps[tsID].taxis.vtime;
  streamptr->record->gridID   = gridID;
  streamptr->record->prec     = vlistInqVarDatatype(vlistID, varID);

  switch (streamptr->filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      grbDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      srvDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      extDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      iegDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);
      cdfDefRecord(streamptr);
      break;
#endif
    default:
      Error("%s support not compiled in!", strfiletype(streamptr->filetype));
      break;
    }
}


void streamReadRecord(int streamID, double *data, int *nmiss)
{
  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);
  stream_check_ptr(__func__, streamptr);

  *nmiss = 0;

  switch (streamptr->filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      grbReadRecord(streamptr, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      srvReadRecord(streamptr, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      extReadRecord(streamptr, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      iegReadRecord(streamptr, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      cdfReadRecord(streamptr, data, nmiss);
      break;
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(streamptr->filetype));
	break;
      }
    }
}

static void
stream_write_record(int streamID, int memtype, const void *data, int nmiss)
{
  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  stream_check_ptr(__func__, streamptr);

  switch (streamptr->filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      grb_write_record(streamptr, memtype, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      if ( memtype == MEMTYPE_FLOAT ) Error("srvWriteRecord not implemented for memtype float!");
      srvWriteRecord(streamptr, (const double *)data);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      if ( memtype == MEMTYPE_FLOAT ) Error("extWriteRecord not implemented for memtype float!");
      extWriteRecord(streamptr, (const double *)data);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      if ( memtype == MEMTYPE_FLOAT ) Error("iegWriteRecord not implemented for memtype float!");
      iegWriteRecord(streamptr, (const double *)data);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	cdf_write_record(streamptr, memtype, data, nmiss);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(streamptr->filetype));
	break;
      }
    }
}

/*
@Function  streamWriteRecord
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteRecord(int streamID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteRecord writes the values of a horizontal slice (record) of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteRecord(int streamID, const double *data, int nmiss)
{
  stream_write_record(streamID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

void streamWriteRecordF(int streamID, const float *data, int nmiss)
{
  stream_write_record(streamID, MEMTYPE_FLOAT, (const void *) data, nmiss);
}


void streamCopyRecord(int streamID2, int streamID1)
{
  stream_t *streamptr1 = stream_to_pointer(streamID1);
  stream_t *streamptr2 = stream_to_pointer(streamID2);

  stream_check_ptr(__func__, streamptr1);
  stream_check_ptr(__func__, streamptr2);

  int filetype1 = streamptr1->filetype;
  int filetype2 = streamptr2->filetype;
  int filetype  = FILETYPE_UNDEF;

  if ( filetype1 == filetype2 ) filetype = filetype2;
  else
    {
      switch (filetype1)
        {
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
          switch (filetype2)
            {
            case FILETYPE_NC:
            case FILETYPE_NC2:
            case FILETYPE_NC4:
            case FILETYPE_NC4C:
              Warning("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));
              filetype = filetype2;
              break;
            }
          break;
        }
    }

  if ( filetype == FILETYPE_UNDEF )
    Error("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      grbCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      srvCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      extCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      iegCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      cdfCopyRecord(streamptr2, streamptr1);
      break;
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}


void cdi_create_records(stream_t *streamptr, int tsID)
{
  unsigned nrecords, maxrecords;
  record_t *records;

  tsteps_t* sourceTstep = streamptr->tsteps;
  tsteps_t* destTstep = sourceTstep + tsID;

  if ( destTstep->records ) return;

  int vlistID = streamptr->vlistID;

  if ( tsID == 0 )
    {
      maxrecords = 0;
      int nvars = streamptr->nvars;
      for ( int varID = 0; varID < nvars; varID++)
        for (int isub=0; isub<streamptr->vars[varID].subtypeSize; isub++)
          maxrecords += (unsigned)streamptr->vars[varID].recordTable[isub].nlevs;
    }
  else
    {
      maxrecords = (unsigned)sourceTstep->recordSize;
    }

  if ( tsID == 0 )
    {
      nrecords = maxrecords;
    }
  else if ( tsID == 1 )
    {
      nrecords = 0;
      maxrecords = (unsigned)sourceTstep->recordSize;
      for ( unsigned recID = 0; recID < maxrecords; recID++ )
	{
	  int varID = sourceTstep->records[recID].varID;
	  nrecords += (varID == CDI_UNDEFID /* varID = CDI_UNDEFID for write mode !!! */
                       || vlistInqVarTsteptype(vlistID, varID) != TSTEP_CONSTANT);
          //    printf("varID nrecords %d %d %d \n", varID, nrecords, vlistInqVarTsteptype(vlistID, varID));
	}
    }
  else
    {
      nrecords = (unsigned)streamptr->tsteps[1].nallrecs;
    }
  //  printf("tsID, nrecords %d %d\n", tsID, nrecords);

  if ( maxrecords > 0 )
    records = (record_t *) malloc(maxrecords * sizeof(record_t));
  else
    records = NULL;

  destTstep->records    = records;
  destTstep->recordSize = (int)maxrecords;
  destTstep->nallrecs   = (int)nrecords;

  if ( tsID == 0 )
    {
      for ( unsigned recID = 0; recID < maxrecords; recID++ )
        recordInitEntry(&destTstep->records[recID]);
    }
  else
    {
      memcpy(destTstep->records, sourceTstep->records, (size_t)maxrecords*sizeof(record_t));

      for ( unsigned recID = 0; recID < maxrecords; recID++ )
	{
          record_t* curRecord = &sourceTstep->records[recID];
          destTstep->records[recID].used = curRecord->used;
          if ( curRecord->used != CDI_UNDEFID && curRecord->varID != -1 ) /* curRecord->varID = -1 for write mode !!! */
            {
              if ( vlistInqVarTsteptype(vlistID, curRecord->varID) != TSTEP_CONSTANT )
                {
                  destTstep->records[recID].position = CDI_UNDEFID;
                  destTstep->records[recID].size     = 0;
                  destTstep->records[recID].used     = FALSE;
                }
            }
	}
    }
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
