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


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
const char *filetypestr(int filetype)
{
  switch ( filetype )
    {
    case FILETYPE_GRB:  return ("GRIB");            break;
    case FILETYPE_GRB2: return ("GRIB2");           break;
    case FILETYPE_NC:   return ("netCDF");          break;
    case FILETYPE_NC2:  return ("netCDF2");         break;
    case FILETYPE_NC4:  return ("netCDF4");         break;
    case FILETYPE_NC4C: return ("netCDF4 classic"); break;
    case FILETYPE_SRV:  return ("SERVICE");         break;
    case FILETYPE_EXT:  return ("EXTRA");           break;
    case FILETYPE_IEG:  return ("IEG");             break;
    default:            return ("");
    }
}

static
const char *datatypestr(int datatype)
{
  static char str[20];

  str[0] = 0;
  sprintf(str, "%d bit packed", datatype);

  if      ( datatype == DATATYPE_PACK   ) return ("P0");
  else if ( datatype > 0 && datatype <= 32 ) return (str);
  else if ( datatype == DATATYPE_CPX32  ) return ("C32");
  else if ( datatype == DATATYPE_CPX64  ) return ("C64");
  else if ( datatype == DATATYPE_FLT32  ) return ("32 bit floats");
  else if ( datatype == DATATYPE_FLT64  ) return ("64 bit floats");
  else if ( datatype == DATATYPE_INT8   ) return ("I8");
  else if ( datatype == DATATYPE_INT16  ) return ("I16");
  else if ( datatype == DATATYPE_INT32  ) return ("I32");
  else if ( datatype == DATATYPE_UINT8  ) return ("U8");
  else if ( datatype == DATATYPE_UINT16 ) return ("U16");
  else if ( datatype == DATATYPE_UINT32 ) return ("U32");
  else                                    return ("");
}

static
void print_stat(const char *sinfo, int memtype, int datatype, int filetype, off_t nvalues, double data_size, double file_size, double tw)
{
  nvalues /= 1000000;
  data_size /= 1024.*1024.*1024.;
  if ( memtype == MEMTYPE_FLOAT )
    cdoPrint("%s Read %.1f GB of 32 bit floats from %s %s, %.1f MVal/s", sinfo, data_size, datatypestr(datatype), filetypestr(filetype), nvalues/tw);
  else
    cdoPrint("%s Read %.1f GB of 64 bit floats from %s %s, %.1f MVal/s", sinfo, data_size, datatypestr(datatype), filetypestr(filetype), nvalues/tw);

  file_size /= 1024.*1024.*1024.;
  cdoPrint("%s Read %.1f GB in %.1f seconds, total %.1f MB/s", sinfo, file_size, tw, 1024*file_size/tw);
}


void *CDIread(void *argument)
{
  int memtype = MEMTYPE_DOUBLE;
  int streamID;
  int tsID, varID, levelID;
  int gridsize, i, nmiss;
  int recID, nrecs;
  int vlistID;
  int filetype = -1, datatype = -1;
  int irun, nruns = 1;
  char sinfo[64];
  off_t nvalues = 0;
  double file_size = 0, data_size = 0;
  double tw, tw0, t0, twsum = 0;
  float *farray = NULL;
  double *darray = NULL;

  sinfo[0] = 0;

  cdoInitialize(argument);

  char *envstr = getenv("MEMTYPE");
  if ( envstr )
    {
      if      ( strcmp(envstr, "float")  == 0 ) memtype = MEMTYPE_FLOAT;
      else if ( strcmp(envstr, "double") == 0 ) memtype = MEMTYPE_DOUBLE;
    }

  if ( cdoVerbose ) cdoPrint("parameter: <nruns>");

  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  if ( operatorArgc() == 1 ) nruns = parameter2int(operatorArgv()[0]);

  if ( nruns <  0 ) nruns = 0;
  if ( nruns > 99 ) nruns = 99;

  if ( cdoVerbose ) cdoPrint("nruns      : %d", nruns);

  // vlistDefNtsteps(vlistID, 1);

  for ( irun = 0; irun < nruns; ++irun )
    {
      tw0 = timer_val(timer_read);
      data_size = 0;
      nvalues = 0;

      streamID = streamOpenRead(cdoStreamName(0));

      vlistID = streamInqVlist(streamID);

      filetype = streamInqFiletype(streamID);
      datatype = vlistInqVarDatatype(vlistID, 0);
	  
      gridsize = vlistGridsizeMax(vlistID);
      
      if ( darray == NULL ) darray = (double*) malloc(gridsize*sizeof(double));
      if ( farray == NULL && memtype == MEMTYPE_FLOAT ) farray = (float*) malloc(gridsize*sizeof(float));

      t0 = timer_val(timer_read);

      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID, &varID, &levelID);

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
	      nvalues += gridsize;

	      if ( memtype == MEMTYPE_FLOAT )
		{
		  cdoAbort("streamReadRecordF not implemented!");
		  // streamReadRecordF(streamID, farray, &nmiss);
		  for ( i = 0; i < gridsize; ++i ) darray[i] = farray[i];
		  data_size += gridsize*4;
		}
	      else
		{
		  streamReadRecord(streamID, darray, &nmiss);
		  data_size += gridsize*8;
		}
	    }

	  if ( cdoVerbose )
	    {
	      tw = timer_val(timer_read) - t0;
	      t0 = timer_val(timer_read);
	      cdoPrint("Timestep %d: %.2f seconds", tsID+1, tw);
	    }

	  tsID++;
	}

      streamClose(streamID);

      tw = timer_val(timer_read) - tw0;
      twsum += tw;

      file_size = (double) fileSize(cdoStreamName(0)->args);

      if ( nruns > 1 ) sprintf(sinfo, "(run %d)", irun+1);

      print_stat(sinfo, memtype, datatype, filetype, nvalues, data_size, file_size, tw);
    }

  if ( nruns > 1 )
    print_stat("(mean)", memtype, datatype, filetype, nvalues, data_size, file_size, twsum/nruns);

  if ( darray ) free(darray);
  if ( farray ) free(farray);

  cdoFinish();

  return (0);
}
