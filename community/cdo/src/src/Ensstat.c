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

      Ensstat    ensmin          Ensemble minimum
      Ensstat    ensmax          Ensemble maximum
      Ensstat    enssum          Ensemble sum
      Ensstat    ensmean         Ensemble mean
      Ensstat    ensavg          Ensemble average
      Ensstat    ensstd          Ensemble standard deviation
      Ensstat    ensstd1         Ensemble standard deviation
      Ensstat    ensvar          Ensemble variance
      Ensstat    ensvar1         Ensemble variance
      Ensstat    enspctl         Ensemble percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Ensstat(void *argument)
{
  int i;
  int varID, recID;
  int gridID;
  int nrecs, nrecs0;
  int levelID;
  int streamID = 0;
  int nmiss;
  int fileID;
  double missval;
  typedef struct
  {
    int streamID;
    int vlistID;
    double missval;
    double *array;
  } ens_file_t;
  int pn = 0;

  cdoInitialize(argument);

  cdoOperatorAdd("ensmin",  func_min,  0, NULL);
  cdoOperatorAdd("ensmax",  func_max,  0, NULL);
  cdoOperatorAdd("enssum",  func_sum,  0, NULL);
  cdoOperatorAdd("ensmean", func_mean, 0, NULL);
  cdoOperatorAdd("ensavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ensstd",  func_std,  0, NULL);
  cdoOperatorAdd("ensstd1", func_std1, 0, NULL);
  cdoOperatorAdd("ensvar",  func_var,  0, NULL);
  cdoOperatorAdd("ensvar1", func_var1, 0, NULL);
  cdoOperatorAdd("enspctl", func_pctl, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int argc = operatorArgc();
  int nargc = argc;
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = parameter2int(operatorArgv()[0]);
      
      if ( pn < 1 || pn > 99 )
        cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
      argc--;
    }

  int count_data = FALSE;
  if ( argc == 1 )
    {
      if ( strcmp("count", operatorArgv()[nargc-1]) == 0 ) count_data = TRUE;
      else cdoAbort("Unknown parameter: >%s<", operatorArgv()[nargc-1]); 
    }
    
  int nfiles = cdoStreamCnt() - 1;

  if ( cdoVerbose ) cdoPrint("Ensemble over %d files.", nfiles);

  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoSilentMode && !cdoOverwriteMode )
    if ( fileExists(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exists!", ofilename);

  ens_file_t *ef = (ens_file_t *) malloc(nfiles*sizeof(ens_file_t));

  field_t *field = (field_t *) malloc(ompNumThreads*sizeof(field_t));
  for ( i = 0; i < ompNumThreads; i++ )
    {
      field_init(&field[i]);
      field[i].size   = nfiles;
      field[i].ptr    = (double*) malloc(nfiles*sizeof(double));
      field[i].weight = (double*) malloc(nfiles*sizeof(double));
      for ( fileID = 0; fileID < nfiles; fileID++ ) field[i].weight[fileID] = 1;
    }

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      ef[fileID].streamID = streamOpenRead(cdoStreamName(fileID));
      ef[fileID].vlistID  = streamInqVlist(ef[fileID].streamID);
    }

  /* check that the contents is always the same */
  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  int vlistID1 = ef[0].vlistID;
  int vlistID2 = vlistDuplicate(vlistID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double*) malloc(gridsize*sizeof(double));

  double *array2 = (double *) malloc(gridsize*sizeof(double));

  int nvars = vlistNvars(vlistID2);
  double *count2 = NULL;
  if ( count_data )
    {
      count2 = (double *) malloc(gridsize*sizeof(double));
      for ( varID = 0; varID < nvars; ++varID )
	{
	  char name[CDI_MAX_NAME];
	  vlistInqVarName(vlistID2, varID, name);
	  strcat(name, "_count");
	  gridID = vlistInqVarGrid(vlistID2, varID);
	  int zaxisID = vlistInqVarZaxis(vlistID2, varID);
	  int tsteptype = vlistInqVarTsteptype(vlistID2, varID);
	  int cvarID = vlistDefVar(vlistID2, gridID, zaxisID, tsteptype);
	  vlistDefVarName(vlistID2, cvarID, name);
	  vlistDefVarDatatype(vlistID2, cvarID, DATATYPE_INT16);
	  if ( cvarID != (varID+nvars) ) cdoAbort("Internal error, varIDs do not match!");
	}
    }

  int streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  streamID = ef[fileID].streamID;
	  nrecs = streamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    {
	      if ( nrecs == 0 )
		cdoAbort("Inconsistent ensemble file, too few time steps in %s!", cdoStreamName(fileID)->args);
	      else
		cdoAbort("Inconsistent ensemble file, number of records at time step %d of %s and %s differ!",
			 tsID+1, cdoStreamName(0)->args, cdoStreamName(fileID)->args);
	    }
	}

      if ( nrecs0 > 0 )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);
	  streamDefTimestep(streamID2, tsID);
	}

      for ( recID = 0; recID < nrecs0; recID++ )
	{
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(fileID, streamID, nmiss) \
                                     lastprivate(varID, levelID)
#endif
	  for ( fileID = 0; fileID < nfiles; fileID++ )
	    {
	      streamID = ef[fileID].streamID;
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, ef[fileID].array, &nmiss);
              ef[fileID].missval = vlistInqVarMissval(ef[fileID].vlistID, varID);
	    }

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  missval  = vlistInqVarMissval(vlistID1, varID);

	  nmiss = 0;
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, fileID)
#endif
	  for ( i = 0; i < gridsize; ++i )
	    {
	      int ompthID = cdo_omp_get_thread_num();

	      field[ompthID].missval = missval;
	      field[ompthID].nmiss = 0;
	      for ( fileID = 0; fileID < nfiles; fileID++ )
		{
		  field[ompthID].ptr[fileID] = ef[fileID].array[i];
		  if ( DBL_IS_EQUAL(field[ompthID].ptr[fileID], ef[fileID].missval) )
                    {
                      field[ompthID].ptr[fileID] = missval;
                      field[ompthID].nmiss++;
                    }
                }

	      if ( operfunc == func_pctl )
	        array2[i] = fldpctl(field[ompthID], pn);
	      else  
	        array2[i] = fldfun(field[ompthID], operfunc);

	      if ( DBL_IS_EQUAL(array2[i], field[ompthID].missval) )
		{
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
		  nmiss++;
		}

	      if ( count_data ) count2[i] = nfiles - field[ompthID].nmiss;
	    }

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);

	  if ( count_data )
	    {
	      streamDefRecord(streamID2, varID+nvars, levelID);
	      streamWriteRecord(streamID2, count2, 0);
	    }
	}

      tsID++;
    }
  while ( nrecs0 > 0 );

  for ( fileID = 0; fileID < nfiles; fileID++ )
    streamClose(ef[fileID].streamID);

  streamClose(streamID2);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);

  if ( ef ) free(ef);
  if ( array2 ) free(array2);
  if ( count2 ) free(count2);

  for ( i = 0; i < ompNumThreads; i++ )
    {
      if ( field[i].ptr    ) free(field[i].ptr);
      if ( field[i].weight ) free(field[i].weight);
    }

  if ( field ) free(field);

  cdoFinish();

  return (0);
}
