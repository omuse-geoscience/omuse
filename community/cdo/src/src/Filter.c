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
 
      Filter    highpass
      Filter    lowpass
      Filter    bandpass
*/

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "statistic.h"
#include "pstream.h"

#if defined(HAVE_LIBFFTW3) 
#include <fftw3.h>
#endif


#define  NALLOC_INC  1024


/* include from Tinfo.c */
void getTimeInc(double jdelta, int vdate0, int vdate1, int *incperiod, int *incunit);

static
void create_fmasc(int nts, double fdata, double fmin, double fmax, int *fmasc)
{  
  double dimin = nts*fmin / fdata;
  double dimax = nts*fmax / fdata;

  int imin = dimin<0 ? 0 : (int)floor(dimin);  
  int imax = ceil(dimax)>nts/2 ? nts/2 : (int) ceil(dimax);

  if ( imin < 0 || imin >= nts ) cdoAbort("Parameter fmin=%g out of bounds (1-%d)!", nts);
  if ( imax < 0 || imax >= nts ) cdoAbort("Parameter fmax=%g out of bounds (1-%d)!", nts);

  fmasc[imin] = 1;
  for ( int i = imin+1; i <= imax; i++ )  
    fmasc[i] = fmasc[nts-i] = 1; 
}

#if defined(HAVE_LIBFFTW3) 
static
void filter_fftw(int nts, const int *fmasc, fftw_complex *fft_out, fftw_plan *p_T2S, fftw_plan *p_S2T)
{  
  int i;

  fftw_execute(*p_T2S);

  for ( i = 0; i < nts; i++ )
    if ( ! fmasc[i] )
      {
        fft_out[i][0] = 0;
        fft_out[i][1] = 0;
      }
  
  fftw_execute(*p_S2T);
  
  return;
}
#endif

static
void filter_intrinsic(int nts, const int *fmasc, double *array1, double *array2)
{  
  int i;
  int lpower2 = FALSE;
  double *work_r = NULL;
  double *work_i = NULL;

  if ( (nts&(nts-1)) == 0 ) lpower2 = TRUE;

  if ( !lpower2 )
    {
      work_r = (double*) malloc(nts*sizeof(double));
      work_i = (double*) malloc(nts*sizeof(double));
    }

  if ( lpower2 )
    fft(array1, array2, nts, 1);
  else
    ft_r(array1, array2, nts, 1, work_r, work_i);

  for ( i = 0; i < nts; i++ )
    if ( ! fmasc[i] )
      array1[i] = array2[i] = 0;

  if ( lpower2 )
    fft(array1, array2, nts, -1);
  else
    ft_r(array1, array2, nts, -1, work_r, work_i);

  if ( work_r ) free(work_r);
  if ( work_i ) free(work_i);
  
  return;
}


void *Filter(void *argument)
{
  enum {BANDPASS, HIGHPASS, LOWPASS};
  char *tunits[] = {"second", "minute", "hour", "day", "month", "year"};
  int iunits[] = {31536000, 525600, 8760, 365, 12, 1};
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int i;
  int nts;
  int nalloc = 0;
  int nmiss;
  int nvars, nlevel;
  int incperiod0, incunit0, incunit, calendar;
  int year0, month0, day0;
  double fdata = 0;
  field_t ***vars = NULL;
  double fmin = 0, fmax = 0;
  int use_fftw = FALSE;
  dtlist_type *dtlist = dtlist_new();
  typedef struct
  {
    double *array1;
    double *array2;
#if defined(HAVE_LIBFFTW3) 
    fftw_complex *in_fft;
    fftw_complex *out_fft;
    fftw_plan p_T2S;
    fftw_plan p_S2T;
#endif
  } memory_t;
  memory_t *ompmem = NULL;
  
  cdoInitialize(argument);

  cdoOperatorAdd("bandpass",  BANDPASS,  0, NULL);
  cdoOperatorAdd("highpass",  HIGHPASS,  0, NULL);
  cdoOperatorAdd("lowpass" ,  LOWPASS,   0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  if ( CDO_Use_FFTW )
    {
#if defined(HAVE_LIBFFTW3) 
      if ( cdoVerbose ) cdoPrint("Using fftw3 lib");
      use_fftw = TRUE;
#else
      if ( cdoVerbose ) cdoPrint("LIBFFTW3 support not compiled in!");
#endif
    }
      
  if ( cdoVerbose && use_fftw  == FALSE ) cdoPrint("Using intrinsic FFT function!");
  
  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  calendar = taxisInqCalendar(taxisID1);  
 
  nvars = vlistNvars(vlistID1);
  
  tsID = 0;    
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
        {
          nalloc += NALLOC_INC;
          vars   = (field_t ***) realloc(vars, nalloc*sizeof(field_t **));
        }
                       
      dtlist_taxisInqTimestep(dtlist, taxisID1, tsID);
   
      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);
           
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          gridID   = vlistInqVarGrid(vlistID1, varID);
          gridsize = gridInqSize(gridID);
          vars[tsID][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
          streamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;
          if ( nmiss ) cdoAbort("Missing value support for operators in module Filter not added yet!");
        }

      /* get and check time increment */                   
      if ( tsID > 0 )
        {    
          int incperiod = 0;
          int year, month, day;
	  int vdate0 = dtlist_get_vdate(dtlist, tsID-1);
	  int vdate  = dtlist_get_vdate(dtlist, tsID);
	  int vtime0 = dtlist_get_vtime(dtlist, tsID-1);
	  int vtime  = dtlist_get_vtime(dtlist, tsID);

          cdiDecodeDate(vdate0, &year0, &month0, &day0);               
          cdiDecodeDate(vdate,  &year,  &month,  &day);

          juldate_t juldate0 = juldate_encode(calendar, vdate0, vtime0);        
          juldate_t juldate  = juldate_encode(calendar, vdate,  vtime);         
          double jdelta   = juldate_to_seconds(juldate_sub(juldate, juldate0));
          
          if ( tsID == 1 ) 
            {           
              getTimeInc(jdelta, vdate0, vdate, &incperiod0, &incunit0);
              incperiod = incperiod0; 
              if ( incperiod == 0 ) cdoAbort("Time step must be different from zero!");
              incunit = incunit0;
              if ( cdoVerbose ) cdoPrint("Time step %i %s", incperiod, tunits[incunit]);
              fdata = 1.*iunits[incunit]/incperiod;
            }
          else 
            getTimeInc(jdelta, vdate0, vdate, &incperiod, &incunit);        

          if ( calendar != CALENDAR_360DAYS && calendar != CALENDAR_365DAYS && calendar != CALENDAR_366DAYS &&
	       incunit0 < 4 && month == 2 && day == 29 && 
               ( day0 != day || month0 != month || year0 != year ) )
            {
              cdoWarning("Filtering of multi-year times series doesn't works properly with a standard calendar.");
              cdoWarning("  Please delete the day %i-02-29 (cdo del29feb)", year);
            }

          if ( ! ( incperiod == incperiod0 && incunit == incunit0 ) )
            cdoWarning("Time increment in step %i (%d%s) differs from step 1 (%d%s)!",
                       tsID, incperiod, tunits[incunit], incperiod0, tunits[incunit0]);        
        }
      tsID++;
    }
  
  nts = tsID;
  if ( nts <= 1 ) cdoAbort("Number of time steps <= 1!");

  if ( use_fftw )
    {
#if defined(HAVE_LIBFFTW3) 
      ompmem = (memory_t*) malloc(ompNumThreads*sizeof(memory_t));
      for ( i = 0; i < ompNumThreads; i++ )
	{
	  ompmem[i].in_fft  = (fftw_complex*) malloc(nts*sizeof(fftw_complex));
	  ompmem[i].out_fft = (fftw_complex*) malloc(nts*sizeof(fftw_complex));
	  ompmem[i].p_T2S = fftw_plan_dft_1d(nts, ompmem[i].in_fft, ompmem[i].out_fft,  1, FFTW_ESTIMATE);
	  ompmem[i].p_S2T = fftw_plan_dft_1d(nts, ompmem[i].out_fft, ompmem[i].in_fft, -1, FFTW_ESTIMATE);
	}
#endif
    }
  else
    {
      ompmem = (memory_t*) malloc(ompNumThreads*sizeof(memory_t));
      for ( i = 0; i < ompNumThreads; i++ )
	{
	  ompmem[i].array1 = (double*) malloc(nts*sizeof(double));
	  ompmem[i].array2 = (double*) malloc(nts*sizeof(double));
	}
    }

  int *fmasc = (int*) calloc(nts, sizeof(int));

  switch(operfunc)
    {
    case BANDPASS: 
      {
        operatorInputArg("lower and upper bound of frequency band");
        operatorCheckArgc(2);
        fmin = parameter2double(operatorArgv()[0]);
        fmax = parameter2double(operatorArgv()[1]);
        break;
      }
    case HIGHPASS:
      {              
        operatorInputArg("lower bound of frequency pass");
        operatorCheckArgc(1);
        fmin = parameter2double(operatorArgv()[0]);
        fmax = fdata;
        break;
      }
    case LOWPASS: 
      { 
        operatorInputArg("upper bound of frequency pass");
        operatorCheckArgc(1);
        fmin = 0;
        fmax = parameter2double(operatorArgv()[0]);
        break;
      }
    }

  if ( cdoVerbose ) cdoPrint("fmin=%g  fmax=%g", fmin, fmax);
  
  create_fmasc(nts, fdata, fmin, fmax, fmasc);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      
      for ( levelID = 0; levelID < nlevel; levelID++ )
        {
          if ( use_fftw )
            {
#if defined(HAVE_LIBFFTW3) 
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, tsID)
#endif
              for ( i = 0; i < gridsize; i++ )
                {
            	  int ompthID = cdo_omp_get_thread_num();

                  for ( tsID = 0; tsID < nts; tsID++ )                              
                    {
                      ompmem[ompthID].in_fft[tsID][0] = vars[tsID][varID][levelID].ptr[i];
                      ompmem[ompthID].in_fft[tsID][1] = 0;
                    }

                  filter_fftw(nts, fmasc, ompmem[ompthID].out_fft, &ompmem[ompthID].p_T2S, &ompmem[ompthID].p_S2T);
                  
                  for ( tsID = 0; tsID < nts; tsID++ )
                    {
                      vars[tsID][varID][levelID].ptr[i] = ompmem[ompthID].in_fft[tsID][0] / nts;  
                    }
                }
#endif
            }
          else
            {
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, tsID)
#endif
              for ( i = 0; i < gridsize; i++ )  
                {
            	  int ompthID = cdo_omp_get_thread_num();

                  for ( tsID = 0; tsID < nts; tsID++ )
                    ompmem[ompthID].array1[tsID] = vars[tsID][varID][levelID].ptr[i];

                  memset(ompmem[ompthID].array2, 0, nts*sizeof(double));

                  filter_intrinsic(nts, fmasc, ompmem[ompthID].array1, ompmem[ompthID].array2);

                  for ( tsID = 0; tsID < nts; tsID++ )
                    vars[tsID][varID][levelID].ptr[i] = ompmem[ompthID].array1[tsID];
                }
            }
        }
    }

  if ( use_fftw )
    {
#if defined(HAVE_LIBFFTW3) 
      for ( i = 0; i < ompNumThreads; i++ )
	{
	  free(ompmem[i].in_fft);
	  free(ompmem[i].out_fft);
	}
      free(ompmem);
#endif
    }
  else
    {
      for ( i = 0; i < ompNumThreads; i++ )
	{
	  free(ompmem[i].array1);
	  free(ompmem[i].array2);
	}
      free(ompmem);
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  
  streamDefVlist(streamID2, vlistID2);
 
  for ( tsID = 0; tsID < nts; tsID++ )
    {
      dtlist_taxisDefTimestep(dtlist, taxisID2, tsID);
      streamDefTimestep(streamID2, tsID);
    
      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID = 0; levelID < nlevel; levelID++ )
            {
              if ( vars[tsID][varID][levelID].ptr )
                {
                  nmiss = vars[tsID][varID][levelID].nmiss;
                  streamDefRecord(streamID2, varID, levelID);
                  streamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);

                  free(vars[tsID][varID][levelID].ptr);
                  vars[tsID][varID][levelID].ptr = NULL;
                }
            }
        }

      field_free(vars[tsID], vlistID1);
    }

  if ( vars   ) free(vars);

  dtlist_delete(dtlist);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();
  
  return (0);
}
