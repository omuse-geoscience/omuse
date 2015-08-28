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

      Adisit      adisit          compute insitu from potential temperature
      Adisit      adipot          compute potential from insitu temperature
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


/*
!>
!! transformation from potential to in situ temperature
!! according to Bryden, 1973, "New polynomials for thermal expansion,
!! adiabatic temperature gradient and potential temperature of sea
!! water". Deep Sea Research and Oceanographic Abstracts. 20, 401-408
!! (GILL P.602), which gives the inverse transformation for an
!! approximate value, all terms linear in t are taken after that one
!! newton step.  for the check value 8.4678516 the accuracy is 0.2
!! mikrokelvin.
!!
*/

/* compute insitu temperature from potential temperature */
static
double adisit_1(double tpot, double sal, double p)
{
  double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9,
         a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
         a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10,
         a_d = 4.1057E-9,
         a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double dc, dv, dvs, fne, fst, qc, qn3, qnq, qv, qvs, t, tpo;

  qc = p * (a_a1 + p * (a_c1 - a_e1 * p));
  qv = p * (a_b1 - a_d * p);
  dc = 1. + p * (-a_a2 + p * (a_c2 - a_e2 * p));
  dv = a_b2 * p;
  qnq  = -p * (-a_a3 + p * a_c3);
  qn3  = -p * a_a4;

    {
      tpo = tpot;
      qvs = qv*(sal - 35.) + qc;
      dvs = dv*(sal - 35.) + dc;
      t   = (tpo + qvs)/dvs;
      fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpo;
      fst = dvs + t*(2.*qnq + 3.*qn3*t);
      t = t - fne/fst;
    }

  return (t);
}

/* compute potential temperature from insitu temperature */
/* Ref: Gill, p. 602, Section A3.5:Potential Temperature */
static
double adipot(double t, double s, double p)
{
  double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9,
         a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
         a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10,
         a_d = 4.1057E-9,
         a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double aa,bb,cc,cc1,dd,tpot,s_rel;

  s_rel = s - 35.0;

   aa = (a_a1+ t*(a_a2 - t*(a_a3 - a_a4*t)));
   bb = s_rel*(a_b1 -a_b2*t)     ;
   cc = (a_c1 + t*(-a_c2 + a_c3*t));
  cc1 = a_d*s_rel;
   dd = (-a_e1 + a_e2*t);

   tpot=t-p*(aa + bb + p*(cc - cc1 + p*dd));

  return (tpot);
}

static
void calc_adisit(long gridsize, long nlevel, double *pressure, field_t tho, field_t sao, field_t tis)
{
  /* pressure units: hPa     */
  /* tho units:      Celsius */
  /* sao units:      psu     */

  long i, levelID, offset;
  double *tisptr, *thoptr, *saoptr;

  for ( levelID = 0; levelID < nlevel; ++levelID )
    {
      offset = gridsize*levelID;
      thoptr = tho.ptr + offset;
      saoptr = sao.ptr + offset;
      tisptr = tis.ptr + offset;

      for ( i = 0; i < gridsize; ++i )
	{
	  if ( DBL_IS_EQUAL(thoptr[i], tho.missval) ||
	       DBL_IS_EQUAL(saoptr[i], sao.missval) )
	    {
	      tisptr[i] = tis.missval;
	    }
	  else
	    {
	      tisptr[i] = adisit_1(thoptr[i], saoptr[i], pressure[levelID]); 
	    }
	}
    }
}

static
void calc_adipot(long gridsize, long nlevel, double *pressure, field_t t, field_t s, field_t tpot)
{
  /* pressure units: hPa     */
  /* t units:      Celsius */
  /* s units:      psu     */

  long i, levelID, offset;
  double *tpotptr, *tptr, *sptr;

  for ( levelID = 0; levelID < nlevel; ++levelID )
    {
      offset = gridsize*levelID;
      tptr = t.ptr + offset;
      sptr = s.ptr + offset;
      tpotptr = tpot.ptr + offset;

      for ( i = 0; i < gridsize; ++i )
	{
	  if ( DBL_IS_EQUAL(tptr[i], t.missval) ||
	       DBL_IS_EQUAL(sptr[i], s.missval) )
	    {
	      tpotptr[i] = tpot.missval;
	    }
	  else
	    {
	      tpotptr[i] = adipot(tptr[i], sptr[i], pressure[levelID]); 
	    }
	}
    }
}


void *Adisit(void *argument)
{
  int operatorID, ADISIT, ADIPOT;
  int streamID1, streamID2;
  int nrecs;
  int tsID, recID, varID, levelID;
  int nlevel1, nlevel2;
  int gridsize;
  int nvars, code, gridID, zaxisID;
  int vlistID1, vlistID2;
  int offset;
  int nlevel;
  int i;
  int nmiss;
  int thoID = -1, saoID = -1;
  int tisID2, saoID2;
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  int taxisID1, taxisID2;
  double pin = -1;
  double *pressure;
  double *single;
  field_t tho, sao, tis;

  cdoInitialize(argument);
  ADISIT = cdoOperatorAdd("adisit", 1, 1, "");
  ADIPOT = cdoOperatorAdd("adipot", 1, 1, "");

  UNUSED(ADIPOT);

  operatorID = cdoOperatorID();

  if ( operatorArgc() == 1 ) pin = parameter2double(operatorArgv()[0]);
  
  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID1, varID);

      code = vlistInqVarCode(vlistID1, varID);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  vlistInqVarStdname(vlistID1,varID, stdname);
	  strtolower(varname);

               if ( strcmp(varname, "s")     == 0 ) code = 5;
          else if ( strcmp(varname, "t")     == 0 ) code = 2;
	  else if ( strcmp(stdname, "sea_water_salinity") == 0 ) code = 5;
          
          if ( operatorID == ADISIT )
          {
	   if ( strcmp(stdname, "sea_water_potential_temperature") == 0 ) code = 2;
          }
          else {
	   if ( strcmp(stdname, "sea_water_temperature") == 0 ) code = 2;
          }
	}

      if      ( code == 2 ) thoID = varID;
      else if ( code == 5 ) saoID = varID;
    }

  if ( saoID == -1 ) cdoAbort("Sea water salinity not found!");
  if ( thoID == -1 ) cdoAbort("Potential or Insitu temperature not found!");

  gridID = vlistGrid(vlistID1, 0);
  gridsize = vlist_check_gridsize(vlistID1);

  zaxisID = vlistInqVarZaxis(vlistID1, saoID);
  nlevel1 = zaxisInqSize(zaxisID);
  zaxisID = vlistInqVarZaxis(vlistID1, thoID);
  nlevel2 = zaxisInqSize(zaxisID);

  if ( nlevel1 != nlevel2 ) cdoAbort("temperature and salinity have different number of levels!");
  nlevel = nlevel1;

  pressure = (double*) malloc(nlevel*sizeof(double));
  zaxisInqLevels(zaxisID, pressure);

  if ( pin >= 0 ) 
    for ( i = 0; i < nlevel; ++i ) pressure[i] = pin;
  else
    for ( i = 0; i < nlevel; ++i ) pressure[i] /= 10;

  if ( cdoVerbose )
    {
      cdoPrint("Level Pressure");
      for ( i = 0; i < nlevel; ++i )
	cdoPrint("%5d  %g", i+1, pressure[i]);
    }

  field_init(&tho);
  field_init(&sao);
  field_init(&tis);
  tho.ptr = (double*) malloc(gridsize*nlevel*sizeof(double));
  sao.ptr = (double*) malloc(gridsize*nlevel*sizeof(double));
  tis.ptr = (double*) malloc(gridsize*nlevel*sizeof(double));

  tho.nmiss = 0;
  sao.nmiss = 0;
  tis.nmiss = 0;
  
  tho.missval = vlistInqVarMissval(vlistID1, thoID);
  sao.missval = vlistInqVarMissval(vlistID1, saoID);
  tis.missval = tho.missval;


  vlistID2 = vlistCreate();

  tisID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARIABLE);
  if (operatorID == ADISIT)
    {
      vlistDefVarParam(vlistID2, tisID2, cdiEncodeParam(20, 255, 255));
      vlistDefVarName(vlistID2, tisID2, "to");
      vlistDefVarLongname(vlistID2, tisID2, "Sea water temperature");
      vlistDefVarStdname(vlistID2, tisID2, "sea_water_temperature");
    }
  else
    {
      vlistDefVarParam(vlistID2, tisID2, cdiEncodeParam(2, 255, 255));
      vlistDefVarName(vlistID2, tisID2, "tho");
      vlistDefVarLongname(vlistID2, tisID2, "Sea water potential temperature");
      vlistDefVarStdname(vlistID2, tisID2, "sea_water_potential_temperature");
    }
  vlistDefVarUnits(vlistID2, tisID2, "K");
  vlistDefVarMissval(vlistID2, tisID2, tis.missval);

  saoID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARIABLE);
  vlistDefVarParam(vlistID2, saoID2, cdiEncodeParam(5, 255, 255));
  vlistDefVarName(vlistID2, saoID2, "s");
  vlistDefVarLongname(vlistID2, saoID2, "Sea water salinity");
  vlistDefVarStdname(vlistID2, saoID2, "sea_water_salinity");
  vlistDefVarUnits(vlistID2, saoID2, "psu");
  vlistDefVarMissval(vlistID2, saoID2, sao.missval);


  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; ++recID )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  offset = gridsize*levelID;

	  if ( varID == thoID ) streamReadRecord(streamID1, tho.ptr+offset, &(tho.nmiss));
	  if ( varID == saoID ) streamReadRecord(streamID1, sao.ptr+offset, &(sao.nmiss));
	}

      if (operatorID == ADISIT )
        {
          calc_adisit(gridsize, nlevel, pressure, tho, sao, tis); 
        }
      else
        {
          calc_adipot(gridsize, nlevel, pressure, tho, sao, tis); 
        }


      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  offset = gridsize*levelID;

	  single = tis.ptr+offset;

	  nmiss = 0;
	  for ( i = 0; i < gridsize; ++i )
	    if ( DBL_IS_EQUAL(single[i], tis.missval) ) nmiss++;
 
	  streamDefRecord(streamID2, tisID2, levelID);
	  streamWriteRecord(streamID2, single, nmiss);     

	  single = sao.ptr+offset;

	  nmiss = 0;
	  for ( i = 0; i < gridsize; ++i )
	    if ( DBL_IS_EQUAL(single[i], sao.missval) ) nmiss++;
 
	  streamDefRecord(streamID2, saoID2, levelID);
	  streamWriteRecord(streamID2, single, nmiss);     
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  free(pressure);
  free(tis.ptr);
  free(tho.ptr);
  free(sao.ptr);

  cdoFinish();

  return (0);
}
