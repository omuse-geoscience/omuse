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

      Rhopot      rhopot          potential density
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

/* compute density from insitu temperature */
static
double potrho_1(double t, double sal, double p)
{
  double r_a0 = 999.842594, r_a1 = 6.793952e-2, r_a2 = -9.095290e-3,
         r_a3 = 1.001685e-4, r_a4 = -1.120083e-6, r_a5 = 6.536332e-9,
         r_b0 = 8.24493e-1, r_b1 = -4.0899e-3, r_b2 = 7.6438e-5,
         r_b3 = -8.2467e-7, r_b4 = 5.3875e-9,
         r_c0 = -5.72466e-3, r_c1 = 1.0227e-4, r_c2 = -1.6546e-6,
         r_d0 = 4.8314e-4,
         r_e0 = 19652.21, r_e1 = 148.4206, r_e2 = -2.327105,
         r_e3 = 1.360477e-2, r_e4 = -5.155288e-5,
         r_f0 = 54.6746, r_f1 = -0.603459, r_f2 = 1.09987e-2,
         r_f3 = -6.1670e-5,
         r_g0 = 7.944e-2, r_g1 = 1.6483e-2, r_g2 = -5.3009e-4,
         r_h0 = 3.239908, r_h1 = 1.43713e-3, r_h2 = 1.16092e-4,
         r_h3 = -5.77905e-7,
         r_ai0 = 2.2838e-3, r_ai1 = -1.0981e-5, r_ai2 = -1.6078e-6,
         r_aj0 = 1.91075e-4,
         r_ak0 = 8.50935e-5, r_ak1 = -6.12293e-6, r_ak2 = 5.2787e-8,
         r_am0 = -9.9348e-7, r_am1 = 2.0816e-8, r_am2 = 9.1697e-10;

  double s, s3h; 
  double rho;

    {
      s = MAX(sal, 0.0);
      s3h = sqrt(s*s*s);

      rho = (r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))
            + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))    
            + r_d0 * s*s                 
            + s3h * (r_c0 + t * (r_c1 + r_c2 * t)))                           
           / (1.                                                            
             - p / (p * (r_h0 + t * (r_h1 + t * (r_h2 + t * r_h3))            
                         + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))              
                         + r_aj0 * s3h                                        
                         + (r_ak0 + t * (r_ak1 + t * r_ak2)                   
                         + s * (r_am0 + t * (r_am1 + t * r_am2))) * p)        
                    + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))  
                    + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))         
                    + s3h * (r_g0 + t * (r_g1 + r_g2 * t))));
    }

  return (rho);
}

/*
#define N 4
int main (int argc, char *argv[])
{
  int i;
  {
    double p    = 0;
    double t[N] = {22, 25, 28, 31};
    double s[N] = {35, 35, 35, 35};
    double x[N] = {24.219, 23.343, 22.397, 21.384};
    double r[N];
  
    potrho_1d(t, s, p, r, N);

    for ( i = 0; i < N; ++i )
      printf("%d %5g %3g %8g %8g %8g %10.3f\n", i, p, s[i], t[i], x[i], r[i], r[i]-x[i]);
  }

  {
    double p    = 300;
    double t[N] = {-2.140, -0.186, 1.771, 3.728};
    double s[N] = {35, 35, 35, 35};
    double x[N] = {42.191, 41.941, 41.649, 41.319};
    double r[N];
  
    potrho_1d(t, s, p, r, N);

    for ( i = 0; i < N; ++i )
      printf("%d %5g %3g %8g %8g %8g %10.3f\n", i, p, s[i], t[i], x[i], r[i], r[i]-x[i]);
  }

  return (0);
}
*/

static
void calc_rhopot(long gridsize, long nlevel, double *pressure, field_t to, field_t sao, field_t rho)
{
  /* pressure units: hPa     */
  /* to units:       Celsius */
  /* sao units:      psu     */

  long i, levelID, offset;
  double *rhoptr, *toptr, *saoptr;

  for ( levelID = 0; levelID < nlevel; ++levelID )
    {
      offset = gridsize*levelID;
      toptr = to.ptr + offset;
      saoptr = sao.ptr + offset;
      rhoptr = rho.ptr + offset;

      for ( i = 0; i < gridsize; ++i )
	{
	  if ( DBL_IS_EQUAL(toptr[i], to.missval) ||
	       DBL_IS_EQUAL(saoptr[i], sao.missval) )
	    {
	      rhoptr[i] = rho.missval;
	    }
	  else
	    {
	      rhoptr[i] = potrho_1(toptr[i], saoptr[i], pressure[levelID]); 
	    }
	}
    }
}


void *Rhopot(void *argument)
{
  int streamID1, streamID2;
  int nrecs;
  int tsID, recID, varID, levelID;
  int nlevel1, nlevel2;
  int nvars, code, zaxisID;
  int vlistID1, vlistID2;
  int offset;
  int nlevel;
  int i;
  int nmiss;
  int toID = -1, saoID = -1, thoID = -1;
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  int taxisID1, taxisID2;
  double pin = -1;
  double *pressure;
  double *single;
  field_t to, sao, rho;

  cdoInitialize(argument);

  if ( operatorArgc() == 1 ) pin = parameter2double(operatorArgv()[0]);
  
  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      code = vlistInqVarCode(vlistID1, varID);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  vlistInqVarStdname(vlistID1,varID, stdname);
	  strtolower(varname);

	  if      ( strcmp(varname, "to")    == 0 ) code = 20;
	  else if ( strcmp(varname, "sao")   == 0 ) code =  5;
	  else if ( strcmp(varname, "tho")   == 0 ) code =  2;

	  else if ( strcmp(varname, "s")     == 0 ) code = 5;
	  else if ( strcmp(varname, "t")     == 0 ) code = 2;

	  else if ( strcmp(stdname, "sea_water_salinity")              == 0 ) code = 5;
	  else if ( strcmp(stdname, "sea_water_potential_temperature") == 0 ) code = 2;
	}

      if      ( code == 20 ) toID  = varID;
      else if ( code ==  5 ) saoID = varID;
      else if ( code ==  2 ) thoID = varID;
    }

  if ( saoID == -1 ) cdoAbort("Sea water salinity not found!");
  if ( toID  == -1 && thoID != -1 )
    {   
      cdoPrint("Use the CDO operator 'adisit' to convert potential temperature to In-situ temperature.");
      cdoPrint("Here is an example:");
      cdoPrint("   cdo rhopot -adisit %s %s", cdoStreamName(0)->args, cdoStreamName(1)->args);
    }
  if ( toID  == -1 ) cdoAbort("In-situ temperature not found!");

  int gridID = vlistGrid(vlistID1, 0);
  int gridsize = vlist_check_gridsize(vlistID1);

  zaxisID = vlistInqVarZaxis(vlistID1, saoID);
  nlevel1 = zaxisInqSize(zaxisID);
  zaxisID = vlistInqVarZaxis(vlistID1, toID);
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

  field_init(&to);
  field_init(&sao);
  field_init(&rho);
  to.ptr = (double*) malloc(gridsize*nlevel*sizeof(double));
  sao.ptr = (double*) malloc(gridsize*nlevel*sizeof(double));
  rho.ptr = (double*) malloc(gridsize*nlevel*sizeof(double));

  to.nmiss = 0;
  sao.nmiss = 0;
  rho.nmiss = 0;
  
  to.missval = vlistInqVarMissval(vlistID1, toID);
  sao.missval = vlistInqVarMissval(vlistID1, saoID);
  rho.missval = to.missval;


  vlistID2 = vlistCreate();
  varID = vlistDefVar(vlistID2, gridID, zaxisID, TSTEP_INSTANT);
  vlistDefVarParam(vlistID2, varID, cdiEncodeParam(18, 255, 255));
  vlistDefVarName(vlistID2, varID, "rhopoto");
  vlistDefVarLongname(vlistID2, varID, "Sea water potential density");
  vlistDefVarStdname(vlistID2, varID, "sea_water_potential_density");
  vlistDefVarUnits(vlistID2, varID, "kg m-3");
  vlistDefVarMissval(vlistID2, varID, rho.missval);


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

	  if ( varID == toID ) streamReadRecord(streamID1, to.ptr+offset, &(to.nmiss));
	  if ( varID == saoID ) streamReadRecord(streamID1, sao.ptr+offset, &(sao.nmiss));
	}

      calc_rhopot(gridsize, nlevel, pressure, to, sao, rho); 

      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  offset = gridsize*levelID;
	  single = rho.ptr+offset;

	  nmiss = 0;
	  for ( i = 0; i < gridsize; ++i )
	    if ( DBL_IS_EQUAL(single[i], rho.missval) ) nmiss++;
 
	  streamDefRecord(streamID2, 0, levelID);
	  streamWriteRecord(streamID2, single, nmiss);     
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  free(pressure);
  free(rho.ptr);
  free(to.ptr);
  free(sao.ptr);

  cdoFinish();

  return (0);
}
