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

      Mrotuv      mrotuv          Forward rotation for MPIOM data
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void rotate_uv(double *u_i, double *v_j, int ix, int iy,
	       double *lon, double *lat, double *u_lon, double *v_lat)
{
  /*
 real,intent(in)      :: u_i(ix,iy,iz),v_j(ix,iy,iz)                  ! vector components in i-j-direction
 real,intent(out)     :: u_lon(ix,iy,iz),v_lat(ix,iy,iz)              ! vector components in lon-lat direction
 real,intent(in)      :: lat(ix,iy),lon(ix,iy)                        ! latitudes and longitudes
  */
  double dlat_i, dlat_j,dlon_i,dlon_j,dist_i,dist_j;
  double lat_factor;
  double absold, absnew;  /* velocity vector lengths */
  int  i, j, ip1, im1, jp1, jm1;
  int  change_sign_u, change_sign_v;
  double pi = 3.14159265359;


  /* specification whether change in sign is needed for the input arrays */
  change_sign_u = FALSE;
  change_sign_v = TRUE;

  /* initialization */
  for ( i = 0; i < ix*iy; i++ )
    {
      v_lat[i] = 0;
      u_lon[i] = 0;
    }

  /* rotation */
  for ( j = 0; j < iy; j++ )
    for ( i = 0; i < ix; i++ )
      {
	ip1 = i + 1;
	im1 = i - 1;
	jp1 = j + 1;
	jm1 = j - 1;
	if ( ip1 >= ix ) ip1 = 0;   /* the 0-meridian */
	if ( im1 <   0 ) im1 = ix-1;
	if ( jp1 >= iy ) jp1 = j;   /* treatment of the last..  */
	if ( jm1 <   0 ) jm1 = j;   /* .. and the fist grid-row */

	/* difference in latitudes */
	dlat_i = lat[IX2D(j,ip1,ix)] - lat[IX2D(j,im1,ix)];
	dlat_j = lat[IX2D(jp1,i,ix)] - lat[IX2D(jm1,i,ix)];

	/* difference in longitudes */
	dlon_i = lon[IX2D(j,ip1,ix)] - lon[IX2D(j,im1,ix)];
	if ( dlon_i >   pi  ) dlon_i -= 2*pi;
	if ( dlon_i < (-pi) ) dlon_i += 2*pi;
	dlon_j = lon[IX2D(jp1,i,ix)] - lon[IX2D(jm1,i,ix)];
	if ( dlon_j >   pi  ) dlon_j -= 2*pi;
	if ( dlon_j < (-pi) ) dlon_j += 2*pi;

	lat_factor = cos(lat[IX2D(j,i,ix)]);
	dlon_i = dlon_i * lat_factor;
	dlon_j = dlon_j * lat_factor;

	/* projection by scalar product */
	u_lon[IX2D(j,i,ix)] = u_i[IX2D(j,i,ix)]*dlon_i + v_j[IX2D(j,i,ix)]*dlat_i;
	v_lat[IX2D(j,i,ix)] = u_i[IX2D(j,i,ix)]*dlon_j + v_j[IX2D(j,i,ix)]*dlat_j;

	dist_i = sqrt(dlon_i*dlon_i + dlat_i*dlat_i);
	dist_j = sqrt(dlon_j*dlon_j + dlat_j*dlat_j);

	if ( fabs(dist_i) > 0 && fabs(dist_j) > 0 )
	  {
	    u_lon[IX2D(j,i,ix)] /= dist_i;
	    v_lat[IX2D(j,i,ix)] /= dist_j;
	  }
	else
	  {
	    u_lon[IX2D(j,i,ix)] = 0.0;
	    v_lat[IX2D(j,i,ix)] = 0.0;
	  }

	absold = sqrt(u_i[IX2D(j,i,ix)]*u_i[IX2D(j,i,ix)] + v_j[IX2D(j,i,ix)]*v_j[IX2D(j,i,ix)]);
	absnew = sqrt(u_lon[IX2D(j,i,ix)]*u_lon[IX2D(j,i,ix)] + v_lat[IX2D(j,i,ix)]*v_lat[IX2D(j,i,ix)]);

	u_lon[IX2D(j,i,ix)] *= absold;
	v_lat[IX2D(j,i,ix)] *= absold;

	if ( absnew > 0 )
	  {
	    u_lon[IX2D(j,i,ix)] /= absnew;
	    v_lat[IX2D(j,i,ix)] /= absnew;
	  }
	else
	  {
	    u_lon[IX2D(j,i,ix)] = 0.0;
	    v_lat[IX2D(j,i,ix)] = 0.0;
	  }

	/* change sign */
	if ( change_sign_u ) u_lon[IX2D(j,i,ix)] *= -1;
	if ( change_sign_v ) v_lat[IX2D(j,i,ix)] *= -1;

	if ( cdoVerbose )
	  {
	    absold = sqrt(u_i[IX2D(j,i,ix)]*u_i[IX2D(j,i,ix)] + v_j[IX2D(j,i,ix)]*v_j[IX2D(j,i,ix)]);
	    absnew = sqrt(u_lon[IX2D(j,i,ix)]*u_lon[IX2D(j,i,ix)] + v_lat[IX2D(j,i,ix)]*v_lat[IX2D(j,i,ix)]);

	    if ( i%20 == 0 && j%20 == 0 && absold > 0 )
	      {
		printf("(absold,absnew) %d %d %g %g %g %g %g %g\n",
		       j+1, i+1, absold, absnew, u_i[IX2D(j,i,ix)], v_j[IX2D(j,i,ix)], u_lon[IX2D(j,i,ix)], v_lat[IX2D(j,i,ix)]);
		
		/* test orthogonality */
		if ( (dlon_i*dlon_j + dlat_j*dlat_i) > 0.1 )            
		  fprintf(stderr, "orthogonal? %d %d %g\n", j+1, i+1, (dlon_i*dlon_j + dlat_j*dlat_i));
	      }
	  }
      }
}


void p_to_uv_grid(int nlon, int nlat, double *grid1x, double *grid1y,
		  double *gridux, double *griduy, double *gridvx, double *gridvy)
{
  int i, j, jp1, ip1;

  
  /* interpolate scalar to u points */
  for ( j = 0; j < nlat; j++ )
    for ( i = 0; i < nlon; i++ )
      {
	ip1 = i + 1;
	if ( ip1 > nlon-1 ) ip1 = 0;

	gridux[IX2D(j,i,nlon)] = (grid1x[IX2D(j,i,nlon)]+grid1x[IX2D(j,ip1,nlon)])*0.5;
	if ( (grid1x[IX2D(j,i,nlon)] > 340 && grid1x[IX2D(j,ip1,nlon)] <  20) ||
             (grid1x[IX2D(j,i,nlon)] < 20  && grid1x[IX2D(j,ip1,nlon)] > 340) )
	  {
	    if ( gridux[IX2D(j,i,nlon)] < 180 )
	      gridux[IX2D(j,i,nlon)] += 180;
	    else
	      gridux[IX2D(j,i,nlon)] -= 180;
	  }
	    
	griduy[IX2D(j,i,nlon)] = (grid1y[IX2D(j,i,nlon)]+grid1y[IX2D(j,ip1,nlon)])*0.5;
      }

  
  /* interpolate scalar to v points */
  for ( j = 0; j < nlat; j++ )
    for ( i = 0; i < nlon; i++ )
      {
	jp1 = j + 1;
	if ( jp1 > nlat-1 ) jp1 = nlat-1;

	gridvx[IX2D(j,i,nlon)] = (grid1x[IX2D(j,i,nlon)]+grid1x[IX2D(jp1,i,nlon)])*0.5;
	if ( (grid1x[IX2D(j,i,nlon)] > 340 && grid1x[IX2D(jp1,i,nlon)] <  20) ||
             (grid1x[IX2D(j,i,nlon)] < 20  && grid1x[IX2D(jp1,i,nlon)] > 340) )
	  {
	    if ( gridvx[IX2D(j,i,nlon)] < 180 )
	      gridvx[IX2D(j,i,nlon)] += 180;
	    else
	      gridvx[IX2D(j,i,nlon)] -= 180;
	  }
	    
	gridvy[IX2D(j,i,nlon)] = (grid1y[IX2D(j,i,nlon)]+grid1y[IX2D(jp1,i,nlon)])*0.5;
      }
}


void *Mrotuv(void *argument)
{
  int streamID1, streamID2, streamID3;
  int nrecs;
  int tsID, recID, levelID;
  int varID, varid;
  int lid, nlevs, code;
  int nvars;
  int gridID1, gridID2, gridIDu, gridIDv;
  int gridsize, gridsizex;
  int nlon, nlat;
  int vlistID1, vlistID2, vlistID3;
  int i, j;
  int taxisID1, taxisID2, taxisID3;
  int nmiss1, nmiss2;
  int uid = -1, vid = -1;
  double missval1, missval2;
  double *ufield = NULL, *vfield = NULL;
  double **urfield = NULL, **vrfield = NULL;
  double *uhelp = NULL, *vhelp = NULL;
  double *grid1x = NULL, *gridux = NULL, *gridvx = NULL;
  double *grid1y = NULL, *griduy = NULL, *gridvy = NULL;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars = vlistNvars(vlistID1);
  for ( varid = 0; varid < nvars; varid++ )
    {
      code = vlistInqVarCode(vlistID1, varid);
      if ( code == 3 || code == 131 ) uid = varid;
      if ( code == 4 || code == 132 ) vid = varid;
    }

  if ( uid == -1 || vid == -1 )
    {
      if ( nvars == 2 )
	{
	  uid = 0;
	  vid = 1;
	}
      else
	cdoAbort("U and V not found in %s",  cdoStreamName(0)->args);
    }

  nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, uid));
  if ( nlevs != zaxisInqSize(vlistInqVarZaxis(vlistID1, vid)) )
    cdoAbort("U and V have different number of levels!");

  gridID1 = vlistInqVarGrid(vlistID1, uid);
  gridID2 = vlistInqVarGrid(vlistID1, vid);
  gridsize = gridInqSize(gridID1);
  if ( gridID1 != gridID2 ) cdoAbort("Input grids differ!");

  if ( gridInqType(gridID1) != GRID_LONLAT      &&
       gridInqType(gridID1) != GRID_GAUSSIAN    &&
       gridInqType(gridID1) != GRID_CURVILINEAR )
    cdoAbort("Grid %s unsupported!", gridNamePtr(gridInqType(gridID1)));
  
  if ( gridInqType(gridID1) != GRID_CURVILINEAR )
    gridID1 = gridToCurvilinear(gridID1, 0);

  if ( gridsize != gridInqSize(gridID1) ) cdoAbort("Internal problem: gridsize changed!");

  nlon    = gridInqXsize(gridID1);
  nlat    = gridInqYsize(gridID1);

  grid1x  = (double*) malloc(gridsize*sizeof(double));
  grid1y  = (double*) malloc(gridsize*sizeof(double));
  gridux  = (double*) malloc(gridsize*sizeof(double));
  griduy  = (double*) malloc(gridsize*sizeof(double));
  gridvx  = (double*) malloc(gridsize*sizeof(double));
  gridvy  = (double*) malloc(gridsize*sizeof(double));

  gridsizex = (nlon+2)*nlat;

  gridInqXvals(gridID1, grid1x);
  gridInqYvals(gridID1, grid1y);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID1, units);
    grid_to_degree(units, gridsize, grid1x, "grid center lon");
    gridInqYunits(gridID1, units);
    grid_to_degree(units, gridsize, grid1y, "grid center lat");
  }

  p_to_uv_grid(nlon, nlat, grid1x, grid1y, gridux, griduy, gridvx, gridvy);

  gridIDu = gridCreate(GRID_CURVILINEAR, nlon*nlat);
  gridDefPrec(gridIDu, gridInqPrec(gridID1));
  gridDefXsize(gridIDu, nlon);
  gridDefYsize(gridIDu, nlat);
  gridDefXvals(gridIDu, gridux);
  gridDefYvals(gridIDu, griduy);

  gridIDv = gridCreate(GRID_CURVILINEAR, nlon*nlat);
  gridDefPrec(gridIDv, gridInqPrec(gridID1));
  gridDefXsize(gridIDv, nlon);
  gridDefYsize(gridIDv, nlat);
  gridDefXvals(gridIDv, gridvx);
  gridDefYvals(gridIDv, gridvy);

  for ( i = 0; i < gridsize; i++ )
    {
      grid1x[i] *= DEG2RAD;
      grid1y[i] *= DEG2RAD;
    }

  vlistClearFlag(vlistID1);
  for ( lid = 0; lid < nlevs; lid++ ) vlistDefFlag(vlistID1, uid, lid, TRUE);
  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);
  vlistChangeVarGrid(vlistID2, 0, gridIDu);

  vlistClearFlag(vlistID1);
  for ( lid = 0; lid < nlevs; lid++ ) vlistDefFlag(vlistID1, vid, lid, TRUE);
  vlistID3 = vlistCreate();
  vlistCopyFlag(vlistID3, vlistID1);
  vlistChangeVarGrid(vlistID3, 0, gridIDv);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  missval1 = vlistInqVarMissval(vlistID1, uid);
  missval2 = vlistInqVarMissval(vlistID1, vid);

  ufield  = (double*) malloc(gridsize*sizeof(double));
  vfield  = (double*) malloc(gridsize*sizeof(double));

  urfield  = (double**) malloc(nlevs*sizeof(double*));
  vrfield  = (double**) malloc(nlevs*sizeof(double*));
  for ( lid = 0; lid < nlevs; lid++ )
    {
      urfield[lid] = (double*) malloc(gridsize*sizeof(double));
      vrfield[lid] = (double*) malloc(gridsize*sizeof(double));
    }

  uhelp   = (double*) malloc(gridsizex*sizeof(double));
  vhelp   = (double*) malloc(gridsizex*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsID);
      taxisCopyTimestep(taxisID3, taxisID1);
      streamDefTimestep(streamID3, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( varID == uid ) streamReadRecord(streamID1, urfield[levelID], &nmiss1);
	  if ( varID == vid ) streamReadRecord(streamID1, vrfield[levelID], &nmiss2);
	}

      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  /* remove missing values */
	  if ( nmiss1 || nmiss2 )
	    {
	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( DBL_IS_EQUAL(urfield[levelID][i], missval1) ) urfield[levelID][i] = 0;
		  if ( DBL_IS_EQUAL(vrfield[levelID][i], missval2) ) vrfield[levelID][i] = 0;
		}
	    }

	  /* rotate*/
	  
	  rotate_uv(urfield[levelID], vrfield[levelID], nlon, nlat, grid1x, grid1y, ufield, vfield);

	  /* load to a help field */
	  for ( j = 0; j < nlat; j++ )
	    for ( i = 0; i < nlon; i++ )
	      {
		uhelp[IX2D(j,i+1,nlon+2)] = ufield[IX2D(j,i,nlon)];
		vhelp[IX2D(j,i+1,nlon+2)] = vfield[IX2D(j,i,nlon)];
	      }

	  /* make help field cyclic */
	  for ( j = 0; j < nlat; j++ )
	    {
	      uhelp[IX2D(j,0,nlon+2)]      = uhelp[IX2D(j,nlon,nlon+2)];
	      uhelp[IX2D(j,nlon+1,nlon+2)] = uhelp[IX2D(j,1,nlon+2)];
	      vhelp[IX2D(j,0,nlon+2)]      = vhelp[IX2D(j,nlon,nlon+2)];
	      vhelp[IX2D(j,nlon+1,nlon+2)] = vhelp[IX2D(j,1,nlon+2)];
	    }

	  /* interpolate on u/v points */
	  for ( j = 0; j < nlat; j++ )
	    for ( i = 0; i < nlon; i++ )
	      {
		ufield[IX2D(j,i,nlon)] = (uhelp[IX2D(j,i+1,nlon+2)]+uhelp[IX2D(j,i+2,nlon+2)])*0.5;
	      }

	  for ( j = 0; j < nlat-1; j++ )
	    for ( i = 0; i < nlon; i++ )
	      {
		vfield[IX2D(j,i,nlon)] = (vhelp[IX2D(j,i+1,nlon+2)]+vhelp[IX2D(j+1,i+1,nlon+2)])*0.5;
	      }

	  for ( i = 0; i < nlon; i++ )
	    {
	      vfield[IX2D(nlat-1,i,nlon)] = vhelp[IX2D(nlat-1,i+1,nlon+2)];
	    }

	  streamDefRecord(streamID2, 0, levelID);
	  streamWriteRecord(streamID2, ufield, nmiss1);     
	  streamDefRecord(streamID3, 0, levelID);
	  streamWriteRecord(streamID3, vfield, nmiss2);     
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( ufield  ) free(ufield);
  if ( vfield  ) free(vfield);
  if ( urfield ) free(urfield);
  if ( vrfield ) free(vrfield);
  if ( uhelp   ) free(uhelp);
  if ( vhelp   ) free(vhelp);
  if ( gridux  ) free(gridux);
  if ( griduy  ) free(griduy);
  if ( gridvx  ) free(gridvx);
  if ( gridvy  ) free(gridvy);

  cdoFinish();

  return (0);
}
