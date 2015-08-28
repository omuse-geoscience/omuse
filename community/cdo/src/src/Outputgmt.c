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

*/
/*
  Output center or bounderies for GMT plotting

  Plotting example:

    - outputcenter
    - outputbounds
    - outputboundscpt
    - outputvector
*/

#if defined(HAVE_CONFIG_H)
#  include "config.h" /* VERSION */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"
#include "color.h"

double intlin(double x, double y1, double x1, double y2, double x2);

static
int pnpoly(int npol, double *xp, double *yp, double x, double y)
{
  int i, j, c = 0;

  for (i = 0, j = npol-1; i < npol; j = i++) {
    if ((((yp[i]<=y) && (y<yp[j])) ||
	 ((yp[j]<=y) && (y<yp[i]))) &&
	(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
      
      c = !c;
  }
  return c;
}


static
double PolygonArea_old(int np, double *xp, double *yp)
{
  int i, j;
  double area = 0;

  for ( i = 0; i < np; i++ )
    {
      j = (i + 1) % np;
      area += xp[i] * yp[j];
      area -= yp[i] * xp[j];
    }

  area /= 2;
  /* return(area < 0 ? -area : area); */
  return (area);
}


static
double PolygonArea(int np, double *xp, double *yp, double yc)
{
  int i, j;
  double area = 0.;

  /* Process area in Radians */
   
  for ( i = 0; i < np; i++ )
    {
      j = (i + 1) % np;
      area += DEG2RAD*xp[i] * DEG2RAD*yp[j];
      area -= DEG2RAD*yp[i] * DEG2RAD*xp[j];
    }
  area *= 0.5 * cos(DEG2RAD*yc);
  return (area);
}

static
int ccw(double p0x, double p0y, double p1x, double p1y, double p2x, double p2y)
{
  /*
    This function says wether the point are orientated clockwise
    +1 positive orientation
    -1 negative orientation
     0 points are on a line --> no orientation
    
    This is done by a comparision of the gradient of
    dy1/dx1 = p1 - p0 vs.
    dy2/dx2 = p2 - p0
    To avoid singularities at dx1=0 OR dx2 = 0 we multiply with dx1*dx2
  */
  double dx1, dx2, dy1, dy2;

  dx1 = p1x - p0x; dy1 = p1y - p0y;
  dx2 = p2x - p0x; dy2 = p2y - p0y;
  if ( dx1*dy2 > dy1*dx2 ) return +1;
  if ( dx1*dy2 < dy1*dx2 ) return -1;
  if ( (dx1*dx2 < 0 ) || (dy1*dy2 < 0)) return -1;
  if ( (dx1*dx1 + dy1*dy1) < (dx2*dx2 + dy2*dy2)) return +1;

  return 0;
}

static
int intersect(double pix, double piy, double pjx, double pjy,
              double pkx, double pky, double plx, double ply)
{
  /*This function returns if there is an intersection between the lines 
    line1 between pi and pj and
    line2 between pk and pl,
    whereas pi = (pix, piy).
      
    This can done by means of ccw since the product of ccw(pi,pj,pk)*ccw(pi,pj,pl)
    shows if pk and pl are on different or the same side(s) of the line1 (They must
    have different signums to be on different sides).
      
    Consequently if and ONLY IF pk as well as pl are on different sides of line1
    AND pi as well as pj are on different sides of line2 there HAS TO be an intersection.
  */
    
  return ( ( ccw(pix, piy, pjx, pjy, pkx, pky) *
	     ccw(pix, piy, pjx, pjy, plx, ply) <= 0 ) &&
	   ( ccw(pkx, pky, plx, ply, pix, piy) *
	     ccw(pkx, pky, plx, ply, pjx, pjy) <= 0 ) );
}

static
int check_ncorner(int ncorner, const double *lon_bounds, const double *lat_bounds)
{
  int ncorner_new = ncorner;
  int k;

  for ( k=ncorner-1; k>0; --k )
    if ( IS_NOT_EQUAL(lon_bounds[k], lon_bounds[k-1]) ||
	 IS_NOT_EQUAL(lat_bounds[k], lat_bounds[k-1]) ) break;

  if ( k < ncorner-1 ) ncorner_new = k+1;

  return ncorner_new;
}

static
void verify_grid(int gridtype, int gridsize, int ncorner,
		double *grid_center_lon, double *grid_center_lat,
		double *grid_corner_lon, double *grid_corner_lat)
{
  int i0, i, j, k, l;
  int l0;
  int nout;
  int isinside, convex, alone, isnegative;
  const int mnv = ncorner+1;
  int cuts[mnv][mnv];  
  int *alone_cell;          
  int check_corners;
  double lon, lat = 0;
  double lon_bounds[mnv], lat_bounds[mnv];
  double area, sumarea;

  alone_cell = (int*) malloc(gridsize*ncorner*sizeof(int));

  check_corners = 0; /* don't execute corner checking (last loop) */
  nout = 0;
  sumarea = 0;
  /*
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];
      for ( k = 0; k < ncorner; ++k )
        {
          lon_bounds[k] = grid_corner_lon[i*ncorner+k];
          lat_bounds[k] = grid_corner_lat[i*ncorner+k];
          if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
          if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
        }      
      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];
      fprintf(stdout, " %6i %6i %9.4f %9.4f :",  nout, i+1, lon, lat);
      for ( k = 0; k < ncorner; k++ )
	fprintf(stdout, " %9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
      fprintf(stdout, "\n");
    }
  */

  /* Check if center is inside bounds of cell */
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
        {
          lon_bounds[k] = grid_corner_lon[i*ncorner+k];
          lat_bounds[k] = grid_corner_lat[i*ncorner+k];
          if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
          if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
        }      
      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];
      
      isinside = pnpoly(ncorner+1, lon_bounds, lat_bounds, lon, lat);

      if ( !isinside ) nout++;
      if ( !isinside && cdoVerbose )
        {
          if ( nout == 1 )
            {
              fprintf(stdout,"\n CENTER IS OUT OF BOUNDS");
              fprintf(stdout,"\n                                               :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "          Corner %2i : ", k+1);
              fprintf(stdout,"\n Number  Index center_lon center_lat area*10^6 :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "   lon_%2.2i    lat_%2.2i : ", k+1, k+1);
              fprintf(stdout, "\n");
            }
          area = PolygonArea(ncorner+1, lon_bounds, lat_bounds,lat);
          fprintf(stdout, " %6i %6i  %9.4f  %9.4f %9.5f :", 
		  nout, i+1, lon, lat, area*pow(10,6));

	  int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

          for ( k = 0; k < ncorner_new; k++ )
	    fprintf(stdout, "%9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
           for ( k = ncorner_new; k < ncorner; k++ )
	     fprintf(stdout, "     ----      ---- : ");
          fprintf(stdout, "\n");
        }
    }

  if ( nout )
    cdoWarning("%d of %d points out of bounds!", nout, gridsize);
  
  /* check that all cell bounds have the same orientation */
  
  nout = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];
      
      for ( k = 0; k < ncorner; ++k )
	{
          lon_bounds[k] = grid_corner_lon[i*ncorner+k];
          lat_bounds[k] = grid_corner_lat[i*ncorner+k];
          if ( (grid_center_lon[i] - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
          if ( (lon_bounds[k] - grid_center_lon[i]) > 270 ) lon_bounds[k] -= 360;
	}
      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];
      
      area = PolygonArea(ncorner+1, lon_bounds, lat_bounds, lat);
      
      isnegative = area < 0 ? 1 : 0;
      sumarea += area < 0 ? -area : area;
      
      if ( isnegative ) nout++;
      
      if ( isnegative && cdoVerbose )
        {
          if ( nout == 1 )
            {
              fprintf(stdout,"\n                                     :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "          Corner %2i : ", k+1);
              fprintf(stdout,"\n Number  Index center_lon center_lat :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "   lon_%2.2i    lat_%2.2i : ", k+1, k+1);
              fprintf(stdout, "\n");
            }
          fprintf(stdout, " %6i %6i  %9.4f  %9.4f :", nout, i+1, lon, lat);

	  int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

          for ( k = 0; k < ncorner_new; k++ )
	    fprintf(stdout, "%9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
           for ( k = ncorner_new; k < ncorner; k++ )
	     fprintf(stdout, "     ----      ---- : ");

          fprintf(stdout, "\n");
        }
    }

  if ( nout )
    cdoWarning("%d of %d grid cells have wrong orientation!", nout, gridsize);

  if ( cdoVerbose ) 
    fprintf(stdout, "area-error: %9.5f%%\n", 100.*(sumarea - 4.*M_PI)/4.*M_PI );

  if ( fabs(100.*(sumarea - 4.*M_PI)/4.*M_PI) > 0.1)
    cdoWarning("area-error: %9.5f%%", 100.*(sumarea - 4.*M_PI)/4.*M_PI );
  
  /* check that all cells are convex */
  
  nout = 0;
  for ( i0 = 0; i0 < gridsize; i0++ )
    {
      lon = grid_center_lon[i0];
      lat = grid_center_lat[i0];

      for ( k = 0; k < ncorner; k++ )
	{
	  lon_bounds[k] = grid_corner_lon[i0*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i0*ncorner+k];
	  /* Find cells that cover left and right border of the grid and adjust
	     coordinates --> they become closed polygons on theta-phi plane! */
	  if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360; 
	  if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
	}
      
      /* Reset found cuts for the current cell before starting the search */
      for ( i = 0; i < ncorner; i++ )
	for ( j = 0; j < ncorner; j++ )
	  cuts[i][j] = 0;
      
      /* Loops cover all combinations between inner lines of the Polygon
	 Check whether each inner line is cut by an other (inner) one at least once. 
	 - Only if there is a cut every inner line the Polygon is convex
	 - We assume: Points are in either cyclic or anticyclic order
      */
      for ( i = 0; i < ncorner-1; i++ )
	{
          /* j = i+2 excludes lines from one corner to an other (j=i+1) and
	     from one point to itself (j=i)*/
          for ( j = i+2 ; j < ncorner; j++ )
	    {
              /* Exclude the line between the last and first corner */
              if ( i == 0 && j == ncorner-1 ) continue;

	      /* k = i+1: if starting point is in common lines to different corners
		 do not intersect */
              for ( k = i+1; k < ncorner - 1; k++ )
		{                  
                  if ( i == k ) l0 = j+1;
                  else          l0 = k+2;

                  for ( l = l0; l < ncorner; l++ )
		    {
                      if ( cuts[k][l] && cuts[i][j] ) continue;
		      /* Exlude the line between the last and first corner 
			 Exlude the line itself (l!=i, k!=j)
			 Check if line ij and kl intersect each other.
			 If so increment respective counters for intersections. 
			 It is not relevant by which line a line is intersected - 
			 it is only relevant if they is itersected! */
                      if ( ! ( k==0 && l == ncorner-1 ) && ( l != j ) && ( k != j )  )
			{
                          if ( intersect(lon_bounds[i], lat_bounds[i], lon_bounds[j], lat_bounds[j],
                                         lon_bounds[k], lat_bounds[k], lon_bounds[l], lat_bounds[l]) )
			    {
			      cuts[i][j]++; cuts[k][l]++; cuts[j][i]++; cuts[l][k]++;
			    }
			}
		    }
		}                  
	    }
	}

      convex = 1;
      /* The following loop covers all inner lines of the Polygon 
	 (The assumption applies that the points are in cyclic order) */
      for ( i = 0; i < ncorner-1; i++ )
	for ( j = i+2; j < ncorner; j++)
	  {
	    if ( i == 0 && j == ncorner-1 ) continue;	   
	    if ( ! cuts[i][j] ) convex = 0;
	  }
      if ( !convex ) nout++;        
      if ( cdoVerbose && ( !convex ) )
	{
          if ( nout == 1 )
	    {
              fprintf(stdout,"\n NO CONVEX POLYGON");
              fprintf(stdout,"\n                                       :");
              for ( k = 0; k < ncorner; k++ )
		fprintf(stdout, "            Corner %2i : ", k);
              fprintf(stdout,"\n Number  Index  center_lon  center_lat :");
              for ( k = 0; k < ncorner; k++ )
		fprintf(stdout, "    lon_%2.2i     lat_%2.2i : ", k, k);
              fprintf(stdout, "\n");
	    }
          
          fprintf(stdout, " %6i %6i   %9.4f   %9.4f :", nout, i0+1, lon, lat);
          for ( k = 0; k < ncorner; k++ )
	    fprintf(stdout, "  %9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
          fprintf(stdout, "\n");         
	}     
    }

  if ( nout )
    cdoWarning("%d of %d cells are not Convex!", nout, gridsize);

  if ( check_corners )
    {
      /* 
	 Check if there is a corner at the same point of 
	 an other cell foreach corner of each cell 
      */
      nout = 0;
      for ( i = 0; i < gridsize*ncorner; i++ )
	alone_cell[i] = 1;
      
      for ( i = 0; i < gridsize*ncorner; i++ )
	{
	  if ( ! alone_cell[i] ) continue;
	  alone = 1;
	  lon = grid_corner_lon[i];
	  lat = grid_corner_lat[i];			
	  for ( j = 0; j < gridsize*ncorner; j++ )
	    if ( j != i && 
		 IS_EQUAL(grid_corner_lat[j], lat) && 
		 IS_EQUAL(grid_corner_lon[j], lon) )
	      { alone = 0; alone_cell[i] = alone_cell[j] = 1; break; }
	  if ( alone )
	    {
	      if      ( lon >= 180. ) lon -= 360.;
	      else if ( lon  < 180. ) lon += 360.;
	      for ( j = i+1; j < gridsize*ncorner; j++ )
		if (j != i  && 
		    IS_EQUAL(grid_corner_lat[j], lat) && 
		    IS_EQUAL(grid_corner_lon[j], lon) )
		  { alone = 0; alone_cell[i] = alone_cell[j] = 0; break; }
	    }
	  if ( alone )
	    { 
	      nout++;
	      if ( cdoVerbose )
		{
		  if ( nout == 1 )
		    {
		      fprintf(stdout,"\n VERTEX ALONE ON GRID\n");
		      fprintf(stdout," number cell-Index  Vert-Index :        lon        lat\n");
		    }							
		  fprintf(stdout, " %6i     %6i      %6i : %10.4f %10.4f\n", 
			  nout, i/ncorner, i, grid_corner_lon[i], grid_corner_lat[i]);
		}					
	    }
	}

      if ( nout )
	cdoWarning("%d of %d corners are lonely on the grid!", nout, gridsize*ncorner);
    }

  free(alone_cell);
}


void verify_grid_old(int gridtype, int gridsize, int ncorner,
		double *grid_center_lon, double *grid_center_lat,
		double *grid_corner_lon, double *grid_corner_lat)
{
  int i, k;
  int nout;
  int isinside;
  int isnegative;
  double area;
  double lon, lat;
  double lon_bounds[9], lat_bounds[9];

  /* check that all centers are inside the bounds */

  nout = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
	{
	  lon_bounds[k] = grid_corner_lon[i*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i*ncorner+k];
	}

      for ( k = 0; k < ncorner; ++k )
	{
	  if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
	  if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
	}

      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];

      isinside = pnpoly(ncorner+1, lon_bounds, lat_bounds, lon, lat);

      if ( !isinside ) nout++;

      if ( !isinside && cdoVerbose )
	printf("center: %d %d %g %g %g %g %g %g %g %g %g %g\n", nout, i, lon, lat, lon_bounds[0], lat_bounds[0],
	       lon_bounds[1], lat_bounds[1], lon_bounds[2], lat_bounds[2], lon_bounds[3], lat_bounds[3]);
    }

  if ( nout > 0 )
    cdoWarning("%d of %d points out of bounds!", nout, gridsize);


  /* check that all cell bounds have the same orientation */

  nout = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
	{
	  lon_bounds[k] = grid_corner_lon[i*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i*ncorner+k];
	}

      for ( k = 0; k < ncorner; ++k )
	{
	  if ( (grid_center_lon[i] - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
	  if ( (lon_bounds[k] - grid_center_lon[i]) > 270 ) lon_bounds[k] -= 360;
	}

      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];

      area = PolygonArea_old(ncorner+1, lon_bounds, lat_bounds);

      if ( area < 0 ) isnegative = 1;
      else            isnegative = 0;

      if ( isnegative ) nout++;

      if ( isnegative && cdoVerbose )
	printf("bounds: %d %d %g %g %g %g %g %g %g %g %g %g\n", nout, i, lon, lat, lon_bounds[0], lat_bounds[0],
	       lon_bounds[1], lat_bounds[1], lon_bounds[2], lat_bounds[2], lon_bounds[3], lat_bounds[3]);
    }

  if ( nout > 0 )
    cdoWarning("%d of %d grid cells have wrong orientation!", nout, gridsize);
}


void make_cyclic(double *array1, double *array2, int nlon, int nlat)
{
  int i, j;
  int ij1, ij2;

  for ( j = 0; j < nlat; ++j )
    {
      for ( i = 0; i < nlon; ++i )
	{
	  ij1 = j*nlon+i;
	  ij2 = j*(nlon+1)+i;
	  array2[ij2] = array1[ij1];
	}
    }

  for ( j = 0; j < nlat; ++j )
    {
      ij2 = j*(nlon+1);
      array2[ij2+nlon] = array2[ij2];
    }
}


void *Outputgmt(void *argument)
{
  int GRIDVERIFY, OUTPUTCENTER, OUTPUTCENTER2, OUTPUTCENTERCPT, OUTPUTBOUNDS;
  int OUTPUTBOUNDSCPT, OUTPUTVECTOR, OUTPUTTRI, OUTPUTVRML;
  int operatorID;
  int process_data = TRUE;
  int i, j;
  int varID0, varID, recID;
  int nvals;
  int gridsize = 0;
  int gridsize2 = 0;
  int gridID, code;
  int nrecs;
  int levelID;
  int tsID;
  int streamID = 0;
  int vlistID;
  int nmiss;
  int nlon, nlat, nalloc;
  int nlev, lzon = FALSE, lmer = FALSE, lhov = FALSE;
  int ncorner = 0, ic;
  int status;
  int lgrid_gen_bounds = FALSE, luse_grid_corner = FALSE;
  int zaxisID, taxisID;
  int ninc = 1;
  int vdate, vtime;
  char varname[CDI_MAX_NAME];
  double level;
  double missval;
  double *array = NULL;
  double *array2 = NULL;
  double *parray;
  double *uf = NULL, *vf = NULL, *alpha = NULL, *auv = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  double *grid_center_lat2 = NULL, *grid_center_lon2 = NULL;
  double *grid_corner_lat = NULL, *grid_corner_lon = NULL;
  double *plat, *plon;
  double *zaxis_center_lev, *zaxis_lower_lev, *zaxis_upper_lev;
  int *grid_mask = NULL;
  FILE *cpt_fp;
  CPT cpt;
  int grid_is_circular;
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32];	  

  cdoInitialize(argument);

  GRIDVERIFY      = cdoOperatorAdd("gridverify",      0, 0, NULL);
  OUTPUTCENTER    = cdoOperatorAdd("outputcenter",    0, 0, NULL);
  OUTPUTCENTER2   = cdoOperatorAdd("outputcenter2",   0, 0, NULL);
  OUTPUTCENTERCPT = cdoOperatorAdd("outputcentercpt", 0, 0, NULL);
  OUTPUTBOUNDS    = cdoOperatorAdd("outputbounds",    0, 0, NULL);
  OUTPUTBOUNDSCPT = cdoOperatorAdd("outputboundscpt", 0, 0, NULL);
  OUTPUTVECTOR    = cdoOperatorAdd("outputvector",    0, 0, NULL);
  OUTPUTTRI       = cdoOperatorAdd("outputtri",       0, 0, NULL);
  OUTPUTVRML      = cdoOperatorAdd("outputvrml",      0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTVECTOR )
    {
      operatorInputArg("increment");
      operatorCheckArgc(1);
      ninc = parameter2int(operatorArgv()[0]);
      if ( ninc < 1 ) cdoAbort("Increment must be greater than 0!");
    }

  if ( operatorID == GRIDVERIFY  )
    {
      process_data = FALSE;
      luse_grid_corner = TRUE;
    }

  if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
    luse_grid_corner = TRUE;

  if ( operatorID == OUTPUTCENTERCPT || operatorID == OUTPUTBOUNDSCPT || operatorID == OUTPUTVRML )
    {
      char *cpt_file;

      cpt_file = operatorArgv()[0];

      if ( (cpt_fp = fopen (cpt_file, "r")) == NULL )
	cdoAbort("Open failed on color palette table %s", cpt_file);

      status = cptRead(cpt_fp, &cpt);
      if ( status != 0 )
	cdoAbort("Error during read of color palette table %s", cpt_file);
      
      if ( cdoVerbose ) cptWrite(stderr, cpt);
    }

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  varID = 0;
  vlistInqVarName(vlistID, varID, varname);
  code    = vlistInqVarCode(vlistID, varID);
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  missval = vlistInqVarMissval(vlistID, varID);

  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 1);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    {
      gridID = gridToCurvilinear(gridID, 1);
      lgrid_gen_bounds = TRUE;
    }

  gridsize = gridInqSize(gridID);
  nlon     = gridInqXsize(gridID);
  nlat     = gridInqYsize(gridID);
  nlev     = zaxisInqSize(zaxisID);

  if ( gridInqMaskGME(gridID, NULL) )
    {
      grid_mask = (int*) malloc(gridsize*sizeof(int));
      gridInqMaskGME(gridID, grid_mask);
    }

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED )
    {
      if ( nlon == 1 && nlat  > 1 && nlev == 1 ) lhov = TRUE;
      if ( nlon == 1 && nlat  > 1 && nlev  > 1 ) lzon = TRUE;
      if ( nlon  > 1 && nlat == 1 && nlev  > 1 ) lmer = TRUE;
    }
  else
    {
      nlat = 1;
    }

  if ( cdoVerbose && lhov ) cdoPrint("Process hovmoeller data");
  if ( cdoVerbose && lzon ) cdoPrint("Process zonal data");
  if ( cdoVerbose && lmer ) cdoPrint("Process meridional data");
  /*
  if ( lzon || lmer ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available for zonal/meridional data!");
    }
  */
  if ( lhov ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available hovmoeller data!");
    }

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    ncorner = gridInqNvertex(gridID);
  else
    ncorner = 4;

  grid_is_circular = gridIsCircular(gridID);

  grid_center_lat = (double*) malloc(gridsize*sizeof(double));
  grid_center_lon = (double*) malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");

  nvals = gridsize;
  plon = grid_center_lon;
  plat = grid_center_lat;

  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
    {
      int ij2;

      gridsize2 = nlat*(nlon+1);

      grid_center_lat2 = (double*) malloc(gridsize2*sizeof(double));
      grid_center_lon2 = (double*) malloc(gridsize2*sizeof(double));

      make_cyclic(grid_center_lat, grid_center_lat2, nlon, nlat);
      make_cyclic(grid_center_lon, grid_center_lon2, nlon, nlat);

      for ( j = 0; j < nlat; ++j )
	{
	  ij2 = j*(nlon+1);
	  grid_center_lon2[ij2+nlon] += 360;
	}

      nvals = gridsize2;
      plon = grid_center_lon2;
      plat = grid_center_lat2;
    }

  zaxis_center_lev = (double*) malloc(nlev*sizeof(double));
  zaxis_lower_lev  = (double*) malloc(nlev*sizeof(double));
  zaxis_upper_lev  = (double*) malloc(nlev*sizeof(double));

  zaxisInqLevels(zaxisID, zaxis_center_lev);

  if ( luse_grid_corner )
    {
      if ( ncorner == 0 ) cdoAbort("grid corner missing!");
      nalloc = ncorner*gridsize;
      grid_corner_lat = (double*) realloc(grid_corner_lat, nalloc*sizeof(double));
      grid_corner_lon = (double*) realloc(grid_corner_lon, nalloc*sizeof(double));

      if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
	{
	  gridInqYbounds(gridID, grid_corner_lat);
	  gridInqXbounds(gridID, grid_corner_lon);
	}
      else
	{
	  if ( lgrid_gen_bounds )
	    {
	      char xunitstr[CDI_MAX_NAME];
	      char yunitstr[CDI_MAX_NAME];
	      gridInqXunits(gridID, xunitstr);
	      gridInqYunits(gridID, yunitstr);
	      if ( ! lzon ) grid_cell_center_to_bounds_X2D(xunitstr, nlon, nlat, grid_center_lon, grid_corner_lon, 0);
	      if ( ! lmer ) grid_cell_center_to_bounds_Y2D(yunitstr, nlon, nlat, grid_center_lat, grid_corner_lat);
	    }
	  else
	    cdoAbort("Grid corner missing!");
	}


      /* Note: using units from latitude instead from bounds */
      grid_to_degree(units, ncorner*gridsize, grid_corner_lon, "grid corner lon");
      grid_to_degree(units, ncorner*gridsize, grid_corner_lat, "grid corner lat");

      if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	{
	  zaxisInqLbounds(zaxisID, zaxis_lower_lev);
	  zaxisInqUbounds(zaxisID, zaxis_upper_lev);
	}
      else
	{
	  zaxis_lower_lev[0] = zaxis_center_lev[0];
	  for ( i = 1; i < nlev; ++i )
	    zaxis_lower_lev[i] = 0.5*(zaxis_center_lev[i] + zaxis_center_lev[i-1]);

	  zaxis_upper_lev[nlev-1] = zaxis_center_lev[nlev-1];
	  for ( i = 0; i < nlev-1; ++i )
	    zaxis_upper_lev[i] = zaxis_lower_lev[i+1];

	  if ( cdoVerbose )
	    for ( i = 0; i < nlev; ++i )
	      printf("level: %d %g %g %g\n",
		     i+1, zaxis_lower_lev[i], zaxis_center_lev[i], zaxis_upper_lev[i]);
	}
    }

  array = (double*) malloc(gridsize*sizeof(double));
  parray = array;
						
  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
    {
      array2 = (double*) malloc(nlat*(nlon+1)*sizeof(double));
      parray = array2;
    }

  if ( operatorID == OUTPUTVECTOR )
    {
      uf    = (double*) malloc(gridsize*sizeof(double));
      vf    = (double*) malloc(gridsize*sizeof(double));
      alpha = (double*) malloc(gridsize*sizeof(double));
      auv   = (double*) malloc(gridsize*sizeof(double));
    }

  if ( operatorID == GRIDVERIFY )
    verify_grid(gridInqType(gridID), gridsize, ncorner,
		grid_center_lon, grid_center_lat,
		grid_corner_lon, grid_corner_lat);

  tsID = 0;
  if ( process_data )
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));

      if ( tsID == 0 && operatorID != OUTPUTTRI )
	{
	  if ( operatorID == OUTPUTVRML )
	    printf("#VRML V2.0 utf8\n\n");
#if defined(VERSION)
	  fprintf(stdout, "# Generated by CDO version %s\n", VERSION);
	  fprintf(stdout, "#\n");
#endif
	  fprintf(stdout, "# Operator = %s\n", cdoOperatorName(operatorID));
	  if      ( lhov )  fprintf(stdout, "# Mode     = hovmoeller\n");
	  else if ( lzon )  fprintf(stdout, "# Mode     = zonal\n");
	  else if ( lmer )  fprintf(stdout, "# Mode     = meridional\n");
	  else              fprintf(stdout, "# Mode     = horizonal\n");

	  if ( operatorID == OUTPUTVECTOR )
	    fprintf(stdout, "# Increment = %d\n", ninc);
	  fprintf(stdout, "#\n");
	  fprintf(stdout, "# File  = %s\n", cdoStreamName(0)->args);
	  fprintf(stdout, "# Date  = %s\n", vdatestr);
	  fprintf(stdout, "# Time  = %s\n", vtimestr);
	  fprintf(stdout, "# Name  = %s\n", varname);
	  fprintf(stdout, "# Code  = %d\n", code);
	}

      varID0 = varID;

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID, &varID, &levelID);

	  if ( varID != varID0 ) continue;
	  if ( recID > 0 && !lzon && !lmer ) continue;

	  streamReadRecord(streamID, array, &nmiss);

	  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
	    make_cyclic(array, array2, nlon, nlat);

	  level = zaxis_center_lev[levelID];

	  if ( (tsID == 0 || lzon || lmer) && operatorID != OUTPUTTRI )
	    fprintf(stdout, "# Level = %g\n", level);
	  if ( lhov )
	    fprintf(stdout, "# Timestep = %d\n", tsID+1);

	  if ( operatorID != OUTPUTTRI ) fprintf(stdout, "#\n");

	  if ( operatorID == OUTPUTCENTER || operatorID == OUTPUTCENTER2 || operatorID == OUTPUTCENTERCPT )
	    {
	      for ( i = 0; i < nvals; i++ )
		{
		  if ( grid_mask )
		    if ( grid_mask[i] == 0 ) continue;

		  if ( operatorID == OUTPUTCENTERCPT )
		    {
		      int r = 0, g = 0, b = 0, n;

		      if ( !DBL_IS_EQUAL(array[i], missval) )
			{
			  for ( n = 0; n < cpt.ncolors; n++ )
			    if ( array[i] > cpt.lut[n].z_low && array[i] <= cpt.lut[n].z_high ) break;

			  if ( n == cpt.ncolors )
			    {
			      r = cpt.bfn[0].rgb[0];  g = cpt.bfn[0].rgb[1];  b = cpt.bfn[0].rgb[2];
			    }
			  else
			    {
			      r = cpt.lut[n].rgb_high[0];  g = cpt.lut[n].rgb_high[1];  b = cpt.lut[n].rgb_high[2];
			    }
			}
		      else
			{
			  r = cpt.bfn[2].rgb[0];  g = cpt.bfn[2].rgb[1];  b = cpt.bfn[2].rgb[2]; 
			}
		    }

		  if ( operatorID == OUTPUTCENTER )
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else if ( lhov )
			fprintf(stdout, " %d  %g  %g\n", tsID+1, grid_center_lat[i], array[i]);
		      else
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i]);
		    }
		  else if ( operatorID == OUTPUTCENTER2 )
		    {
		      fprintf(stdout, " %g  %g  %g\n", plon[i], plat[i], parray[i]);
		    }
		  else
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else
			fprintf(stdout, " %g  %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i], array[i]);
		    }
		}
	      fprintf(stdout, "#\n");
	    }
	  else if ( operatorID == OUTPUTTRI )
	    {
	      int ij, c1, c2, c3;
	      int mlon, ip1;
	      if ( gridInqType(gridID) != GRID_CURVILINEAR ) cdoAbort("Unsupported grid!");

	      mlon = nlon-1;
	      /* if ( gridIsCircular(gridID) ) mlon = nlon; */
	      for ( j = 0; j < nlat-1; ++j )
		{
		  for ( i = 0; i < mlon; ++i )
		    {
		      ip1 = i+1;
		      if ( i == nlon-1 ) ip1 = 0;
		      ij = j*nlon+i;
		      c1 = (j)*nlon+ip1;
		      c2 = (j)*nlon+i;
		      c3 = (j+1)*nlon+i;
		      fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
		      c1 = (j)*nlon+i+1;
		      c2 = (j+1)*nlon+i;
		      c3 = (j+1)*nlon+ip1;
		      fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
		    }
		}
	    }
	  else if ( operatorID == OUTPUTVECTOR )
	    {
	      if ( nrecs < 2 ) cdoAbort("Too few fields!");

	      memcpy(uf, array, gridsize*sizeof(double));
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, vf, &nmiss);

	      for ( j = 0; j < nlat; j += ninc )
		for ( i = 0; i < nlon; i += ninc )
		  {
		    /* compute length of velocity vector */
		    auv[IX2D(j,i,nlon)] = sqrt(uf[IX2D(j,i,nlon)]*uf[IX2D(j,i,nlon)] + 
					       vf[IX2D(j,i,nlon)]*vf[IX2D(j,i,nlon)]);

		    alpha[IX2D(j,i,nlon)] = atan2(vf[IX2D(j,i,nlon)],uf[IX2D(j,i,nlon)]);
		    alpha[IX2D(j,i,nlon)] = 90. - alpha[IX2D(j,i,nlon)]*RAD2DEG;

		    if ( alpha[IX2D(j,i,nlon)] <   0 ) alpha[IX2D(j,i,nlon)] += 360;
		    if ( alpha[IX2D(j,i,nlon)] > 360 ) alpha[IX2D(j,i,nlon)] -= 360;

		    if ( fabs(auv[IX2D(j,i,nlon)]) > 0 )
		      fprintf(stdout, " %g  %g  %g  %g\n",
			      grid_center_lon[IX2D(j,i,nlon)], grid_center_lat[IX2D(j,i,nlon)],
			      alpha[IX2D(j,i,nlon)], auv[IX2D(j,i,nlon)]);
		  }

	      fprintf(stdout, "#\n");
	      break;
	    }
	  else if ( operatorID == OUTPUTVRML )
	    {
	      double minval = 1e33;
	      double maxval = -1e33;
	      double meanval = 0;
	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( array[i] < minval ) minval = array[i];
		  if ( array[i] > maxval ) maxval = array[i];
		  meanval += array[i];
		}
	      meanval /= gridsize;

	      double dx = 10./nlon;

	      printf("Viewpoint {\n");
	      printf("  description \"viewpoint1\"\n");
	      printf("  orientation 0 0 1 0\n");
	      printf("  position 0.0 0.0 10.0\n");
	      printf("}\n");
	      printf("\n");
	      printf("Background {\n");
	      printf("  skyColor [\n");
	      printf("    0.0 0.1 0.8,\n");
	      printf("    0.0 0.5 1.0,\n");
	      printf("    1.0 1.0 1.0\n");
	      printf("  ]\n");
	      printf("  skyAngle [0.785, 1.571]\n");
	      printf("\n");
	      printf("  groundColor [\n");
	      printf("    0.0 0.0 0.0,\n");
	      printf("    0.3 0.3 0.3,\n");
	      printf("    0.5 0.5 0.5\n");
	      printf("  ]\n");
	      printf("  groundAngle [0.785, 1.571]\n");
	      printf("}\n");
	      printf("\n");
	      printf("Transform {\n");
	      printf("  children [\n");
	      printf("    Shape {\n");
	      printf("      appearance Appearance {\n");
	      printf("        material Material {}\n");
	      printf("      }\n");
	      printf("      geometry ElevationGrid {\n");
    	      printf("        colorPerVertex TRUE\n");
    	      printf("        solid FALSE\n");
    	      printf("        xDimension %d\n", nlon);
    	      printf("        zDimension %d\n", nlat);
    	      printf("        xSpacing %g\n", dx);
    	      printf("        zSpacing %g\n", dx);
	      printf("        color Color {\n");
	      printf("          color [\n");
	      for ( j = nlat-1; j >= 0 ; --j )
		for ( i = 0; i < nlon; ++i )
		{
		  int r = 0, g = 0, b = 0, n;
		  double val = array[j*nlon+i];

		  if ( !DBL_IS_EQUAL(val, missval) )
		    {
		      for ( n = 0; n < cpt.ncolors; n++ )
			if ( val > cpt.lut[n].z_low && val <= cpt.lut[n].z_high ) break;
		      
		      if ( n == cpt.ncolors )
			{
			  r = cpt.bfn[0].rgb[0];  g = cpt.bfn[0].rgb[1];  b = cpt.bfn[0].rgb[2];
			}
		      else
			{
			  //  r = cpt.lut[n].rgb_high[0];  g = cpt.lut[n].rgb_high[1];  b = cpt.lut[n].rgb_high[2];
			  r = intlin(val, cpt.lut[n].rgb_low[0], cpt.lut[n].z_low, cpt.lut[n].rgb_high[0], cpt.lut[n].z_high);
			  g = intlin(val, cpt.lut[n].rgb_low[1], cpt.lut[n].z_low, cpt.lut[n].rgb_high[1], cpt.lut[n].z_high);
			  b = intlin(val, cpt.lut[n].rgb_low[2], cpt.lut[n].z_low, cpt.lut[n].rgb_high[2], cpt.lut[n].z_high);
			}
		    }
		  else
		    {
		      r = cpt.bfn[2].rgb[0];  g = cpt.bfn[2].rgb[1];  b = cpt.bfn[2].rgb[2]; 
		    }
		  printf(" %.3g %.3g %.3g,\n", r/255., g/255., b/255.);
		}
	      printf("          ]\n");
	      printf("        }\n");
    	      printf("        height [\n");
	      //	      for ( j = 0; j < nlat; ++j )
	      for ( j = nlat-1; j >= 0 ; --j )
		{
		  for ( i = 0; i < nlon; ++i )
		    {
		      printf("%g,\n", array[j*nlon+i]);
		    }
		}
    	      printf("        ]\n");
	      printf("      }\n");
	      printf("    }\n");
	      printf("  ]\n");
	      printf("  translation -5 0 %g\n", -5.*nlat/nlon);
	      printf("  rotation 0.0 0.0 0.0 0.0\n");
	      printf("  scale 1.0 %g 1.0\n", 0.5/(maxval-minval));
	      printf("}\n");
            }
	  else if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( grid_mask )
		    if ( grid_mask[i] == 0 ) continue;

		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    fprintf(stdout, "> -Z%g", array[i]);
		  else
		    fprintf(stdout, "> -ZNaN");

		  if ( operatorID == OUTPUTBOUNDSCPT )
		    {
		      int r = 0, g = 0, b = 0, n;

		      if ( !DBL_IS_EQUAL(array[i], missval) )
			{
			  for ( n = 0; n < cpt.ncolors; n++ )
			    if ( array[i] > cpt.lut[n].z_low && array[i] <= cpt.lut[n].z_high ) break;

			  if ( n == cpt.ncolors )
			    {
			      r = cpt.bfn[0].rgb[0];  g = cpt.bfn[0].rgb[1];  b = cpt.bfn[0].rgb[2];
			    }
			  else
			    {
			      r = cpt.lut[n].rgb_high[0];  g = cpt.lut[n].rgb_high[1];  b = cpt.lut[n].rgb_high[2];
			    }
			}
		      else
			{
			  r = cpt.bfn[2].rgb[0];  g = cpt.bfn[2].rgb[1];  b = cpt.bfn[2].rgb[2]; 
			}

		      fprintf(stdout, " -G%d/%d/%d", r, g, b);
		    }

		  fprintf(stdout, "\n");

		  if ( lzon )
		    {
		      double xlev[4], xlat[4];
		      double levmin = zaxis_lower_lev[levelID];
		      double levmax = zaxis_upper_lev[levelID];
		      double latmin = grid_corner_lat[i*4+0];
		      double latmax = grid_corner_lat[i*4+0];
		      for ( ic = 1; ic < 4; ic++ )
			{
			  if ( grid_corner_lat[i*4+ic] < latmin ) latmin = grid_corner_lat[i*4+ic];
			  if ( grid_corner_lat[i*4+ic] > latmax ) latmax = grid_corner_lat[i*4+ic];
			}
		      xlev[0] = levmin;
		      xlev[1] = levmax;
		      xlev[2] = levmax;
		      xlev[3] = levmin;
		      xlat[0] = latmin;
		      xlat[1] = latmin;
		      xlat[2] = latmax;
		      xlat[3] = latmax;
		      for ( ic = 0; ic < 4; ic++ )
			fprintf(stdout, "   %g  %g\n", xlat[ic], xlev[ic]);
		      fprintf(stdout, "   %g  %g\n", xlat[0], xlev[0]);
		    }
		  else if ( lmer )
		    {
		      double xlev[4], xlon[4];
		      double levmin = zaxis_lower_lev[levelID];
		      double levmax = zaxis_upper_lev[levelID];
		      double lonmin = grid_corner_lon[i*4+0];
		      double lonmax = grid_corner_lon[i*4+0];
		      for ( ic = 1; ic < 4; ic++ )
			{
			  if ( grid_corner_lon[i*4+ic] < lonmin ) lonmin = grid_corner_lon[i*4+ic];
			  if ( grid_corner_lon[i*4+ic] > lonmax ) lonmax = grid_corner_lon[i*4+ic];
			}
		      xlev[0] = levmin;
		      xlev[1] = levmin;
		      xlev[2] = levmax;
		      xlev[3] = levmax;
		      xlon[0] = lonmin;
		      xlon[1] = lonmax;
		      xlon[2] = lonmax;
		      xlon[3] = lonmin;
		      for ( ic = 0; ic < 4; ic++ )
			fprintf(stdout, "   %g  %g\n", xlon[ic], xlev[ic]);
		      fprintf(stdout, "   %g  %g\n", xlon[0], xlev[0]);
		    }
		  else if ( lhov )
		    {
		      cdoAbort("Implementation for hovmoeller data missing!\n");
		    }
		  else
		    {
		      const double *lon_bounds = grid_corner_lon+i*ncorner;
		      const double *lat_bounds = grid_corner_lat+i*ncorner;
		      int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

		      for ( ic = 0; ic < ncorner_new; ic++ )
			fprintf(stdout, "   %g  %g\n", lon_bounds[ic], lat_bounds[ic]);
		      fprintf(stdout, "   %g  %g\n", lon_bounds[0], lat_bounds[0]);
		    }
		}
	      fprintf(stdout, "\n");
	    }
	}

      if ( ! lhov ) break;

      tsID++;
    }

  streamClose(streamID);

  if ( array  ) free(array);
  if ( array2 ) free(array2);
  if ( grid_mask ) free(grid_mask);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);
  if ( grid_center_lon2 ) free(grid_center_lon2);
  if ( grid_center_lat2 ) free(grid_center_lat2);
  if ( grid_corner_lon ) free(grid_corner_lon);
  if ( grid_corner_lat ) free(grid_corner_lat);

  free(zaxis_center_lev);
  free(zaxis_lower_lev);
  free(zaxis_upper_lev);

  cdoFinish();

  return (0);
}
