#include "cdo.h"
#include "remap.h"


int rect_grid_search(long *ii, long *jj, double x, double y, long nxm, long nym, const double *restrict xm, const double *restrict ym);


int grid_search_reg2d_nn(long nx, long ny, int *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
			 const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  long n, srch_add;
  long i;
  long ii, jj;
  long jjskip;
  double coslat, sinlat;
  double dist_min, distance; /* For computing dist-weighted avg */
  double *sincoslon;
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  dist_min = BIGNUM;
  for ( n = 0; n < 4; ++n ) nbr_dist[n] = BIGNUM;  

  long jjf = 0, jjl = ny-1;
  if ( plon >= src_center_lon[0] && plon <= src_center_lon[nx-1] )
    {
      if ( src_center_lat[0] < src_center_lat[ny-1] )
	{
	  if ( plat <= src_center_lat[0] )
	    { jjf = 0; jjl = 1; }
	  else
	    { jjf = ny-2; jjl = ny-1; }
	}
      else
	{
	  if ( plat >= src_center_lat[0] )
	    { jjf = 0; jjl = 1; }
	  else
	    { jjf = ny-2; jjl = ny-1; }
	}
    }

  sincoslon = (double*) malloc(nx*sizeof(double));

  for ( ii = 0; ii < nx; ++ii )
    sincoslon[ii] = coslon_dst*cos(src_center_lon[ii]) + sinlon_dst*sin(src_center_lon[ii]);

  for ( jj = jjf; jj <= jjl; ++jj )
    {
      coslat = coslat_dst*cos(src_center_lat[jj]);
      sinlat = sinlat_dst*sin(src_center_lat[jj]);

      jjskip = jj > 1 && jj < (ny-2);

      for ( ii = 0; ii < nx; ++ii )
	{
	  if ( jjskip && ii > 1 && ii < (nx-2) ) continue;

	  srch_add = jj*nx + ii;

	  distance = acos(coslat*sincoslon[ii] + sinlat);

	  if ( distance < dist_min )
	    {
	      for ( n = 0; n < 4; ++n )
		{
		  if ( distance < nbr_dist[n] )
		    {
		      for ( i = 3; i > n; --i )
			{
			  nbr_add [i] = nbr_add [i-1];
			  nbr_dist[i] = nbr_dist[i-1];
			}
		      search_result = -1;
		      nbr_add [n] = srch_add;
		      nbr_dist[n] = distance;
		      dist_min = nbr_dist[3];
		      break;
		    }
		}
	    }
	}
    }

  free(sincoslon);

  for ( n = 0; n < 4; ++n ) nbr_dist[n] = ONE/(nbr_dist[n] + TINY);
  distance = 0.0;
  for ( n = 0; n < 4; ++n ) distance += nbr_dist[n];
  for ( n = 0; n < 4; ++n ) nbr_dist[n] /= distance;

  return (search_result);
}


int grid_search_reg2d(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		      double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		      const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  /*
    Output variables:

    int    src_add[4]              ! address of each corner point enclosing P
    double src_lats[4]             ! latitudes  of the four corner points
    double src_lons[4]             ! longitudes of the four corner points

    Input variables:

    double plat                    ! latitude  of the search point
    double plon                    ! longitude of the search point

    int src_grid_dims[2]           ! size of each src grid dimension

    double src_center_lat[]        ! latitude  of each src grid center 
    double src_center_lon[]        ! longitude of each src grid center
  */
  /*  Local variables */
  int search_result = 0;
  long n;
  long ii, iix, jj;

  for ( n = 0; n < 4; ++n ) src_add[n] = 0;

  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];

  long nxm = nx;
  if ( src_grid->is_cyclic ) nxm++;

  if ( /*plon < 0   &&*/ plon < src_center_lon[0]     ) plon += PI2;
  if ( /*plon > PI2 &&*/ plon > src_center_lon[nxm-1] ) plon -= PI2;

  int lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);

  if ( lfound )
    {
      iix = ii;
      if ( src_grid->is_cyclic && iix == (nxm-1) ) iix = 0;
      src_add[0] = (jj-1)*nx+(ii-1);
      src_add[1] = (jj-1)*nx+(iix);
      src_add[2] = (jj)*nx+(iix);
      src_add[3] = (jj)*nx+(ii-1);

      src_lons[0] = src_center_lon[ii-1];
      src_lons[1] = src_center_lon[iix];
      /* For consistency, we must make sure all lons are in same 2pi interval */
      if ( src_lons[0] > PI2 ) src_lons[0] -= PI2;
      if ( src_lons[0] < 0   ) src_lons[0] += PI2;
      if ( src_lons[1] > PI2 ) src_lons[1] -= PI2;
      if ( src_lons[1] < 0   ) src_lons[1] += PI2;
      src_lons[2] = src_lons[1];
      src_lons[3] = src_lons[0];

      src_lats[0] = src_center_lat[jj-1];
      src_lats[1] = src_lats[0];
      src_lats[2] = src_center_lat[jj];
      src_lats[3] = src_lats[2];

      search_result = 1;
      
      return (search_result);
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside 
    the grid. Fall back to a distance-weighted average of the four closest points.
    Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
    parent routine from computing bilinear weights.
  */
  if ( !src_grid->lextrapolate ) return (search_result);

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_reg2d_nn(nx, ny, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return (search_result);
}  /* grid_search_reg2d */
