#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BILINEAR INTERPOLATION                                             */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */


int find_ij_weights(double plon, double plat, double* restrict src_lats, double* restrict src_lons, double *ig, double *jg)
{
  int lfound = 0;
  long iter;                     /*  iteration counters   */
  double iguess, jguess;         /*  current guess for bilinear coordinate  */
  double deli, delj;             /*  corrections to iw,jw                   */
  double dth1, dth2, dth3;       /*  some latitude  differences             */
  double dph1, dph2, dph3;       /*  some longitude differences             */
  double dthp, dphp;             /*  difference between point and sw corner */
  double mat1, mat2, mat3, mat4; /*  matrix elements                        */
  double determinant;            /*  matrix determinant                     */
  double converge = 1.e-10;      /* Convergence criterion                   */
  extern long remap_max_iter;

  /* Iterate to find iw,jw for bilinear approximation  */

  dth1 = src_lats[1] - src_lats[0];
  dth2 = src_lats[3] - src_lats[0];
  dth3 = src_lats[2] - src_lats[1] - dth2;

  dph1 = src_lons[1] - src_lons[0];
  dph2 = src_lons[3] - src_lons[0];
  dph3 = src_lons[2] - src_lons[1];

  if ( dph1 >  THREE*PIH ) dph1 -= PI2;
  if ( dph2 >  THREE*PIH ) dph2 -= PI2;
  if ( dph3 >  THREE*PIH ) dph3 -= PI2;
  if ( dph1 < -THREE*PIH ) dph1 += PI2;
  if ( dph2 < -THREE*PIH ) dph2 += PI2;
  if ( dph3 < -THREE*PIH ) dph3 += PI2;

  dph3 = dph3 - dph2;

  iguess = HALF;
  jguess = HALF;

  for ( iter = 0; iter < remap_max_iter; ++iter )
    {
      dthp = plat - src_lats[0] - dth1*iguess - dth2*jguess - dth3*iguess*jguess;
      dphp = plon - src_lons[0];
      
      if ( dphp >  THREE*PIH ) dphp -= PI2;
      if ( dphp < -THREE*PIH ) dphp += PI2;

      dphp = dphp - dph1*iguess - dph2*jguess - dph3*iguess*jguess;

      mat1 = dth1 + dth3*jguess;
      mat2 = dth2 + dth3*iguess;
      mat3 = dph1 + dph3*jguess;
      mat4 = dph2 + dph3*iguess;

      determinant = mat1*mat4 - mat2*mat3;

      deli = (dthp*mat4 - dphp*mat2)/determinant;
      delj = (dphp*mat1 - dthp*mat3)/determinant;

      if ( fabs(deli) < converge && fabs(delj) < converge ) break;

      iguess += deli;
      jguess += delj;
    }

  *ig = iguess;
  *jg = jguess;

  if ( iter < remap_max_iter ) lfound = 1;

  return (lfound);
}

static
void set_bilinear_weights(double iw, double jw, double wgts[4])
{
  wgts[0] = (1.-iw) * (1.-jw);
  wgts[1] =     iw  * (1.-jw);
  wgts[2] =     iw  *     jw;
  wgts[3] = (1.-iw) *     jw;
}


int num_src_points(const int* restrict mask, const int src_add[4], double src_lats[4])
{
  int icount = 0;

  for ( int n = 0; n < 4; ++n )
    {
      if ( mask[src_add[n]] )
	icount++;
      else
	src_lats[n] = 0.;
    }

  return (icount);
}

static
void renormalize_weights(const double src_lats[4], double wgts[4])
{
  int n;
  double sum_wgts = 0.0; /* sum of weights for normalization */
  /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
  for ( n = 0; n < 4; ++n ) sum_wgts += fabs(src_lats[n]);
  for ( n = 0; n < 4; ++n ) wgts[n] = fabs(src_lats[n])/sum_wgts;
}

static
void bilinear_warning(double plon, double plat, double iw, double jw, int* src_add, double* src_lons, double* src_lats, remapgrid_t* src_grid)
{
  static int lwarn = TRUE;

  if ( cdoVerbose )
    {
      cdoPrint("Point coords: %g %g", plat, plon);
      cdoPrint("Src grid lats: %g %g %g %g", src_lats[0], src_lats[1], src_lats[2], src_lats[3]);
      cdoPrint("Src grid lons: %g %g %g %g", src_lons[0], src_lons[1], src_lons[2], src_lons[3]);
      cdoPrint("Src grid addresses: %d %d %d %d", src_add[0], src_add[1], src_add[2], src_add[3]);
      cdoPrint("Src grid lats: %g %g %g %g",
	       src_grid->cell_center_lat[src_add[0]], src_grid->cell_center_lat[src_add[1]],
	       src_grid->cell_center_lat[src_add[2]], src_grid->cell_center_lat[src_add[3]]);
      cdoPrint("Src grid lons: %g %g %g %g",
	       src_grid->cell_center_lon[src_add[0]], src_grid->cell_center_lon[src_add[1]],
	       src_grid->cell_center_lon[src_add[2]], src_grid->cell_center_lon[src_add[3]]);
      cdoPrint("Current iw,jw : %g %g", iw, jw);
    }

  if ( cdoVerbose || lwarn )
    {
      lwarn = FALSE;
      //  cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
      cdoWarning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

static
void bilinear_remap(double* restrict tgt_point, const double* restrict src_array, const double wgts[4], const int src_add[4])
{
  // *tgt_point = 0.;
  // for ( int n = 0; n < 4; ++n ) *tgt_point += src_array[src_add[n]]*wgts[n];
  *tgt_point = src_array[src_add[0]]*wgts[0] + src_array[src_add[1]]*wgts[1]
             + src_array[src_add[2]]*wgts[2] + src_array[src_add[3]]*wgts[3];
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_weights_bilinear(remapgrid_t* src_grid, remapgrid_t* tgt_grid, remapvars_t* rv)
{
  /*   Local variables */
  int  search_result;
  long tgt_cell_add;             /*  destination addresss                   */
  int src_add[4];                /*  address for the four source points     */
  double src_lats[4];            /*  latitudes  of four bilinear corners    */
  double src_lons[4];            /*  longitudes of four bilinear corners    */
  double wgts[4];                /*  bilinear weights for four corners      */
  double plat, plon;             /*  lat/lon coords of destination point    */
  extern int timer_remap_bil;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  if ( cdoTimer ) timer_start(timer_remap_bil);

  progressInit();

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bilinear interpolation when source grid rank != 2"); 

  long tgt_grid_size = tgt_grid->size;

  weightlinks_t *weightlinks = (weightlinks_t *) malloc(tgt_grid_size*sizeof(weightlinks_t));

  double findex = 0;

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoVerbose, weightlinks, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, findex) \
  private(tgt_cell_add, src_add, src_lats, src_lons, wgts, plat, plon, search_result)    \
  schedule(static)
#endif
  for ( tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      int lprogress = 1;
      if ( cdo_omp_get_thread_num() != 0 ) lprogress = 0;

#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      weightlinks[tgt_cell_add].nlinks = 0;	

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      plat = tgt_grid->cell_center_lat[tgt_cell_add];
      plon = tgt_grid->cell_center_lon[tgt_cell_add];

      /* Find nearest square of grid points on source grid  */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);

      /* Check to see if points are land points */
      if ( search_result > 0 )
	{
	  for ( int n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[tgt_cell_add] = 1.;

          if ( find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw) )
	    {
	      /* Successfully found iw,jw - compute weights */
	      set_bilinear_weights(iw, jw, wgts);

	      store_weightlinks(4, src_add, wgts, tgt_cell_add, weightlinks);
	    }
          else
	    {
	      bilinear_warning(plon, plat, iw, jw, src_add, src_lons, src_lats, src_grid);

	      search_result = -1;
	    }
	}

      /*
	Search for bilinear failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
          if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[tgt_cell_add] = 1.;

	      store_weightlinks(4, src_add, wgts, tgt_cell_add, weightlinks);
	    }
        }
    }

  weightlinks2remaplinks(tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) free(weightlinks);

  if ( cdoTimer ) timer_stop(timer_remap_bil);
} /* scrip_remap_weights_bilinear */

/*
  -----------------------------------------------------------------------

  This routine computes and apply the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_bilinear(remapgrid_t* src_grid, remapgrid_t* tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
  /*   Local variables */
  int  search_result;
  long tgt_cell_add;             /*  destination addresss                   */
  int src_add[4];                /*  address for the four source points     */
  double src_lats[4];            /*  latitudes  of four bilinear corners    */
  double src_lons[4];            /*  longitudes of four bilinear corners    */
  double wgts[4];                /*  bilinear weights for four corners      */
  double plat, plon;             /*  lat/lon coords of destination point    */
  extern int timer_remap_bil;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  if ( cdoTimer ) timer_start(timer_remap_bil);

  progressInit();

  long tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bilinear interpolation when source grid rank != 2"); 

  double findex = 0;

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoVerbose, cdoSilentMode, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, src_array, tgt_array, missval, findex) \
  private(tgt_cell_add, src_add, src_lats, src_lons, wgts, plat, plon, search_result)    \
  schedule(static)
#endif
  for ( tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      int lprogress = 1;
      if ( cdo_omp_get_thread_num() != 0 ) lprogress = 0;

      if ( !cdoSilentMode )
	{
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
	  findex++;
	  if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);
	}

      tgt_array[tgt_cell_add] = missval;

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      if ( tgt_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
        {
          long nx = tgt_grid->dims[0];
	  long iy = tgt_cell_add/nx;
	  long ix = tgt_cell_add - iy*nx;

          plat = tgt_grid->reg2d_center_lat[iy];
          plon = tgt_grid->reg2d_center_lon[ix];
        }
      else
        {
          plat = tgt_grid->cell_center_lat[tgt_cell_add];
          plon = tgt_grid->cell_center_lon[tgt_cell_add];
        }

      /* Find nearest square of grid points on source grid  */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);

      /* Check to see if points are land points */
      if ( search_result > 0 )
	{
	  for ( int n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[tgt_cell_add] = 1.;

          if ( find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw) )
	    {
	      /* Successfully found iw,jw - compute weights */
	      set_bilinear_weights(iw, jw, wgts);

	      sort_add_and_wgts(4, src_add, wgts);

	      bilinear_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add);
	    }
          else
	    {
	      bilinear_warning(plon, plat, iw, jw, src_add, src_lons, src_lats, src_grid);

	      search_result = -1;
	    }
	}

      /*
	Search for bilinear failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
          if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[tgt_cell_add] = 1.;

	      sort_add_and_wgts(4, src_add, wgts);

	      bilinear_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add);
	    }
        }
    }

  if ( cdoTimer ) timer_stop(timer_remap_bil);
} /* scrip_remap_bilinear */
