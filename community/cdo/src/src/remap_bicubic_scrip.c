#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BICUBIC INTERPOLATION                                              */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static
void set_bicubic_weights(double iw, double jw, double wgts[4][4])
{
  wgts[0][0] = (1.-jw*jw*(3.-2.*jw))  * (1.-iw*iw*(3.-2.*iw));
  wgts[1][0] = (1.-jw*jw*(3.-2.*jw))  *     iw*iw*(3.-2.*iw);
  wgts[2][0] =     jw*jw*(3.-2.*jw)   *     iw*iw*(3.-2.*iw);
  wgts[3][0] =     jw*jw*(3.-2.*jw)   * (1.-iw*iw*(3.-2.*iw));
  wgts[0][1] = (1.-jw*jw*(3.-2.*jw))  *     iw*(iw-1.)*(iw-1.);
  wgts[1][1] = (1.-jw*jw*(3.-2.*jw))  *     iw*iw*(iw-1.);
  wgts[2][1] =     jw*jw*(3.-2.*jw)   *     iw*iw*(iw-1.);
  wgts[3][1] =     jw*jw*(3.-2.*jw)   *     iw*(iw-1.)*(iw-1.);
  wgts[0][2] =     jw*(jw-1.)*(jw-1.) * (1.-iw*iw*(3.-2.*iw));
  wgts[1][2] =     jw*(jw-1.)*(jw-1.) *     iw*iw*(3.-2.*iw);
  wgts[2][2] =     jw*jw*(jw-1.)      *     iw*iw*(3.-2.*iw);
  wgts[3][2] =     jw*jw*(jw-1.)      * (1.-iw*iw*(3.-2.*iw));
  wgts[0][3] =     iw*(iw-1.)*(iw-1.) *     jw*(jw-1.)*(jw-1.);
  wgts[1][3] =     iw*iw*(iw-1.)      *     jw*(jw-1.)*(jw-1.);
  wgts[2][3] =     iw*iw*(iw-1.)      *     jw*jw*(jw-1.);
  wgts[3][3] =     iw*(iw-1.)*(iw-1.) *     jw*jw*(jw-1.);
}

int num_src_points(const int* restrict mask, const int src_add[4], double src_lats[4]);

static
void renormalize_weights(const double src_lats[4], double wgts[4][4])
{
  int n;
  double sum_wgts = 0.0; /* sum of weights for normalization */
  /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
  for ( n = 0; n < 4; ++n ) sum_wgts  += fabs(src_lats[n]);
  for ( n = 0; n < 4; ++n ) wgts[n][0] = fabs(src_lats[n])/sum_wgts;
  for ( n = 0; n < 4; ++n ) wgts[n][1] = 0.;
  for ( n = 0; n < 4; ++n ) wgts[n][2] = 0.;
  for ( n = 0; n < 4; ++n ) wgts[n][3] = 0.;
}

static
void bicubic_warning(void)
{
  static int lwarn = TRUE;

  if ( cdoVerbose || lwarn )
    {
      lwarn = FALSE;
      // cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
      cdoWarning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

static
void bicubic_remap(double* restrict tgt_point, const double* restrict src_array, double wgts[4][4], const int src_add[4],
		   const double* restrict grad1, const double* restrict grad2, const double* restrict grad3)
{
  *tgt_point = 0.;
  for ( int n = 0; n < 4; ++n )
    *tgt_point += src_array[src_add[n]]*wgts[n][0] +
                      grad1[src_add[n]]*wgts[n][1] +
                      grad2[src_add[n]]*wgts[n][2] +
                      grad3[src_add[n]]*wgts[n][3];  
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_weights_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*   Local variables */
  int  search_result;
  long tgt_cell_add;        /*  destination addresss                   */
  int src_add[4];      /*  address for the four source points     */
  double src_lats[4];  /*  latitudes  of four bilinear corners    */
  double src_lons[4];  /*  longitudes of four bilinear corners    */
  double wgts[4][4];   /*  bicubic weights for four corners       */
  double plat, plon;   /*  lat/lon coords of destination point    */
  extern int timer_remap_bic;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  if ( cdoTimer ) timer_start(timer_remap_bic);

  progressInit();

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bicubic interpolation when source grid rank != 2"); 

  long tgt_grid_size = tgt_grid->size;

  weightlinks4_t *weightlinks = (weightlinks4_t *) malloc(tgt_grid_size*sizeof(weightlinks4_t));

  /* Loop over destination grid */

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoVerbose, weightlinks, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, findex) \
  private(tgt_cell_add, src_add, src_lats, src_lons, wgts, plat, plon, search_result)
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
	      set_bicubic_weights(iw, jw, wgts);

	      store_weightlinks4(4, src_add, wgts, tgt_cell_add, weightlinks);
	    }
          else
	    {
	      bicubic_warning();

	      search_result = -1;
	    }
	}
	  
      /*
	Search for bicubic failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
	  if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[tgt_cell_add] = 1.;

	      store_weightlinks4(4, src_add, wgts, tgt_cell_add, weightlinks);
	    }
        }
    }

  if ( cdoTimer ) timer_stop(timer_remap_bic);

  weightlinks2remaplinks4(tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) free(weightlinks);

} /* scrip_remap_weights_bicubic */

/*
  -----------------------------------------------------------------------

  This routine computes ans apply the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
  /*   Local variables */
  int  search_result;
  long tgt_cell_add;        /*  destination addresss                 */
  int src_add[4];      /*  address for the four source points   */
  double src_lats[4];  /*  latitudes  of four bilinear corners  */
  double src_lons[4];  /*  longitudes of four bilinear corners  */
  double wgts[4][4];   /*  bicubic weights for four corners     */
  double plat, plon;   /*  lat/lon coords of destination point  */
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  long tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bicubic interpolation when source grid rank != 2"); 

  double *grad1_lat    = (double*) malloc(src_grid->size*sizeof(double));
  double *grad1_lon    = (double*) malloc(src_grid->size*sizeof(double));
  double *grad1_latlon = (double*) malloc(src_grid->size*sizeof(double));

  remap_gradients(*src_grid, src_array, grad1_lat, grad1_lon, grad1_latlon);

  /* Loop over destination grid */

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, src_array, tgt_array, missval, grad1_lat, grad1_lon, grad1_latlon, findex) \
  private(tgt_cell_add, src_add, src_lats, src_lons, wgts, plat, plon, search_result)
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

      tgt_array[tgt_cell_add] = missval;

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
	      set_bicubic_weights(iw, jw, wgts);

	      sort_add_and_wgts4(4, src_add, wgts);

	      bicubic_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add, grad1_lat, grad1_lon, grad1_latlon);
	    }
          else
	    {
	      bicubic_warning();

	      search_result = -1;
	    }
	}
	  
      /*
	Search for bicubic failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
	  if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[tgt_cell_add] = 1.;

	      sort_add_and_wgts4(4, src_add, wgts);

	      bicubic_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add, grad1_lat, grad1_lon, grad1_latlon);
	    }
        }
    }

  free(grad1_lat);
  free(grad1_lon);
  free(grad1_latlon);

} /* scrip_remap_bicubic */
