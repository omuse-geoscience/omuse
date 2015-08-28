#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

#include <time.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"


void remapgrid_init(remapgrid_t *grid);
void remapGridRealloc(int map_type, remapgrid_t *grid);


#if defined(HAVE_LIBNETCDF)
static
void nce(int istat)
{
  /*
    This routine provides a simple interface to netCDF error message routine.
  */

  if ( istat != NC_NOERR ) cdoAbort(nc_strerror(istat));
}
#endif


void write_remap_scrip(const char *interp_file, int map_type, int submap_type, int num_neighbors,
		       int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv)
{
  /*
    Writes remap data to a netCDF file using SCRIP conventions
  */
  /*
    Input variables:

    interp_file  ! filename for remap data
  */

#if defined(HAVE_LIBNETCDF)

  /* Local variables */

  int nc_file_id;           /* id for netCDF file                       */
  int nc_srcgrdsize_id;     /* id for source grid size                  */
  int nc_dstgrdsize_id;     /* id for destination grid size             */
  int nc_srcgrdcorn_id = 0; /* id for number of source grid corners     */
  int nc_dstgrdcorn_id = 0; /* id for number of dest grid corners       */
  int nc_srcgrdrank_id;     /* id for source grid rank                  */
  int nc_dstgrdrank_id;     /* id for dest grid rank                    */
  int nc_numlinks_id;       /* id for number of links in mapping        */
  int nc_numwgts_id;        /* id for number of weights for mapping     */
  int nc_srcgrddims_id;     /* id for source grid dimensions            */
  int nc_dstgrddims_id;     /* id for dest grid dimensions              */
  int nc_srcgrdcntrlat_id;  /* id for source grid center latitude       */
  int nc_dstgrdcntrlat_id;  /* id for dest grid center latitude         */
  int nc_srcgrdcntrlon_id;  /* id for source grid center longitude      */
  int nc_dstgrdcntrlon_id;  /* id for dest grid center longitude        */
  int nc_srcgrdimask_id;    /* id for source grid mask                  */
  int nc_dstgrdimask_id;    /* id for dest grid mask                    */
  int nc_srcgrdcrnrlat_id;  /* id for latitude of source grid corners   */
  int nc_srcgrdcrnrlon_id;  /* id for longitude of source grid corners  */
  int nc_dstgrdcrnrlat_id;  /* id for latitude of dest grid corners     */
  int nc_dstgrdcrnrlon_id;  /* id for longitude of dest grid corners    */
  int nc_srcgrdarea_id;     /* id for area of source grid cells         */
  int nc_dstgrdarea_id;     /* id for area of dest grid cells           */
  int nc_srcgrdfrac_id;     /* id for area fraction on source grid      */
  int nc_dstgrdfrac_id;     /* id for area fraction on dest grid        */
  int nc_srcadd_id;         /* id for map source address                */
  int nc_dstadd_id;         /* id for map destination address           */
  int nc_rmpmatrix_id;      /* id for remapping matrix                  */

  int nc_dims2_id[2];       /* netCDF ids for 2d array dims             */

  char *map_name = "SCRIP remapping with CDO";
  char normalize_opt[64] = "unknown";
  char map_method[64] = "unknown";
  char tmp_string[64] = "unknown";
  char history[1024] = "date and time";
  char src_grid_name[64] = "source grid";
  char tgt_grid_name[64] = "dest grid";
  char *src_grid_units = "radians";
  char *tgt_grid_units = "radians";
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  long i;
  int lgridarea = FALSE;
  int writemode = NC_CLOBBER;

  switch ( rv.norm_opt )
    {
    case NORM_OPT_NONE:
      strcpy(normalize_opt, "none");
      break;
    case NORM_OPT_FRACAREA:
      strcpy(normalize_opt, "fracarea");
      break;
    case NORM_OPT_DESTAREA:
      strcpy(normalize_opt, "destarea");
      break;
    }

  switch ( map_type )
    {
    case MAP_TYPE_CONSERV:
      lgridarea = TRUE;
      if ( submap_type == SUBMAP_TYPE_LAF )
	{
	  strcpy(map_method, "Largest area fraction");
	  break;
	}
      else
	{
	  strcpy(map_method, "Conservative remapping");
	  break;
	}
    case MAP_TYPE_CONSERV_YAC:
      lgridarea = TRUE;
      /*
      if ( submap_type == SUBMAP_TYPE_LAF )
	{
	  strcpy(map_method, "Largest area fraction");
	  break;
	}
      else
      */
	{
	  strcpy(map_method, "Conservative remapping using clipping on sphere");
	  break;
	}
    case MAP_TYPE_BILINEAR:
      strcpy(map_method, "Bilinear remapping");
      break;
    case MAP_TYPE_BICUBIC:
      strcpy(map_method, "Bicubic remapping");
      break;
    case MAP_TYPE_DISTWGT:
      if ( num_neighbors == 1 )
	strcpy(map_method, "Nearest neighbor");
      else
	strcpy(map_method, "Distance weighted avg of nearest neighbors");
      break;
    }

  if ( rv.num_links == 0 )
    cdoAbort("Number of remap links is 0, no remap weights found!");

  {
    size_t nele1 = 4*8 + 4;
    size_t nele2 = 4*8 + 4;
    if ( src_grid.lneed_cell_corners ) nele1 += src_grid.num_cell_corners*2*8;
    if ( tgt_grid.lneed_cell_corners ) nele2 += tgt_grid.num_cell_corners*2*8;
    size_t filesize = src_grid.size*(nele1) +
                      tgt_grid.size*(nele2) +
                      rv.num_links*(4 + 4 + rv.num_wts*8);

    if ( cdoVerbose )
      cdoPrint("Filesize for remap weights: ~%lu", (unsigned long) filesize);
    
    if ( filesize > 0x7FFFFC00 ) /* 2**31 - 1024 (<2GB) */
      {
#if defined(NC_64BIT_OFFSET)
	writemode = NC_CLOBBER | NC_64BIT_OFFSET;
#else
	cdoAbort("Filesize for remap weights maybe too large!");
#endif
      }
  }

  /* Create netCDF file for mapping and define some global attributes */
  nce(nc_create(interp_file, writemode, &nc_file_id));

  /* Map name */
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "title", strlen(map_name), map_name));

  /* Normalization option */
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "normalization", strlen(normalize_opt), normalize_opt));

  /* Map method */
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "map_method", strlen(map_method), map_method));

  /* Remap order */
  if ( map_type == MAP_TYPE_CONSERV && submap_type == SUBMAP_TYPE_NONE )
    nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "remap_order", NC_INT, 1L, &remap_order));

  /* History */
  date_and_time_in_sec = time(NULL);

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(history, 1024, "%d %b %Y : ", date_and_time);
      strcat(history, commandLine());
    }

  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "history", strlen(history), history));

  /* File convention */
  strcpy(tmp_string, "SCRIP");
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "conventions", strlen(tmp_string), tmp_string));

  /* Source and destination grid names */
  gridName(gridInqType(src_grid.gridID), src_grid_name);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "source_grid", strlen(src_grid_name), src_grid_name));

  gridName(gridInqType(tgt_grid.gridID), tgt_grid_name);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "dest_grid", strlen(tgt_grid_name), tgt_grid_name));

  /* Prepare netCDF dimension info */

  /* Define grid size dimensions */
  nce(nc_def_dim(nc_file_id, "src_grid_size", src_grid.size, &nc_srcgrdsize_id));
  nce(nc_def_dim(nc_file_id, "dst_grid_size", tgt_grid.size, &nc_dstgrdsize_id));

  /* Define grid corner dimension */
  if ( src_grid.lneed_cell_corners )
    nce(nc_def_dim(nc_file_id, "src_grid_corners", src_grid.num_cell_corners, &nc_srcgrdcorn_id));
  if ( tgt_grid.lneed_cell_corners )
    nce(nc_def_dim(nc_file_id, "dst_grid_corners", tgt_grid.num_cell_corners, &nc_dstgrdcorn_id));

  /* Define grid rank dimension */
  nce(nc_def_dim(nc_file_id, "src_grid_rank", src_grid.rank, &nc_srcgrdrank_id));
  nce(nc_def_dim(nc_file_id, "dst_grid_rank", tgt_grid.rank, &nc_dstgrdrank_id));

  /* Define map size dimensions */
  nce(nc_def_dim(nc_file_id, "num_links", rv.num_links, &nc_numlinks_id));
  nce(nc_def_dim(nc_file_id, "num_wgts", rv.num_wts, &nc_numwgts_id));
       
  /* Define grid dimensions */
  nce(nc_def_var(nc_file_id, "src_grid_dims", NC_INT, 1, &nc_srcgrdrank_id, &nc_srcgrddims_id));
  nce(nc_def_var(nc_file_id, "dst_grid_dims", NC_INT, 1, &nc_dstgrdrank_id, &nc_dstgrddims_id));

  /* Define all arrays for netCDF descriptors */

  /* Define grid center latitude array */
  nce(nc_def_var(nc_file_id, "src_grid_center_lat", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlat_id));
  nce(nc_def_var(nc_file_id, "dst_grid_center_lat", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlat_id));

  /* Define grid center longitude array */
  nce(nc_def_var(nc_file_id, "src_grid_center_lon", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlon_id));
  nce(nc_def_var(nc_file_id, "dst_grid_center_lon", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlon_id));

  /* Define grid corner lat/lon arrays */

  nc_dims2_id[0] = nc_srcgrdsize_id;
  nc_dims2_id[1] = nc_srcgrdcorn_id;

  if ( src_grid.lneed_cell_corners )
    {
      nce(nc_def_var(nc_file_id, "src_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id, &nc_srcgrdcrnrlat_id));
      nce(nc_def_var(nc_file_id, "src_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id, &nc_srcgrdcrnrlon_id));
    }

  nc_dims2_id[0] = nc_dstgrdsize_id;
  nc_dims2_id[1] = nc_dstgrdcorn_id;

  if ( tgt_grid.lneed_cell_corners )
    {
      nce(nc_def_var(nc_file_id, "dst_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id, &nc_dstgrdcrnrlat_id));
      nce(nc_def_var(nc_file_id, "dst_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id, &nc_dstgrdcrnrlon_id));
    }

  /* Define units for all coordinate arrays */
  nce(nc_put_att_text(nc_file_id, nc_srcgrdcntrlat_id, "units", strlen(src_grid_units), src_grid_units));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdcntrlat_id, "units", strlen(tgt_grid_units), tgt_grid_units));
  nce(nc_put_att_text(nc_file_id, nc_srcgrdcntrlon_id, "units", strlen(src_grid_units), src_grid_units));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdcntrlon_id, "units", strlen(tgt_grid_units), tgt_grid_units));
  if ( src_grid.lneed_cell_corners )
    {
      nce(nc_put_att_text(nc_file_id, nc_srcgrdcrnrlat_id, "units", strlen(src_grid_units), src_grid_units));
      nce(nc_put_att_text(nc_file_id, nc_srcgrdcrnrlon_id, "units", strlen(src_grid_units), src_grid_units));
    }
  if ( tgt_grid.lneed_cell_corners )
    {
      nce(nc_put_att_text(nc_file_id, nc_dstgrdcrnrlat_id, "units", strlen(tgt_grid_units), tgt_grid_units));
      nce(nc_put_att_text(nc_file_id, nc_dstgrdcrnrlon_id, "units", strlen(tgt_grid_units), tgt_grid_units));
    }

  /* Define grid mask */

  nce(nc_def_var(nc_file_id, "src_grid_imask", NC_INT, 1, &nc_srcgrdsize_id, &nc_srcgrdimask_id));
  nce(nc_put_att_text(nc_file_id, nc_srcgrdimask_id, "units", 8, "unitless"));

  nce(nc_def_var(nc_file_id, "dst_grid_imask", NC_INT, 1, &nc_dstgrdsize_id, &nc_dstgrdimask_id));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdimask_id, "units", 8, "unitless"));

  /* Define grid area arrays */

  if ( lgridarea )
    {
      nce(nc_def_var(nc_file_id, "src_grid_area", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdarea_id));
      nce(nc_put_att_text(nc_file_id, nc_srcgrdarea_id, "units", 14, "square radians"));

      nce(nc_def_var(nc_file_id, "dst_grid_area", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdarea_id));   
      nce(nc_put_att_text(nc_file_id, nc_dstgrdarea_id, "units", 14, "square radians"));
    }

  /* Define grid fraction arrays */

  nce(nc_def_var(nc_file_id, "src_grid_frac", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdfrac_id));
  nce(nc_put_att_text(nc_file_id, nc_srcgrdfrac_id, "units", 8, "unitless"));

  nce(nc_def_var(nc_file_id, "dst_grid_frac", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdfrac_id));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdfrac_id, "units", 8, "unitless"));

  /* Define mapping arrays */

  nce(nc_def_var(nc_file_id, "src_address", NC_INT, 1, &nc_numlinks_id, &nc_srcadd_id));      
  nce(nc_def_var(nc_file_id, "dst_address", NC_INT, 1, &nc_numlinks_id, &nc_dstadd_id));

  nc_dims2_id[0] = nc_numlinks_id;
  nc_dims2_id[1] = nc_numwgts_id;

  nce(nc_def_var(nc_file_id, "remap_matrix", NC_DOUBLE, 2, nc_dims2_id, &nc_rmpmatrix_id));

  /* End definition stage */
  
  nce(nc_enddef(nc_file_id));


  /* Write mapping data */

  nce(nc_put_var_int(nc_file_id, nc_srcgrddims_id, src_grid.dims));
  nce(nc_put_var_int(nc_file_id, nc_dstgrddims_id, tgt_grid.dims));

  nce(nc_put_var_int(nc_file_id, nc_srcgrdimask_id, src_grid.mask));
  nce(nc_put_var_int(nc_file_id, nc_dstgrdimask_id, tgt_grid.mask));

  if ( src_grid.cell_center_lat ) nce(nc_put_var_double(nc_file_id, nc_srcgrdcntrlat_id, src_grid.cell_center_lat)); 
  if ( src_grid.cell_center_lon ) nce(nc_put_var_double(nc_file_id, nc_srcgrdcntrlon_id, src_grid.cell_center_lon));

  if ( src_grid.lneed_cell_corners )
    {
      nce(nc_put_var_double(nc_file_id, nc_srcgrdcrnrlat_id, src_grid.cell_corner_lat));
      nce(nc_put_var_double(nc_file_id, nc_srcgrdcrnrlon_id, src_grid.cell_corner_lon));
    }

  if ( tgt_grid.cell_center_lat ) nce(nc_put_var_double(nc_file_id, nc_dstgrdcntrlat_id, tgt_grid.cell_center_lat));
  if ( tgt_grid.cell_center_lon ) nce(nc_put_var_double(nc_file_id, nc_dstgrdcntrlon_id, tgt_grid.cell_center_lon));

  if ( tgt_grid.lneed_cell_corners )
    {
      nce(nc_put_var_double(nc_file_id, nc_dstgrdcrnrlat_id, tgt_grid.cell_corner_lat));
      nce(nc_put_var_double(nc_file_id, nc_dstgrdcrnrlon_id, tgt_grid.cell_corner_lon));
    }

  if ( lgridarea )
    nce(nc_put_var_double(nc_file_id, nc_srcgrdarea_id, src_grid.cell_area));

  nce(nc_put_var_double(nc_file_id, nc_srcgrdfrac_id, src_grid.cell_frac));

  /*
  if ( luse_cell_area )
    nce(nc_put_var_double(nc_file_id, nc_dstgrdarea_id, tgt_grid.cell_area_in));
  else
  */
  if ( lgridarea )
    nce(nc_put_var_double(nc_file_id, nc_dstgrdarea_id, tgt_grid.cell_area));

  nce(nc_put_var_double(nc_file_id, nc_dstgrdfrac_id, tgt_grid.cell_frac));

  for ( i = 0; i < rv.num_links; i++ )
    {
      rv.src_cell_add[i]++;
      rv.tgt_cell_add[i]++;
    }

  nce(nc_put_var_int(nc_file_id, nc_srcadd_id, rv.src_cell_add));
  nce(nc_put_var_int(nc_file_id, nc_dstadd_id, rv.tgt_cell_add));

  nce(nc_put_var_double(nc_file_id, nc_rmpmatrix_id, rv.wts));

  nce(nc_close(nc_file_id));

#else
  cdoAbort("netCDF support not compiled in!");
#endif

}  /* write_remap_scrip */

/*****************************************************************************/

void read_remap_scrip(const char *interp_file, int gridID1, int gridID2, int *map_type, int *submap_type, int *num_neighbors,
		      int *remap_order, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*
    The routine reads a netCDF file to extract remapping info in SCRIP format
  */
  /*
    Input variables

    interp_file        ! filename for remap data
  */
#if defined(HAVE_LIBNETCDF)

  /* Local variables */

  int lgridarea = FALSE;
  int status;
  int nc_file_id;           /* id for netCDF file                       */
  int nc_srcgrdsize_id;     /* id for source grid size                  */
  int nc_dstgrdsize_id;     /* id for destination grid size             */
  int nc_srcgrdcorn_id;     /* id for number of source grid corners     */
  int nc_dstgrdcorn_id;     /* id for number of dest grid corners       */
  int nc_srcgrdrank_id;     /* id for source grid rank                  */
  int nc_dstgrdrank_id;     /* id for dest grid rank                    */
  int nc_numlinks_id;       /* id for number of links in mapping        */
  int nc_numwgts_id;        /* id for number of weights for mapping     */
  int nc_srcgrddims_id;     /* id for source grid dimensions            */
  int nc_dstgrddims_id;     /* id for dest grid dimensions              */
  int nc_srcgrdcntrlat_id;  /* id for source grid center latitude       */
  int nc_dstgrdcntrlat_id;  /* id for dest grid center latitude         */
  int nc_srcgrdcntrlon_id;  /* id for source grid center longitude      */
  int nc_dstgrdcntrlon_id;  /* id for dest grid center longitude        */
  int nc_srcgrdimask_id;    /* id for source grid mask                  */
  int nc_dstgrdimask_id;    /* id for dest grid mask                    */
  int nc_srcgrdcrnrlat_id;  /* id for latitude of source grid corners   */
  int nc_srcgrdcrnrlon_id;  /* id for longitude of source grid corners  */
  int nc_dstgrdcrnrlat_id;  /* id for latitude of dest grid corners     */
  int nc_dstgrdcrnrlon_id;  /* id for longitude of dest grid corners    */
  int nc_srcgrdarea_id;     /* id for area of source grid cells         */
  int nc_dstgrdarea_id;     /* id for area of dest grid cells           */
  int nc_srcgrdfrac_id;     /* id for area fraction on source grid      */
  int nc_dstgrdfrac_id;     /* id for area fraction on dest grid        */
  int nc_srcadd_id;         /* id for map source address                */
  int nc_dstadd_id;         /* id for map destination address           */
  int nc_rmpmatrix_id;      /* id for remapping matrix                  */

  long i;                   /* dummy index */

  char map_name[1024];
  char map_method[64];      /* character string for map_type             */
  char normalize_opt[64];   /* character string for normalization option */
  char convention[64];      /* character string for output convention    */
  char src_grid_name[64];      /* grid name for source grid                 */
  char tgt_grid_name[64];      /* grid name for dest   grid                 */
  char src_grid_units[64];
  char tgt_grid_units[64];
  size_t attlen, dimlen;

  int gridID1_gme_c = -1;


  /* Open file and read some global information */

  /* nce(nc_open(interp_file, NC_NOWRITE, &nc_file_id)); */
  nc_file_id = cdf_openread(interp_file);

  /* Map name */

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "title", map_name));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "title", &attlen));
  map_name[attlen] = 0;

  if ( cdoVerbose )
    {
      cdoPrint("Reading remapping: %s", map_name);
      cdoPrint("From file: %s", interp_file);
    }

  /* Normalization option */

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "normalization", normalize_opt));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "normalization", &attlen));
  normalize_opt[attlen] = 0;

  rv->sort_add = FALSE;

  if ( strcmp(normalize_opt, "none") == 0 )
    rv->norm_opt = NORM_OPT_NONE;
  else if ( strcmp(normalize_opt, "fracarea") == 0 )
    rv->norm_opt = NORM_OPT_FRACAREA;
  else if ( strcmp(normalize_opt, "destarea") == 0 )
    rv->norm_opt = NORM_OPT_DESTAREA;
  else
    {
      cdoPrint("normalize_opt = %s", normalize_opt);
      cdoAbort("Invalid normalization option");
    }

  if ( cdoVerbose )
    cdoPrint("normalize_opt = %s", normalize_opt);

  /* Map method */

  nce(nc_get_att_text (nc_file_id, NC_GLOBAL, "map_method", map_method));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "map_method", &attlen));
  map_method[attlen] = 0;

  *submap_type = SUBMAP_TYPE_NONE;
  *remap_order = 1;

  if ( cmpstr(map_method, "Conservative") == 0 )
    {
      int iatt;
      if ( cmpstr(map_method, "Conservative remapping using clipping on sphere") == 0 )
	rv->map_type = MAP_TYPE_CONSERV_YAC;
      else
	rv->map_type = MAP_TYPE_CONSERV;

      status = nc_get_att_int(nc_file_id, NC_GLOBAL, "remap_order", &iatt);
      if ( status == NC_NOERR ) *remap_order = iatt;
    }
  else if ( cmpstr(map_method, "Bilinear") == 0 ) rv->map_type = MAP_TYPE_BILINEAR;
  else if ( cmpstr(map_method, "Bicubic")  == 0 ) rv->map_type = MAP_TYPE_BICUBIC;
  else if ( cmpstr(map_method, "Distance") == 0 )
    {
      rv->map_type = MAP_TYPE_DISTWGT;
      *num_neighbors = 4;
    }
  else if ( cmpstr(map_method, "Nearest") == 0 )
    {
      rv->map_type = MAP_TYPE_DISTWGT;
      *num_neighbors = 1;
    }
  else if ( cmpstr(map_method, "Largest") == 0 )
    {
      rv->map_type = MAP_TYPE_CONSERV;
      *submap_type = SUBMAP_TYPE_LAF;
    }
  else
    {
      cdoPrint("map_type = %s", map_method);
      cdoAbort("Invalid Map Type");
    }

  if ( cdoVerbose )
    cdoPrint("map_type = %s", map_method);

  if ( rv->map_type == MAP_TYPE_CONSERV ) lgridarea = TRUE;

  *map_type = rv->map_type;

  /* File convention */

  nce(nc_get_att_text (nc_file_id, NC_GLOBAL, "conventions", convention));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "conventions", &attlen));
  convention[attlen] = 0;

  if ( strcmp(convention, "SCRIP") != 0 )
    {
      cdoPrint("convention = %s", convention);
      if ( strcmp(convention, "NCAR-CSM") == 0 )
        cdoAbort("Unsupported file convention!");
      else
        cdoAbort("Unknown file convention!");
    }

  /* Read some additional global attributes */

  /* Source and destination grid names */

  nce(nc_get_att_text (nc_file_id, NC_GLOBAL, "source_grid", src_grid_name));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "source_grid", &attlen));
  src_grid_name[attlen] = 0;
  
  nce(nc_get_att_text (nc_file_id, NC_GLOBAL, "dest_grid", tgt_grid_name));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "dest_grid", &attlen));
  tgt_grid_name[attlen] = 0;
 
  if ( cdoVerbose )
    cdoPrint("Remapping between: %s and %s", src_grid_name, tgt_grid_name);

  /* Initialize remapgrid structure */
  remapgrid_init(src_grid);
  remapgrid_init(tgt_grid);

  /* Read dimension information */

  nce(nc_inq_dimid(nc_file_id, "src_grid_size", &nc_srcgrdsize_id));
  nce(nc_inq_dimlen(nc_file_id, nc_srcgrdsize_id, &dimlen));
  src_grid->size = dimlen;
  /*
  if (  src_grid->size != gridInqSize(gridID1) )
    cdoAbort("Source grids have different size!");
  */
  nce(nc_inq_dimid(nc_file_id, "dst_grid_size", &nc_dstgrdsize_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dstgrdsize_id, &dimlen));
  tgt_grid->size = dimlen;
  /*
  if ( tgt_grid->size != gridInqSize(gridID2) )
    cdoAbort("Target grids have different size!");
  */
  status = nc_inq_dimid(nc_file_id, "src_grid_corners", &nc_srcgrdcorn_id);
  if ( status == NC_NOERR )
    {
      nce(nc_inq_dimlen(nc_file_id, nc_srcgrdcorn_id, &dimlen));
      src_grid->num_cell_corners = dimlen;
      src_grid->luse_cell_corners = TRUE;
      src_grid->lneed_cell_corners = TRUE;
    }

  status = nc_inq_dimid(nc_file_id, "dst_grid_corners", &nc_dstgrdcorn_id);
  if ( status == NC_NOERR )
    {
      nce(nc_inq_dimlen(nc_file_id, nc_dstgrdcorn_id, &dimlen));
      tgt_grid->num_cell_corners = dimlen;
      tgt_grid->luse_cell_corners = TRUE;
      tgt_grid->lneed_cell_corners = TRUE;
    }

  nce(nc_inq_dimid(nc_file_id, "src_grid_rank", &nc_srcgrdrank_id));
  nce(nc_inq_dimlen(nc_file_id, nc_srcgrdrank_id, &dimlen));
  src_grid->rank = dimlen;

  nce(nc_inq_dimid(nc_file_id, "dst_grid_rank", &nc_dstgrdrank_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dstgrdrank_id, &dimlen));
  tgt_grid->rank = dimlen;

  nce(nc_inq_dimid(nc_file_id, "num_links", &nc_numlinks_id));
  nce(nc_inq_dimlen(nc_file_id, nc_numlinks_id, &dimlen));
  rv->num_links = dimlen;

  if ( rv->num_links == 0 )
    cdoAbort("Number of remap links is 0, no remap weights found!");

  nce(nc_inq_dimid(nc_file_id, "num_wgts", &nc_numwgts_id));
  nce(nc_inq_dimlen(nc_file_id, nc_numwgts_id, &dimlen));
  rv->num_wts = dimlen;

  src_grid->gridID = gridID1;
  tgt_grid->gridID = gridID2;

  if ( gridID1 != -1 && gridInqType(gridID1) == GRID_GME )
    {
      src_grid->nvgp = gridInqSize(gridID1);
      gridID1_gme_c = gridToUnstructured(gridID1, 1);
    }

  remapGridRealloc(rv->map_type, src_grid);
  remapGridRealloc(rv->map_type, tgt_grid);

  if ( gridID1 != -1 && gridInqType(gridID1) == GRID_GME ) gridInqMaskGME(gridID1_gme_c, src_grid->vgpm);    

  rv->pinit = TRUE;
  rv->wts = NULL;

  rv->max_links = rv->num_links;

  rv->resize_increment = (int) (0.1 * MAX(src_grid->size, tgt_grid->size));

  /* Allocate address and weight arrays for mapping 1 */

  rv->src_cell_add = (int*) malloc(rv->num_links*sizeof(int));
  rv->tgt_cell_add = (int*) malloc(rv->num_links*sizeof(int));

  rv->wts = (double*) malloc(rv->num_wts*rv->num_links*sizeof(double));

  /* Get variable ids */

  nce(nc_inq_varid(nc_file_id, "src_grid_dims", &nc_srcgrddims_id));
  nce(nc_inq_varid(nc_file_id, "src_grid_imask", &nc_srcgrdimask_id));
  nce(nc_inq_varid(nc_file_id, "src_grid_center_lat", &nc_srcgrdcntrlat_id));
  nce(nc_inq_varid(nc_file_id, "src_grid_center_lon", &nc_srcgrdcntrlon_id));
  if ( src_grid->num_cell_corners )
    {
      nce(nc_inq_varid(nc_file_id, "src_grid_corner_lat", &nc_srcgrdcrnrlat_id));
      nce(nc_inq_varid(nc_file_id, "src_grid_corner_lon", &nc_srcgrdcrnrlon_id));
    }
  if ( lgridarea )
    {
      nce(nc_inq_varid(nc_file_id, "src_grid_area", &nc_srcgrdarea_id));
    }
  nce(nc_inq_varid(nc_file_id, "src_grid_frac", &nc_srcgrdfrac_id));

  nce(nc_inq_varid(nc_file_id, "dst_grid_dims", &nc_dstgrddims_id));
  nce(nc_inq_varid(nc_file_id, "dst_grid_imask", &nc_dstgrdimask_id));
  nce(nc_inq_varid(nc_file_id, "dst_grid_center_lat", &nc_dstgrdcntrlat_id));
  nce(nc_inq_varid(nc_file_id, "dst_grid_center_lon", &nc_dstgrdcntrlon_id));
  if ( tgt_grid->num_cell_corners )
    {
      nce(nc_inq_varid(nc_file_id, "dst_grid_corner_lat", &nc_dstgrdcrnrlat_id));
      nce(nc_inq_varid(nc_file_id, "dst_grid_corner_lon", &nc_dstgrdcrnrlon_id));
    }
  if ( lgridarea )
    {
      nce(nc_inq_varid(nc_file_id, "dst_grid_area", &nc_dstgrdarea_id));
    }
  nce(nc_inq_varid(nc_file_id, "dst_grid_frac", &nc_dstgrdfrac_id));

  nce(nc_inq_varid(nc_file_id, "src_address", &nc_srcadd_id));
  nce(nc_inq_varid(nc_file_id, "dst_address", &nc_dstadd_id));
  nce(nc_inq_varid(nc_file_id, "remap_matrix", &nc_rmpmatrix_id));

  /* Read all variables */

  nce(nc_get_var_int(nc_file_id, nc_srcgrddims_id, src_grid->dims));

  nce(nc_get_var_int(nc_file_id, nc_srcgrdimask_id, src_grid->mask));

  nce(nc_get_var_double(nc_file_id, nc_srcgrdcntrlat_id, src_grid->cell_center_lat));
  nce(nc_get_var_double(nc_file_id, nc_srcgrdcntrlon_id, src_grid->cell_center_lon));

  nce(nc_get_att_text(nc_file_id, nc_srcgrdcntrlat_id, "units", src_grid_units));
  nce(nc_inq_attlen(nc_file_id, nc_srcgrdcntrlat_id, "units", &attlen));
  src_grid_units[attlen] = 0;

  grid_to_radian(src_grid_units, src_grid->size, src_grid->cell_center_lon, "source grid center lon"); 
  grid_to_radian(src_grid_units, src_grid->size, src_grid->cell_center_lat, "source grid center lat"); 

  if ( src_grid->num_cell_corners )
    {
      nce(nc_get_var_double(nc_file_id, nc_srcgrdcrnrlat_id, src_grid->cell_corner_lat));
      nce(nc_get_var_double(nc_file_id, nc_srcgrdcrnrlon_id, src_grid->cell_corner_lon));

      nce(nc_get_att_text(nc_file_id, nc_srcgrdcrnrlat_id, "units", src_grid_units));
      nce(nc_inq_attlen(nc_file_id, nc_srcgrdcrnrlat_id, "units", &attlen));
      src_grid_units[attlen] = 0;

      grid_to_radian(src_grid_units, src_grid->num_cell_corners*src_grid->size, src_grid->cell_corner_lon, "source grid corner lon"); 
      grid_to_radian(src_grid_units, src_grid->num_cell_corners*src_grid->size, src_grid->cell_corner_lat, "source grid corner lat"); 
    }

  if ( lgridarea )
    nce(nc_get_var_double(nc_file_id, nc_srcgrdarea_id, src_grid->cell_area));

  nce(nc_get_var_double(nc_file_id, nc_srcgrdfrac_id, src_grid->cell_frac));

  nce(nc_get_var_int(nc_file_id, nc_dstgrddims_id, tgt_grid->dims));

  nce(nc_get_var_int(nc_file_id, nc_dstgrdimask_id, tgt_grid->mask));

  nce(nc_get_var_double(nc_file_id, nc_dstgrdcntrlat_id, tgt_grid->cell_center_lat));
  nce(nc_get_var_double(nc_file_id, nc_dstgrdcntrlon_id, tgt_grid->cell_center_lon));

  nce(nc_get_att_text(nc_file_id, nc_dstgrdcntrlat_id, "units", tgt_grid_units));
  nce(nc_inq_attlen(nc_file_id, nc_dstgrdcntrlat_id, "units", &attlen));
  tgt_grid_units[attlen] = 0;

  grid_to_radian(tgt_grid_units, tgt_grid->size, tgt_grid->cell_center_lon, "target grid center lon"); 
  grid_to_radian(tgt_grid_units, tgt_grid->size, tgt_grid->cell_center_lat, "target grid center lat"); 

  if ( tgt_grid->num_cell_corners )
    {
      nce(nc_get_var_double(nc_file_id, nc_dstgrdcrnrlat_id, tgt_grid->cell_corner_lat));
      nce(nc_get_var_double(nc_file_id, nc_dstgrdcrnrlon_id, tgt_grid->cell_corner_lon));

      nce(nc_get_att_text(nc_file_id, nc_dstgrdcrnrlat_id, "units", tgt_grid_units));
      nce(nc_inq_attlen(nc_file_id, nc_dstgrdcrnrlat_id, "units", &attlen));
      tgt_grid_units[attlen] = 0;
      
      grid_to_radian(tgt_grid_units, tgt_grid->num_cell_corners*tgt_grid->size, tgt_grid->cell_corner_lon, "target grid corner lon"); 
      grid_to_radian(tgt_grid_units, tgt_grid->num_cell_corners*tgt_grid->size, tgt_grid->cell_corner_lat, "target grid corner lat"); 
    }

  if ( lgridarea )
    nce(nc_get_var_double(nc_file_id, nc_dstgrdarea_id, tgt_grid->cell_area));

  nce(nc_get_var_double(nc_file_id, nc_dstgrdfrac_id, tgt_grid->cell_frac));

  nce(nc_get_var_int(nc_file_id, nc_srcadd_id, rv->src_cell_add));
  nce(nc_get_var_int(nc_file_id, nc_dstadd_id, rv->tgt_cell_add));

  for ( i = 0; i < rv->num_links; i++ )
    {
      rv->src_cell_add[i]--;
      rv->tgt_cell_add[i]--;
    }

  nce(nc_get_var_double(nc_file_id, nc_rmpmatrix_id, rv->wts));

  /* Close input file */

  nce(nc_close(nc_file_id));

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  rv->links.option    = FALSE;
  rv->links.max_links = 0;
  rv->links.num_blks  = 0;
  rv->links.num_links = NULL;
  rv->links.src_add   = NULL;
  rv->links.dst_add   = NULL;
  rv->links.w_index   = NULL;
}  /* read_remap_scrip */
