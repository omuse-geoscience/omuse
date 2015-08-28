#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#define H5_USE_16_API

#if defined(HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "griddes.h"
#include "error.h"


#if defined(HAVE_LIBHDF5)
static herr_t
obj_info(hid_t loc_id, const char *name, void *objname)
{
  H5G_obj_t obj_type;
  H5G_stat_t statbuf;
  herr_t lexist = 0;

  H5Gget_objinfo(loc_id, name, FALSE, &statbuf);

  if ( strcmp(name, (char *) objname) == 0 )
    {
      lexist = 1;

      obj_type = statbuf.type;

      switch (obj_type) {
      case H5G_GROUP:
	if ( cdoVerbose ) cdoPrint(" Object with name %s is a group", name);
	break;
      case H5G_DATASET: 
	if ( cdoVerbose ) cdoPrint(" Object with name %s is a dataset", name);
	break;
      case H5G_TYPE: 
	if ( cdoVerbose ) cdoPrint(" Object with name %s is a named datatype", name);
	break;
      default:
	/*cdoAbort(" Unable to identify an object %s", name);*/
	break;
      }
  }

  return lexist;
}
#endif


#if defined(HAVE_LIBHDF5)
static
int h5find_object(hid_t file_id, char *name)
{
  int lexist = 0;

  lexist = (int) H5Giterate(file_id, "/", NULL, obj_info, (void *) name);

  return lexist;
}
#endif

static
void fill_gridvals(int xsize, int ysize, double *xvals, double *yvals)
{
  int i, j, ii, jj;
  int index, index2;
  double xmin, xmax, ymin, ymax;

  xmin = -180;
  xmax =  180;
  ymin = -90;
  ymax =  90;

  for ( ii = 0; ii < xsize/2; ++ii )
    {
      index2 = ysize/2*xsize + ii;
      if ( xvals[index2] > -180 && xvals[index2] < 360 )
	{
	  xmin = xvals[index2];
	  break;
	}
    }

  for ( ii = xsize-1; ii > xsize/2; --ii )
    {
      index2 = ysize/2*xsize + ii;
      if ( xvals[index2] > -180 && xvals[index2] < 360 )
	{
	  xmax = xvals[index2];
	  break;
	}
    }
  /*
  for ( jj = 0; jj < ysize; ++jj )
    {
      index2 = jj*xsize + xsize/2;
      if ( xvals[index2] < -180 || xvals[index2] > 360 ) xvals[index2] = 0;
      index2 = jj*xsize + xsize/2-1;
      if ( xvals[index2] < -180 || xvals[index2] > 360 ) xvals[index2] = 0;
    }
  */
  for ( jj = 0; jj < ysize/2; ++jj )
    {
      index2 = jj*xsize + xsize/2;
      if ( yvals[index2] > -90 && yvals[index2] < 90 )
	{
	  ymax = yvals[index2];
	  break;
	}
    }

  for ( jj = ysize-1; jj > ysize/2; --jj )
    {
      index2 = jj*xsize + xsize/2;
      if ( yvals[index2] > -90 && yvals[index2] < 90 )
	{
	  ymin = yvals[index2];
	  break;
	}
    }

  /* printf("xmin %g, xmax %g, ymin %g, ymax %g\n", xmin, xmax, ymin, ymax); */
  
  for ( i = 0; i < xsize*ysize; ++i )
    {
      if ( xvals[i] > -180 && xvals[i] < 360 )
	{
	  if ( xvals[i] < xmin ) xmin = xvals[i];
	  if ( xvals[i] > xmax ) xmax = xvals[i];
	}

      if ( yvals[i] > -90 && yvals[i] < 90 )
	{
	  if ( yvals[i] < ymin ) ymin = yvals[i];
	  if ( yvals[i] > ymax ) ymax = yvals[i];
	}
    }

  for ( j = 0; j < ysize; ++j )
    for ( i = 0; i < xsize; ++i )
      {
	index = j*xsize + i;

	if ( xvals[index] < -180 || xvals[index] > 360 )
	  {
	    if ( i < xsize/2 )
	      xvals[index] = xmin;
	    else
	      xvals[index] = xmax;
	    /*
	    if ( j < ysize/2 )
	      for ( jj = j+1; jj < ysize/2; ++jj )
		{
		  index2 = jj*xsize + i;
		  if ( xvals[index2] > -180 && xvals[index2] < 360 )
		    {
		      xvals[index] = xvals[index2];
		      break;
		    }
		}
	    else
	      for ( jj = j-1; jj > ysize/2; --jj )
		{
		  index2 = jj*xsize + i;
		  if ( xvals[index2] > -180 && xvals[index2] < 360 )
		    {
		      xvals[index] = xvals[index2];
		      break;
		    }
		}
	    */
	    /*
	    if ( i < xsize/2 )
	      {
		xvals[index] = xmin;
		for ( ii = i+1; ii < xsize/2; ++ii )
		  {
		    index2 = j*xsize + ii;
		    if ( xvals[index2] > -180 && xvals[index2] < 360 )
		      {
			xvals[index] = (xmin*(ii-i) + xvals[index2]*(i))/ii;
			break;
		      }
		  }
	      }
	    else
	      {
		for ( ii = i-1; ii >= xsize/2; --ii )
		  {
		    index2 = j*xsize + ii;
		    if ( xvals[index2] > -180 && xvals[index2] < 360 )
		      {
			xvals[index] = (xmax*(i-ii) + xvals[index2]*((xsize-1)-i))/(xsize-1-ii);
			break;
		      }
		  }
	      }
	    */
	  }

	if ( yvals[index] < -90 || yvals[index] > 90 )
	  {
	    if ( j < ysize/2 )
	      yvals[index] = ymax;
	    else
	      yvals[index] = ymin;

	    if ( i < xsize/2 )
	      for ( ii = i+1; ii < xsize/2; ++ii )
		{
		  index2 = j*xsize + ii;
		  if ( yvals[index2] > -90 && yvals[index2] < 90 )
		    {
		      yvals[index] = yvals[index2];
		      break;
		    }
		}
	    else
	      for ( ii = i-1; ii > xsize/2; --ii )
		{
		  index2 = j*xsize + ii;
		  if ( yvals[index2] > -90 && yvals[index2] < 90 )
		    {
		      yvals[index] = yvals[index2];
		      break;
		    }
		}
	  }
      }
}


void correct_sinxvals(int xsize, int ysize, double *xvals)
{
  long i, j, istart, index;
  double xmin, xmax;

  xmin = -180;
  xmax =  180;

  for ( j = 0; j < ysize; ++j )
    {
      istart = xsize/2-1;
      xmin = xvals[j*xsize+istart];
      for ( i = istart-1; i >= 0; i-- )
	{
	  index = j*xsize+i;
	  if ( xvals[index] > xmin ) break;
	  xmin = xvals[index];
	}

      if ( i >= 0 )
	{
	  istart = i;
	  // printf("%d %d %g\n",j,i, xmin);
	  for ( i = 0; i <= istart; ++i )
	    {
	      index = j*xsize+i;
	      xvals[index] = xmin;
	    }
	}

      istart = xsize/2;
      xmax = xvals[j*xsize+istart];
      for ( i = istart+1; i < xsize; i++ )
	{
	  index = j*xsize+i;
	  if ( xvals[index] < xmax ) break;
	  xmax = xvals[index];
	}

      if ( i < xsize )
	{
	  istart = i;
	  // printf("%d %d %g\n",j,i, xmax);
	  for ( i = istart; i < xsize; ++i )
	    {
	      index = j*xsize+i;
	      xvals[index] = xmax;
	    }
	}
    }
}


int gridFromH5file(const char *gridfile)
{
  int       gridID = -1;
#if defined(HAVE_LIBHDF5)
  hid_t     fapl_id = H5P_DEFAULT;
  hid_t	    file_id;	    /* HDF5 File ID	        	*/
  hid_t	    lon_id = -1;    /* Dataset ID	        	*/
  hid_t	    lat_id = -1;    /* Dataset ID	        	*/
  hid_t     att_id;
  hid_t     dataspace;   
  hsize_t   dims_out[9];    /* dataset dimensions               */
  herr_t    status;	    /* Generic return value		*/
  int       rank;
  griddes_t grid;


  gridInit(&grid);

  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);

  /* Open an existing file. */
  file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, fapl_id);

  H5Pclose(fapl_id);

  if ( file_id < 0 ) return(gridID);

  if ( h5find_object(file_id, "lon") > 0 && 
       h5find_object(file_id, "lat") > 0 )
    {
      lon_id = H5Dopen(file_id, "/lon");
      lat_id = H5Dopen(file_id, "/lat");
    }
  else if ( h5find_object(file_id, "Longitude") > 0 && 
	    h5find_object(file_id, "Latitude") > 0 )
    {
      lon_id = H5Dopen(file_id, "/Longitude");
      lat_id = H5Dopen(file_id, "/Latitude");
    }
  
  if ( lon_id >= 0 && lat_id >= 0 )
    {
      hid_t type_id;
      hid_t native_type;
      int ftype = 0;

      dataspace = H5Dget_space(lon_id);    /* dataspace handle */
      rank      = H5Sget_simple_extent_ndims(dataspace);
      status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

      if ( rank != 2 )
	{
	  //if ( cdoVerbose ) cdoWarning("Unexpected rank = %d!", rank);
	  goto RETURN;
	}

      att_id = H5Aopen_name(lon_id, "bounds");
      if ( att_id >= 0 )
	{
	  H5Aclose(att_id);
	  goto RETURN;
	}

      att_id = H5Aopen_name(lat_id, "bounds");
      if ( att_id >= 0 )
	{
	  H5Aclose(att_id);
	  goto RETURN;
	}

      /*
      printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
	     (unsigned long)(dims_out[1]), (unsigned long)(dims_out[0]));
      */

      type_id = H5Dget_type(lon_id);  /* get datatype*/

      native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
      if      ( H5Tequal(native_type, H5T_NATIVE_SCHAR)  > 0 ) {ftype=0;}
      else if ( H5Tequal(native_type, H5T_NATIVE_UCHAR)  > 0 ) {ftype=0;}
      else if ( H5Tequal(native_type, H5T_NATIVE_SHORT)  > 0 ) {ftype=0;}
      else if ( H5Tequal(native_type, H5T_NATIVE_USHORT) > 0 ) {ftype=0;}
      else if ( H5Tequal(native_type, H5T_NATIVE_INT)    > 0 ) {ftype=0;}
      else if ( H5Tequal(native_type, H5T_NATIVE_UINT)   > 0 ) {ftype=0;}
      else if ( H5Tequal(native_type, H5T_NATIVE_FLOAT)  > 0 ) {ftype=1;}
      else if ( H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0 ) {ftype=1;}
      else
	{
	  cdoWarning("Grid has unsupported native datatype!");
	  goto RETURN;
	}
      H5Tclose(native_type);

      grid.xsize = (int)dims_out[1];
      grid.ysize = (int)dims_out[0];
      grid.size  = grid.xsize*grid.ysize;

      grid.xvals = (double*) malloc(grid.size*sizeof(double));
      grid.yvals = (double*) malloc(grid.size*sizeof(double));

      if ( ftype )
	{
	  status = H5Dread(lon_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.xvals);
	  status = H5Dread(lat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.yvals);
	}
      else
	{
	  int *iarray, i;
	  iarray = (int*) malloc(grid.size*sizeof(int));
	  status = H5Dread(lon_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray);
	  for ( i = 0; i < grid.size; ++i ) grid.xvals[i] = iarray[i];
	  status = H5Dread(lat_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray);
	  for ( i = 0; i < grid.size; ++i ) grid.yvals[i] = iarray[i];
	  free(iarray);
	}

      status = H5Sclose(dataspace);

      /* Close the dataset. */
      status = H5Dclose(lon_id);
      status = H5Dclose(lat_id);

      fill_gridvals(grid.xsize, grid.ysize, grid.xvals, grid.yvals);

      grid.type = GRID_CURVILINEAR;
      grid.prec = DATATYPE_FLT32;

      gridID = gridDefine(grid);
    }
  else  if ( h5find_object(file_id, "where") > 0 )
    {
      double xscale = 1, yscale = 1;
      double xoffset = 0, yoffset = 0;
      hid_t grp_id;
      int i;

      grp_id = H5Gopen(file_id, "/where/lon/what");
      if ( grp_id >= 0 )
	{
	  att_id = H5Aopen_name(grp_id, "gain");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xscale);
	      H5Aclose(att_id);
	    }

	  att_id = H5Aopen_name(grp_id, "offset");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xoffset);
	      H5Aclose(att_id);
	    }
	  
	  H5Gclose(grp_id);
	}

      grp_id = H5Gopen(file_id, "/where/lat/what");
      if ( grp_id >= 0 )
	{
	  att_id = H5Aopen_name(grp_id, "gain");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yscale);
	      H5Aclose(att_id);
	    }

	  att_id = H5Aopen_name(grp_id, "offset");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yoffset);
	      H5Aclose(att_id);
	    }
	  
	  H5Gclose(grp_id);
	}

      /* Open an existing dataset. */
      lon_id = H5Dopen(file_id, "/where/lon/data");
      if ( lon_id >= 0 )
	lat_id = H5Dopen(file_id, "/where/lat/data");

      if ( lon_id >= 0 && lat_id >= 0 )
	{
	  hid_t type_id;
	  hid_t native_type;
	  int ftype = 0;

	  dataspace = H5Dget_space(lon_id);    /* dataspace handle */
	  rank      = H5Sget_simple_extent_ndims(dataspace);
	  status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

	  if ( rank != 2 )
	    {
	      //if ( cdoVerbose ) cdoWarning("Unexpected rank = %d!", rank);
	      goto RETURN;
	    }
	  /*
	  printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
		 (unsigned long)(dims_out[1]), (unsigned long)(dims_out[0]));
	  */

	  type_id = H5Dget_type(lon_id);  /* get datatype*/

	  native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
	  if      ( H5Tequal(native_type, H5T_NATIVE_SCHAR)  > 0 ) {ftype=0;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_UCHAR)  > 0 ) {ftype=0;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_SHORT)  > 0 ) {ftype=0;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_USHORT) > 0 ) {ftype=0;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_INT)    > 0 ) {ftype=0;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_UINT)   > 0 ) {ftype=0;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_FLOAT)  > 0 ) {ftype=1;}
	  else if ( H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0 ) {ftype=1;}
	  else
	    {
	      cdoWarning("Grid has unsupported native datatype!");
	      goto RETURN;
	    }
	  H5Tclose(native_type);

	  grid.xsize = (int)dims_out[1];
	  grid.ysize = (int)dims_out[0];
	  grid.size  = grid.xsize*grid.ysize;

	  grid.xvals = (double*) malloc(grid.size*sizeof(double));
	  grid.yvals = (double*) malloc(grid.size*sizeof(double));

	  if ( ftype )
	    {
	      status = H5Dread(lon_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.xvals);
	      status = H5Dread(lat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.yvals);
	    }
	  else
	    {
	      int *iarray, i;
	      iarray = (int*) malloc(grid.size*sizeof(int));
	      status = H5Dread(lon_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray);
	      for ( i = 0; i < grid.size; ++i ) grid.xvals[i] = iarray[i];
	      status = H5Dread(lat_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray);
	      for ( i = 0; i < grid.size; ++i ) grid.yvals[i] = iarray[i];
	      free(iarray);
	    }

	  status = H5Sclose(dataspace);
	  
	  /* Close the dataset. */
	  status = H5Dclose(lon_id);
	  status = H5Dclose(lat_id);

	  for ( i = 0; i < grid.size; ++i ) grid.xvals[i] = grid.xvals[i]*xscale + xoffset;
	  for ( i = 0; i < grid.size; ++i ) grid.yvals[i] = grid.yvals[i]*yscale + yoffset;

	  grid.type = GRID_CURVILINEAR;
	  grid.prec = DATATYPE_FLT32;

	  gridID = gridDefine(grid);
	}
    }

 RETURN:

  /* Close file */
  if ( file_id >= 0 )  status = H5Fclose(file_id);

#else
  cdoWarning("HDF5 support not compiled in!");
#endif

  return (gridID);
}
