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

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <time.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#if defined(HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

static int nvars_ml = 4;
static const char strfiletype_ml[]  = "Initial file spectral";
static const char strfiletype_res[] = "Restart history file";


typedef struct {
  int  gridtype;
  int  zaxistype;
  int  code;
  char *name;
  char *longname;
  char *units;
  int gridID;
  int zaxisID;
  int gridsize;
  int nlev;
  double *ptr;
} VAR;


typedef struct {
  int    naint;
  int    naflt;
  int    natxt;
  char   *aintname[1024];
  int    *aintentry[1024];
  char   *afltname[1024];
  double *afltentry[1024];
  char   *atxtname[1024];
  char   *atxtentry[1024];
} ATTS;

static
void iniatts(ATTS *atts)
{
  atts->naint = 0;
  atts->naflt = 0;
  atts->natxt = 0;
}

static
void inivar(VAR *var, int gridtype, int zaxistype, int code, const char *name,
	       const char *longname, const char *units)
{  
  var->gridtype  = gridtype;
  var->zaxistype = zaxistype;
  var->code      = code;
  var->name      = NULL;
  var->longname  = NULL;
  var->units     = NULL;
  if ( name )     var->name      = strdup(name);
  if ( longname ) var->longname  = strdup(longname);
  if ( units )    var->units     = strdup(units);
}

static
void inivars_ml(VAR **vars)
{
  *vars = (VAR*) malloc((nvars_ml+1)*sizeof(VAR));

  inivar(&(*vars)[0], GRID_GAUSSIAN, ZAXIS_HYBRID,  133, "Q",   "specific humidity", "kg/kg");
  inivar(&(*vars)[1], GRID_SPECTRAL, ZAXIS_HYBRID,  138, "SVO", "vorticity", "1/s");
  inivar(&(*vars)[2], GRID_SPECTRAL, ZAXIS_HYBRID,  155, "SD",  "divergence", "1/s");
  inivar(&(*vars)[3], GRID_SPECTRAL, ZAXIS_HYBRID,  130, "STP", "temperature", "K");
  /* Don't change the order (lsp must be the last one)! */
  inivar(&(*vars)[4], GRID_SPECTRAL, ZAXIS_SURFACE, 152, "LSP", "log surface pressure", "");
}

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

static
int import_e5ml(const char *filename, VAR **vars)
{
  int nvars = 0;
#if defined(HAVE_LIBNETCDF)
  int nc_dim_id, nc_var_id;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nlevp1, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  int nc_file_id;
  char filetype[256];
  size_t attlen;
  double *xvals, *yvals, *vct, *levs;

  /* open file and check file type */
  /* nce(nc_open(filename, NC_NOWRITE, &nc_file_id)); */
  nc_file_id = cdf_openread(filename);

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( strcmp(filetype, strfiletype_ml) != 0 ) return (0);

  inivars_ml(vars);

  /* read dimensions */

  nce(nc_inq_dimid(nc_file_id, "lon", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nlon = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "lat", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nlat = (int) dimlen;

  gridIDgp = gridCreate(GRID_GAUSSIAN, nlon*nlat);
  gridDefXsize(gridIDgp, nlon);
  gridDefYsize(gridIDgp, nlat);

  nce(nc_inq_dimid(nc_file_id, "nsp", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nsp = (int) dimlen;

  gridIDsp = gridCreate(GRID_SPECTRAL, nsp*2);
  gridDefComplexPacking(gridIDsp, 1);

  nce(nc_inq_dimid(nc_file_id, "nlev", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nlev = (int) dimlen;
  nlevp1 = nlev + 1;
  nvct = nlevp1*2;

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisIDml  = zaxisCreate(ZAXIS_HYBRID, nlev);

  levs = (double*) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = i+1;
  zaxisDefLevels(zaxisIDml, levs);
  free(levs);

  /* read variables */

  xvals = (double*) malloc(nlon*sizeof(double));
  yvals = (double*) malloc(nlat*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "lon", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, xvals));

  nce(nc_inq_varid(nc_file_id, "lat", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, yvals));

  gridDefXvals(gridIDgp, xvals);
  gridDefYvals(gridIDgp, yvals);

  free(xvals);
  free(yvals);

  vct   = (double*) malloc(nvct*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "vct_a", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_inq_varid(nc_file_id, "vct_b", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  zaxisDefVct(zaxisIDml, 2*nlevp1, vct);
  free(vct);

  for ( iv = 0; iv < nvars_ml; iv++ )
    {
      nvals = 0;

      gridtype  = (*vars)[iv].gridtype;
      zaxistype = (*vars)[iv].zaxistype;

      UNUSED(zaxistype);

      if ( gridtype == GRID_GAUSSIAN )
	{
	  (*vars)[iv].gridID = gridIDgp;
	  nvals += nlon*nlat;
	}
      else
	{
	  (*vars)[iv].gridID = gridIDsp;
	  nvals += nsp*2;
	}

      (*vars)[iv].zaxisID   = zaxisIDml;
      (*vars)[iv].gridsize  = nvals;
      (*vars)[iv].nlev      = nlev;

      (*vars)[iv].ptr = (double*) malloc(nlev*nvals*sizeof(double));
      
      for ( i = 0; i < nlev; i++ )
	{
	  if ( gridtype == GRID_GAUSSIAN )
	    {
	      start[0] = 0;     start[1] = i;  start[2] = 0;
	      count[0] = nlat;  count[1] = 1;  count[2] = nlon;     
	    }
	  else
	    {
	      start[0] = 0;     start[1] = 0;  start[2] = i;
	      count[0] = nsp;   count[1] = 2;  count[2] = 1;
	    }

	  nce(nc_inq_varid(nc_file_id, (*vars)[iv].name, &nc_var_id));
	  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[iv].ptr+i*nvals));
	}
    }

  /* read lsp */

  (*vars)[nvars_ml].gridID    = gridIDsp;
  (*vars)[nvars_ml].zaxisID   = zaxisIDsfc;
  (*vars)[nvars_ml].gridsize  = nsp*2;
  (*vars)[nvars_ml].nlev      = 1;

  start[0] = 0;    start[1] = 0;  start[2] = nlev;
  count[0] = nsp;  count[1] = 2;  count[2] = 1;

  (*vars)[nvars_ml].ptr = (double*) malloc(nsp*2*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "STP", &nc_var_id));
  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[nvars_ml].ptr));

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = nvars_ml + 1;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}

static
void export_e5ml(const char *filename, VAR *vars, int nvars, int vdate, int vtime, int ntr)
{
#if defined(HAVE_LIBNETCDF)
  int nc_var_id;
  size_t nvals;
  size_t start[3], count[3];
  int dimidsp[9];
  int varid, code;
  int ilev;
  int lon, lat;
  int nlon, nlat, nlev, nlevp1/*, nvct*/, nsp, n2, i, nvclev;
  int lat_dimid, lon_dimid, nlev_dimid, nlevp1_dimid, nsp_dimid, nvclev_dimid, n2_dimid;
  int gridIDgp = -1, gridIDsp, zaxisIDml = -1;
  int gridtype, zaxistype;
  int nc_file_id;
  int nc_stpid, lspid;
  double *xvals, *yvals;
  const double *vct;
  char atttext[1024];
  size_t attlen;
  //int attint;
  char *username;
  char timestr[30];
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  int writemode = NC_CLOBBER;
  unsigned long  data_size;
  
  date_and_time_in_sec = time(NULL);
  timestr[0] = 0;

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%d/%m/%Y %H:%M", date_and_time);
    }

  username = getenv("LOGNAME");
  if ( username == NULL )
    {
      username = getenv("USER");
      if ( username == NULL ) username = "unknown";
    }


  n2 = 2;

  lon = 0; lat = 0; nsp = 0; nlev = 0; nlevp1 = 0; nvclev = 0;
  for ( varid = 0; varid < nvars; ++varid )
    {
      gridtype  = vars[varid].gridtype;
      zaxistype = vars[varid].zaxistype;

      if ( gridtype == GRID_GAUSSIAN && lat == 0 )
	{
	  gridIDgp = vars[varid].gridID;
	  lon = gridInqXsize(vars[varid].gridID);
	  lat = gridInqYsize(vars[varid].gridID);
	}
      else if ( gridtype == GRID_SPECTRAL && nsp == 0 )
	{
	  gridIDsp = vars[varid].gridID;
	  UNUSED(gridIDsp);
	  nsp = gridInqSize(vars[varid].gridID);
	  nsp = nsp/2;
	}

      if ( zaxistype == ZAXIS_HYBRID && nlev == 0 )
	{
	  zaxisIDml = vars[varid].zaxisID;
	  nlev = zaxisInqSize(vars[varid].zaxisID);
	  nlevp1 = nlev + 1;
	  nvclev = nlev + 1;
	}
    }

  if ( lat  == 0 ) cdoAbort("Gaussian grid not found!");
  if ( nsp  == 0 ) cdoAbort("Spectral data not found!");
  if ( nlev == 0 ) cdoAbort("Hybrid level not found!");

  nlon = lon;
  nlat = lat;


  data_size = nlon+nlat + 2*nvclev + 2*nsp*2*nlev + nsp*2*nlevp1 + nlon*nlat*nlev;

  if ( data_size*8 > 2147000000 )
    {
#if defined(NC_64BIT_OFFSET)
      writemode = NC_CLOBBER | NC_64BIT_OFFSET;
#else
      cdoWarning("Datasize > 2GB and NC_64BIT_OFFSET not available!");
#endif
    }

  /* create file */
  nce(nc_create(filename, writemode, &nc_file_id));

  strcpy(atttext, "IEEE");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "source_type", attlen, atttext));

  strcpy(atttext, commandLine());
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "history", attlen, atttext));

  strcpy(atttext, username);
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "user", attlen, atttext));

  strcpy(atttext, timestr);
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "created", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_1", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_2", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_3", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_4", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_5", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_6", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_7", attlen, atttext));

  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "label_8", attlen, atttext));

  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "fdate", NC_INT, 1, &vdate));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "ftime", NC_INT, 1, &vtime));

  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "vdate", NC_INT, 1, &vdate));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "vtime", NC_INT, 1, &vtime));

  //attint = 31;
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "spherical_truncation_n", NC_INT, 1, &ntr));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "spherical_truncation_m", NC_INT, 1, &ntr));
  nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "spherical_truncation_k", NC_INT, 1, &ntr));

  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "file_type", strlen(strfiletype_ml), strfiletype_ml));

  
  nce(nc_def_dim(nc_file_id, "lat", lat, &lat_dimid));

  nce(nc_def_dim(nc_file_id, "lon", lon, &lon_dimid));

  nce(nc_def_dim(nc_file_id, "nlev", nlev, &nlev_dimid));
  nce(nc_def_dim(nc_file_id, "nlevp1", nlevp1, &nlevp1_dimid));

  nce(nc_def_dim(nc_file_id, "nsp", nsp, &nsp_dimid));

  nce(nc_def_dim(nc_file_id, "nvclev", nvclev, &nvclev_dimid));

  nce(nc_def_dim(nc_file_id, "n2", n2, &n2_dimid));

  nce(nc_enddef(nc_file_id));

  /* define gaussian grid */

  xvals = (double*) malloc(nlon*sizeof(double));
  yvals = (double*) malloc(nlat*sizeof(double));

  gridInqXvals(gridIDgp, xvals);
  gridInqYvals(gridIDgp, yvals);

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "lat", NC_DOUBLE, 1, &lat_dimid, &nc_var_id));
  strcpy(atttext, "Gaussian latitude");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  strcpy(atttext, "degrees_N");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, yvals));
   
  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "lon", NC_DOUBLE, 1, &lon_dimid, &nc_var_id));
  strcpy(atttext, "longitude");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  strcpy(atttext, "degrees_E");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, xvals));

  free(xvals);
  free(yvals);

  /* define model level */

  // nvct = nvclev*2;

  /* vct   = (double*) malloc(nvct*sizeof(double)); */

  vct = zaxisInqVctPtr(zaxisIDml);

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "vct_a", NC_DOUBLE, 1, &nvclev_dimid, &nc_var_id));
  strcpy(atttext, "vertical-coordinate parameter set A");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "vct_b", NC_DOUBLE, 1, &nvclev_dimid, &nc_var_id));
  strcpy(atttext, "vertical-coordinate parameter set B");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", attlen, atttext));
  strcpy(atttext, "");
  attlen = strlen(atttext);
  nce(nc_put_att_text(nc_file_id, nc_var_id, "units", attlen, atttext));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  /* free(vct); */

  lspid = -1;
  nc_stpid = -1;

  for ( varid = 0; varid < nvars; varid++ )
    {
      nvals = 0;

      code      = vars[varid].code;
      gridtype  = vars[varid].gridtype;
      zaxistype = vars[varid].zaxistype;

      ilev = zaxisInqSize(vars[varid].zaxisID);

      if ( ilev == 1 )
	{
	  if ( code == 152 )
	    {
	      lspid = varid;
	      if ( gridtype != GRID_SPECTRAL )
		cdoAbort("%s has wrong gridtype!", vars[varid].name);
	    }
	  continue;
	}

      if ( nlev != ilev )
	cdoAbort("Unexpected number of level %d!", ilev);

      if ( gridtype == GRID_GAUSSIAN )
	{
	  nvals = nlon*nlat;

	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = nlev_dimid;
	  dimidsp[2] = lon_dimid;
	}
      else if ( gridtype == GRID_SPECTRAL )
	{
	  nvals = nsp*2;

	  dimidsp[0] = nsp_dimid;
	  dimidsp[1] = n2_dimid;

	  if ( strcmp(vars[varid].name, "STP") == 0 || strcmp(vars[varid].name, "T") == 0 )
	    dimidsp[2] = nlevp1_dimid;
	  else
	    dimidsp[2] = nlev_dimid;
	}
      else
	cdoAbort("Unsupported grid!");

      nce(nc_redef(nc_file_id));
      nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 3, dimidsp, &nc_var_id));
      if ( vars[varid].longname && *vars[varid].longname)
	nce(nc_put_att_text(nc_file_id, nc_var_id, "long_name", strlen(vars[varid].longname), vars[varid].longname));
      if ( vars[varid].units && *vars[varid].units)
	nce(nc_put_att_text(nc_file_id, nc_var_id, "units", strlen(vars[varid].units), vars[varid].units));
      nce(nc_enddef(nc_file_id));

      if (  dimidsp[2] == nlevp1_dimid ) nc_stpid = nc_var_id;

      for ( i = 0; i < nlev; i++ )
	{
	  if ( gridtype == GRID_GAUSSIAN )
	    {
	      start[0] = 0;     start[1] = i;  start[2] = 0;
	      count[0] = nlat;  count[1] = 1;  count[2] = nlon;     
	    }
	  else
	    {
	      start[0] = 0;     start[1] = 0;  start[2] = i;
	      count[0] = nsp;   count[1] = 2;  count[2] = 1;
	    }
	  
	  nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	}
    }

  if ( lspid == -1 ) cdoAbort("LSP not found!");
  if ( nc_stpid == -1 ) cdoAbort("STP not found!");

  /* write lsp */
  start[0] = 0;    start[1] = 0;  start[2] = nlev;
  count[0] = nsp;  count[1] = 2;  count[2] = 1;

  nce(nc_put_vara_double(nc_file_id, nc_stpid, start, count, vars[lspid].ptr));

  /*close input file */
  nce(nc_close(nc_file_id));

#else
  cdoAbort("netCDF support not compiled in!");
#endif
}

#if defined(HAVE_LIBNETCDF)
static
void read_gg3d(int nc_file_id, const char *name, VAR *var, int gridID, int zaxisID)
{
  int nlev, nlat, nlon, gridsize, i;
  int gridtype, zaxistype;
  int nc_var_id;
  size_t start[3], count[3];

  gridtype = gridInqType(gridID);
  zaxistype = zaxisInqType(zaxisID);

  inivar(var, gridtype, zaxistype,  0, name, NULL, NULL);

  nlon = gridInqXsize(gridID);
  nlat = gridInqYsize(gridID);
  nlev = zaxisInqSize(zaxisID);

  gridsize = nlon*nlat;

  var->gridID    = gridID;
  var->zaxisID   = zaxisID;
  var->gridsize  = gridsize;
  var->nlev      = nlev;

  var->ptr = (double*) malloc(nlev*gridsize*sizeof(double));
  
  for ( i = 0; i < nlev; i++ )
    {
      start[0] = 0;     start[1] = i;  start[2] = 0;
      count[0] = nlat;  count[1] = 1;  count[2] = nlon;     

      nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
      nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, var->ptr+i*gridsize));
    }
}
#endif

#if defined(HAVE_LIBNETCDF)
static
void read_fc4d(int nc_file_id, const char *name, VAR *var, int gridID, int zaxisID, int nhgl, int nmp1)
{
  int nlev, nfc, i;
  int gridtype, zaxistype;
  int nc_var_id;
  size_t start[4], count[4];

  gridtype = gridInqType(gridID);
  zaxistype = zaxisInqType(zaxisID);

  inivar(var, gridtype, zaxistype,  0, name, NULL, NULL);

  nfc  = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);

  if ( nfc != nhgl*nmp1*2 ) cdoAbort("%s: inconsistent dimension length!", __func__);

  var->gridID    = gridID;
  var->zaxisID   = zaxisID;
  var->gridsize  = nfc;
  var->nlev      = nlev;

  var->ptr = (double*) malloc(nlev*nfc*sizeof(double));
  
  for ( i = 0; i < nlev; i++ )
    {
      start[0] = 0;     start[1] = 0;    start[2] = 0;  start[3] = i;
      count[0] = nhgl;  count[1] = nmp1; count[2] = 2;  count[3] = 1;     

      nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
      nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, var->ptr+i*nfc));
    }
}
#endif

static
int import_e5res(const char *filename, VAR **vars, ATTS *atts)
{
  int nvars = 0;
#if defined(HAVE_LIBNETCDF)
  int nc_var_id;
  int varid;
  size_t nvals;
  size_t start[3], count[3];
  int nlon, nlat, nlev, nvct, nfc, i;
  int gridIDgp, gridIDfc, gridIDhgl, zaxisIDml, zaxisIDmlh, zaxisIDsfc, zaxisIDbsfc, zaxisIDn2;
  int nc_file_id;
  char filetype[256];
  double *xvals, *yvals, *vct, *levs;
  int ncvarid;
  int ndims, ngatts, unlimdimid;
  int nvdims, nvatts;
  int dimidsp[9];
  int max_vars = 4096;
  char name[CDI_MAX_NAME];
  int lon_dimid, lat_dimid, nhgl_dimid, nlevp1_dimid, spc_dimid, nvclev_dimid;
  int complex_dimid, nmp1_dimid, belowsurface_dimid, lev_dimid, ilev_dimid;
  int /* surface_dimid, height2m_dimid, height10m_dimid,*/ n2_dimid;
  int lon, lat, nhgl, nlevp1, spc, nvclev;
  int complex, nmp1, belowsurface, lev, ilev;
  int /* surface, height2m, height10m,*/ n2;
  int iatt;
  int attint;
  double attflt;
  size_t attlen, dimlen;
  nc_type xtype;
  char attname[CDI_MAX_NAME];
  const int attstringlen = 8192; char attstring[8192];

  /* open file and check file type */
  /* nce(nc_open(filename, NC_NOWRITE, &nc_file_id)); */
  nc_file_id = cdf_openread(filename);

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( memcmp(filetype, strfiletype_res, strlen(strfiletype_res)) != 0 ) return (0);

  /* printf("%s\n", filetype); */

  nce(nc_inq_dimid(nc_file_id, "lon", &lon_dimid));
  nce(nc_inq_dimlen(nc_file_id, lon_dimid, &dimlen));
  lon = dimlen;

  nce(nc_inq_dimid(nc_file_id, "lat", &lat_dimid));
  nce(nc_inq_dimlen(nc_file_id, lat_dimid, &dimlen));
  lat = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "nhgl", &nhgl_dimid));
  nce(nc_inq_dimlen(nc_file_id, nhgl_dimid, &dimlen));
  nhgl = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "nlevp1", &nlevp1_dimid));
  nce(nc_inq_dimlen(nc_file_id, nlevp1_dimid, &dimlen));
  nlevp1 = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "spc", &spc_dimid));
  nce(nc_inq_dimlen(nc_file_id, spc_dimid, &dimlen));
  spc = (int) dimlen;

  UNUSED(spc);

  nce(nc_inq_dimid(nc_file_id, "nvclev", &nvclev_dimid));
  nce(nc_inq_dimlen(nc_file_id, nvclev_dimid, &dimlen));
  nvclev = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "complex", &complex_dimid));
  nce(nc_inq_dimlen(nc_file_id, complex_dimid, &dimlen));
  complex = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "nmp1", &nmp1_dimid));
  nce(nc_inq_dimlen(nc_file_id, nmp1_dimid, &dimlen));
  nmp1 = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "belowsurface", &belowsurface_dimid));
  nce(nc_inq_dimlen(nc_file_id, belowsurface_dimid, &dimlen));
  belowsurface = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "lev", &lev_dimid));
  nce(nc_inq_dimlen(nc_file_id, lev_dimid, &dimlen));
  lev = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "ilev", &ilev_dimid));
  nce(nc_inq_dimlen(nc_file_id, ilev_dimid, &dimlen));
  ilev = (int) dimlen;

  UNUSED(ilev);
  /*
  nce(nc_inq_dimid(nc_file_id, "surface", &surface_dimid));
  nce(nc_inq_dimlen(nc_file_id, surface_dimid, &surface));

  nce(nc_inq_dimid(nc_file_id, "height2m", &height2m_dimid));
  nce(nc_inq_dimlen(nc_file_id, height2m_dimid, &height2m));

  nce(nc_inq_dimid(nc_file_id, "height10m", &height10m_dimid));
  nce(nc_inq_dimlen(nc_file_id, height10m_dimid, &height10m));
  */
  nce(nc_inq_dimid(nc_file_id, "n2", &n2_dimid));
  nce(nc_inq_dimlen(nc_file_id, n2_dimid, &dimlen));
  n2 = (int) dimlen;

  /* define gaussian grid */

  nlon = lon;
  nlat = lat;

  gridIDgp = gridCreate(GRID_GAUSSIAN, nlon*nlat);
  gridDefXsize(gridIDgp, nlon);
  gridDefYsize(gridIDgp, nlat);

  xvals = (double*) malloc(nlon*sizeof(double));
  yvals = (double*) malloc(nlat*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "lon", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, xvals));

  nce(nc_inq_varid(nc_file_id, "lat", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, yvals));

  gridDefXvals(gridIDgp, xvals);
  gridDefYvals(gridIDgp, yvals);

  free(xvals);
  free(yvals);

  /* define fourier grid */

  nfc = nhgl*nmp1*complex;

  gridIDfc = gridCreate(GRID_FOURIER, nfc);

  /* define nhgl grid */

  gridIDhgl = gridCreate(GRID_GENERIC, nhgl);


  /* define surface level */

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);

  /* define below surface level */

  nlev = belowsurface;
  zaxisIDbsfc = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, nlev);

  levs = (double*) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = 0;
  zaxisDefLevels(zaxisIDbsfc, levs);
  free(levs);


  /* define n2 level */

  nlev = n2;
  zaxisIDn2 = zaxisCreate(ZAXIS_GENERIC, nlev);

  levs = (double*) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = 0;
  zaxisDefLevels(zaxisIDn2, levs);
  free(levs);

  /* define model level */

  nlev = lev;
  nvct = nvclev*2;

  vct   = (double*) malloc(nvct*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "vct_a", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_inq_varid(nc_file_id, "vct_b", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  /* ZAXIS_HYBRID */ 

  zaxisIDml  = zaxisCreate(ZAXIS_HYBRID, nlev);

  levs = (double*) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = i+1;
  zaxisDefLevels(zaxisIDml, levs);
  free(levs);

  zaxisDefVct(zaxisIDml, 2*nlevp1, vct);

  /* ZAXIS_HYBRID_HALF */ 

  zaxisIDmlh  = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1);

  levs = (double*) malloc(nlevp1*sizeof(double));
  for ( i = 0; i < nlevp1; i++ ) levs[i] = i;
  zaxisDefLevels(zaxisIDmlh, levs);
  free(levs);

  zaxisDefVct(zaxisIDmlh, 2*nlevp1, vct);

  free(vct);


  nce(nc_inq(nc_file_id, &ndims, &nvars, &ngatts, &unlimdimid));

  /* read global attributtes*/
  for ( iatt = 0; iatt < ngatts; iatt++ )
    {
      nce(nc_inq_attname(nc_file_id, NC_GLOBAL, iatt, attname));
      nce(nc_inq_atttype(nc_file_id, NC_GLOBAL, attname, &xtype));
      nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, attname, &attlen));

      if ( xtype == NC_CHAR )
	{
	  if ( (int)attlen > attstringlen )
	    {
	      fprintf(stderr, "Attribute %s too large, skipped!\n", attname);
	      continue;
	    }
	  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, attname, attstring));
	  attstring[attlen] = 0;

	  atts->atxtname[atts->natxt]  = strdup(attname);
	  atts->atxtentry[atts->natxt] = strdup(attstring);
	  atts->natxt++;
	}
      else if ( xtype == NC_INT )
	{
	  if ( attlen > 1 )
	    {
	      fprintf(stderr, "Attribute %s too large, skipped!\n", attname);
	      continue;
	    }
	  nce(nc_get_att_int(nc_file_id, NC_GLOBAL, attname, &attint));
	}
      else if ( xtype == NC_DOUBLE )
	{
	  if ( attlen > 1 )
	    {
	      fprintf(stderr, "Attribute %s too large, skipped!\n", attname);
	      continue;
	    }
	  nce(nc_get_att_double(nc_file_id, NC_GLOBAL, attname, &attflt));
	}
      else
	{
	  fprintf(stderr, "Unsupported attribute type for %s\n", attname);
	}
    }

  *vars = (VAR*) malloc(max_vars*sizeof(VAR));

  varid = 0;
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {      
      nce(nc_inq_var(nc_file_id, ncvarid, name, &xtype, &nvdims, dimidsp, &nvatts));

      if ( varid >= max_vars ) cdoAbort("Too many variables (max = %d)!", max_vars);

      if ( nvdims == 4 )
	{
	  if ( dimidsp[0] == nhgl_dimid    && dimidsp[1] == nmp1_dimid &&
	       dimidsp[2] == complex_dimid && dimidsp[3] == lev_dimid )
	    {
	      read_fc4d(nc_file_id, name, &(*vars)[varid], gridIDfc, zaxisIDml, nhgl, nmp1);
	      varid++;
	    }
	  else if ( dimidsp[0] == nhgl_dimid    && dimidsp[1] == nmp1_dimid &&
	            dimidsp[2] == complex_dimid && dimidsp[3] == ilev_dimid )
	    {
	      read_fc4d(nc_file_id, name, &(*vars)[varid], gridIDfc, zaxisIDmlh, nhgl, nmp1);
	      varid++;
	    }
	  else
	    {
	      fprintf(stderr, "4D structure of %s unsupported!\n", name);
	    }
	}
      else if ( nvdims == 3 )
	{
	  if ( dimidsp[0] == lat_dimid && dimidsp[1] == lev_dimid &&
	       dimidsp[2] == lon_dimid)
	    {
	      read_gg3d(nc_file_id, name, &(*vars)[varid], gridIDgp, zaxisIDml);
	      varid++;
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == ilev_dimid &&
		    dimidsp[2] == lon_dimid)
	    {  
	      read_gg3d(nc_file_id, name, &(*vars)[varid], gridIDgp, zaxisIDmlh);
	      varid++;
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == belowsurface_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      read_gg3d(nc_file_id, name, &(*vars)[varid], gridIDgp, zaxisIDbsfc);
	      varid++;
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == n2_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      read_gg3d(nc_file_id, name, &(*vars)[varid], gridIDgp, zaxisIDn2);
	      varid++;
	    }
	  else
	    {
	      fprintf(stderr, "3D structure of %s unsupported!\n", name);
	    }
	}
      else if ( nvdims == 2 )
	{
	  if ( dimidsp[0] == lat_dimid && dimidsp[1] == lon_dimid)
	    {
	      nvals = nlon*nlat;
  
	      inivar(&(*vars)[varid], GRID_GAUSSIAN, ZAXIS_SURFACE,  0, name, NULL, NULL);

	      (*vars)[varid].gridID    = gridIDgp;
	      (*vars)[varid].zaxisID   = zaxisIDsfc;
	      (*vars)[varid].gridsize  = nvals;
	      (*vars)[varid].nlev      = 1;

	      (*vars)[varid].ptr = (double*) malloc(nvals*sizeof(double));

	      nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
	      nce(nc_get_var_double(nc_file_id, nc_var_id, (*vars)[varid].ptr));

	      varid++;
	    }
	  else if ( dimidsp[0] == nhgl_dimid && dimidsp[1] == lev_dimid)
	    {
	      nvals = nhgl;
  
	      inivar(&(*vars)[varid], GRID_GENERIC, ZAXIS_HYBRID,  0, name, NULL, NULL);

	      (*vars)[varid].gridID    = gridIDhgl;
	      (*vars)[varid].zaxisID   = zaxisIDml;
	      (*vars)[varid].gridsize  = nvals;
	      (*vars)[varid].nlev      = nlev;

	      (*vars)[varid].ptr = (double*) malloc(nvals*nlev*sizeof(double));

	      for ( i = 0; i < nlev; i++ )
		{
		  start[0] = 0;      start[1] = i;
		  count[0] = nvals;  count[1] = 1;  

		  nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
		  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[varid].ptr+i*nvals));
		}

	      varid++;
	    }
	  else
	    {
	      fprintf(stderr, "2D structure of %s unsupported!\n", name);
	    }
	}
      else if ( nvdims == 1 )
	{
	  if ( dimidsp[0] == lat_dimid && strcmp(name, "lat") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == lon_dimid && strcmp(name, "lon") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == nvclev_dimid && strcmp(name, "vct_a") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == nvclev_dimid && strcmp(name, "vct_b") == 0 )
	    {
	    }
	  else
	    {
	      fprintf(stderr, "1D structure of %s unsupported!\n", name);
	    }
	}
      else
	{
	  fprintf(stderr, "structure of %s unsupported!\n", name);
	}
    }

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = varid;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}

static
void export_e5res(const char *filename, VAR *vars, int nvars)
{
#if defined(HAVE_LIBNETCDF)
  int nc_var_id;
  int varid;
  size_t nvals;
  size_t start[4], count[4];
  int nlon, nlat, nlev/*, nvct*/, nfc, i;
  int gridIDgp = -1, zaxisIDml = -1;
  int nc_file_id;
  double *xvals, *yvals;
  const double *vct;
  int dimidsp[9];
  int lon_dimid, lat_dimid, nhgl_dimid, nlevp1_dimid, spc_dimid, nvclev_dimid;
  int complex_dimid, nmp1_dimid, belowsurface_dimid, lev_dimid, ilev_dimid;
  int /* surface_dimid, height2m_dimid, height10m_dimid,*/ n2_dimid;
  int lon = 0, lat = 0, nhgl = 0, nlevp1 = 0, spc = 0, nvclev = 0;
  int complex, nmp1, belowsurface, lev = 0, ilev = 0;
  int /* surface, height2m, height10m,*/ n2;
  int gridtype, zaxistype;

  /* create file */
  nce(nc_create(filename, NC_CLOBBER, &nc_file_id));

  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "file_type", strlen(strfiletype_res), strfiletype_res));

  n2 = 2;
  complex = 2;

  lon = 0; lat = 0; nfc = 0; lev = 0; belowsurface = 0;
  for ( varid = 0; varid < nvars; ++varid )
    {
      gridtype  = vars[varid].gridtype;
      zaxistype = vars[varid].zaxistype;

      if ( gridtype == GRID_GAUSSIAN && lat == 0 )
	{
	  gridIDgp = vars[varid].gridID;
	  lon = gridInqXsize(vars[varid].gridID);
	  lat = gridInqYsize(vars[varid].gridID);
	  nhgl = lat/2;
	}
      else if ( gridtype == GRID_FOURIER && nfc == 0 )
	{
	  nfc = gridInqSize(vars[varid].gridID);
	}

      if ( zaxistype == ZAXIS_HYBRID && lev == 0 )
	{
	  zaxisIDml = vars[varid].zaxisID;
	  lev = zaxisInqSize(vars[varid].zaxisID);
	  nlevp1 = lev + 1;
	  ilev = lev + 1;
	  nvclev = ilev;
	}
      else if ( zaxistype == ZAXIS_DEPTH_BELOW_LAND && belowsurface == 0 )
	{
	  belowsurface = zaxisInqSize(vars[varid].zaxisID);
	}
    }

  if ( lat == 0 ) cdoAbort("Gaussian grid not found!");
  if ( nfc == 0 ) cdoAbort("Fourier data not found!");
  if ( lev == 0 ) cdoAbort("Hybrid level not found!");

  nmp1 = (nfc/nhgl)/2;
  spc  = nmp1*(nmp1+1)/2;

  nlon = lon;
  nlat = lat;
  nlev = lev;
  
  nce(nc_def_dim(nc_file_id, "lon", lon, &lon_dimid));

  nce(nc_def_dim(nc_file_id, "lat", lat, &lat_dimid));

  nce(nc_def_dim(nc_file_id, "nhgl", nhgl, &nhgl_dimid));

  nce(nc_def_dim(nc_file_id, "nlevp1", nlevp1, &nlevp1_dimid));

  nce(nc_def_dim(nc_file_id, "spc", spc, &spc_dimid));

  nce(nc_def_dim(nc_file_id, "nvclev", nvclev, &nvclev_dimid));

  nce(nc_def_dim(nc_file_id, "complex", complex, &complex_dimid));

  nce(nc_def_dim(nc_file_id, "nmp1", nmp1, &nmp1_dimid));

  nce(nc_def_dim(nc_file_id, "belowsurface", belowsurface, &belowsurface_dimid));

  nce(nc_def_dim(nc_file_id, "lev", lev, &lev_dimid));

  nce(nc_def_dim(nc_file_id, "ilev", ilev, &ilev_dimid));

  /*
  nce(nc_def_dim(nc_file_id, "surface", surface, &surface_dimid));

  nce(nc_def_dim(nc_file_id, "height2m", height2m, &height2m_dimid));

  nce(nc_def_dim(nc_file_id, "height10m", height10m, &height10m_dimid));
  */
  nce(nc_def_dim(nc_file_id, "n2", n2, &n2_dimid));

  nce(nc_enddef(nc_file_id));


  for ( varid = 0; varid < nvars; varid++ )
    {
      gridtype  = vars[varid].gridtype;
      zaxistype = vars[varid].zaxistype;

      if ( gridtype == GRID_FOURIER && zaxistype == ZAXIS_HYBRID )
	{
	  nvals = nfc;

	  dimidsp[0] = nhgl_dimid;
	  dimidsp[1] = nmp1_dimid;
	  dimidsp[2] = complex_dimid;
	  dimidsp[3] = lev_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 4, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < nlev; i++ )
	    {
	      start[0] = 0;     start[1] = 0;    start[2] = 0;  start[3] = i;
	      count[0] = nhgl;  count[1] = nmp1; count[2] = 2;  count[3] = 1; 

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_FOURIER && zaxistype == ZAXIS_HYBRID_HALF )
	{
	  nvals = nfc;

	  dimidsp[0] = nhgl_dimid;
	  dimidsp[1] = nmp1_dimid;
	  dimidsp[2] = complex_dimid;
	  dimidsp[3] = ilev_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 4, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < ilev; i++ )
	    {
	      start[0] = 0;     start[1] = 0;    start[2] = 0;  start[3] = i;
	      count[0] = nhgl;  count[1] = nmp1; count[2] = 2;  count[3] = 1; 

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_GAUSSIAN && zaxistype == ZAXIS_HYBRID )
	{
	  nvals = lat*lon;

	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = lev_dimid;
	  dimidsp[2] = lon_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 3, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < nlev; i++ )
	    {
	      start[0] = 0;    start[1] = i;  start[2] = 0;
	      count[0] = lat;  count[1] = 1;  count[2] = lon;

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_GAUSSIAN && zaxistype == ZAXIS_HYBRID_HALF )
	{
	  nvals = lat*lon;

	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = ilev_dimid;
	  dimidsp[2] = lon_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 3, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < ilev; i++ )
	    {
	      start[0] = 0;    start[1] = i;  start[2] = 0;
	      count[0] = lat;  count[1] = 1;  count[2] = lon;

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_GAUSSIAN && zaxistype == ZAXIS_DEPTH_BELOW_LAND )
	{
	  nvals = lat*lon;

	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = belowsurface_dimid;
	  dimidsp[2] = lon_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 3, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < belowsurface; i++ )
	    {
	      start[0] = 0;    start[1] = i;  start[2] = 0;
	      count[0] = lat;  count[1] = 1;  count[2] = lon;

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_GAUSSIAN && zaxistype == ZAXIS_GENERIC )
	{
	  nvals = lat*lon;

	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = n2_dimid;
	  dimidsp[2] = lon_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 3, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < n2; i++ )
	    {
	      start[0] = 0;    start[1] = i;  start[2] = 0;
	      count[0] = lat;  count[1] = 1;  count[2] = lon;

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_GAUSSIAN && zaxistype == ZAXIS_HYBRID_HALF )
	{
	  nvals = lat*lon;

	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = ilev_dimid;
	  dimidsp[2] = lon_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 3, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < ilev; i++ )
	    {
	      start[0] = 0;    start[1] = i;  start[2] = 0;
	      count[0] = lat;  count[1] = 1;  count[2] = lon;

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else if ( gridtype == GRID_GAUSSIAN && zaxistype == ZAXIS_SURFACE )
	{
	  dimidsp[0] = lat_dimid;
	  dimidsp[1] = lon_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 2, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));
	  nce(nc_put_var_double(nc_file_id, nc_var_id, vars[varid].ptr));
	}
      else if ( gridtype == GRID_GENERIC && zaxistype == ZAXIS_HYBRID )
	{
	  nvals = nhgl;

	  dimidsp[0] = nhgl_dimid;
	  dimidsp[1] = lev_dimid;

	  nce(nc_redef(nc_file_id));
	  nce(nc_def_var(nc_file_id, vars[varid].name, NC_DOUBLE, 2, dimidsp, &nc_var_id));
	  nce(nc_enddef(nc_file_id));

	  for ( i = 0; i < nlev; i++ )
	    {
	      start[0] = 0;      start[1] = i;
	      count[0] = nvals;  count[1] = 1;

	      nce(nc_put_vara_double(nc_file_id, nc_var_id, start, count, vars[varid].ptr+i*nvals));
	    }  
	}
      else
	{
	  fprintf(stderr, "structure of %s unsupported!\n", vars[varid].name);
	}
    }

  /* define gaussian grid */

  nlon = lon;
  nlat = lat;

  xvals = (double*) malloc(nlon*sizeof(double));
  yvals = (double*) malloc(nlat*sizeof(double));

  gridInqXvals(gridIDgp, xvals);
  gridInqYvals(gridIDgp, yvals);


  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "lat", NC_DOUBLE, 1, &lat_dimid, &nc_var_id));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, yvals));
   
  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "lon", NC_DOUBLE, 1, &lon_dimid, &nc_var_id));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, xvals));

  free(xvals);
  free(yvals);

  /* define model level */

  nlev = lev;
  //nvct = nvclev*2;

  /* vct   = (double*) malloc(nvct*sizeof(double)); */

  vct = zaxisInqVctPtr(zaxisIDml);

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "vct_a", NC_DOUBLE, 1, &nvclev_dimid, &nc_var_id));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_redef(nc_file_id));
  nce(nc_def_var(nc_file_id, "vct_b", NC_DOUBLE, 1, &nvclev_dimid, &nc_var_id));
  nce(nc_enddef(nc_file_id));
  nce(nc_put_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  /* free(vct); */

  /*close input file */
  nce(nc_close(nc_file_id));

#else
  cdoAbort("netCDF support not compiled in!");
#endif
}


void *Echam5ini(void *argument)
{
  int operatorID;
  int operfunc;
  int IMPORT_E5ML, IMPORT_E5RES;
  int EXPORT_E5ML, EXPORT_E5RES;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs = 0;
  int recID, varID, levelID;
  int vlistID1, vlistID2;
  int nvars = 0;
  int iv, nlev;
  int gridsize, nmiss;
  int taxisID, tsID;

  cdoInitialize(argument);

  IMPORT_E5ML  = cdoOperatorAdd("import_e5ml",  func_read,  0, NULL);
  IMPORT_E5RES = cdoOperatorAdd("import_e5res", func_read,  0, NULL);
  EXPORT_E5ML  = cdoOperatorAdd("export_e5ml",  func_write, 0, NULL);
  EXPORT_E5RES = cdoOperatorAdd("export_e5res", func_write, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  if ( operatorID == EXPORT_E5ML && processSelf() != 0 )
    cdoAbort("This operator can't be linked with other operators!");

  if ( operfunc == func_read )
    {
      VAR *vars = NULL;
      ATTS atts;
      int iatt;

      iniatts(&atts);

      if ( operatorID == IMPORT_E5ML )
	nvars = import_e5ml(cdoStreamName(0)->args, &vars);
      else if ( operatorID == IMPORT_E5RES )
	nvars = import_e5res(cdoStreamName(0)->args, &vars, &atts);
      else
	cdoAbort("Operator not implemented!");

      if ( nvars == 0 ) cdoAbort("Unsupported file type!");
      
      vlistID2 = vlistCreate();
      vlistDefNtsteps(vlistID2, 0);

      for ( iv = 0; iv < nvars; iv++ )
	{/*
	  fprintf(stderr, "%d %s %d %d %d %d\n", iv, vars[iv].name, vars[iv].gridID, vars[iv].zaxisID, gridInqSize(vars[iv].gridID), zaxisInqSize(vars[iv].zaxisID));
	 */
	  varID = vlistDefVar(vlistID2, vars[iv].gridID, vars[iv].zaxisID, TSTEP_CONSTANT);
	  if ( vars[iv].code > 0 ) vlistDefVarCode(vlistID2, varID, vars[iv].code);
	  if ( vars[iv].name )     vlistDefVarName(vlistID2, varID, vars[iv].name);
	  if ( vars[iv].longname ) vlistDefVarLongname(vlistID2, varID, vars[iv].longname);
          if ( vars[iv].units )    vlistDefVarUnits(vlistID2, varID, vars[iv].units);
	  vlistDefVarDatatype(vlistID2, varID, DATATYPE_FLT64);
	}

      for ( iatt = 0; iatt < atts.natxt; ++iatt )
	{
	  /* printf("%s: %s\n", atts.atxtname[iatt], atts.atxtentry[iatt]); */
	  vlistDefAttTxt(vlistID2, CDI_GLOBAL, atts.atxtname[iatt],
			 (int)strlen(atts.atxtentry[iatt])+1, atts.atxtentry[iatt]);
	}

      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      vlistDefTaxis(vlistID2, taxisID);

      if ( cdoDefaultFileType == CDI_UNDEFID )
	cdoDefaultFileType = FILETYPE_NC;

      streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

      streamDefVlist(streamID2, vlistID2);

      tsID = 0;
      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = vars[varID].gridsize;
	  nlev     = vars[varID].nlev;

	  for ( levelID = 0; levelID < nlev; levelID++ )
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, vars[varID].ptr+levelID*gridsize, 0);
	    }
	}

      streamClose(streamID2);

      vlistDestroy(vlistID2);
    }
  else if ( operfunc == func_write )
    {
      VAR *vars = NULL;
      int code, gridID, zaxisID, gridtype, zaxistype, gridsize, nlev;
      char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
      int taxisID, vdate, vtime;
      int ntr = 0;

      streamID1 = streamOpenRead(cdoStreamName(0));

      vlistID1 = streamInqVlist(streamID1);
      taxisID = vlistInqTaxis(vlistID1);

      nvars = vlistNvars(vlistID1);

      vars = (VAR*) malloc(nvars*sizeof(VAR));

      for ( varID = 0; varID < nvars; ++varID )
	{
	  code = vlistInqVarCode(vlistID1, varID);
	  vlistInqVarName(vlistID1, varID, name);
	  vlistInqVarLongname(vlistID1, varID, longname);
	  vlistInqVarUnits(vlistID1, varID, units);

	  if ( code < 0 ) code = 0;
	  if ( strncmp(name, "var", 3) == 0 )
	    {
	      if ( code > 0 )
		{
		  if ( code == 133 )
		    { strcpy(name, "Q"); strcpy(longname, "specific humidity"); strcpy(units, "kg/kg"); }
		  if ( code == 138 )
		    { strcpy(name, "SVO"); strcpy(longname, "vorticity"); strcpy(units, "1/s"); }
		  if ( code == 155 )
		    { strcpy(name, "SD"); strcpy(longname, "divergence"); strcpy(units, "1/s"); }
		  if ( code == 130 )
		    { strcpy(name, "STP"); strcpy(longname, "temperature"); strcpy(units, "K"); }
		  if ( code == 152 )
		    { strcpy(name, "LSP"); strcpy(longname, "log surface pressure"); }
		}
	    }
	  else if ( strncmp(name, "LSP", 3) == 0 ) code = 152;

	  gridID  = vlistInqVarGrid(vlistID1, varID);
	  zaxisID = vlistInqVarZaxis(vlistID1, varID);

	  gridtype  = gridInqType(gridID);
	  zaxistype = zaxisInqType(zaxisID);

	  if ( gridtype == GRID_SPECTRAL && ntr == 0 )
	    {
	      ntr = gridInqTrunc(gridID);
	    }

	  gridsize = gridInqSize(gridID);
	  nlev     = zaxisInqSize(zaxisID);

	  if ( zaxistype == ZAXIS_HYBRID && nlev == 1 ) zaxistype = ZAXIS_SURFACE;

	  inivar(&vars[varID], gridtype, zaxistype, code, name, longname, units);
	  
	  vars[varID].gridID    = gridID;
	  vars[varID].zaxisID   = zaxisID;
	  vars[varID].gridsize  = gridsize;
	  vars[varID].nlev      = nlev;

	  vars[varID].ptr = (double*) malloc(nlev*gridsize*sizeof(double));
	}

      nrecs = streamInqTimestep(streamID1, 0);
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      if ( vdate == 0 )
	{
	  vdate = 19890101;
	  vtime = 120000;
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  gridID = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);

	  streamReadRecord(streamID1, vars[varID].ptr+levelID*gridsize, &nmiss);
	}

      streamClose(streamID1);

      if ( operatorID == EXPORT_E5ML )
	export_e5ml(cdoStreamName(1)->args, vars, nvars, vdate, vtime, ntr);
      else if ( operatorID == EXPORT_E5RES )
	export_e5res(cdoStreamName(1)->args, vars, nvars);
      else
	cdoAbort("Operator not implemented!");
    }
  else
    {
      cdoAbort("Internal error!");
    }

  /*
  vlistDestroy(vlistID2);

  */
  cdoFinish();

  return (0);
}
