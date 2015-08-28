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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"

#define UNDEFID -1

#define MAX_LINE_LEN 65536


typedef struct {
  double *vals;
  double *lbounds;
  double *ubounds;
  double *vct;
  int     vctsize;
  int     type;
  int     size;
  char    name[CDI_MAX_NAME];
  char    longname[CDI_MAX_NAME];
  char    units[CDI_MAX_NAME];
}
zaxis_t;


void zaxisInit(zaxis_t *zaxis)
{
  zaxis->vals        = NULL;
  zaxis->lbounds     = NULL;
  zaxis->ubounds     = NULL;
  zaxis->vct         = NULL;
  zaxis->type        = UNDEFID;
  zaxis->vctsize     = 0;
  zaxis->size        = 0;
  zaxis->name[0]     = 0;
  zaxis->longname[0] = 0;
  zaxis->units[0]    = 0;
}


static int getoptname(char *optname, const char *optstring, int nopt)
{
  int i, nerr = 0;
  size_t namelen;
  const char *pname;
  const char *pend;

  pname = optstring;
  pend  = optstring;

  for ( i = 0; i < nopt; i++ )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL )
	break;
      else
	pname = pend + 1;
    }

  if ( pend )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL )
	namelen = strlen(pname);
      else
	namelen = pend - pname;

      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return (nerr);
}


int zaxisDefine(zaxis_t zaxis)
{
  int zaxisID = UNDEFID;

  if ( zaxis.type == -1 ) Error("zaxistype undefined!");

  if ( zaxis.size == 0 ) Error("zaxis size undefined!");

  zaxisID = zaxisCreate(zaxis.type, zaxis.size);

  if ( zaxis.vals )
    {
      zaxisDefLevels(zaxisID, zaxis.vals);
      free(zaxis.vals);
    }
  if ( zaxis.lbounds )
    {
      zaxisDefLbounds(zaxisID, zaxis.lbounds);
      free(zaxis.lbounds);
    }
  if ( zaxis.ubounds )
    {
      zaxisDefUbounds(zaxisID, zaxis.ubounds);
      free(zaxis.ubounds);
    }

  if ( zaxis.name[0] )     zaxisDefName(zaxisID, zaxis.name);
  if ( zaxis.longname[0] ) zaxisDefLongname(zaxisID, zaxis.longname);
  if ( zaxis.units[0] )    zaxisDefUnits(zaxisID, zaxis.units);

  if ( zaxis.type == ZAXIS_HYBRID || zaxis.type == ZAXIS_HYBRID_HALF )
    {
      if ( zaxis.vctsize && zaxis.vct )
	zaxisDefVct(zaxisID, zaxis.vctsize, zaxis.vct);
      else
	Warning("vct undefined!");	    
    }

  return (zaxisID);
}


static char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return (pline);
}


int zaxisFromFile(FILE *gfp, const char *dname)
{
  char line[MAX_LINE_LEN], *pline;
  int zaxisID;
  int lerror;
  size_t i, len;
  zaxis_t zaxis;

  zaxisInit(&zaxis);

  while ( readline(gfp, line, MAX_LINE_LEN) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      len = strlen(line);

      lerror = FALSE;
      for ( i = 0; i < len; ++i )
	if ( !(line[i] == 9 || (line[i] > 31 && line[i] < 127)) )
	  {
	    lerror = TRUE;
	    line[i] = '#';
	  }
      if ( lerror ) cdoAbort("Zaxis description file >%s< contains illegal characters (line: %s)!", dname, line);

      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      if ( cmpstrlen(pline, "zaxistype", len) == 0 || 
	   cmpstrlen(pline, "type", len) == 0 )
	{
	  if ( *pline == 'z' )
	    pline = skipSeparator(pline + 9);
	  else
	    pline = skipSeparator(pline + 4);

	  if ( cmpstrlen(pline, "pressure", len) == 0 )
	    zaxis.type = ZAXIS_PRESSURE;
	  else if ( cmpstrlen(pline, "hybrid_half", len)  == 0 )
	    zaxis.type = ZAXIS_HYBRID_HALF;
	  else if ( cmpstrlen(pline, "hybrid", len)  == 0 )
	    zaxis.type = ZAXIS_HYBRID;
	  else if ( cmpstrlen(pline, "height", len) == 0 )
	    zaxis.type = ZAXIS_HEIGHT;
	  else if ( cmpstrlen(pline, "depth below sea", len) == 0 ||
		    cmpstrlen(pline, "depth_below_sea", len) == 0 )
	    zaxis.type = ZAXIS_DEPTH_BELOW_SEA;
	  else if ( cmpstrlen(pline, "depth below land", len) == 0 ||
		    cmpstrlen(pline, "depth_below_land", len) == 0 )
	    zaxis.type = ZAXIS_DEPTH_BELOW_LAND;
	  else if ( cmpstrlen(pline, "isentropic", len)  == 0 )
	    zaxis.type = ZAXIS_ISENTROPIC;
	  else if ( cmpstrlen(pline, "surface", len)  == 0 )
	    zaxis.type = ZAXIS_SURFACE;
	  else if ( cmpstrlen(pline, "generic", len)  == 0 )
	    zaxis.type = ZAXIS_GENERIC;
	  else
	    cdoAbort("Invalid zaxisname : %s (zaxis description file: %s)", pline, dname);
	}
      else if ( cmpstrlen(pline, "size", len)  == 0 )
	{
	  zaxis.size = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "vctsize", len)  == 0 )
	{
	  zaxis.vctsize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "name", len)  == 0 )
	{
	  strcpy(zaxis.name, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "longname", len)  == 0 )
	{
	  strcpy(zaxis.longname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "units", len)  == 0 )
	{
	  strcpy(zaxis.units, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "levels", len)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + len);
	  
	      zaxis.vals = (double*) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			cdoAbort("Incomplete command: >levels< (zaxis description file: %s)", dname);

		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.vals[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      cdoAbort("size undefined (zaxis description file: %s)!", dname);
	    }
	}
      else if ( cmpstrlen(pline, "vct", len)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.vctsize > 0 )
	    {
	      pline = skipSeparator(pline + len);
	  
	      zaxis.vct = (double*) malloc(zaxis.vctsize*sizeof(double));
	      for ( i = 0; i < zaxis.vctsize; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			cdoAbort("Incomplete command: >vct< (zaxis description file: %s)", dname);

		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.vct[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      cdoAbort("vctsize undefined (zaxis description file: %s)!", dname);
	    }
	}
      else if ( cmpstrlen(pline, "lbounds", len)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + len);
	  
	      zaxis.lbounds = (double*) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			cdoAbort("Incomplete command: >lbounds< (zaxis description file: %s)", dname);

		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.lbounds[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      cdoAbort("size undefined (zaxis description file: %s)!", dname);
	    }
	}
      else if ( cmpstrlen(pline, "ubounds", len)  == 0 )
	{
	  int i;
	  double flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + len);
	  
	      zaxis.ubounds = (double*) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			cdoAbort("Incomplete command: >ubounds< (zaxis description file: %s)", dname);

		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  flev = 0;
		  sscanf(pline, "%lg", &flev);
		  zaxis.ubounds[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      cdoAbort("size undefined (zaxis description file: %s)!", dname);
	    }
	}
      else
	cdoAbort("Invalid zaxis command : >%s< (zaxis description file: %s)", pline, dname);
    }

  zaxisID = zaxisDefine(zaxis);

  return (zaxisID);
}


int zaxisFromName(const char *zaxisname)
{
  const char *pline;
  int zaxisID = UNDEFID;
  zaxis_t zaxis;

  zaxisInit(&zaxis);

  pline = zaxisname;
  if ( cmpstr(pline, "surface") == 0 ) /* surface */
    {
      zaxis.type = ZAXIS_SURFACE;
      zaxis.size = 1;
      zaxis.vals = (double*) malloc(zaxis.size*sizeof(double));
      zaxis.vals[0] = 0;
    }

  if ( zaxis.type != -1 ) zaxisID = zaxisDefine(zaxis);

  return (zaxisID);
}


int cdoDefineZaxis(const char *zaxisfile)
{
  FILE *zfp;
  int zaxisID = -1;

  zfp = fopen(zaxisfile, "r");
  if ( zfp == NULL )
    {
      zaxisID = zaxisFromName(zaxisfile);

      if ( zaxisID == -1 ) cdoAbort("Open failed on %s!", zaxisfile);
    }
  else
    {
      zaxisID = zaxisFromFile(zfp, zaxisfile);
      fclose(zfp);
    }

  if ( zaxisID == -1 ) cdoAbort("Invalid zaxis description file %s!", zaxisfile);

  return (zaxisID);
}


void defineZaxis(const char *zaxisarg)
{
  char zaxisfile[4096];
  int nfile = 0;

  while ( getoptname(zaxisfile, zaxisarg, nfile++) == 0 )
    {      
      (void) cdoDefineZaxis(zaxisfile);
    }
}

static
int ztype2ltype(int zaxistype)
{
  int ltype = -1;

  if      ( zaxistype == ZAXIS_SURFACE           )  ltype =   1;
  else if ( zaxistype == ZAXIS_PRESSURE          )  ltype = 100;
  else if ( zaxistype == ZAXIS_ALTITUDE          )  ltype = 103;
  else if ( zaxistype == ZAXIS_HEIGHT            )  ltype = 105;
  else if ( zaxistype == ZAXIS_SIGMA             )  ltype = 107;
  else if ( zaxistype == ZAXIS_HYBRID            )  ltype = 109;
  else if ( zaxistype == ZAXIS_HYBRID_HALF       )  ltype = 109;
  else if ( zaxistype == ZAXIS_DEPTH_BELOW_LAND  )  ltype = 111;
  else if ( zaxistype == ZAXIS_ISENTROPIC        )  ltype = 113;
  else if ( zaxistype == ZAXIS_DEPTH_BELOW_SEA   )  ltype = 160;

  return (ltype);
}


int zaxis2ltype(int zaxisID)
{
  int ltype;
  int zaxistype;

  zaxistype = zaxisInqType(zaxisID);

  ltype = zaxisInqLtype(zaxisID);

  if ( ltype <= 0 ) ltype = ztype2ltype(zaxistype);

  return (ltype);
}
