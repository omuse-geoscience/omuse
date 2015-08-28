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

      Filedes    pardes          Parameter description
      Filedes    griddes         Grid description
      Filedes    vct             Vertical coordinate table
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


static
void printAtts(FILE *fp, int vlistID, int varID)
{
#define MAXATT 8192
  int natts, ia;
  char attname[1024];
  int atttype, attlen;
  char atttxt[MAXATT];
  int attint[MAXATT];
  double attflt[MAXATT];
  int i;

  vlistInqNatts(vlistID, varID, &natts);

  for ( ia = 0; ia < natts; ++ia )
    {
      vlistInqAtt(vlistID, varID, ia, attname, &atttype, &attlen);
      if ( atttype == DATATYPE_INT )
	{
	  if ( attlen > MAXATT ) attlen = MAXATT;
	  vlistInqAttInt(vlistID, varID, attname, attlen, attint);
	  fprintf(fp, "  %s=", attname);
	  for ( i = 0; i < attlen; ++i)
	    {
	      if ( i > 0 ) fprintf(fp, ", ");
	      fprintf(fp, "%d", attint[i]);
	    }
	  fprintf(fp, "\n");
	}
      else if ( atttype == DATATYPE_FLT )
	{
	  if ( attlen > MAXATT ) attlen = MAXATT;
	  vlistInqAttFlt(vlistID, varID, attname, MAXATT, attflt);
	  fprintf(fp, "  %s=", attname);
	  for ( i = 0; i < attlen; ++i)
	    {
	      if ( i > 0 ) fprintf(fp, ", ");
	      fprintf(fp, "%g", attflt[i]);
	    }
	  fprintf(fp, "\n");
	}
      else if ( atttype == DATATYPE_TXT )
	{
	  vlistInqAttTxt(vlistID, varID, attname, sizeof(atttxt), atttxt);
	  atttxt[attlen] = 0;
	  fprintf(fp, "  %s=\"%s\"\n", attname, atttxt);
	}
    }
}

static
void partab(FILE *fp, int vlistID, int option)
{
  int varID, datatype = -1;
  int param;
  char pstr[32];
  char paramstr[32];
  char varname[CDI_MAX_NAME], varlongname[CDI_MAX_NAME], varstdname[CDI_MAX_NAME], varunits[CDI_MAX_NAME];
  int natts;
  int nvars;
  int chunktype;
  int linebreak = 1;
  double missval;
      
  nvars  = vlistNvars(vlistID);
  if ( option == 4 ) linebreak = 0;

  if ( option == 2 )
    {
      vlistInqNatts(vlistID, CDI_GLOBAL, &natts);
      if ( natts > 0 )
	{
	  fprintf(fp, "&parameter\n");
	  fprintf(fp, "  name=_GLOBAL_\n");
	  printAtts(fp, vlistID, CDI_GLOBAL);
	  fprintf(fp, "/\n");
	}
    }

  if ( nvars > 1 )
    {
      datatype = vlistInqVarDatatype(vlistID, 0);
      for ( varID = 1; varID < nvars; varID++ )
	{
	  if ( datatype != vlistInqVarDatatype(vlistID, varID) )
	    {
	      datatype = -1;
	      break;
	    }
	}

      if ( datatype != -1 )
	{
	  fprintf(fp, "&parameter");
	  if ( linebreak ) fprintf(fp, "\n");
	  fprintf(fp, "  name=_default_");
	  if ( linebreak ) fprintf(fp, "\n");
	  if ( datatype2str(datatype, pstr) == 0 )
	    {
	      fprintf(fp, "  datatype=%s", pstr);
	      if ( linebreak ) fprintf(fp, "\n");
	    }
	  fprintf(fp, "/\n");
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      fprintf(fp, "&parameter");
      if ( linebreak ) fprintf(fp, "\n");
      
      varname[0]     = 0;
      varlongname[0] = 0;
      varunits[0]    = 0;
      param    = vlistInqVarParam(vlistID, varID);
      missval  = vlistInqVarMissval(vlistID, varID);
      vlistInqVarName(vlistID, varID, varname);
      /* printf("1>%s<\n", varname); */
      vlistInqVarStdname(vlistID, varID, varstdname);
      /* printf("2>%s<\n", varname); */
      vlistInqVarLongname(vlistID, varID, varlongname);
      /* printf("3>%s<\n", varname); */
      vlistInqVarUnits(vlistID, varID, varunits);
            
      fprintf(fp, "  name=%s", varname);
      if ( linebreak ) fprintf(fp, "\n");
      // if ( code   > 0 ) fprintf(fp, "  code=%d\n", code);
      // if ( tabnum > 0 ) fprintf(fp, "  table=%d\n", tabnum);
      if ( param >= 0 )
	{
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  fprintf(fp, "  param=%s", paramstr);
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( strlen(varstdname) )
	{
	  fprintf(fp, "  standard_name=%s", varstdname);
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( strlen(varlongname) )
	{
	  fprintf(fp, "  long_name=\"%s\"", varlongname);
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( strlen(varunits) )
	{
	  fprintf(fp, "  units=\"%s\"", varunits);
	  if ( linebreak ) fprintf(fp, "\n");
	}

      if ( datatype == -1 )
	if ( datatype2str(vlistInqVarDatatype(vlistID, varID), pstr) == 0 )
	  {
	    fprintf(fp, "  datatype=%s", pstr);
	    if ( linebreak ) fprintf(fp, "\n");
	  }

      chunktype = vlistInqVarChunkType(vlistID, varID);
      if ( chunktype == CHUNK_AUTO )
	{
	  fprintf(fp, "  chunktype=auto");
	  if ( linebreak ) fprintf(fp, "\n");
	}
      else if ( chunktype == CHUNK_GRID )
	{
	  fprintf(fp, "  chunktype=grid");
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( chunktype == CHUNK_LINES )
	{
	  fprintf(fp, "  chunktype=lines");
	  if ( linebreak ) fprintf(fp, "\n");
	}
      
      if ( option == 2 ) printAtts(fp, vlistID, varID);
      if ( option == 2 ) 
	fprintf(fp, "  missing_value=%g\n", missval);
      
      if ( !linebreak ) fprintf(fp, "  ");
      fprintf(fp, "/\n");
    }   
}

static
void filedes(int streamID)
{
  int filetype;

  printf("\n");
  filetype = streamInqFiletype(streamID);
  switch ( filetype )
    {
    case FILETYPE_GRB:
      printf("  GRIB data\n");
      break;
    case FILETYPE_GRB2:
      printf("  GRIB2 data\n");
      break;
    case FILETYPE_NC:
      printf("  netCDF data\n");
      break;
    case FILETYPE_NC2:
      printf("  netCDF2 data\n");
      break;
    case FILETYPE_NC4:
      printf("  netCDF4 data\n");
      break;
    case FILETYPE_NC4C:
      printf("  netCDF4 classic data\n");
      break;
    case FILETYPE_SRV:
      printf("  SERVICE data\n");
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  byteorder is BIGENDIAN\n"); break;
	case CDI_LITTLEENDIAN:
	  printf("  byteorder is LITTLEENDIAN\n"); break;
	default:
	  printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	}
      break;
    case FILETYPE_EXT:
      printf("  EXTRA data\n");
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  byteorder is BIGENDIAN\n"); break;
	case CDI_LITTLEENDIAN:
	  printf("  byteorder is LITTLEENDIAN\n"); break;
	default:
	  printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	}
      break;
    case FILETYPE_IEG:
      printf("  IEG data\n");
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  byteorder is BIGENDIAN\n"); break;
	case CDI_LITTLEENDIAN:
	  printf("  byteorder is LITTLEENDIAN\n"); break;
	default:
	  printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	}
      break;
    default:
      printf("  unsupported filetype %d\n" , filetype);
      break;
    }
  
  printf("\n");
}


void *Filedes(void *argument)
{
  int GRIDDES, GRIDDES2, ZAXISDES, VCT, VCT2, PARDES, FILEDES, VLIST, SPARTAB, PARTAB, PARTAB2;
  int operatorID;
  int streamID = 0;
  int zaxisID;
  int nvars, ngrids, nzaxis;
  int type, index;
  int vlistID;

  cdoInitialize(argument);

  GRIDDES  = cdoOperatorAdd("griddes",   0, 0, NULL);
  GRIDDES2 = cdoOperatorAdd("griddes2",  0, 0, NULL);
  ZAXISDES = cdoOperatorAdd("zaxisdes",  0, 0, NULL);
  VCT      = cdoOperatorAdd("vct",       0, 0, NULL);
  VCT2     = cdoOperatorAdd("vct2",      0, 0, NULL);
  PARDES   = cdoOperatorAdd("pardes",    0, 0, NULL);
  FILEDES  = cdoOperatorAdd("filedes",   0, 0, NULL);
  VLIST    = cdoOperatorAdd("vlist",     0, 0, NULL);
  SPARTAB  = cdoOperatorAdd("spartab",   0, 0, NULL);
  PARTAB   = cdoOperatorAdd("partab",    0, 0, NULL);
  PARTAB2  = cdoOperatorAdd("partab2",   0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  nvars  = vlistNvars(vlistID);
  ngrids = vlistNgrids(vlistID);
  nzaxis = vlistNzaxis(vlistID);

  if ( operatorID == GRIDDES || operatorID == GRIDDES2 )
    {
      int opt = 0;
      if ( operatorID == GRIDDES ) opt = 1;
      for ( index = 0; index < ngrids; index++ )
	gridPrint(vlistGrid(vlistID, index), index+1, opt);
    }
  else if ( operatorID == ZAXISDES )
    {
      for ( index = 0; index < nzaxis; index++ )
	zaxisPrint(vlistZaxis(vlistID, index), index+1);
    }
  else if ( operatorID == VCT || operatorID == VCT2 )
    {
      for ( index = 0; index < nzaxis; index++)
	{
	  zaxisID = vlistZaxis(vlistID, index);
	  type = zaxisInqType(zaxisID);
	  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
	    {
	      int i, vctsize;
	      const double *vct;

	      vctsize = zaxisInqVctSize(zaxisID);
	      vct     = zaxisInqVctPtr(zaxisID);
		
	      if ( vctsize%2 == 0 )
		{
		  if ( operatorID == VCT )
		    {
		      fprintf(stdout, "#   k         vct_a(k) [Pa]             vct_b(k) []\n");
		      for ( i = 0; i < vctsize/2; i++ )
			fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
		    }
		  else
		    {
		      int nbyte0, nbyte;
		      fprintf(stdout, "vctsize   = %d\n", vctsize);
		      nbyte0 = fprintf(stdout, "vct       = ");
		      nbyte = nbyte0;
		      for ( i = 0; i < vctsize; i++ )
			{
			  if ( nbyte > 70 || i == vctsize/2 )
			    {
			      fprintf(stdout, "\n%*s", nbyte0, "");
			      nbyte = nbyte0;
			    }
			  nbyte += fprintf(stdout, "%.9g ", vct[i]);
			}
		      fprintf(stdout, "\n");
		    }
		}
	      else
		for ( i = 0; i < vctsize; i++ )
		  fprintf(stdout, "%5d %25.17f\n", i, vct[i]);

	      break;
	    }
	}
    }
  else if ( operatorID == VLIST )
    {
      vlistPrint(vlistID);
    }
  else if ( operatorID == PARDES )
    {
      int varID, code;
      char varname[CDI_MAX_NAME], varlongname[CDI_MAX_NAME], varunits[CDI_MAX_NAME];

      for ( varID = 0; varID < nvars; varID++ )
	{
	  varname[0]     = 0;
	  varlongname[0] = 0;
	  varunits[0]    = 0;
	  code     = vlistInqVarCode(vlistID, varID);
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarLongname(vlistID, varID, varlongname);
	  vlistInqVarUnits(vlistID, varID, varunits);
	  fprintf(stdout, "%4d  %-12s", code, varname);
	  if ( strlen(varlongname) )
	    {
	      fprintf(stdout, "  %s", varlongname);
	      if ( strlen(varunits) )
		fprintf(stdout, " [%s]", varunits);
	    }
	  fprintf(stdout, "\n");
	}   
    }
  else if ( operatorID == PARTAB || operatorID == SPARTAB || operatorID == PARTAB2 )
    {
      int option = 1;

      if ( operatorID == SPARTAB ) option = 4;
      if ( operatorID == PARTAB2 ) option = 2;
      
      partab(stdout, vlistID, option);
    }
  else if ( operatorID == FILEDES )
    {
      filedes(streamID);
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}
