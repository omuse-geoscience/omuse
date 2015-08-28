#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "cdi_int.h"
#include "varscan.h"
#include "datetime.h"
#include "ieg.h"
#include "stream_fcommon.h"
#include "stream_ieg.h"
#include "vlist.h"


#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID

#define SINGLE_PRECISION  4
#define DOUBLE_PRECISION  8

#if defined (HAVE_LIBIEG)


typedef struct {
  int param;
  int level;
} IEGCOMPVAR;


int iegInqDatatype(int prec)
{
  int datatype;

  if ( prec == DOUBLE_PRECISION ) datatype = DATATYPE_FLT64;
  else                            datatype = DATATYPE_FLT32;

  return (datatype);
}


int iegDefDatatype(int datatype)
{
  int prec;

  if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    Error("CDI/IEG library does not support complex numbers!");

  if ( datatype != DATATYPE_FLT32 && datatype != DATATYPE_FLT64 )
    datatype = DATATYPE_FLT32;

  if ( datatype == DATATYPE_FLT64 ) prec = DOUBLE_PRECISION;
  else                              prec = SINGLE_PRECISION;

  return (prec);
}

/* not used
int iegInqRecord(stream_t *streamptr, int *varID, int *levelID)
{
  int status;
  int fileID;
  int icode, ilevel;
  int zaxisID = -1;
  int vlistID;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  *varID   = -1;
  *levelID = -1;

  status = iegRead(fileID, iegp);
  if ( status != 0 ) return (0);

  icode  = IEG_P_Parameter(iegp->ipdb);
  if ( IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER )
    ilevel = IEG_P_Level1(iegp->ipdb);
  else
    ilevel = IEG_P_Level2(iegp->ipdb);

  *varID = vlistInqVarID(vlistID, icode);

  if ( *varID == UNDEFID ) Error("Code %d undefined", icode);

  zaxisID = vlistInqVarZaxis(vlistID, *varID);

  *levelID = zaxisInqLevelID(zaxisID, (double) ilevel);

  return (1);
}
*/

void iegReadRecord(stream_t *streamptr, double *data, int *nmiss)
{
  int vlistID, fileID;
  int status;
  int recID, vrecID, tsID;
  off_t recpos;
  int varID, gridID;
  int i, size;
  double missval;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;
  tsID    = streamptr->curTsID;
  vrecID  = streamptr->tsteps[tsID].curRecID;
  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  recpos  = streamptr->tsteps[tsID].records[recID].position;
  varID   = streamptr->tsteps[tsID].records[recID].varID;

  fileSetPos(fileID, recpos, SEEK_SET);

  status = iegRead(fileID, iegp);
  if ( status != 0 )
    Error("Could not read IEG record!");

  iegInqDataDP(iegp, data);

  missval = vlistInqVarMissval(vlistID, varID);
  gridID  = vlistInqVarGrid(vlistID, varID);
  size    = gridInqSize(gridID);

  streamptr->numvals += size;

  *nmiss = 0;
  for ( i = 0; i < size; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}

static
int iegGetZaxisType(int iegleveltype)
{
  int leveltype = 0;

  switch ( iegleveltype )
    {
    case IEG_LTYPE_SURFACE:
      {
	leveltype = ZAXIS_SURFACE;
	break;
      }
    case IEG_LTYPE_99:
    case IEG_LTYPE_ISOBARIC:
      {
	leveltype = ZAXIS_PRESSURE;
	break;
      }
    case IEG_LTYPE_HEIGHT:
      {
	leveltype = ZAXIS_HEIGHT;
	break;
      }
    case IEG_LTYPE_ALTITUDE:
      {
	leveltype = ZAXIS_ALTITUDE;
	break;
      }
    case IEG_LTYPE_HYBRID:
    case IEG_LTYPE_HYBRID_LAYER:
      {
	leveltype = ZAXIS_HYBRID;
	break;
      }
    case IEG_LTYPE_LANDDEPTH:
    case IEG_LTYPE_LANDDEPTH_LAYER:
      {
	leveltype = ZAXIS_DEPTH_BELOW_LAND;
	break;
      }
    case IEG_LTYPE_SEADEPTH:
      {
	leveltype = ZAXIS_DEPTH_BELOW_SEA;
	break;
      }
    default:
      {
	leveltype = ZAXIS_GENERIC;
	break;
      }
    }

  return (leveltype);
}


void iegDefTime(int *pdb, int date, int time, int taxisID)
{
  int year, month, day, hour, minute, second;
  int timetype = -1;

  if ( taxisID != -1 ) timetype = taxisInqType(taxisID);

  if ( timetype == TAXIS_ABSOLUTE || timetype == TAXIS_RELATIVE )
    {
      cdiDecodeDate(date, &year, &month, &day);
      cdiDecodeTime(time, &hour, &minute, &second);

      IEG_P_Year(pdb)     = year;
      IEG_P_Month(pdb)    = month;
      IEG_P_Day(pdb)      = day;
      IEG_P_Hour(pdb)     = hour;
      IEG_P_Minute(pdb)   = minute;

      pdb[15] = 1;
      pdb[16] = 0;
      pdb[17] = 0;
      pdb[18] = 10;
      pdb[36] = 1;
    }

  pdb[5] = 128;
}

static
int calc_resfac(double xfirst, double xlast, double xinc, double yfirst, double ylast, double yinc)
{
  int i, j;
  int iresfac = 1000;
  int ifact;
  int ifacarr[5] = {1000, 10000, 100000, 1000000, 10000000};
  double vals[6] = {xfirst, xlast, xinc, yfirst, ylast, yinc};

  for ( j = 0; j < 5; ++j )
    {
      ifact = ifacarr[j];
      for ( i = 0; i < 6; ++i )
        {
          if ( fabs(vals[i]*ifact - round(vals[i]*ifact)) > FLT_EPSILON ) break;
        }
      if ( i == 6 )
        {
          iresfac = ifact;
          break;
        }
    }

  return (iresfac);
}

static
void iegDefGrid(int *gdb, int gridID)
{
  int gridtype;

  gridtype = gridInqType(gridID);

  if ( gridtype == GRID_GENERIC )
    {
      int xsize, ysize;

      xsize = gridInqXsize(gridID);
      ysize = gridInqYsize(gridID);

      if ( (ysize == 32  || ysize == 48 || ysize == 64 ||
	    ysize == 96  || ysize == 160) &&
	   (xsize == 2*ysize || xsize == 1) )
	{
	  gridtype = GRID_GAUSSIAN;
	  gridChangeType(gridID, gridtype);
	}
      else if ( (xsize == 1 && ysize == 1) || (xsize == 0 && ysize == 0) )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      gridtype = GRID_LONLAT;
    }

  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
    {
      double xfirst = 0, xlast = 0, xinc = 0;
      double yfirst = 0, ylast = 0, yinc = 0;

      int nlon = gridInqXsize(gridID),
        nlat = gridInqYsize(gridID);

      if ( nlon == 0 )
	{
	  nlon = 1;
	}
      else
	{
	  xfirst = gridInqXval(gridID,      0);
	  xlast  = gridInqXval(gridID, nlon-1);
	  xinc   = gridInqXinc(gridID);
	}

      if ( nlat == 0 )
	{
	  nlat = 1;
	}
      else
	{
	  yfirst = gridInqYval(gridID,      0);
	  ylast  = gridInqYval(gridID, nlat-1);
	  yinc   = gridInqYinc(gridID);
	}

      if ( gridtype == GRID_GAUSSIAN )
	IEG_G_GridType(gdb) = 4;
      else if ( gridtype == GRID_LONLAT && gridIsRotated(gridID) )
	IEG_G_GridType(gdb) = 10;
      else
	IEG_G_GridType(gdb) = 0;

      int iresfac = calc_resfac(xfirst, xlast, xinc, yfirst, ylast, yinc);
      double resfac = (double) iresfac;
      if ( iresfac == 1000 ) iresfac = 0;

      IEG_G_ResFac(gdb)   = iresfac;

      IEG_G_NumLon(gdb)   = nlon;
      IEG_G_NumLat(gdb)   = nlat;
      IEG_G_FirstLat(gdb) = (int)lround(yfirst*resfac);
      IEG_G_LastLat(gdb)  = (int)lround(ylast*resfac);
      IEG_G_FirstLon(gdb) = (int)lround(xfirst*resfac);
      IEG_G_LastLon(gdb)  = (int)lround(xlast*resfac);
      IEG_G_LonIncr(gdb)  = (int)lround(xinc*resfac);
      if ( fabs(xinc*resfac - IEG_G_LonIncr(gdb)) > FLT_EPSILON )
	IEG_G_LonIncr(gdb) = 0;

      if ( gridtype == GRID_GAUSSIAN )
	IEG_G_LatIncr(gdb) = nlat/2;
      else
	{
	  IEG_G_LatIncr(gdb) = (int)lround(yinc*resfac);
	  if ( fabs(yinc*resfac - IEG_G_LatIncr(gdb)) > FLT_EPSILON )
	    IEG_G_LatIncr(gdb) = 0;

	  if ( IEG_G_LatIncr(gdb) < 0 ) IEG_G_LatIncr(gdb) = -IEG_G_LatIncr(gdb);
	}

      if ( IEG_G_NumLon(gdb) > 1 && IEG_G_NumLat(gdb) == 1 )
	if ( IEG_G_LonIncr(gdb) != 0 && IEG_G_LatIncr(gdb) == 0 ) IEG_G_LatIncr(gdb) = IEG_G_LonIncr(gdb);

      if ( IEG_G_NumLon(gdb) == 1 && IEG_G_NumLat(gdb) > 1 )
	if ( IEG_G_LonIncr(gdb) == 0 && IEG_G_LatIncr(gdb) != 0 ) IEG_G_LonIncr(gdb) = IEG_G_LatIncr(gdb);

      if ( IEG_G_LatIncr(gdb) == 0 || IEG_G_LonIncr(gdb) == 0 )
	IEG_G_ResFlag(gdb) = 0;
      else
	IEG_G_ResFlag(gdb) = 128;

      if ( gridIsRotated(gridID) )
	{
	  IEG_G_LatSP(gdb) = - (int)lround(gridInqYpole(gridID) * resfac);
	  IEG_G_LonSP(gdb) =   (int)lround((gridInqXpole(gridID) + 180) * resfac);
	  IEG_G_Size(gdb)  = 42;
	}
      else
	{
	  IEG_G_Size(gdb)  = 32;
	}
    }
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }

  IEG_G_ScanFlag(gdb) = 64;
}

static
void iegDefLevel(int *pdb, int *gdb, double *vct, int zaxisID, int levelID)
{
  double level;
  int ilevel, leveltype;
  static int warning = 1;
  static int vct_warning = 1;

  leveltype = zaxisInqType(zaxisID);

  if ( leveltype == ZAXIS_GENERIC )
    {
      Message("Changed zaxis type from %s to %s",
	      zaxisNamePtr(leveltype),
	      zaxisNamePtr(ZAXIS_PRESSURE));
      leveltype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, leveltype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  /*  IEG_G_NumVCP(gdb) = 0; */

  switch (leveltype)
    {
    case ZAXIS_SURFACE:
      {
	IEG_P_LevelType(pdb) = IEG_LTYPE_SURFACE;
	IEG_P_Level1(pdb)    = 0;
	IEG_P_Level2(pdb)    = (int) zaxisInqLevel(zaxisID, levelID);
	break;
      }
    case ZAXIS_HYBRID:
      {
	int vctsize;

	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    IEG_P_LevelType(pdb) = IEG_LTYPE_HYBRID_LAYER;
	    IEG_P_Level1(pdb)    = (int) zaxisInqLbound(zaxisID, levelID);
	    IEG_P_Level2(pdb)    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
	    IEG_P_LevelType(pdb) = IEG_LTYPE_HYBRID;
	    IEG_P_Level1(pdb)    = 0;
	    IEG_P_Level2(pdb)    = (int) zaxisInqLevel(zaxisID, levelID);
	  }

	vctsize = zaxisInqVctSize(zaxisID);
	if ( vctsize == 0 && warning )
	  {
	    Warning("VCT missing. ( code = %d, zaxisID = %d )",
		    IEG_P_Parameter(pdb), zaxisID);
	    warning = 0;
	  }
	if ( vctsize > 100 )
	  {
	    /*	    IEG_G_NumVCP(gdb) = 0; */
	    if ( vct_warning )
	      {
		Warning("VCT size of %d is too large (maximum is 100). Set to 0!", vctsize);
		vct_warning = 0;
	      }
	  }
	else
	  {
	    IEG_G_Size(gdb) += (vctsize*4);
	    memcpy(vct, zaxisInqVctPtr(zaxisID), (size_t)vctsize/2*sizeof(double));
	    memcpy(vct+50, zaxisInqVctPtr(zaxisID)+vctsize/2, (size_t)vctsize/2*sizeof(double));
	  }
	break;
      }
    case ZAXIS_PRESSURE:
      {
	double dum;
	char units[128];

	level = zaxisInqLevel(zaxisID, levelID);
	if ( level < 0 )
	  Warning("pressure level of %f Pa is below 0.", level);

	zaxisInqUnits(zaxisID, units);
	if ( memcmp(units, "hPa", 3) == 0 || memcmp(units, "mb",2 ) == 0 )
	  level = level*100;

	ilevel = (int) level;
	if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
	  {
	    IEG_P_LevelType(pdb) = IEG_LTYPE_99;
	    IEG_P_Level1(pdb)    = 0;
	    IEG_P_Level2(pdb)    = ilevel;
	  }
	else
	  {
	    IEG_P_LevelType(pdb) = IEG_LTYPE_ISOBARIC;
	    IEG_P_Level1(pdb)    = 0;
	    IEG_P_Level2(pdb)    = ilevel/100;
	  }
	break;
      }
    case ZAXIS_HEIGHT:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	IEG_P_LevelType(pdb) = IEG_LTYPE_HEIGHT;
	IEG_P_Level1(pdb)    = 0;
	IEG_P_Level2(pdb)    = ilevel;

	break;
      }
    case ZAXIS_ALTITUDE:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	IEG_P_LevelType(pdb) = IEG_LTYPE_ALTITUDE;
	IEG_P_Level1(pdb)    = 0;
	IEG_P_Level2(pdb)    = ilevel;

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    IEG_P_LevelType(pdb) = IEG_LTYPE_LANDDEPTH_LAYER;
	    IEG_P_Level1(pdb)    = (int) zaxisInqLbound(zaxisID, levelID);
	    IEG_P_Level2(pdb)    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
	    level = zaxisInqLevel(zaxisID, levelID);

	    ilevel = (int) level;
	    IEG_P_LevelType(pdb) = IEG_LTYPE_LANDDEPTH;
	    IEG_P_Level1(pdb)    = 0;
	    IEG_P_Level2(pdb)    = ilevel;
	  }

	break;
      }
    case ZAXIS_DEPTH_BELOW_SEA:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	IEG_P_LevelType(pdb) = IEG_LTYPE_SEADEPTH;
	IEG_P_Level1(pdb)    = 0;
	IEG_P_Level2(pdb)    = ilevel;

	break;
      }
    case ZAXIS_ISENTROPIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	IEG_P_LevelType(pdb) = 113;
	IEG_P_Level1(pdb)    = 0;
	IEG_P_Level2(pdb)    = ilevel;

	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(leveltype));
	break;
      }
    }
}


void iegCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  streamFCopyRecord(streamptr2, streamptr1, "IEG");
}


void iegDefRecord(stream_t *streamptr)
{
  int vlistID;
  int gridID;
  int date, time;
  int datatype;
  int i;
  int param, pdis, pcat, pnum;
  int varID, levelID, tsID, zaxisID;
  int byteorder;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  byteorder = streamptr->byteorder;

  varID   = streamptr->record->varID;
  levelID = streamptr->record->levelID;
  tsID    = streamptr->curTsID;

  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);

  iegInitMem(iegp);
  for ( i = 0; i < 37; i++ ) iegp->ipdb[i] = -1;

  iegp->byteswap = getByteswap(byteorder);

  param =  vlistInqVarParam(vlistID, varID);
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  IEG_P_Parameter(iegp->ipdb) = pnum;
  if ( pdis == 255 ) IEG_P_CodeTable(iegp->ipdb) = pcat;
  date     = streamptr->tsteps[tsID].taxis.vdate;
  time     = streamptr->tsteps[tsID].taxis.vtime;

  iegDefTime(iegp->ipdb, date, time, vlistInqTaxis(vlistID));
  iegDefGrid(iegp->igdb, gridID);
  iegDefLevel(iegp->ipdb, iegp->igdb, iegp->vct, zaxisID, levelID);

  datatype = streamptr->record->prec;

  iegp->dprec = iegDefDatatype(datatype);
}


void iegWriteRecord(stream_t *streamptr, const double *data)
{
  int fileID;
  int i, gridsize, gridID;
  double refval;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  fileID = streamptr->fileID;
  gridID = streamptr->record->gridID;

  gridsize = gridInqSize(gridID);

  refval = data[0];
  for ( i = 1; i < gridsize; i++ )
    if ( data[i] < refval ) refval = data[i];

  iegp->refval = refval;

  iegDefDataDP(iegp, data);

  iegWrite(fileID, iegp);
}

static
void iegAddRecord(stream_t *streamptr, int param, int *pdb, int *gdb, double *vct,
		  size_t recsize, off_t position, int prec)
{
  int leveltype;
  int gridID = UNDEFID;
  int levelID = 0;
  int tsID, recID, varID;
  int datatype;
  int level1, level2;
  int gridtype;
  int lbounds = 0;
  record_t *record;
  grid_t grid;
  int vlistID;

  vlistID = streamptr->vlistID;
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamptr, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  if ( IEG_P_LevelType(pdb) == IEG_LTYPE_HYBRID_LAYER )
    {
      level1 = IEG_P_Level1(pdb);
      level2 = IEG_P_Level2(pdb);
    }
  else
    {
      level1 = IEG_P_Level2(pdb);
      level2 = 0;
      if ( IEG_P_LevelType(pdb) == 100 ) level1 *= 100;
    }

  record->size     = recsize;
  record->position = position;
  record->param    = param;
  record->ilevel   = level1;
  record->ilevel2  = level2;
  record->ltype    = IEG_P_LevelType(pdb);

  if ( IEG_G_GridType(gdb) == 0 || IEG_G_GridType(gdb) == 10 )
    gridtype = GRID_LONLAT;
  else if ( IEG_G_GridType(gdb) == 4 )
    gridtype = GRID_GAUSSIAN;
  else
    gridtype = GRID_GENERIC;

  memset(&grid, 0, sizeof(grid_t));
  grid.type  = gridtype;
  grid.size  = IEG_G_NumLon(gdb)*IEG_G_NumLat(gdb);
  grid.xsize = IEG_G_NumLon(gdb);
  grid.ysize = IEG_G_NumLat(gdb);
  grid.xinc  = 0;
  grid.yinc  = 0;
  grid.xdef  = 0;

  int iresfac = IEG_G_ResFac(gdb);
  if ( iresfac == 0 ) iresfac = 1000;
  double resfac = 1./(double) iresfac;

  /* if ( IEG_G_FirstLon != 0 || IEG_G_LastLon != 0 ) */
  {
    if ( grid.xsize > 1 )
      {
	if ( IEG_G_ResFlag(gdb) && IEG_G_LonIncr(gdb) > 0 )
	  grid.xinc = IEG_G_LonIncr(gdb) * resfac;
	else
	  grid.xinc = (IEG_G_LastLon(gdb) - IEG_G_FirstLon(gdb)) * resfac / (grid.xsize - 1);

	/* correct xinc if necessary */
	if ( IEG_G_FirstLon(gdb) == 0 && IEG_G_LastLon(gdb) > 354000 )
	  {
	    double xinc = 360. / grid.xsize;
            /* FIXME: why not use grid.xinc != xinc as condition? */
	    if ( fabs(grid.xinc-xinc) > 0.0 )
	      {
		grid.xinc = xinc;
		if ( CDI_Debug ) Message("set xinc to %g", grid.xinc);
	      }
	  }
      }
    grid.xfirst = IEG_G_FirstLon(gdb) * resfac;
    grid.xlast  = IEG_G_LastLon(gdb)  * resfac;
    grid.xdef   = 2;
  }
  grid.ydef  = 0;
  /* if ( IEG_G_FirstLat != 0 || IEG_G_LastLat != 0 ) */
  {
    if ( grid.ysize > 1 )
      {
	if ( IEG_G_ResFlag(gdb) && IEG_G_LatIncr(gdb) > 0 )
	  grid.yinc = IEG_G_LatIncr(gdb) * resfac;
	else
	  grid.yinc = (IEG_G_LastLat(gdb) - IEG_G_FirstLat(gdb)) * resfac / (grid.ysize - 1);
      }
    grid.yfirst = IEG_G_FirstLat(gdb) * resfac;
    grid.ylast  = IEG_G_LastLat(gdb)  * resfac;
    grid.ydef   = 2;
  }
  /*
  grid.xfirst= IEG_G_FirstLon(gdb) * resfac;
  grid.xlast = IEG_G_LastLon(gdb) * resfac;
  grid.xinc  = IEG_G_LonIncr(gdb) * resfac;
  grid.xdef  = 2;
  grid.yfirst= IEG_G_FirstLat(gdb) * resfac;
  grid.ylast = IEG_G_LastLat(gdb) * resfac;
  grid.yinc  = IEG_G_LatIncr(gdb) * resfac;
  grid.ydef  = 2;
  */
  grid.xvals = NULL;
  grid.yvals = NULL;

  grid.isRotated = FALSE;
  if ( IEG_G_GridType(gdb) == 10 )
    {
      grid.isRotated = TRUE;
      grid.ypole     = - IEG_G_LatSP(gdb) * resfac;
      grid.xpole     =   IEG_G_LonSP(gdb) * resfac - 180;
      grid.angle     = 0;
    }

  gridID = varDefGrid(vlistID, &grid, 0);

  leveltype = iegGetZaxisType(IEG_P_LevelType(pdb));

  if ( leveltype == ZAXIS_HYBRID )
    {
      double tmpvct[100];
      size_t vctsize = (size_t)IEG_G_NumVCP(gdb);

      for (size_t i = 0; i < vctsize/2; i++ ) tmpvct[i] = vct[i];
      for (size_t i = 0; i < vctsize/2; i++ ) tmpvct[i+vctsize/2] = vct[i+50];

      varDefVCT(vctsize, tmpvct);
    }

  if ( IEG_P_LevelType(pdb) == IEG_LTYPE_HYBRID_LAYER ) lbounds = 1;

  datatype = iegInqDatatype(prec);

  varAddRecord(recID, param, gridID, leveltype, lbounds, level1, level2, 0, 0,
	       datatype, &varID, &levelID, TSTEP_INSTANT, 0, 0, -1,
               NULL, NULL, NULL, NULL, NULL, NULL);

  record->varID   = (short)varID;
  record->levelID = (short)levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d gridID = %d levelID = %d",
	    varID, gridID, levelID);
}

#if 0
static
void iegCmpRecord(stream_t *streamptr, int tsID, int recID, off_t position, int param,
		  int level, int xsize, int ysize)
{
  int varID = 0;
  int levelID = 0;
  record_t *record;

  record  = &streamptr->tsteps[tsID].records[recID];

  if ( param != (*record).param || level != (*record).ilevel )
    Error("inconsistent timestep");

  (*record).position = position;
  /*
  varID   = (*record).varID;
  levelID = (*record).levelID;

  streamptr->vars[varID].level[levelID] = recID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;
  */
  if ( CDI_Debug )
    Message("varID = %d levelID = %d", varID, levelID);
}
#endif

void iegDateTime(int *pdb, int *date, int *time)
{
  int ryear, rmonth, rday, rhour, rminute;

  ryear   = IEG_P_Year(pdb);

  rmonth  = IEG_P_Month(pdb);
  rday    = IEG_P_Day(pdb);

  rhour   = IEG_P_Hour(pdb);
  rminute = IEG_P_Minute(pdb);

  if ( rminute == -1 ) rminute = 0;

  *date = cdiEncodeDate(ryear, rmonth, rday);
  *time = cdiEncodeTime(rhour, rminute, 0);
}

static
void iegScanTimestep1(stream_t *streamptr)
{
  int prec = 0;
  int status;
  int fileID;
  int tabnum;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID;
  size_t recsize;
  off_t recpos;
  int nrecords, nrecs, recID;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  IEGCOMPVAR compVar, compVar0;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  streamptr->curTsID = 0;

  tsID  = tstepsNewEntry(streamptr);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamptr->fileID;

  nrecs = 0;
  while ( TRUE )
    {
      recpos = fileGetPos(fileID);
      status = iegRead(fileID, iegp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 1;
	  break;
	}
      recsize = (size_t)(fileGetPos(fileID) - recpos);

      prec   = iegp->dprec;
      rcode  = IEG_P_Parameter(iegp->ipdb);
      tabnum = IEG_P_CodeTable(iegp->ipdb);
      param  = cdiEncodeParam(rcode, tabnum, 255);

      if ( IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER )
	rlevel = IEG_P_Level1(iegp->ipdb);
      else
	rlevel = IEG_P_Level2(iegp->ipdb);

      if ( IEG_P_LevelType(iegp->ipdb) == 100 ) rlevel *= 100;

      iegDateTime(iegp->ipdb, &vdate, &vtime);

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}
      else
	{
	  datetime.date = vdate;
	  datetime.time = vtime;
	  compVar.param = param;
          compVar.level = rlevel;
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      compVar0.param = streamptr->tsteps[0].records[recID].param;
	      compVar0.level = streamptr->tsteps[0].records[recID].ilevel;

	      if ( memcmp(&compVar0, &compVar, sizeof(IEGCOMPVAR)) == 0 ) break;
	    }
	  if ( recID < nrecs ) break;
	  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) )
	    Warning("Inconsistent verification time for param %d level %d", param, rlevel);
	}

      nrecs++;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", nrecs, (int)recpos, param, rlevel, vdate, vtime);

      iegAddRecord(streamptr, param, iegp->ipdb, iegp->igdb, iegp->vct, recsize, recpos, prec);
    }

  streamptr->rtsteps = 1;

  cdi_generate_vars(streamptr);

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  taxis->type  = TAXIS_ABSOLUTE;
  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  vlist_check_contents(vlistID);

  nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
	(record_t *)xrealloc(streamptr->tsteps[0].records,
                             (size_t)nrecords * sizeof (record_t));
    }

  streamptr->tsteps[0].recIDs = (int *)xmalloc((size_t)nrecords * sizeof (int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for ( varID = 0; varID < streamptr->nvars; varID++ )
	    {
	      vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	    }
	}
    }
}

static
int iegScanTimestep2(stream_t *streamptr)
{
  int status;
  int fileID;
  int tabnum;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  int tsID;
  int varID;
  size_t recsize;
  off_t recpos = 0;
  int nrecords, nrecs, recID, rindex;
  int nextstep;
  taxis_t *taxis;
  int vlistID;
  IEGCOMPVAR compVar, compVar0;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  streamptr->curTsID = 1;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpected timestep %d", tsID+1);

  taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  nrecords = streamptr->tsteps[0].nallrecs;
  streamptr->tsteps[1].recIDs = (int *)xmalloc((size_t)nrecords * sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position =
	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     =
	streamptr->tsteps[0].records[recID].size;
    }

  for ( rindex = 0; rindex <= nrecords; rindex++ )
    {
      recpos = fileGetPos(fileID);
      status = iegRead(fileID, iegp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      recsize = (size_t)(fileGetPos(fileID) - recpos);

      rcode  = IEG_P_Parameter(iegp->ipdb);
      tabnum = IEG_P_CodeTable(iegp->ipdb);
      param  = cdiEncodeParam(rcode, tabnum, 255);

      if ( IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER )
	rlevel = IEG_P_Level1(iegp->ipdb);
      else
	rlevel = IEG_P_Level2(iegp->ipdb);

      if ( IEG_P_LevelType(iegp->ipdb) == 100 ) rlevel *= 100;

      iegDateTime(iegp->ipdb, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  taxis->type  = TAXIS_ABSOLUTE;
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;
	}

      compVar.param = param;
      compVar.level = rlevel;
      nextstep = FALSE;
      for ( recID = 0; recID < nrecords; recID++ )
	{
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(IEGCOMPVAR)) == 0 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  nextstep = TRUE;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = TRUE;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}
	      break;
	    }
	}
      if ( recID == nrecords )
	{
	  char paramstr[32];
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  Warning("param %s level %d not defined at timestep 1", paramstr, rlevel);
	  return (CDI_EUFSTRUCT);
	}

      if ( nextstep ) break;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", rindex+1, (int)recpos, param, rlevel, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      compVar0.param = streamptr->tsteps[tsID].records[recID].param;
      compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

      if ( memcmp(&compVar0, &compVar, sizeof(IEGCOMPVAR)) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
    }

  nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  varID = streamptr->tsteps[tsID].records[recID].varID;
          vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  return (0);
}


int iegInqContents(stream_t *streamptr)
{
  int fileID;
  int status = 0;

  fileID = streamptr->fileID;

  streamptr->curTsID = 0;

  iegScanTimestep1(streamptr);

  if ( streamptr->ntsteps == -1 ) status = iegScanTimestep2(streamptr);

  fileSetPos(fileID, 0, SEEK_SET);

  return (status);
}

static
long iegScanTimestep(stream_t *streamptr)
{
  int status;
  int fileID;
  int tsID;
  int tabnum;
  int param = 0;
  int rcode = 0, rlevel = 0, vdate = 0, vtime = 0;
  size_t recsize = 0;
  off_t recpos = 0;
  int recID;
  taxis_t *taxis;
  int rindex, nrecs = 0;
  IEGCOMPVAR compVar, compVar0;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }

  if ( streamptr->rtsteps == 0 )
    Error("Internal problem! Missing contents.");

  tsID  = streamptr->rtsteps;
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs
        = (int *)xmalloc((size_t)nrecs * sizeof (int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      for ( rindex = 0; rindex <= nrecs; rindex++ )
	{
	  recpos = fileGetPos(fileID);
	  status = iegRead(fileID, iegp);
	  if ( status != 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  recsize = (size_t)(fileGetPos(fileID) - recpos);

	  rcode  = IEG_P_Parameter(iegp->ipdb);
	  tabnum = IEG_P_CodeTable(iegp->ipdb);
	  param  = cdiEncodeParam(rcode, tabnum, 255);

	  if ( IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER )
	    rlevel = IEG_P_Level1(iegp->ipdb);
	  else
	    rlevel = IEG_P_Level2(iegp->ipdb);

	  if ( IEG_P_LevelType(iegp->ipdb) == 100 ) rlevel *= 100;

	  iegDateTime(iegp->ipdb, &vdate, &vtime);

	  // if ( rindex == nrecs ) break; gcc-4.5 internal compiler error
	  if ( rindex == nrecs ) continue;
	  recID = streamptr->tsteps[tsID].recIDs[rindex];

	  if ( rindex == 0 )
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;
	    }

	  compVar.param = param;
          compVar.level = rlevel;
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(IEGCOMPVAR)) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d%8d%4d%8d%8d%6d", rindex, (int)recpos, param, rlevel, vdate, vtime);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamptr);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = 1;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return (streamptr->ntsteps);
}


int iegInqTimestep(stream_t *streamptr, int tsID)
{
  int nrecs;

  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);

  long ntsteps = UNDEFID;
  while ( ( tsID + 1 ) > streamptr->rtsteps && ntsteps == UNDEFID )
    ntsteps = iegScanTimestep(streamptr);

  if ( tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID )
    {
      nrecs = 0;
    }
  else
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return (nrecs);
}


void iegReadVarDP(stream_t *streamptr, int varID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int levID, nlevs, gridID, gridsize;
  off_t recpos, currentfilepos;
  int tsid;
  int recID;
  int i;
  double missval;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  /* NOTE: tiles are not supported here! */
  nlevs    = streamptr->vars[varID].recordTable[0].nlevs;
  missval  = vlistInqVarMissval(vlistID, varID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsid     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  for (levID = 0; levID < nlevs; levID++)
    {
      /* NOTE: tiles are not supported here! */
      recID = streamptr->vars[varID].recordTable[0].recordID[levID];
      recpos = streamptr->tsteps[tsid].records[recID].position;
      fileSetPos(fileID, recpos, SEEK_SET);
      iegRead(fileID, iegp);
      iegInqDataDP(iegp, &data[levID*gridsize]);
    }
  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  for ( i = 0; i < nlevs*gridsize; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}


void iegReadVarSliceDP(stream_t *streamptr, int varID, int levID, double *data, int *nmiss)
{
  int vlistID, fileID;
  int nlevs, gridID, gridsize;
  off_t recpos, currentfilepos;
  int tsid;
  int recID;
  int i;
  double missval;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  /* NOTE: tiles are not supported here! */
  nlevs    = streamptr->vars[varID].recordTable[0].nlevs;
  missval  = vlistInqVarMissval(vlistID, varID);
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  tsid     = streamptr->curTsID;

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d",
	     nlevs, gridID, gridsize);

  currentfilepos = fileGetPos(fileID);

  /* NOTE: tiles are not supported here! */
  recID = streamptr->vars[varID].recordTable[0].recordID[levID];
  recpos = streamptr->tsteps[tsid].records[recID].position;
  fileSetPos(fileID, recpos, SEEK_SET);
  iegRead(fileID, iegp);
  iegInqDataDP(iegp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  for ( i = 0; i < gridsize; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}


void iegWriteVarDP(stream_t *streamptr, int varID, const double *data)
{
  int fileID;
  int levID, nlevs, gridID, gridsize;
  int zaxisID;
  int datatype;
  int tsID;
  int vlistID;
  int i;
  int date, time;
  int param, pdis, pcat, pnum;
  double refval;
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamptr->self, varID);

  iegInitMem(iegp);
  for ( i = 0; i < 37; i++ ) iegp->ipdb[i] = -1;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  tsID     = streamptr->curTsID;
  gridID   = vlistInqVarGrid(vlistID, varID);
  gridsize = gridInqSize(gridID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  nlevs    = zaxisInqSize(zaxisID);

  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);

  param    = vlistInqVarParam(vlistID, varID);
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  IEG_P_Parameter(iegp->ipdb) = pnum;
  if ( pdis == 255 ) IEG_P_CodeTable(iegp->ipdb) = pcat;
  date     = streamptr->tsteps[tsID].taxis.vdate;
  time     = streamptr->tsteps[tsID].taxis.vtime;

  iegDefTime(iegp->ipdb, date, time, vlistInqTaxis(vlistID));
  iegDefGrid(iegp->igdb, gridID);

  datatype = vlistInqVarDatatype(vlistID, varID);

  iegp->dprec = iegDefDatatype(datatype);

  for ( levID = 0;  levID < nlevs; levID++ )
    {
      iegDefLevel(iegp->ipdb, iegp->igdb, iegp->vct, zaxisID, levID);

      refval = data[0];
      for ( i = 1; i < gridsize; i++ )
	if ( data[levID*gridsize+i] < refval ) refval = data[levID*gridsize+i];

      iegp->refval = refval;

      iegDefDataDP(iegp, &data[levID*gridsize]);
      iegWrite(fileID, iegp);
    }
}


void iegWriteVarSliceDP(stream_t *streamptr, int varID, int levID, const double *data)
{
  int fileID;
  int gridID;
  int zaxisID;
  /* double level; */
  int datatype;
  /* int tsID; */
  int vlistID;
  /* int param, date, time, datasize; */
  iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;

  vlistID  = streamptr->vlistID;
  fileID   = streamptr->fileID;
  /* tsID     = streamptr->curTsID; */
  gridID   = vlistInqVarGrid(vlistID, varID);
  zaxisID  = vlistInqVarZaxis(vlistID, varID);
  (void)levID;
  /* level    = zaxisInqLevel(zaxisID, levID); */

  if ( CDI_Debug )
    Message("gridID = %d zaxisID = %d", gridID, zaxisID);

  /* param = vlistInqVarParam(vlistID, varID); */
  /* date = streamptr->tsteps[tsID].taxis.vdate; */
  /* time = streamptr->tsteps[tsID].taxis.vtime; */
  /* datasize = gridInqSize(gridID); */

  datatype = vlistInqVarDatatype(vlistID, varID);

  iegp->dprec = iegDefDatatype(datatype);

  iegDefDataDP(iegp, data);
  iegWrite(fileID, iegp);
}

#endif /* HAVE_LIBIEG */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
