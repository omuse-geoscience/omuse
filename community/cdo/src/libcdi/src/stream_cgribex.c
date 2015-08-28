#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
// #include <float.h>  /* FLT_EPSILON */

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"

#if  defined  (HAVE_LIBCGRIBEX)
#  include "cgribex.h"
#endif

extern int cdiInventoryMode;

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
  int tsteptype;
} compvar_t;


#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetGridType(int *isec2)
{
  int gridtype = GRID_GENERIC;

  switch (ISEC2_GridType)
    {
    case  GRIB1_GTYPE_LATLON:     { if ( ISEC2_Reduced )      break; }
    case  GRIB1_GTYPE_LATLON_ROT: { gridtype = GRID_LONLAT;   break; }
    case  GRIB1_GTYPE_LCC:        { gridtype = GRID_LCC;      break; }
    case  GRIB1_GTYPE_GAUSSIAN:   { if ( ISEC2_Reduced )
	                              gridtype = GRID_GAUSSIAN_REDUCED;
                         	    else
				      gridtype = GRID_GAUSSIAN;
          	                    break;
                                  }
    case  GRIB1_GTYPE_SPECTRAL:   { gridtype = GRID_SPECTRAL; break; }
    case  GRIB1_GTYPE_GME:        { gridtype = GRID_GME;      break; }
    }

  return (gridtype);
}

static
int cgribexGetIsRotated(int *isec2)
{
  int isRotated = 0;

  if ( ISEC2_GridType == GRIB1_GTYPE_LATLON_ROT )
    {
      isRotated = 1;
    }

  return (isRotated);
}

static
int cgribexGetZaxisHasBounds(int grb_ltype)
{
  int lbounds = 0;

  switch (grb_ltype)
    {
    case GRIB1_LTYPE_SIGMA_LAYER:
    case GRIB1_LTYPE_HYBRID_LAYER:
    case GRIB1_LTYPE_LANDDEPTH_LAYER:
      {
	lbounds = 1;
	break;
      }
    }

  return (lbounds);
}

static
int cgribexGetTimeUnit(int *isec1)
{
  int timeunit = TUNIT_HOUR;
  static int lprint = TRUE;

  switch ( ISEC1_TimeUnit )
    {
    case ISEC1_TABLE4_MINUTE:    timeunit = TUNIT_MINUTE;    break;
    case ISEC1_TABLE4_QUARTER:   timeunit = TUNIT_QUARTER;   break;
    case ISEC1_TABLE4_30MINUTES: timeunit = TUNIT_30MINUTES; break;
    case ISEC1_TABLE4_HOUR:      timeunit = TUNIT_HOUR;      break;
    case ISEC1_TABLE4_3HOURS:    timeunit = TUNIT_3HOURS;    break;
    case ISEC1_TABLE4_6HOURS:    timeunit = TUNIT_6HOURS;    break;
    case ISEC1_TABLE4_12HOURS:   timeunit = TUNIT_12HOURS;   break;
    case ISEC1_TABLE4_DAY:       timeunit = TUNIT_DAY;       break;
    default:
      if ( lprint )
	{
	  Message("GRIB time unit %d unsupported!", ISEC1_TimeUnit);
	  lprint = FALSE;
	}
      break;
    }

  return (timeunit);
}

static
int cgribexTimeIsFC(int *isec1)
{
  int isFC = TRUE;

  if ( ISEC1_TimeRange == 10 && ISEC1_TimePeriod1 == 0 && ISEC1_TimePeriod2 == 0 )
    isFC = FALSE;

  return (isFC);
}

static
int cgribexGetTsteptype(int timerange)
{
  int tsteptype = TSTEP_INSTANT;
  static int lprint = TRUE;

  switch ( timerange )
    {
    case  0:  tsteptype = TSTEP_INSTANT;  break;
    case  1:  tsteptype = TSTEP_INSTANT2; break;
    case  2:  tsteptype = TSTEP_RANGE;    break;
    case  3:  tsteptype = TSTEP_AVG;      break;
    case  4:  tsteptype = TSTEP_ACCUM;    break;
    case  5:  tsteptype = TSTEP_DIFF;     break;
    case 10:  tsteptype = TSTEP_INSTANT3; break;
    default:
      if ( lprint )
	{
	  Message("Time range indicator %d unsupported, set to 0!", timerange);
	  lprint = FALSE;
	}
      break;
    }

  return (tsteptype);
}

static
void cgribexGetGrid(stream_t *streamptr, int *isec2, double *fsec2, int *isec4, grid_t *grid, int iret)
{
  int compyinc = TRUE;
  int gridtype = cgribexGetGridType(isec2);

  if ( streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED && iret != -801 )
    {
      int ilat, nlon = 0;
      for ( ilat = 0; ilat < ISEC2_NumLat; ++ilat )
        if ( ISEC2_RowLon(ilat) > nlon ) nlon = ISEC2_RowLon(ilat);
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = nlon;
      ISEC4_NumValues = nlon*ISEC2_NumLat;
      compyinc = FALSE;
    }

  memset(grid, 0, sizeof(grid_t));
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	if ( ISEC4_NumValues != ISEC2_NumLon*ISEC2_NumLat )
	  Error("numberOfPoints (%d) and gridSize (%d) differ!", ISEC4_NumValues, ISEC2_NumLon*ISEC2_NumLat);
	grid->size  = ISEC4_NumValues;
	grid->xsize = ISEC2_NumLon;
	grid->ysize = ISEC2_NumLat;
        if ( gridtype == GRID_GAUSSIAN ) grid->np = ISEC2_NumPar;
	grid->xinc  = 0;
	grid->yinc  = 0;
	grid->xdef  = 0;
	/* if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 ) */
	  {
	    if ( grid->xsize > 1 )
	      {
                int recompinc = TRUE;

                if ( ISEC2_LastLon < ISEC2_FirstLon && ISEC2_LastLon < 0 ) ISEC2_LastLon += 360000;

		if ( ISEC2_ResFlag && ISEC2_LonIncr > 0 )
                  {
                    if ( abs(ISEC2_LastLon - (ISEC2_FirstLon+ISEC2_LonIncr*(grid->xsize-1))) <= 2 )
                      {
                        recompinc = FALSE;
                        grid->xinc = ISEC2_LonIncr * 0.001;
                      }
                  }

		/* recompute xinc if necessary */
                if ( recompinc ) grid->xinc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->xsize-1);

		/* correct xinc if necessary */
		if ( ISEC2_FirstLon == 0 && ISEC2_LastLon > 354000 && ISEC2_LastLon < 360000 )
		  {
		    double xinc = 360. / grid->xsize;

		    if ( fabs(grid->xinc-xinc) > 0.0 )
		      {
			grid->xinc = xinc;
			if ( CDI_Debug ) Message("set xinc to %g", grid->xinc);
		      }
		  }
	      }
	    grid->xfirst = ISEC2_FirstLon * 0.001;
	    grid->xlast  = ISEC2_LastLon  * 0.001;
	    grid->xdef   = 2;
	  }
	grid->ydef  = 0;
	/* if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 ) */
	  {
	    if ( grid->ysize > 1 && compyinc )
	      {
                int recompinc = TRUE;
		if ( ISEC2_ResFlag && ISEC2_LatIncr > 0 )
                  {
                    if ( abs(ISEC2_LastLat - (ISEC2_FirstLat+ISEC2_LatIncr*(grid->ysize-1))) <= 2 )
                      {
                        recompinc = FALSE;
                        grid->yinc = ISEC2_LatIncr * 0.001;
                      }
                  }

		/* recompute yinc if necessary */
                if ( recompinc ) grid->yinc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->ysize - 1);
	      }
	    grid->yfirst = ISEC2_FirstLat * 0.001;
	    grid->ylast  = ISEC2_LastLat  * 0.001;
	    grid->ydef   = 2;
	  }
	break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
        grid->np     = ISEC2_NumPar;
	grid->size   = ISEC4_NumValues;
        grid->rowlon = ISEC2_RowLonPtr;
	grid->ysize  = ISEC2_NumLat;
	grid->xinc   = 0;
	grid->yinc   = 0;
	grid->xdef   = 0;
	/* if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 ) */
	  {
	    if ( grid->xsize > 1 )
	      {
                if ( ISEC2_LastLon < ISEC2_FirstLon && ISEC2_LastLon < 0 ) ISEC2_LastLon += 360000;

		if ( ISEC2_ResFlag && ISEC2_LonIncr > 0 )
		  grid->xinc = ISEC2_LonIncr * 0.001;
		else
		  grid->xinc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->xsize - 1);
	      }
	    grid->xfirst = ISEC2_FirstLon * 0.001;
	    grid->xlast  = ISEC2_LastLon  * 0.001;
	    grid->xdef   = 2;
	  }
	grid->ydef  = 0;
	/* if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 ) */
	  {
	    if ( grid->ysize > 1 )
	      {
		if ( ISEC2_ResFlag && ISEC2_LatIncr > 0 )
		  grid->yinc = ISEC2_LatIncr * 0.001;
		else
		  grid->yinc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->ysize - 1);
	      }
	    grid->yfirst = ISEC2_FirstLat * 0.001;
	    grid->ylast  = ISEC2_LastLat  * 0.001;
	    grid->ydef   = 2;
	  }
	break;
      }
    case GRID_LCC:
      {
	if ( ISEC4_NumValues != ISEC2_NumLon*ISEC2_NumLat )
	  Error("numberOfPoints (%d) and gridSize (%d) differ!",
		ISEC4_NumValues, ISEC2_NumLon*ISEC2_NumLat);

	grid->size  = ISEC4_NumValues;
	grid->xsize = ISEC2_NumLon;
	grid->ysize = ISEC2_NumLat;

	grid->lcc_xinc      = ISEC2_Lambert_dx;
	grid->lcc_yinc      = ISEC2_Lambert_dy;
	grid->lcc_originLon = ISEC2_FirstLon * 0.001;
	grid->lcc_originLat = ISEC2_FirstLat * 0.001;
	grid->lcc_lonParY   = ISEC2_Lambert_Lov * 0.001;
	grid->lcc_lat1      = ISEC2_Lambert_LatS1 * 0.001;
	grid->lcc_lat2      = ISEC2_Lambert_LatS2 * 0.001;
	grid->lcc_projflag  = ISEC2_Lambert_ProjFlag;
	grid->lcc_scanflag  = ISEC2_ScanFlag;

	grid->xdef   = 0;
	grid->ydef   = 0;

	break;
      }
    case GRID_SPECTRAL:
      {
	grid->size  = ISEC4_NumValues;
	grid->trunc = ISEC2_PentaJ;
	if ( ISEC2_RepMode == 2 )
	  grid->lcomplex = 1;
	else
	  grid->lcomplex = 0;

	break;
      }
    case GRID_GME:
      {
	grid->size  = ISEC4_NumValues;
	grid->nd    = ISEC2_GME_ND;
	grid->ni    = ISEC2_GME_NI;
	grid->ni2   = ISEC2_GME_NI2;
	grid->ni3   = ISEC2_GME_NI3;
	break;
      }
    case GRID_GENERIC:
      {
	grid->size  = ISEC4_NumValues;
	grid->xsize = 0;
	grid->ysize = 0;
	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }

  grid->isRotated = FALSE;
  if ( cgribexGetIsRotated(isec2) )
    {
      grid->isRotated = TRUE;
      grid->ypole     = - ISEC2_LatSP*0.001;
      grid->xpole     =   ISEC2_LonSP*0.001 - 180;
      grid->angle     = FSEC2_RotAngle;
    }

  grid->xvals = NULL;
  grid->yvals = NULL;
  grid->type  = gridtype;
}

static
void cgribexAddRecord(stream_t * streamptr, int param, int *isec1, int *isec2, double *fsec2, double *fsec3,
		      int *isec4, long recsize, off_t position, int datatype, int comptype, int lmv, int iret)
{
  int varID;
  int levelID = 0;
  grid_t grid;

  int vlistID = streamptr->vlistID;
  int tsID    = streamptr->curTsID;
  int recID   = recordNewEntry(streamptr, tsID);
  record_t *record  = &streamptr->tsteps[tsID].records[recID];

  int tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);
  int numavg    = ISEC1_AvgNum;

  int level1  = ISEC1_Level1;
  int level2  = ISEC1_Level2;

  /* fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, ISEC1_LevelType); */

  record->size      = (size_t)recsize;
  record->position  = position;
  record->param     = param;
  record->ilevel    = level1;
  record->ilevel2   = level2;
  record->ltype     = ISEC1_LevelType;
  record->tsteptype = tsteptype;

  cgribexGetGrid(streamptr, isec2, fsec2, isec4, &grid, iret);

  int gridID = varDefGrid(vlistID, &grid, 0);

  int zaxistype = grib1ltypeToZaxisType(ISEC1_LevelType);

  if ( zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
    {
      size_t vctsize = (size_t)ISEC2_NumVCP;
      double *vctptr = &fsec2[10];

      varDefVCT(vctsize, vctptr);
    }

  int lbounds = cgribexGetZaxisHasBounds(ISEC1_LevelType);

  if ( datatype > 32 ) datatype = DATATYPE_PACK32;
  if ( datatype <  0 ) datatype = DATATYPE_PACK;

  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, 0, 0,
	       datatype, &varID, &levelID, tsteptype, numavg, ISEC1_LevelType, -1,
               NULL, NULL, NULL, NULL, NULL, NULL);

  (*record).varID   = (short)varID;
  (*record).levelID = (short)levelID;

  varDefCompType(varID, comptype);

  if ( ISEC1_LocalFLag )
    {
      if      ( ISEC1_CenterID == 78  && isec1[36] == 253 ) // DWD local extension
        varDefEnsembleInfo(varID, isec1[54], isec1[53], isec1[52]);
      else if ( ISEC1_CenterID == 252 && isec1[36] ==   1 ) // MPIM local extension
        varDefEnsembleInfo(varID, isec1[38], isec1[39], isec1[37]);
    }

  if ( lmv ) varDefMissval(varID, FSEC3_MissVal);

  if ( varInqInst(varID) == CDI_UNDEFID )
    {
      int center, subcenter, instID;
      center    = ISEC1_CenterID;
      subcenter = ISEC1_SubCenterID;
      instID    = institutInq(center, subcenter, NULL, NULL);
      if ( instID == CDI_UNDEFID )
	instID = institutDef(center, subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if ( varInqModel(varID) == CDI_UNDEFID )
    {
      int modelID;
      modelID = modelInq(varInqInst(varID), ISEC1_ModelID, NULL);
      if ( modelID == CDI_UNDEFID )
	modelID = modelDef(varInqInst(varID), ISEC1_ModelID, NULL);
      varDefModel(varID, modelID);
    }

  if ( varInqTable(varID) == CDI_UNDEFID )
    {
      int tableID;

      tableID = tableInq(varInqModel(varID), ISEC1_CodeTable, NULL);

      if ( tableID == CDI_UNDEFID )
	tableID = tableDef(varInqModel(varID), ISEC1_CodeTable, NULL);
      varDefTable(varID, tableID);
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;
}

static
void MCH_get_undef(int *isec1, double *undef_pds, double *undef_eps)
{
  /* 2010-01-13: Oliver Fuhrer */
  if ( ISEC1_CenterID == 215 ) {
    if (isec1[34] != 0 && isec1[34] != 255) {
      if (isec1[34] & 2) {
        if (isec1[34] & 1) {
          *undef_pds = -0.99*pow(10.0,-isec1[35]);
        } else {
          *undef_pds = +0.99*pow(10.0,-isec1[35]);
        }
        *undef_eps = pow(10.0,-isec1[35]-1);
      } else {
        if (isec1[34] & 1) {
          *undef_pds = -0.99*pow(10.0,+isec1[35]);
        } else {
          *undef_pds = +0.99*pow(10.0,+isec1[35]);
        }
        *undef_eps = pow(10.0,isec1[35]-1);
      }
    }
  }
}

static
void cgribexDecodeHeader(int *isec0, int *isec1, int *isec2, double *fsec2,
			 int *isec3, double *fsec3, int *isec4, double *fsec4,
			 int *gribbuffer, int recsize, int *lmv, int *iret)
{
  int ipunp = 0, iword = 0;

  memset(isec1, 0, 256*sizeof(int));

  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
	   ipunp, (int *) gribbuffer, recsize, &iword, "J", iret);

  *lmv = 0;

  if ( ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255) )
    {
      double undef_pds, undef_eps;

      MCH_get_undef(isec1, &undef_pds, &undef_eps);
      FSEC3_MissVal = undef_pds;
      *lmv = 1;
    }
}

static
compvar_t cgribexVarSet(int param, int level1, int level2, int leveltype, int trange)
{
  compvar_t compVar;
  int tsteptype = cgribexGetTsteptype(trange);

  compVar.param     = param;
  compVar.level1    = level1;
  compVar.level2    = level2;
  compVar.ltype     = leveltype;
  compVar.tsteptype = tsteptype;

  return (compVar);
}

static inline int
cgribexVarCompare(compvar_t compVar, record_t record, int flag)
{
  int tstepDiff = (!((flag == 0) & (((compVar.tsteptype == TSTEP_INSTANT)
                                     & (record.tsteptype == TSTEP_INSTANT3))
                                    |((compVar.tsteptype == TSTEP_INSTANT3)
                                      & (record.tsteptype == TSTEP_INSTANT)))))
    & (compVar.tsteptype != record.tsteptype);
  int rstatus = (compVar.param != record.param)
    |           (compVar.level1 != record.ilevel)
    |           (compVar.level2 != record.ilevel2)
    |           (compVar.ltype != record.ltype)
    |           tstepDiff;
  return (rstatus);
}
#endif

#define gribWarning(text, nrecs, timestep, paramstr, level1, level2) \
            Warning("Record %2d (id=%s lev1=%d lev2=%d) timestep %d: %s", nrecs, paramstr, level1, level2, timestep, text)

#if  defined  (HAVE_LIBCGRIBEX)

static inline void
cgribexScanTsFixNtsteps(stream_t *streamptr, off_t recpos)
{
  if ( streamptr->ntsteps == -1 )
    {
      int tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }
}

static inline void
cgribexScanTsConstAdjust(stream_t *streamptr, taxis_t *taxis)
{
  int vlistID = streamptr->vlistID;
  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for (int varID = 0; varID < streamptr->nvars; varID++ )
	    {
	      vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	    }
	}
    }
}


int cgribexScanTimestep1(stream_t * streamptr)
{
  int *isec0, *isec1, *isec2, *isec3, *isec4;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  size_t buffersize = 0;
  int rstatus;
  int fileID;
  int param = 0;
  int level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  size_t readsize;
  unsigned nrecords, recID;
  int nrecs_scanned = 0;
  int datatype;
  long recsize = 0;
  int warn_time = TRUE;
  int warn_numavg = TRUE;
  int taxisID = -1;
  int rdate = 0, rtime = 0, tunit = 0, fcast = 0;
  taxis_t *taxis;
  int vlistID;
  int comptype;
  long unzipsize;
  char paramstr[32];
  extern int cdiSkipRecords;
  int nskip = cdiSkipRecords;

  streamptr->curTsID = 0;

  isec0 = streamptr->record->sec0;
  isec1 = streamptr->record->sec1;
  isec2 = streamptr->record->sec2;
  isec3 = streamptr->record->sec3;
  isec4 = streamptr->record->sec4;

  tsID  = tstepsNewEntry(streamptr);
  taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  fileID = streamptr->fileID;

  while ( nskip-- > 0 )
    {
      recsize = gribGetSize(fileID);
      if ( recsize == 0 )
	Error("Skipping of %d records failed!", cdiSkipRecords);

      recpos  = fileGetPos(fileID);
      fileSetPos(fileID, (off_t)recsize, SEEK_CUR);
    }

  unsigned nrecs = 0;
  while ( TRUE )
    {
      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);

      if ( recsize == 0 )
	{
	  if ( nrecs == 0 )
	    Error("No GRIB records found!");

	  streamptr->ntsteps = 1;
	  break;
	}
      if ( (size_t)recsize > buffersize )
	{
	  buffersize = (size_t)recsize;
	  gribbuffer = (unsigned char *)xrealloc(gribbuffer, buffersize);
	}

      readsize = (size_t)recsize;
      rstatus = gribRead(fileID, gribbuffer, &readsize);
      if ( rstatus ) break;

      comptype = COMPRESS_NONE;
      if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
	{
	  comptype = COMPRESS_SZIP;
	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	  if ( buffersize < (size_t)unzipsize )
	    {
	      buffersize = (size_t)unzipsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }
	}

      nrecs_scanned++;
      cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			  (int *) gribbuffer, (int)recsize, &lmv, &iret);

      param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      cdiParamToString(param, paramstr, sizeof(paramstr));

      if ( ISEC1_LevelType == 100 ) ISEC1_Level1 *= 100;
      if ( ISEC1_LevelType ==  99 ) ISEC1_LevelType = 100;
      level1   = ISEC1_Level1;
      level2   = ISEC1_Level2;

      gribDateTime(isec1, &vdate, &vtime);

      if ( ISEC4_NumBits > 0 && ISEC4_NumBits <= 32 )
	datatype = ISEC4_NumBits;
      else
        datatype = DATATYPE_PACK;

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	  rdate = gribRefDate(isec1);
	  rtime = gribRefTime(isec1);
	  tunit = cgribexGetTimeUnit(isec1);
	  fcast = cgribexTimeIsFC(isec1);
	}
      else
	{
	  datetime.date  = vdate;
	  datetime.time  = vtime;

	  compvar_t compVar = cgribexVarSet(param, level1, level2, ISEC1_LevelType, ISEC1_TimeRange);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      if ( cgribexVarCompare(compVar, streamptr->tsteps[0].records[recID], 0) == 0 ) break;
	    }

	  if ( cdiInventoryMode == 1 )
	    {
	      if ( recID < nrecs ) break;
	      if ( warn_time )
		if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 )
		  {
                    gribWarning("Inconsistent verification time!", nrecs_scanned, tsID+1, paramstr, level1, level2);
		    warn_time = FALSE;
		  }
	    }
	  else
	    {
	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

	      if ( recID < nrecs )
		{
		  gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, paramstr, level1, level2);
		  continue;
		}
	    }
	}

      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum) )
	    {
	      Warning("Changing numavg from %d to %d not supported!", taxis->numavg, ISEC1_AvgNum);
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}

      nrecs++;

      if ( CDI_Debug )
	Message("Read record %2d (id=%s lev1=%d lev2=%d) %8d %6d", nrecs_scanned, paramstr, level1, level2, vdate, vtime);

      cgribexAddRecord(streamptr, param, isec1, isec2, fsec2, fsec3,
		       isec4, recsize, recpos, datatype, comptype, lmv, iret);
    }

  streamptr->rtsteps = 1;

  if ( nrecs == 0 ) return (CDI_EUFSTRUCT);

  cdi_generate_vars(streamptr);

  if ( fcast )
    {
      taxisID = taxisCreate(TAXIS_RELATIVE);
      taxis->type  = TAXIS_RELATIVE;
      taxis->rdate = rdate;
      taxis->rtime = rtime;
      taxis->unit  = tunit;
    }
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      taxis->type  = TAXIS_ABSOLUTE;
      taxis->unit  = tunit;
    }

  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  nrecords = (unsigned)streamptr->tsteps[0].nallrecs;
  if ( nrecords < (unsigned)streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = (int)nrecords;
      streamptr->tsteps[0].records =
      (record_t *) realloc(streamptr->tsteps[0].records, nrecords*sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) malloc(nrecords*sizeof(int));
  streamptr->tsteps[0].nrecs = (int)nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = (int)recID;

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = (size_t)buffersize;

  cgribexScanTsFixNtsteps(streamptr, recpos);
  cgribexScanTsConstAdjust(streamptr, taxis);

  return (0);
}


int cgribexScanTimestep2(stream_t * streamptr)
{
  int rstatus = 0;
  int *isec0, *isec1, *isec2, *isec3, *isec4;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  size_t buffersize = 0;
  int fileID;
  int param = 0;
  int level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int tsID;
  int varID, gridID;
  size_t readsize;
  int nrecords, nrecs, recID, rindex;
  int nrecs_scanned = 0;
  long recsize = 0;
  int warn_numavg = TRUE;
  int tsteptype;
  int taxisID = -1;
  taxis_t *taxis;
  int vlistID;
  long unzipsize;
  char paramstr[32];

  streamptr->curTsID = 1;

  isec0 = streamptr->record->sec0;
  isec1 = streamptr->record->sec1;
  isec2 = streamptr->record->sec2;
  isec3 = streamptr->record->sec3;
  isec4 = streamptr->record->sec4;

  fileID  = streamptr->fileID;
  vlistID = streamptr->vlistID;
  taxisID = vlistInqTaxis(vlistID);

  gribbuffer = (unsigned char *) streamptr->record->buffer;
  buffersize = streamptr->record->buffersize;

  tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpected timestep %d", tsID+1);

  taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  nrecords = streamptr->tsteps[tsID].nallrecs;
  if ( nrecords ) streamptr->tsteps[1].recIDs = (int *)xmalloc((size_t)nrecords * sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position =	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     =	streamptr->tsteps[0].records[recID].size;
    }

  nrecs_scanned = nrecords;
  rindex = 0;
  while ( TRUE )
    {
      if ( rindex > nrecords ) break;

      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);
      if ( recsize == 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      if ( (size_t)recsize > buffersize )
	{
	  buffersize = (size_t)recsize;
	  gribbuffer = (unsigned char *)xrealloc(gribbuffer, buffersize);
	}

      readsize = (size_t)recsize;
      rstatus = gribRead(fileID, gribbuffer, &readsize);
      if ( rstatus ) break;

      if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
	{
	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	  if ( buffersize < (size_t)unzipsize )
	    {
	      buffersize = (size_t)unzipsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }
	}

      cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			  (int *) gribbuffer, (int)recsize, &lmv, &iret);

      nrecs_scanned++;

      param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      cdiParamToString(param, paramstr, sizeof(paramstr));

      if ( ISEC1_LevelType == 100 ) ISEC1_Level1 *= 100;
      if ( ISEC1_LevelType ==  99 ) ISEC1_LevelType = 100;
      level1    = ISEC1_Level1;
      level2    = ISEC1_Level2;

      gribDateTime(isec1, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;
	      taxis->rdate = gribRefDate(isec1);
	      taxis->rtime = gribRefTime(isec1);
	    }
	  else
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	    }
	  taxis->unit  = cgribexGetTimeUnit(isec1);
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;

	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}

      tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);

      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg &&
        	(taxis->numavg != ISEC1_AvgNum) )
	    {
	  /*
	      Warning("Changing numavg from %d to %d not supported!", taxis->numavg, ISEC1_AvgNum);
	  */
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}

      datetime.date  = vdate;
      datetime.time  = vtime;

      compvar_t compVar = cgribexVarSet(param, level1, level2, ISEC1_LevelType, ISEC1_TimeRange);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;
	}

      if ( recID == nrecords )
	{
	  gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, paramstr, level1, level2);
	  return (CDI_EUFSTRUCT);
	}

      if ( cdiInventoryMode == 1 )
	{
	  if ( streamptr->tsteps[tsID].records[recID].used )
	    {
	      break;
	    }
	  else
	    {
	      streamptr->tsteps[tsID].records[recID].used = TRUE;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	}
      else
	{
	  if ( streamptr->tsteps[tsID].records[recID].used )
	    {
	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

              gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, paramstr, level1, level2);
	      continue;
	    }
	  else
	    {
	      streamptr->tsteps[tsID].records[recID].used = TRUE;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	}

      if ( CDI_Debug )
	Message("Read record %2d (id=%s lev1=%d lev2=%d) %8d %6d", nrecs_scanned, paramstr, level1, level2, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = (size_t)recsize;

      if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, level1);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
      varID = streamptr->tsteps[tsID].records[recID].varID;
      gridID = vlistInqVarGrid(vlistID, varID);
      if ( gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT )
	{
	  if ( IS_NOT_EQUAL(gridInqXval(gridID, 0),ISEC2_FirstLon*0.001) ||
	       IS_NOT_EQUAL(gridInqYval(gridID, 0),ISEC2_FirstLat*0.001) )
	    gridChangeType(gridID, GRID_TRAJECTORY);
	}

      if ( tsteptype != vlistInqVarTsteptype(vlistID, varID) )
	vlistDefVarTsteptype(vlistID, varID, tsteptype);

      rindex++;
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

  cgribexScanTsFixNtsteps(streamptr, recpos);

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  return (rstatus);
}
#endif


#if  defined  (HAVE_LIBCGRIBEX)
int cgribexScanTimestep(stream_t * streamptr)
{
  int rstatus = 0;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0, iret = 0;
  long recsize = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer;
  size_t buffersize = 0;
  int fileID;
  int param = 0;
  int level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0;
  int vrecID, recID;
  int warn_numavg = TRUE;
  size_t readsize;
  int taxisID = -1;
  int rindex, nrecs = 0;
  int nrecs_scanned;
  long unzipsize;
  char paramstr[32];

  /*
  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }
  */
  int *isec0 = streamptr->record->sec0;
  int *isec1 = streamptr->record->sec1;
  int *isec2 = streamptr->record->sec2;
  int *isec3 = streamptr->record->sec3;
  int *isec4 = streamptr->record->sec4;

  int tsID  = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      gribbuffer = (unsigned char *) streamptr->record->buffer;
      buffersize = streamptr->record->buffersize;

      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *)xmalloc((size_t)nrecs * sizeof (int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      nrecs_scanned = streamptr->tsteps[0].nallrecs + streamptr->tsteps[1].nrecs*(tsID-1);
      rindex = 0;
      while ( TRUE )
	{
	  if ( rindex > nrecs ) break;

	  recsize = gribGetSize(fileID);
	  recpos  = fileGetPos(fileID);
	  if ( recsize == 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  if ( recsize > 0 && (size_t)recsize > buffersize )
	    {
	      buffersize = (size_t)recsize;
	      gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
	    }

	  if ( rindex >= nrecs ) break;

	  readsize = (size_t)recsize;
	  rstatus = gribRead(fileID, gribbuffer, &readsize);
	  if ( rstatus )
	    {
	      Warning("Inconsistent timestep %d (GRIB record %d/%d)!", tsID+1, rindex+1,
                      streamptr->tsteps[tsID].recordSize);
	      break;
	    }

	  if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
	    {
	      unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	      if ( buffersize < (size_t)unzipsize )
		{
		  buffersize = (size_t)unzipsize;
		  gribbuffer = (unsigned char *) realloc(gribbuffer, buffersize);
		}
	    }

	  cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			      (int *) gribbuffer, (int)recsize, &lmv, &iret);

          nrecs_scanned++;

	  param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
          cdiParamToString(param, paramstr, sizeof(paramstr));

	  if ( ISEC1_LevelType == 100 ) ISEC1_Level1 *= 100;
	  if ( ISEC1_LevelType ==  99 ) ISEC1_LevelType = 100;
	  level1   = ISEC1_Level1;
	  level2   = ISEC1_Level2;

	  gribDateTime(isec1, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
              int vlistID = streamptr->vlistID;
	      taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;
		  taxis->rdate = gribRefDate(isec1);
		  taxis->rtime = gribRefTime(isec1);
		}
	      else
		{
		  taxis->type  = TAXIS_ABSOLUTE;
		}
	      taxis->unit  = cgribexGetTimeUnit(isec1);
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;

	      datetime0.date = vdate;
	      datetime0.time = vtime;
	    }

	  if ( ISEC1_AvgNum )
	    {
	      if (  taxis->numavg && warn_numavg &&
		   (taxis->numavg != ISEC1_AvgNum) )
		{
	      /*
	          Warning("Changing numavg from %d to %d not supported!", streamptr->tsteps[tsID].taxis.numavg, ISEC1_AvgNum);
	      */
		  warn_numavg = FALSE;
		}
	      else
		{
		  taxis->numavg = ISEC1_AvgNum;
		}
	    }

	  datetime.date  = vdate;
	  datetime.time  = vtime;

	  compvar_t compVar = cgribexVarSet(param, level1, level2, ISEC1_LevelType, ISEC1_TimeRange);

	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;
	    }

	  if ( vrecID == nrecs )
	    {
	      gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, paramstr, level1, level2);

	      if ( cdiInventoryMode == 1 )
		return (CDI_EUFSTRUCT);
	      else
		continue;
	    }

	  if ( cdiInventoryMode == 1 )
	    {
	      streamptr->tsteps[tsID].records[recID].used = TRUE;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	  else
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  char paramstr[32];
		  cdiParamToString(param, paramstr, sizeof(paramstr));

		  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

		  if ( CDI_Debug )
                    gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, paramstr, level1, level2);

		  continue;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = TRUE;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}
	    }

	  if ( CDI_Debug )
            Message("Read record %2d (id=%s lev1=%d lev2=%d) %8d %6d", nrecs_scanned, paramstr, level1, level2, vdate, vtime);

	  if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, level1);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = (size_t)recsize;

	  rindex++;
	}

      for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	{
	  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
	  if ( ! streamptr->tsteps[tsID].records[recID].used ) break;
	}

      if ( vrecID < nrecs )
	{
	  cdiParamToString(streamptr->tsteps[tsID].records[recID].param, paramstr, sizeof(paramstr));
	  gribWarning("Paramameter not found!", nrecs_scanned, tsID+1, paramstr,
                      streamptr->tsteps[tsID].records[recID].ilevel, streamptr->tsteps[tsID].records[recID].ilevel2);
	  return (CDI_EUFSTRUCT);
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

      streamptr->record->buffer     = gribbuffer;
      streamptr->record->buffersize = buffersize;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  rstatus = (int)streamptr->ntsteps;

  return (rstatus);
}
#endif

#ifdef gribWarning
#undef gribWarning
#endif

#if  defined  (HAVE_LIBCGRIBEX)
int cgribexDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, double missval)
{
  int status = 0;
  int iret = 0, iword = 0;
  int isec0[2], isec1[4096], isec2[4096], isec3[2], isec4[512];
  double fsec2[512], fsec3[2];
  char hoper[2];

  if ( unreduced ) strcpy(hoper, "R");
  else             strcpy(hoper, "D");

  FSEC3_MissVal = missval;

  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, data,
	   gridsize, (int *) gribbuffer, gribsize, &iword, hoper, &iret);

  if ( ISEC1_Sec2Or3Flag & 64 )
    *nmiss = ISEC4_NumValues - ISEC4_NumNonMissValues;
  else
    *nmiss = 0;

  if ( ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255) )
    {
      double undef_pds, undef_eps;
      int i;

      MCH_get_undef(isec1, &undef_pds, &undef_eps);

      *nmiss = 0;
      for ( i = 0; i < gridsize; i++ )
        if ( (fabs(data[i]-undef_pds) < undef_eps) || IS_EQUAL(data[i],FSEC3_MissVal) ) {
          data[i] = missval;
          (*nmiss)++;
        }
    }

  return (status);
}
#endif


#if  defined  (HAVE_LIBCGRIBEX)
static
void cgribexDefInstitut(int *isec1, int vlistID, int varID)
{
  int instID;

  if ( vlistInqInstitut(vlistID) != CDI_UNDEFID )
    instID = vlistInqInstitut(vlistID);
  else
    instID = vlistInqVarInstitut(vlistID, varID);

  if ( instID != CDI_UNDEFID )
    {
      int center, subcenter;
      center    = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);
      ISEC1_CenterID    = center;
      ISEC1_SubCenterID = subcenter;
    }
}

static
void cgribexDefModel(int *isec1, int vlistID, int varID)
{
  int modelID;

  if ( vlistInqModel(vlistID) != CDI_UNDEFID )
    modelID = vlistInqModel(vlistID);
  else
    modelID = vlistInqVarModel(vlistID, varID);

  if ( modelID != CDI_UNDEFID )
    ISEC1_ModelID = modelInqGribID(modelID);
}

static
void cgribexDefParam(int *isec1, int param)
{
  int pdis, pcat, pnum;

  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  if ( pnum < 0 ) pnum = -pnum;

  static bool lwarn_pdis = true;
  if ( pdis != 255 && lwarn_pdis )
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
      lwarn_pdis = false;
    }

  static bool lwarn_pnum = true;
  if ( pnum > 255 && lwarn_pnum )
    {
      Warning("Parameter number %d out of range (1-255), set to %d!", pnum, pnum%256);
      lwarn_pnum = false;
      pnum = pnum%256;
    }

  ISEC1_CodeTable = pcat;
  ISEC1_Parameter = pnum;
}

static
int cgribexDefTimerange(int tsteptype, int factor, int calendar,
			int rdate, int rtime, int vdate, int vtime, int *pip1, int *pip2)
{
  int timerange = -1;
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;
  int ip, ip1 = 0, ip2 = 0;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday1, &secofday1);

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

  (void) julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

  if ( !(int) fmod(days*86400.0 + secs, factor) )
    {
      ip = (int) ((days*86400.0 + secs)/factor);

      switch ( tsteptype )
	{
	case TSTEP_INSTANT:  timerange =  0; ip1 = ip; ip2 = 0;  break;
	case TSTEP_INSTANT2: timerange =  1; ip1 = 0;  ip2 = 0;  break;
	case TSTEP_RANGE:    timerange =  2; ip1 = 0;  ip2 = ip; break;
	case TSTEP_AVG:      timerange =  3; ip1 = 0;  ip2 = ip; break;
	case TSTEP_ACCUM:    timerange =  4; ip1 = 0;  ip2 = ip; break;
	case TSTEP_DIFF:     timerange =  5; ip1 = 0;  ip2 = ip; break;
	case TSTEP_INSTANT3:
	default:             timerange = 10; ip1 = ip/256; ip2 = ip%256;  break;
	}
    }

  *pip1 = ip1;
  *pip2 = ip2;

  return (timerange);
}

static
int cgribexDefDateTime(int *isec1, int timeunit, int date, int time)
{
  int year, month, day, hour, minute, second;
  int century = 0;
  int factor = 1;

  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  century =  year / 100;

  ISEC1_Year = year - century*100;

  if ( year < 0 )
    {
      century = -century;
      ISEC1_Year = -ISEC1_Year;
    }

  if ( ISEC1_Year == 0 )
    {
      century -= 1;
      ISEC1_Year = 100;
    }

  century += 1;
  if ( year < 0 ) century = -century;

  ISEC1_Month  = month;
  ISEC1_Day    = day;
  ISEC1_Hour   = hour;
  ISEC1_Minute = minute;

  ISEC1_Century = century;

  switch (timeunit)
    {
    case TUNIT_MINUTE:    factor =    60; ISEC1_TimeUnit = ISEC1_TABLE4_MINUTE;    break;
    case TUNIT_QUARTER:   factor =   900; ISEC1_TimeUnit = ISEC1_TABLE4_QUARTER;   break;
    case TUNIT_30MINUTES: factor =  1800; ISEC1_TimeUnit = ISEC1_TABLE4_30MINUTES; break;
    case TUNIT_HOUR:      factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;      break;
    case TUNIT_3HOURS:    factor = 10800; ISEC1_TimeUnit = ISEC1_TABLE4_3HOURS;    break;
    case TUNIT_6HOURS:    factor = 21600; ISEC1_TimeUnit = ISEC1_TABLE4_6HOURS;    break;
    case TUNIT_12HOURS:   factor = 43200; ISEC1_TimeUnit = ISEC1_TABLE4_12HOURS;   break;
    case TUNIT_DAY:       factor = 86400; ISEC1_TimeUnit = ISEC1_TABLE4_DAY;       break;
    default:              factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;      break;
    }

  return (factor);
}

static
void cgribexDefTime(int *isec1, int vdate, int vtime, int tsteptype, int numavg, int taxisID)
{
  int timetype = TAXIS_ABSOLUTE;
  int timerange = 0;
  int timeunit = TUNIT_HOUR;

  if ( taxisID != -1 )
    {
      timetype = taxisInqType(taxisID);
      timeunit = taxisInqTunit(taxisID);
    }

  if ( timetype == TAXIS_RELATIVE )
    {
      int factor = 1;
      int rdate, rtime;
      int ip1 = 0, ip2 = 0;
      int calendar;

      calendar = taxisInqCalendar(taxisID);
      rdate    = taxisInqRdate(taxisID);
      rtime    = taxisInqRtime(taxisID);

      factor = cgribexDefDateTime(isec1, timeunit, rdate, rtime);

      timerange = cgribexDefTimerange(tsteptype, factor, calendar,
				      rdate, rtime, vdate, vtime, &ip1, &ip2);

      if ( timerange == -1 || timerange == 3 )
	{
	  timetype = TAXIS_ABSOLUTE;
	}
      /*
      else if ( timerange == 10 )
	{
	  if ( ip1 < 0 || ip1 > 0xFFFF ) timetype = TAXIS_ABSOLUTE;
	  if ( ip2 < 0 || ip2 > 0xFFFF ) timetype = TAXIS_ABSOLUTE;
	}
      */
      else
	{
	  if ( ip1 < 0 || ip1 > 0xFF   ) timetype = TAXIS_ABSOLUTE;
	  if ( ip2 < 0 || ip2 > 0xFF   ) timetype = TAXIS_ABSOLUTE;
	}

      ISEC1_TimeRange   = timerange;
      ISEC1_TimePeriod1 = ip1;
      ISEC1_TimePeriod2 = ip2;
    }

  if ( timetype == TAXIS_ABSOLUTE )
    {
      (void) cgribexDefDateTime(isec1, timeunit, vdate, vtime);

      /*
      if ( numavg > 0 )
	ISEC1_TimeRange = 0;
      else
      */
      if ( ISEC1_TimeRange != 3 )
	ISEC1_TimeRange   = 10;

      ISEC1_TimePeriod1 = 0;
      ISEC1_TimePeriod2 = 0;
    }

  ISEC1_AvgNum         = numavg;
  ISEC1_AvgMiss        = 0;
  ISEC1_DecScaleFactor = 0;
}

static
void cgribexDefGrid(int *isec1, int *isec2, double *fsec2, int *isec4, int gridID)
{
  int gridtype;
  bool lcurvi = false;
  static bool lwarning = true;

  memset(isec2, 0, 16*sizeof(int));

  ISEC1_Sec2Or3Flag = 128;

  gridtype = gridInqType(gridID);

  ISEC1_GridDefinition = 255;

  if ( gridtype == GRID_GENERIC )
    {
      int xsize, ysize, gridsize;

      gridsize = gridInqSize(gridID);
      xsize = gridInqXsize(gridID);
      ysize = gridInqYsize(gridID);

      if ( (ysize ==  32 || ysize ==  48 || ysize ==  64 ||
	    ysize ==  96 || ysize == 160 || ysize == 192 ||
	    ysize == 240 || ysize == 320 || ysize == 384 ||
	    ysize == 480 || ysize == 768 ) &&
	   (xsize == 2*ysize || xsize == 1) )
	{
	  gridtype = GRID_GAUSSIAN;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridsize == 1 )
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
      if ( lwarning && gridInqSize(gridID) > 1 )
	{
	  lwarning = false;
	  Warning("Curvilinear grids are unsupported in GRIB1! Created wrong GDS!");
	}
      gridtype = GRID_LONLAT;
      lcurvi = true;
    }

  ISEC2_Reduced  = FALSE;
  ISEC2_ScanFlag = 0;

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
	int nlon = 0, nlat;
	double xfirst = 0, xlast = 0, xinc = 0;
	double yfirst = 0, ylast = 0, yinc = 0;

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          ISEC2_GridType = GRIB1_GTYPE_GAUSSIAN;
        else if ( gridtype == GRID_LONLAT && gridIsRotated(gridID) )
	  ISEC2_GridType = GRIB1_GTYPE_LATLON_ROT;
	else
	  ISEC2_GridType = GRIB1_GTYPE_LATLON;

	nlon = gridInqXsize(gridID);
	nlat = gridInqYsize(gridID);

	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    ISEC2_Reduced = TRUE;
	    nlon = 0;
	    gridInqRowlon(gridID, ISEC2_RowLonPtr);
	  }
	else
	  {
	    if ( nlon == 0 )
	      {
		nlon = 1;
	      }
	    else
	      {
		xfirst = gridInqXval(gridID,      0);
		if ( lcurvi )
		  xlast  = gridInqXval(gridID, nlon*nlat-1);
		else
		  xlast  = gridInqXval(gridID, nlon-1);
		xinc   = gridInqXinc(gridID);
	      }
	  }

	if ( nlat == 0 )
	  {
	    nlat = 1;
	  }
	else
	  {
	    yfirst = gridInqYval(gridID,      0);
	    if ( lcurvi )
	      ylast  = gridInqYval(gridID, nlon*nlat-1);
	    else
	      ylast  = gridInqYval(gridID, nlat-1);
	    yinc   = gridInqYinc(gridID);
	    if ( yinc < 0 ) yinc = -yinc;
	  }

	ISEC2_NumLon   = nlon;
	ISEC2_NumLat   = nlat;
	ISEC2_FirstLat = (int)lround(yfirst*1000);
	ISEC2_LastLat  = (int)lround(ylast*1000);
	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    ISEC2_FirstLon = 0;
	    ISEC2_LastLon  = (int)lround(1000*(360.-360./(nlat*2)));
	    ISEC2_LonIncr  = (int)lround(1000*360./(nlat*2));
	  }
	else
	  {
	    ISEC2_FirstLon = (int)lround(xfirst*1000);
	    ISEC2_LastLon  = (int)lround(xlast*1000);
	    ISEC2_LonIncr  = (int)lround(xinc*1000);
	  }

	// if ( fabs(xinc*1000 - ISEC2_LonIncr) > FLT_EPSILON ) ISEC2_LonIncr = 0;

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          {
            int np = gridInqNP(gridID);
            if ( np == 0 ) np = nlat/2;
            ISEC2_NumPar = np;
          }
	else
	  {
	    ISEC2_LatIncr = (int)lround(yinc*1000);
	    // if ( fabs(yinc*1000 - ISEC2_LatIncr) > FLT_EPSILON ) ISEC2_LatIncr = 0;

	    if ( ISEC2_LatIncr < 0 ) ISEC2_LatIncr = -ISEC2_LatIncr;
	  }

	if ( ISEC2_NumLon > 1 && ISEC2_NumLat == 1 )
	  if ( ISEC2_LonIncr != 0 && ISEC2_LatIncr == 0 ) ISEC2_LatIncr = ISEC2_LonIncr;

	if ( ISEC2_NumLon == 1 && ISEC2_NumLat > 1 )
	  if ( ISEC2_LonIncr == 0 && ISEC2_LatIncr != 0 ) ISEC2_LonIncr = ISEC2_LatIncr;

	if ( ISEC2_LatIncr == 0 || ISEC2_LonIncr == 0 )
	  ISEC2_ResFlag = 0;
	else
	  ISEC2_ResFlag = 128;

	if ( gridIsRotated(gridID) )
	  {
	    ISEC2_LatSP = - (int)lround(gridInqYpole(gridID) * 1000);
	    ISEC2_LonSP =   (int)lround((gridInqXpole(gridID) + 180) * 1000);
            FSEC2_RotAngle = gridInqAngle(gridID);
	  }

	/* East -> West */
	if ( ISEC2_LastLon < ISEC2_FirstLon ) ISEC2_ScanFlag += 128;

	/* South -> North */
	if ( ISEC2_LastLat > ISEC2_FirstLat ) ISEC2_ScanFlag += 64;

	break;
      }
    case GRID_LCC:
      {
	double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	int xsize, ysize;
	int projflag, scanflag;

	xsize = gridInqXsize(gridID);
	ysize = gridInqYsize(gridID);

	gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		   &projflag, &scanflag);

	ISEC2_GridType = GRIB1_GTYPE_LCC;
	ISEC2_NumLon   = xsize;
	ISEC2_NumLat   = ysize;
	ISEC2_FirstLon = (int)lround(originLon * 1000);
	ISEC2_FirstLat = (int)lround(originLat * 1000);
	ISEC2_Lambert_Lov    = (int)lround(lonParY * 1000);
	ISEC2_Lambert_LatS1  = (int)lround(lat1 * 1000);
	ISEC2_Lambert_LatS2  = (int)lround(lat2 * 1000);
	ISEC2_Lambert_dx     = (int)lround(xincm);
	ISEC2_Lambert_dy     = (int)lround(yincm);
	ISEC2_Lambert_LatSP  = 0;
	ISEC2_Lambert_LatSP  = 0;
	ISEC2_Lambert_ProjFlag = projflag;
	ISEC2_ScanFlag = scanflag;

	break;
      }
    case GRID_SPECTRAL:
      {
	ISEC2_GridType = GRIB1_GTYPE_SPECTRAL;
	ISEC2_PentaJ   = gridInqTrunc(gridID);
	ISEC2_PentaK   = ISEC2_PentaJ;
	ISEC2_PentaM   = ISEC2_PentaJ;
	ISEC2_RepType  = 1;
	isec4[2]       = 128;
	if ( gridInqComplexPacking(gridID) && ISEC2_PentaJ >= 21 )
	  {
	    ISEC2_RepMode  = 2;
	    isec4[3]       = 64;
	    isec4[16]      = 0;
	    isec4[17]      = 20;
	    isec4[18]      = 20;
	    isec4[19]      = 20;
	  }
	else
	  {
	    ISEC2_RepMode  = 1;
	    isec4[3]       = 0;
	  }
	break;
      }
    case GRID_GME:
      {
	ISEC2_GridType   = GRIB1_GTYPE_GME;
	ISEC2_GME_ND     = gridInqGMEnd(gridID);
	ISEC2_GME_NI     = gridInqGMEni(gridID);
	ISEC2_GME_NI2    = gridInqGMEni2(gridID);
	ISEC2_GME_NI3    = gridInqGMEni3(gridID);
	ISEC2_GME_AFlag  = 0;
	ISEC2_GME_LatPP  = 90000;
	ISEC2_GME_LonPP  = 0;
	ISEC2_GME_LonMPL = 0;
	ISEC2_GME_BFlag  = 0;
	break;
      }
    default:
      {
	Warning("The CGRIBEX library can not store fields on the used grid!");
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }
}

static
void cgribexDefLevel(int *isec1, int *isec2, double *fsec2, int zaxisID, int levelID)
{
  double level;
  int ilevel, zaxistype, ltype;
  static bool lwarning = true;
  static bool lwarning_vct = true;

  zaxistype = zaxisInqType(zaxisID);
  ltype = zaxisInqLtype(zaxisID);

  if ( zaxistype == ZAXIS_GENERIC && ltype == 0 )
    {
      Message("Changed zaxis type from %s to %s",
	      zaxisNamePtr(zaxistype),
	      zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  ISEC2_NumVCP = 0;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
      {
	ISEC1_LevelType = GRIB1_LTYPE_SURFACE;
	ISEC1_Level1    = (int) zaxisInqLevel(zaxisID, levelID);
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_CLOUD_BASE:
      {
	ISEC1_LevelType = GRIB1_LTYPE_CLOUD_BASE;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_CLOUD_TOP:
      {
	ISEC1_LevelType = GRIB1_LTYPE_CLOUD_TOP;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_ISOTHERM_ZERO:
      {
	ISEC1_LevelType = GRIB1_LTYPE_ISOTHERM0;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_TOA:
      {
	ISEC1_LevelType = GRIB1_LTYPE_TOA;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_SEA_BOTTOM:
      {
	ISEC1_LevelType = GRIB1_LTYPE_SEA_BOTTOM;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_ATMOSPHERE:
      {
	ISEC1_LevelType = GRIB1_LTYPE_ATMOSPHERE;
	ISEC1_Level1    = 0;
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_MEANSEA:
      {
	ISEC1_LevelType = GRIB1_LTYPE_MEANSEA;
	ISEC1_Level1    = (int) zaxisInqLevel(zaxisID, levelID);
	ISEC1_Level2    = 0;
	break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
	int vctsize;

	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_HYBRID_LAYER;
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_HYBRID;
	    ISEC1_Level1    = (int) zaxisInqLevel(zaxisID, levelID);
	    ISEC1_Level2    = 0;
	  }

	vctsize = zaxisInqVctSize(zaxisID);
	if ( vctsize == 0 && lwarning )
	  {
	    Warning("VCT missing. ( param = %d, zaxisID = %d )", ISEC1_Parameter, zaxisID);
	    lwarning = false;
	  }
	if ( vctsize > 255 )
	  {
	    ISEC2_NumVCP = 0;
	    if ( lwarning_vct )
	      {
		Warning("VCT size of %d is too large (maximum is 255). Set to 0!", vctsize);
		lwarning_vct = false;
	      }
	  }
	else
	  {
	    ISEC2_NumVCP = vctsize;
	    zaxisInqVct(zaxisID, &fsec2[10]);
	  }
	break;
      }
    case ZAXIS_PRESSURE:
      {
	double dum;
	char units[128];

	level = zaxisInqLevel(zaxisID, levelID);
	if ( level < 0 )
	  Warning("Pressure level of %f Pa is below zero!", level);

	zaxisInqUnits(zaxisID, units);
	if ( memcmp(units, "Pa", 2) != 0 ) level *= 100;

	ilevel = (int) level;
	if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_99;
	    ISEC1_Level1    = ilevel;
	    ISEC1_Level2    = 0;
	  }
	else
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_ISOBARIC;
	    ISEC1_Level1    = ilevel/100;
	    ISEC1_Level2    = 0;
	  }
	break;
      }
    case ZAXIS_HEIGHT:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_HEIGHT;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_ALTITUDE:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_ALTITUDE;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_SIGMA:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
	    ISEC1_LevelType = GRIB1_LTYPE_SIGMA_LAYER;
	    ISEC1_Level1    = (int) zaxisInqLbound(zaxisID, levelID);
	    ISEC1_Level2    = (int) zaxisInqUbound(zaxisID, levelID);
	  }
	else
	  {
            level = zaxisInqLevel(zaxisID, levelID);

            ilevel = (int) level;
            ISEC1_LevelType = GRIB1_LTYPE_SIGMA;
            ISEC1_Level1    = ilevel;
            ISEC1_Level2    = 0;
          }

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	char units[128];
	double factor;

	zaxisInqUnits(zaxisID, units);

        if      ( memcmp(units, "mm", 2) == 0 ) factor =   0.1;
        else if ( memcmp(units, "cm", 2) == 0 ) factor =   1;
        else if ( memcmp(units, "dm", 2) == 0 ) factor =  10;
        else                                    factor = 100; // meter

	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	  {
            double level1, level2;
            level1 = zaxisInqLbound(zaxisID, levelID);
            level2 = zaxisInqUbound(zaxisID, levelID);
	    ISEC1_LevelType = GRIB1_LTYPE_LANDDEPTH_LAYER;
	    ISEC1_Level1    = (int) (level1*factor);
	    ISEC1_Level2    = (int) (level2*factor);
	  }
	else
	  {
	    level = zaxisInqLevel(zaxisID, levelID);

	    ilevel = (int) (level*factor);
	    ISEC1_LevelType = GRIB1_LTYPE_LANDDEPTH;
	    ISEC1_Level1    = ilevel;
	    ISEC1_Level2    = 0;
	  }

	break;
      }
    case ZAXIS_DEPTH_BELOW_SEA:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_SEADEPTH;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_ISENTROPIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = GRIB1_LTYPE_ISENTROPIC;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    case ZAXIS_GENERIC:
      {
	level = zaxisInqLevel(zaxisID, levelID);

	ilevel = (int) level;
	ISEC1_LevelType = ltype;
	ISEC1_Level1    = ilevel;
	ISEC1_Level2    = 0;

	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
	break;
      }
    }
}

static
void cgribexDefaultSec0(int *isec0)
{
  ISEC0_GRIB_Len     = 0;
  ISEC0_GRIB_Version = 0;
}

static
void cgribexDefaultSec1(int *isec1)
{
  ISEC1_CenterID    = 0;
  ISEC1_SubCenterID = 0;
  ISEC1_LocalFLag   = 0;
}

static
void cgribexDefaultSec4(int *isec4)
{
  long i;

  for ( i = 2; i <= 10; ++i ) isec4[i] = 0;
}

static
void cgribexDefEnsembleVar(int *isec1, int vlistID, int varID)
{
  int ensID, ensCount, forecast_type;

  /* For Ensemble info  */

  //Put1Byte(isec1[36]);        /* MPIM local GRIB use definition identifier  */
                                /*    (extension identifier)                  */
  //Put1Byte(isec1[37]);        /* type of ensemble forecast                  */
  //Put2Byte(isec1[38]);        /* individual ensemble member                 */
  //Put2Byte(isec1[39]);        /* number of forecasts in ensemble            */

  if ( vlistInqVarEnsemble(vlistID, varID, &ensID, &ensCount, &forecast_type) )
    {
      if ( ISEC1_CenterID == 252 )
        {
          ISEC1_LocalFLag = 1;
          isec1[36] = 1;

          isec1[37] =  forecast_type;
          isec1[38] =  ensID;
          isec1[39] =  ensCount;
        }
    }
}
#endif


#if  defined  (HAVE_LIBCGRIBEX)
size_t cgribexEncode(int memtype, int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg,
		     long datasize, const double *data, int nmiss, unsigned char *gribbuffer, size_t gribbuffersize)
{
  size_t nbytes = 0;
  int gribsize;
  int iret = 0, iword = 0;
  int isec0[2], isec1[4096], isec2[4096], isec3[2], isec4[512];
  float fsec2f[512], fsec3f[2];
  double fsec2[512], fsec3[2];
  int datatype;
  int param;

  memset(isec1, 0, 256*sizeof(int));
  fsec2[0] = 0; fsec2[1] = 0;
  fsec2f[0] = 0; fsec2f[1] = 0;

  gribsize = (int)(gribbuffersize / sizeof(int));
  param    = vlistInqVarParam(vlistID, varID);

  cgribexDefaultSec0(isec0);
  cgribexDefaultSec1(isec1);
  cgribexDefaultSec4(isec4);

  cgribexDefInstitut(isec1, vlistID, varID);
  cgribexDefModel(isec1, vlistID, varID);

  datatype = vlistInqVarDatatype(vlistID, varID);

  cgribexDefParam(isec1, param);
  cgribexDefTime(isec1, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID));
  cgribexDefGrid(isec1, isec2, fsec2, isec4, gridID);
  cgribexDefLevel(isec1, isec2, fsec2, zaxisID, levelID);

  cgribexDefEnsembleVar(isec1, vlistID, varID);

  ISEC4_NumValues = gridInqSize(gridID);
  ISEC4_NumBits   = grbBitsPerValue(datatype);

  if ( nmiss > 0 )
    {
      FSEC3_MissVal = vlistInqVarMissval(vlistID, varID);
      ISEC1_Sec2Or3Flag |= 64;
    }

  if ( isec4[2] == 128 && isec4[3] == 64 )
    {
      isec4[16] = (int) (1000*calculate_pfactor(data, ISEC2_PentaJ, isec4[17]));
      if ( isec4[16] < -10000 ) isec4[16] = -10000;
      if ( isec4[16] >  10000 ) isec4[16] =  10000;
    }
  //printf("isec4[16] %d\n", isec4[16]);

  if ( memtype == MEMTYPE_FLOAT )
    {
      size_t numVCP = ISEC2_NumVCP > 0 ? (size_t)ISEC2_NumVCP : (size_t)0;
      for ( size_t i = 0; i < numVCP; ++i ) fsec2f[10+i] = (float)fsec2[10+i];
      fsec3f[ 1] = (float)fsec3[ 1];
    }

  if ( memtype == MEMTYPE_FLOAT )
    gribExSP(isec0, isec1, isec2, fsec2f, isec3, fsec3f, isec4, (float*) data,
             (int)datasize, (int *) gribbuffer, gribsize, &iword, "C", &iret);
  else
    gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, (double*) data,
             (int)datasize, (int *) gribbuffer, gribsize, &iword, "C", &iret);

  if ( iret ) Error("Problem during GRIB encode (errno = %d)!", iret);

  nbytes = (size_t)iword * sizeof (int);
  return (nbytes);
}
#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
