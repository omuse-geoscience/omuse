#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

void uuid2str(const unsigned char uuid[CDI_UUID_SIZE], char *uuidstr);

static inline
int cdiUUIDIsNull(const unsigned char uuid[CDI_UUID_SIZE])
{
  static unsigned char uuid_nil[CDI_UUID_SIZE];
  return !memcmp(uuid, uuid_nil, CDI_UUID_SIZE);
}


void datetime2str(int date, int time, char *datetimestr, int maxlen)
{
  int year, month, day;
  int hour, minute, second;
  int len;

  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  len = sprintf(datetimestr, DATE_FORMAT "T" TIME_FORMAT, year, month, day, hour, minute, second);

  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void date2str(int date, char *datestr, int maxlen)
{
  int year, month, day;
  int len;

  cdiDecodeDate(date, &year, &month, &day);

  len = sprintf(datestr, DATE_FORMAT, year, month, day);

  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void time2str(int time, char *timestr, int maxlen)
{
  int hour, minute, second;
  int len;

  cdiDecodeTime(time, &hour, &minute, &second);

  len = sprintf(timestr, TIME_FORMAT, hour, minute, second);

  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void printFiletype(int streamID, int vlistID)
{
  int filetype;

  filetype = streamInqFiletype(streamID);

  switch ( filetype )
    {
    case FILETYPE_GRB:
      printf("GRIB");
      break;
    case FILETYPE_GRB2:
      printf("GRIB2");
      break;
    case FILETYPE_NC:
      printf("netCDF");
      break;
    case FILETYPE_NC2:
      printf("netCDF2");
      break;
    case FILETYPE_NC4:
      printf("netCDF4");
      break;
    case FILETYPE_NC4C:
      printf("netCDF4 classic");
      break;
    case FILETYPE_SRV:
      printf("SERVICE");
      break;
    case FILETYPE_EXT:
      printf("EXTRA");
      break;
    case FILETYPE_IEG:
      printf("IEG");
      break;
    default:
      printf("  File format: unsupported filetype %d" , filetype);
      break;
    }

  if ( filetype == FILETYPE_SRV || filetype == FILETYPE_EXT || filetype == FILETYPE_IEG )
    {
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
    }

  if ( filetype == FILETYPE_GRB || filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C )
    {
      int nvars, varID;
      int comptype;

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == COMPRESS_SZIP )
		printf(" SZIP");
	      else if ( comptype == COMPRESS_ZIP )
		printf(" ZIP");

	      break;
	    }
	}
    }

  if ( filetype == FILETYPE_GRB2 )
    {
      int nvars, varID;
      int comptype;

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == COMPRESS_JPEG )
		printf(" JPEG");

	      break;
	    }
	}
    }

  printf("\n");
}

static
void printGridInfo(int vlistID)
{
  int ngrids, index;
  int gridID, gridtype, trunc, gridsize, xsize, ysize, xysize;
  char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME], xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];

  ngrids = vlistNgrids(vlistID);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID   = vlistGrid(vlistID, index);
      gridtype = gridInqType(gridID);
      trunc    = gridInqTrunc(gridID);
      gridsize = gridInqSize(gridID);
      xsize    = gridInqXsize(gridID);
      ysize    = gridInqYsize(gridID);
      xysize   = xsize*ysize;
      gridInqXname(gridID, xname);
      gridInqYname(gridID, yname);
      gridInqXunits(gridID, xunits);
      gridInqYunits(gridID, yunits);

      fprintf(stdout, "  %4d : %-24s", index+1, gridNamePtr(gridtype));

      if ( gridtype == GRID_LONLAT   ||
	   gridtype == GRID_LCC2 ||
	   gridtype == GRID_LAEA ||
	   gridtype == GRID_SINUSOIDAL ||
	   gridtype == GRID_GENERIC ||
	   gridtype == GRID_GAUSSIAN ||
	   gridtype == GRID_GAUSSIAN_REDUCED )
	{
          int lxcoord = 1, lycoord = 1;
	  double xfirst = 0.0, xlast = 0.0;
	  double yfirst = 0.0, ylast = 0.0;
	  double xinc = 0.0, yinc = 0.0;

	  yfirst = gridInqYval(gridID, 0);
	  ylast  = gridInqYval(gridID, ysize-1);
	  yinc   = gridInqYinc(gridID);

          fprintf(stdout, " : points=%d", gridsize);
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    fprintf(stdout, "  nlat=%d", ysize);
	  else if ( xysize )
	    fprintf(stdout, " (%dx%d)", xsize, ysize);

	  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
	    fprintf(stdout, "  np=%d", gridInqNP(gridID));

	  fprintf(stdout, "\n");

          if ( gridInqXvals(gridID, NULL) == 0 ) lxcoord = 0;
          if ( gridInqYvals(gridID, NULL) == 0 ) lycoord = 0;

	  if ( xsize > 0 && lxcoord )
	    {
              xfirst = gridInqXval(gridID, 0);
              xlast  = gridInqXval(gridID, xsize-1);
              xinc   = gridInqXinc(gridID);
              fprintf(stdout, "%33s : %g", xname, xfirst);
              if ( xsize > 1 )
                {
                  fprintf(stdout, " to %g", xlast);
                  if ( IS_NOT_EQUAL(xinc, 0) )
                    fprintf(stdout, " by %g", xinc);
                }
              fprintf(stdout, " %s", xunits);
              if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
              fprintf(stdout, "\n");
	    }

	  if ( ysize > 0 && lycoord )
	    {
	      fprintf(stdout, "%33s : %g", yname, yfirst);
	      if ( ysize > 1 )
                {
                  fprintf(stdout, " to %g", ylast);
                  if ( IS_NOT_EQUAL(yinc, 0) && gridtype != GRID_GAUSSIAN && gridtype != GRID_GAUSSIAN_REDUCED )
                    fprintf(stdout, " by %g", yinc);
                }
              fprintf(stdout, " %s", yunits);
	      fprintf(stdout, "\n");
	    }

	  if ( gridIsRotated(gridID) )
	    {
	      double lonpole, latpole, angle;
	      lonpole = gridInqXpole(gridID);
	      latpole = gridInqYpole(gridID);
	      angle   = gridInqAngle(gridID);
	      fprintf(stdout, "%33s : lon=%g  lat=%g", "northpole", lonpole, latpole);
	      if ( IS_NOT_EQUAL(angle, 0) ) fprintf(stdout, "  angle=%g", angle);
	      fprintf(stdout, "\n");
	    }

	  if ( gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
	    {
	      fprintf(stdout, "%33s :", "available");
	      if ( gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ) fprintf(stdout, " cellbounds");
	      if ( gridHasArea(gridID) )          fprintf(stdout, " area");
	      if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
	      fprintf(stdout, "\n");
	    }

	  if ( gridtype == GRID_LAEA )
	    {
	      double a, lon_0, lat_0;
	      gridInqLaea(gridID, &a, &lon_0, &lat_0);
	      fprintf(stdout, "%33s : a=%g  lon_0=%g  lat_0=%g\n", "projpar", a, lon_0, lat_0);
	    }

	  if ( gridtype == GRID_LCC2 )
	    {
	      double a, lon_0, lat_0, lat_1, lat_2;
	      gridInqLcc2(gridID, &a, &lon_0, &lat_0, &lat_1, &lat_2);
	      fprintf(stdout, "%33s : a=%7.0f  lon_0=%g  lat_0=%g  lat_1=%g  lat_2=%g\n",
                      "projpar", a, lon_0, lat_0, lat_1, lat_2);
	    }
	}
      else if ( gridtype == GRID_SPECTRAL )
	{
	  fprintf(stdout, " : points=%d  nsp=%d  truncation=%d", gridsize, gridsize/2, trunc);
          if ( gridInqComplexPacking(gridID) ) fprintf(stdout, "  complexPacking");
          fprintf(stdout, "\n");
	}
      else if ( gridtype == GRID_FOURIER )
	{
	  fprintf(stdout, " : points=%d  nfc=%d  truncation=%d\n", gridsize, gridsize/2, trunc);
	}
      else if ( gridtype == GRID_GME )
	{
	  int ni, nd;
	  ni = gridInqGMEni(gridID);
	  nd = gridInqGMEnd(gridID);
	  fprintf(stdout, " : points=%d  nd=%d  ni=%d\n", gridsize, nd, ni);
	}
      else if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	{
	  if ( gridtype == GRID_CURVILINEAR )
	    fprintf(stdout, " : points=%d (%dx%d)", gridsize, xsize, ysize);
	  else
	    fprintf(stdout, " : points=%d", gridsize);

          if ( gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) > 0 )
	    fprintf(stdout, "  nvertex=%d", gridInqNvertex(gridID));

          fprintf(stdout, "\n");

          if ( gridtype == GRID_UNSTRUCTURED )
            {
              int number   = gridInqNumber(gridID);
              int position = gridInqPosition(gridID);

              if ( number > 0 )
                {
                  fprintf(stdout, "%33s : number=%d  position=%d\n", "grid", number, position);
                }

              if ( gridInqReference(gridID, NULL) )
                {
                  char reference_link[8192];
                  gridInqReference(gridID, reference_link);
                  fprintf(stdout, "%33s : %s\n", "uri", reference_link);
                }
            }

	  if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	    {
	      int i;
	      double *xvals, *yvals;
	      double xfirst, xlast, yfirst, ylast;
	      xvals = (double*) malloc((size_t)gridsize*sizeof(double));
	      yvals = (double*) malloc((size_t)gridsize*sizeof(double));

	      gridInqXvals(gridID, xvals);
	      gridInqYvals(gridID, yvals);

	      xfirst = xvals[0];
	      xlast  = xvals[0];
	      yfirst = yvals[0];
	      ylast  = yvals[0];
	      for ( i = 1; i < gridsize; i++ )
		{
		  if ( xvals[i] < xfirst ) xfirst = xvals[i];
		  if ( xvals[i] > xlast  ) xlast  = xvals[i];
		  if ( yvals[i] < yfirst ) yfirst = yvals[i];
		  if ( yvals[i] > ylast  ) ylast  = yvals[i];
		}

	      fprintf(stdout, "%33s : %g to %g %s", xname, xfirst, xlast, xunits);
	      if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
	      fprintf(stdout, "\n");
	      fprintf(stdout, "%33s : %g to %g %s\n", yname, yfirst, ylast, yunits);

	      free(xvals);
	      free(yvals);
	    }
	}
      else if ( gridtype == GRID_LCC )
	{
	  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	  int projflag, scanflag;

	  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		     &projflag, &scanflag);

	  fprintf(stdout, " : points=%d (%dx%d)  ", gridsize, xsize, ysize);
	  if ( (projflag&128) == 0 )
	    fprintf(stdout, "North Pole\n");
	  else
	    fprintf(stdout, "South Pole\n");

	  fprintf(stdout, "%33s : originLon=%g  originLat=%g  lonParY=%g\n", " ", originLon, originLat, lonParY);
	  fprintf(stdout, "%33s : lat1=%g  lat2=%g  xinc=%g m  yinc=%g m\n", " ", lat1, lat2, xincm, yincm);
	}
      else /* if ( gridtype == GRID_GENERIC ) */
	{
	  if ( ysize == 0 )
	    fprintf(stdout, " : points=%d\n", gridsize);
	  else
            fprintf(stdout, " : points=%d (%dx%d)\n", gridsize, xsize, ysize);
	}

      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || gridtype == GRID_LCC )
	{
	  if ( gridHasArea(gridID) ||
	       gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
	    {
	      fprintf(stdout, "%33s :", "available");
	      if ( gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ) fprintf(stdout, " cellbounds");
	      if ( gridHasArea(gridID) )          fprintf(stdout, " area");
	      if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
	      fprintf(stdout, "\n");
	    }
	}

      gridInqUUID(gridID, uuidOfHGrid);
      if ( !cdiUUIDIsNull(uuidOfHGrid) )
        {
          char uuidOfHGridStr[37];
          uuid2str(uuidOfHGrid, uuidOfHGridStr);
          if ( uuidOfHGridStr[0] != 0  && strlen(uuidOfHGridStr) == 36 )
            {
	      fprintf(stdout, "%33s : %s\n", "uuid", uuidOfHGridStr);
            }
        }
    }
}

static
void printZaxisInfo(int vlistID)
{
  int zaxisID, zaxistype, levelsize, levelID;
  int ltype;
  double *levels = NULL;
  char zaxisname[CDI_MAX_NAME], zname[CDI_MAX_NAME], zunits[CDI_MAX_NAME];

  int nzaxis = vlistNzaxis(vlistID);
  for ( int index = 0; index < nzaxis; index++)
    {
      double zfirst = 0, zlast = 0, zinc = 0;
      zaxisID   = vlistZaxis(vlistID, index);
      zaxistype = zaxisInqType(zaxisID);
      ltype     = zaxisInqLtype(zaxisID);
      levelsize = zaxisInqSize(zaxisID);
      zaxisName(zaxistype, zaxisname);
      zaxisInqName(zaxisID, zname);
      zaxisInqUnits(zaxisID, zunits);
      zunits[12] = 0;

      if ( zaxistype == ZAXIS_GENERIC && ltype != 0 )
        fprintf(stdout, "  %4d : %-12s (ltype=%3d) :", vlistZaxisIndex(vlistID, zaxisID)+1, zaxisname, ltype);
      else
        fprintf(stdout, "  %4d : %-24s :", vlistZaxisIndex(vlistID, zaxisID)+1, zaxisname);

      fprintf(stdout, " levels=%d", levelsize);
      fprintf(stdout, "\n");

      levels = (double*) malloc((size_t)levelsize*sizeof(double));
      zaxisInqLevels(zaxisID, levels);

      if ( !(zaxistype == ZAXIS_SURFACE && levelsize == 1 && !(fabs(levels[0]) > 0)) )
        {
          zfirst = levels[0];
          zlast  = levels[levelsize-1];
          if ( levelsize > 2 )
            {
              zinc = (levels[levelsize-1] - levels[0]) / (levelsize-1);
              for ( levelID = 2; levelID < levelsize; ++levelID )
                if ( fabs(fabs(levels[levelID] - levels[levelID-1]) - zinc) > 0.001*zinc ) break;

              if ( levelID < levelsize ) zinc = 0;
            }

          fprintf(stdout, "%33s : %g", zname, zfirst);
          if ( levelsize > 1 )
            {
              fprintf(stdout, " to %g", zlast);
              if ( IS_NOT_EQUAL(zinc, 0) )
                fprintf(stdout, " by %g", zinc);
            }
          fprintf(stdout, " %s", zunits);
          fprintf(stdout, "\n");
        }

      free(levels);

      if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
        {
          double level1, level2;
          fprintf(stdout, "%33s : ", "bounds");

          level1 = zaxisInqLbound(zaxisID, 0);
          level2 = zaxisInqUbound(zaxisID, 0);
          fprintf(stdout, "%.9g-%.9g", level1, level2);
          if ( levelsize > 1 )
            {
              level1 = zaxisInqLbound(zaxisID, levelsize-1);
              level2 = zaxisInqUbound(zaxisID, levelsize-1);
              fprintf(stdout, " to %.9g-%.9g", level1, level2);
              if ( IS_NOT_EQUAL(zinc, 0) )
                fprintf(stdout, " by %g", zinc);
            }
          fprintf(stdout, " %s", zunits);
          fprintf(stdout, "\n");
        }

      if ( zaxistype == ZAXIS_REFERENCE )
        {
          int number   = zaxisInqNumber(zaxisID);

          if ( number > 0 )
            {
              fprintf(stdout, "%33s : ", "zaxis");
              fprintf(stdout, "number = %d\n", number);
            }

          unsigned char uuidOfVGrid[CDI_UUID_SIZE];
          zaxisInqUUID(zaxisID, uuidOfVGrid);
          if ( !cdiUUIDIsNull(uuidOfVGrid) )
            {
              char uuidOfVGridStr[37];
              uuid2str(uuidOfVGrid, uuidOfVGridStr);
              if ( uuidOfVGridStr[0] != 0  && strlen(uuidOfVGridStr) == 36 )
                {
                  fprintf(stdout, "%33s : ", "uuid");
                  fprintf(stdout, "%s\n", uuidOfVGridStr);
                }
            }
        }
    }
}

static
void printSubtypeInfo(int vlistID)
{
  int nsubtypes = vlistNsubtypes(vlistID);
  for ( int index = 0; index < nsubtypes; index++)
    {
      int subtypeID = vlistSubtype(vlistID, index);
      int subtypesize = subtypeInqSize(subtypeID);
      // subtypePrint(subtypeID);
      fprintf(stdout, "  %4d : %-24s :", vlistSubtypeIndex(vlistID, subtypeID)+1, "tiles");
      fprintf(stdout, " ntiles=%d", subtypesize);
      fprintf(stdout, "\n");

    }
}

static
int printDateTime(int ntimeout, int vdate, int vtime)
{
  char vdatestr[32], vtimestr[32];

  if ( ntimeout == 4 )
    {
      ntimeout = 0;
      fprintf(stdout, "\n");
    }

  date2str(vdate, vdatestr, sizeof(vdatestr));
  time2str(vtime, vtimestr, sizeof(vtimestr));

  fprintf(stdout, " %s %s", vdatestr, vtimestr);

  return (++ntimeout);
}

#define NUM_TIMESTEP 60
#define MAX_DOTS     80

static
int printDot(int ndotout, int *nfact, int *ncout)
{

  //printf("ncout %d %d %d\n",*ncout, (*ncout)%(*nfact), *nfact);
  if ( (*ncout)%(*nfact) == 0 )
    {
      if ( ndotout == MAX_DOTS )
	{
	  *ncout = 0;
	  ndotout = 0;
	  fprintf(stdout, "\n   ");
	  (*nfact) *= 10;
	}

      fprintf(stdout, ".");
      fflush(stdout);
      ndotout++;
    }

  (*ncout)++;

  return (ndotout);
}

static
void printTimesteps(int streamID, int taxisID, int verbose)
{
  int nrecs;
  int vdate, vtime;
  struct datetime {
    int vdate;
    int vtime;
    struct datetime *next;
  };
  struct datetime vdatetime[NUM_TIMESTEP];
  struct datetime *next_vdatetime = vdatetime;

  for ( int i = 0; i < NUM_TIMESTEP-1; ++i ) vdatetime[i].next = &vdatetime[i+1];
  vdatetime[NUM_TIMESTEP-1].next = &vdatetime[0];

  int ntimeout = 0;
  int ndotout = 0;
  int nvdatetime = 0;
  int ncout = 0;
  int nfact = 1;
  int tsID = 0;

  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      if ( verbose || tsID < NUM_TIMESTEP )
	{
	  ntimeout = printDateTime(ntimeout, vdate, vtime);
	}
      else
	{
	  if ( tsID == 2*NUM_TIMESTEP ) fprintf(stdout, "\n   ");
	  if ( tsID >= 2*NUM_TIMESTEP ) ndotout = printDot(ndotout, &nfact, &ncout);

	  if ( nvdatetime < NUM_TIMESTEP )
	    {
	      vdatetime[nvdatetime].vdate = vdate;
	      vdatetime[nvdatetime].vtime = vtime;
	      nvdatetime++;
	    }
	  else
	    {
	      next_vdatetime->vdate = vdate;
	      next_vdatetime->vtime = vtime;
	      next_vdatetime = next_vdatetime->next;
	    }
	}

      tsID++;
    }

  if ( nvdatetime )
    {
      fprintf(stdout, "\n");

      ntimeout = 0;
      int toff = tsID%4;
      if ( toff > 0 ) toff = 4 - toff;
      for ( int i = 0; i < toff; ++i ) next_vdatetime = next_vdatetime->next;

      for ( int i = toff; i < nvdatetime; ++i )
	{
	  vdate = next_vdatetime->vdate;
	  vtime = next_vdatetime->vtime;
	  ntimeout = printDateTime(ntimeout, vdate, vtime);
	  next_vdatetime = next_vdatetime->next;
	}
    }
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
