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

      Gradsdes   gradsdes       GrADS data descriptor file
*/


#if defined(HAVE_CONFIG_H)
#  include "config.h" /* VERSION */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


/*
  Values output into the grib map file:

  Header:

  hipnt info:  0 - version number (1)
               1 - number of times in file
               2 - number of records per time
               3 - Grid type
                 255 - user defined grid.  descriptor
                       describes grid exactly; one record
                       per grid.
                  29 - Predefined grid set 29 and 30.
                       Two records per grid.

  hfpnt info:  None

  Info:

  intpnt info (for each mapped grib record) :
                 0 - position of start of data in file
                 1 - position of start of bit map in file
                 2 - number of bits per data element

  fltpnt info :
                 0 - decimal scale factor for this record
                 1 - binary scale factor
                 2 - reference value

*/
struct gaindx {
  int    type;      /* Indexing file type             */
  int    hinum;     /* Number of ints in header       */
  int    hfnum;     /* Number of floats in header     */
  int    intnum;    /* Number of index ints (long)    */
  int    fltnum;    /* Number of index floats         */
  int   *hipnt;     /* Pointer to header int values   */
  float *hfpnt;     /* Pointer to header float values */
  int   *intpnt;    /* Pointer to int index values    */
  float *fltpnt;    /* Pointer to float index values  */
};
struct gaindxb {
  int    bignum;    /* Number of off_t values */
  off_t *bigpnt;    /* Pointer to off_t values */
};


/* Byte swap requested number of 4 byte elements */

static
void gabswp (void *r, int cnt) {
int i;
char *ch1,*ch2,*ch3,*ch4,cc1,cc2;

  ch1 = (char *)r;
  ch2 = ch1+1;
  ch3 = ch2+1;
  ch4 = ch3+1;
  for (i=0; i<cnt; i++) {
    cc1 = *ch1;
    cc2 = *ch2;
    *ch1 = *ch4;
    *ch2 = *ch3;
    *ch3 = cc2;
    *ch4 = cc1;
    ch1+=4; ch2+=4; ch3+=4; ch4+=4;
  }
}

/*
 * convert an IBM float to single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 */
static
float ibm2flt(unsigned char *ibm) {

  int positive, power;
  unsigned int abspower;
  long int mant;
  double value, exp;

  positive = (ibm[0] & 0x80) == 0;
  mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
  power = (int) (ibm[0] & 0x7f) - 64;
  abspower = power > 0 ? power : -power;


  /* calc exp */
  exp = 16.0;
  value = 1.0;
  while (abspower) {
    if (abspower & 1) {
      value *= exp;
    }
    exp = exp * exp;
    abspower >>= 1;
  }

  if (power < 0) value = 1.0 / value;
  value = value * mant / 16777216.0;
  if (positive == 0) value = -value;
  return (float)value;
}

/*
 * convert a float to an IBM single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 *
 * doesn't handle subnormal numbers
 */
static
int flt2ibm(float x, unsigned char *ibm) {

  int sign, exp, i;
  double mant;

  if ( !(fabs((double)x) > 0) ) {
    ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
    return 0;
  }

  /* sign bit */
  if (x < 0.0) {
    sign = 128;
    x = -x;
  }
  else sign = 0;

  mant = frexp((double) x, &exp);

  /* round up by adding 2**-24 */
  /* mant = mant + 1.0/16777216.0; */

  if (mant >= 1.0) {
    mant = 0.5;
    exp++;
  }
  while (exp & 3) {
    mant *= 0.5;
    exp++;
  }

  exp = exp/4 + 64;

  if (exp < 0) {
    fprintf(stderr,"underflow in flt2ibm\n");
    ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
    return 0;
  }
  if (exp > 127) {
    fprintf(stderr,"overflow in flt2ibm\n");
    ibm[0] = sign | 127;
    ibm[1] = ibm[2] = ibm[3] = 255;
    return -1;
  }

  /* normal number */

  ibm[0] = sign | exp;

  mant = mant * 256.0;
  i = (int) floor(mant);
  mant = mant - i;
  ibm[1] = i;

  mant = mant * 256.0;
  i = (int) floor(mant);
  mant = mant - i;
  ibm[2] = i;

  ibm[3] = (int) floor(mant*256.0);

  return 0;
}

#define  GET_UINT4(a,b,c,d) ((int) ((a << 24) + (b << 16) + (c << 8) + (d)))
#define  Put1Byte(buf, cnt, ival)  (buf[cnt++] = (ival))
#define  Put2Byte(buf, cnt, ival) ((buf[cnt++] = (ival) >>  8), \
                                   (buf[cnt++] = (ival)))
#define  Put4Byte(buf, cnt, ival) ((buf[cnt++] = (ival) >> 24), \
                                   (buf[cnt++] = (ival) >> 16), \
                                   (buf[cnt++] = (ival) >>  8), \
                                   (buf[cnt++] = (ival)))

#define  PutInt(buf, cnt, ival)   (ival < 0 ? Put4Byte(buf, cnt, 0x7fffffff - ival + 1) : Put4Byte(buf, cnt, ival))

static
void dumpmap()
{
  unsigned char urec[4];
  unsigned char vermap;
  unsigned char mrec[512];
  int swpflg = 0;
  int i;
  int nrecords = 0;
  struct gaindx indx;
  struct gaindxb indxb;
  size_t nbytes;
  FILE *mapfp;

  indx.hipnt = NULL;
  indx.hfpnt = NULL;
  indx.intpnt = NULL;
  indx.fltpnt = NULL;
  indxb.bigpnt = NULL;
  indxb.bignum = 0;

  mapfp = fopen(cdoStreamName(0)->args, "r");
  if ( mapfp == NULL ) cdoAbort("Open failed on %s", cdoStreamName(0)->args);

  /* check the version number */

  fseek(mapfp, 1, 0);
  nbytes = fread(&vermap, sizeof(unsigned char), 1, mapfp);

  if ( vermap == 0 ) vermap = 1;

  printf("gribmap version = %d\n", vermap);

  if ( vermap == 2 )
    {
      fseek(mapfp, 2, 0);

      nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
      indx.hinum = GET_UINT4(mrec[0],mrec[1],mrec[2],mrec[3]);

      nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
      indx.hfnum = GET_UINT4(mrec[0],mrec[1],mrec[2],mrec[3]);

      nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
      indx.intnum = GET_UINT4(mrec[0],mrec[1],mrec[2],mrec[3]);

      nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
      indx.fltnum = GET_UINT4(mrec[0],mrec[1],mrec[2],mrec[3]);

      nbytes = fread(mrec, sizeof(unsigned char), 7, mapfp);

      if ( indx.hinum > 0 )
        {
          indx.hipnt = (int*) malloc(sizeof(int)*indx.hinum);
          for ( i = 0; i < indx.hinum; i++ )
            {
              nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
              indx.hipnt[i] = GET_UINT4(mrec[0],mrec[1],mrec[2],mrec[3]);
            }
        }
      if ( indx.hfnum > 0 )
        {
          indx.hfpnt = (float*) malloc(sizeof(float)*indx.hfnum);
          nbytes = fread (indx.hfpnt,sizeof(float),indx.hfnum,mapfp);
        }
      if ( indx.intnum > 0 )
        {
          indx.intpnt = (int*) malloc(sizeof(int)*indx.intnum);
          for ( i = 0; i < indx.intnum; i++ )
            {
              nbytes = fread(mrec, sizeof(unsigned char), 4, mapfp);
              indx.intpnt[i] = GET_UINT4(mrec[0],mrec[1],mrec[2],mrec[3]);
              if ( indx.intpnt[i] < 0 ) indx.intpnt[i] = 0x7fffffff - indx.intpnt[i] + 1;
            }
        }
      if ( indx.fltnum > 0 )
        {
          indx.fltpnt = (float*) malloc(sizeof(float)*indx.fltnum);
          for ( i = 0; i < indx.fltnum; i++ )
            {
              nbytes = fread(urec, sizeof(unsigned char), 4, mapfp);
              indx.fltpnt[i] = ibm2flt(urec);
            }
        }
    }
  else
    {
      fseek(mapfp, 0, 0);
      nbytes = fread (&indx, sizeof(struct gaindx), 1, mapfp);

      if ( indx.type>>24 > 0 ) swpflg = 1;
      if ( swpflg ) printf("swap endian!\n");
      if ( swpflg ) gabswp((float *)&indx.type, 5);

      if ( indx.hinum > 0 )
        {
          indx.hipnt = (int*) malloc(sizeof(int)*indx.hinum);
          nbytes = fread (indx.hipnt, sizeof(int), indx.hinum, mapfp);
          if ( swpflg ) gabswp((float *)(indx.hipnt),indx.hinum);
        }
      if ( indx.hfnum > 0 )
        {
          indx.hfpnt = (float*) malloc(sizeof(float)*indx.hfnum);
          nbytes = fread (indx.hfpnt,sizeof(float),indx.hfnum,mapfp);
          if ( swpflg ) gabswp(indx.hfpnt,indx.hfnum);
        }

      if ( indx.intnum > 0 )
        {
          indx.intpnt = (int*) malloc(sizeof(int)*indx.intnum);
          nbytes = fread (indx.intpnt,sizeof(int),indx.intnum,mapfp);
          if ( swpflg ) gabswp((float *)(indx.intpnt),indx.intnum);
        }
      if ( indx.fltnum > 0 )
        {
          indx.fltpnt = (float*) malloc(sizeof(float)*indx.fltnum);
          nbytes = fread (indx.fltpnt,sizeof(float),indx.fltnum,mapfp);
          if ( swpflg ) gabswp(indx.fltpnt,indx.fltnum);
        }

      if ( indx.hipnt[0] == 4 )
        {
          indxb.bignum = indx.hipnt[4];
          if ( indxb.bignum > 0 )
            {
              indxb.bigpnt = (off_t*) malloc(sizeof(off_t)*indxb.bignum);
              nbytes = fread (indxb.bigpnt,sizeof(off_t),indxb.bignum,mapfp);
              if ( swpflg ) gabswp(indxb.bigpnt,indxb.bignum);
            }
        }
    }

  UNUSED(nbytes);

  fclose(mapfp);

  printf("hinum: %d\n", indx.hinum);
  for ( i = 0; i < indx.hinum; i++ )
    printf("%3d %5d\n", i+1, indx.hipnt[i]);

  printf("\n");
  printf("hfnum: %d\n", indx.hfnum);
  for ( i = 0; i < indx.hfnum; i++ )
    printf("%3d %g\n", i+1, indx.hfpnt[i]);

  printf("\n");

  nrecords = indx.hipnt[1]*indx.hipnt[2];

  if ( indx.intnum == indx.fltnum )
    {
      printf("num: %d\n", indx.intnum);
      for ( i = 0; i < indx.intnum/3; i++ )
        printf("%3d %8d %6d %4d %8g %10g %8g\n", i+1,
               indx.intpnt[i*3], indx.intpnt[i*3+1], indx.intpnt[i*3+2],
               indx.fltpnt[i*3], indx.fltpnt[i*3+1], indx.fltpnt[i*3+2]);
    }
  else if ( indx.intnum == nrecords && indx.fltnum == nrecords*3 && indxb.bignum == nrecords*2 )
    {
      printf("nrecords: %d\n", nrecords);
      for ( i = 0; i < nrecords; i++ )
        printf("%3d %8zd %6zd %4d %8g %10g %8g\n", i+1,
               (size_t)indxb.bigpnt[i*2], (size_t)indxb.bigpnt[i*2+1], indx.intpnt[i],
               indx.fltpnt[i*3], indx.fltpnt[i*3+1], indx.fltpnt[i*3+2]);
    }
  else
    {
      printf("intnum: %d\n", indx.intnum);
      for ( i = 0; i < indx.intnum; i++ )
        printf("%3d %d\n", i+1, indx.intpnt[i]);

      printf("\n");
      printf("fltnum: %d\n", indx.fltnum);
      for ( i = 0; i < indx.fltnum; i++ )
        printf("%3d %g\n", i+1, indx.fltpnt[i]);

      printf("\n");
      printf("bignum: %d\n", indxb.bignum);
      for ( i = 0; i < indxb.bignum; i++ )
        printf("%3d %zd\n", i+1, (size_t)indxb.bigpnt[i]);
    }
}

static
void ctl_xydef(FILE *gdp, int gridID, int *yrev)
{
  int gridtype;
  int i, j;
  int xsize, ysize;
  double xfirst, yfirst, xinc, yinc;
  double *xvals, *yvals;

  *yrev = FALSE;

  xsize  = gridInqXsize(gridID);
  ysize  = gridInqYsize(gridID);

  gridtype = gridInqType(gridID);

  /* XDEF */

  if ( gridtype == GRID_LCC )
    {
      double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
      double xmin = 1.e10, xmax = -1.e10, ymin = 1.e10, ymax = -1.e10;
      double xrange, yrange;
      int projflag, scanflag;
      int nx, ny, ni;
      double inc[] = { 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001 };

      gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
                 &projflag, &scanflag);
      fprintf(gdp, "PDEF %d %d LCCR %g %g 1 1 %g %g %g %g %g\n",
              xsize, ysize, originLat, originLon, lat1, lat2, lonParY, xincm, yincm);

      gridID = gridToCurvilinear(gridID, 0);
      xvals = (double*) malloc(xsize*ysize*sizeof(double));
      yvals = (double*) malloc(xsize*ysize*sizeof(double));
      gridInqXvals(gridID, xvals);
      gridInqYvals(gridID, yvals);
      for ( i = 0; i < xsize*ysize; ++i )
        {
          if ( xvals[i] > 180  ) xvals[i] -= 360;
          if ( xvals[i] < xmin ) xmin = xvals[i];
          if ( xvals[i] > xmax ) xmax = xvals[i];
          if ( yvals[i] < ymin ) ymin = yvals[i];
          if ( yvals[i] > ymax ) ymax = yvals[i];
        }
      free(xvals);
      free(yvals);

      xfirst = ((int)(xmin-0.0));
      yfirst = ((int)(ymin-0.0));
      xrange = ((int)(xmax+1.5)) - xfirst;
      yrange = ((int)(ymax+1.5)) - yfirst;

      ni = sizeof(inc)/sizeof(inc[0]);
      for ( i = 0; i < ni; i++ )
        {
          xinc = yinc = inc[i];
          nx = 1 + (int) (xrange / xinc);
          ny = 1 + (int) (yrange / yinc);

          if ( nx > 1.5*xsize && ny > 1.5*ysize ) break;
        }

      fprintf(gdp, "XDEF %d LINEAR %f %f\n", nx, xfirst, xinc);
      fprintf(gdp, "YDEF %d LINEAR %f %f\n", ny, yfirst, yinc);

      fprintf(gdp, "* XDEF 3600 LINEAR -179.95 0.1\n");
      fprintf(gdp, "* YDEF 1800 LINEAR  -89.95 0.1\n");
    }
  else
    {
      xfirst = gridInqXval(gridID, 0);
      xinc   = gridInqXinc(gridID);
      if ( IS_EQUAL(xinc, 0) && gridInqXvals(gridID, NULL) )
        {
          xvals = (double*) malloc(xsize*sizeof(double));
          gridInqXvals(gridID, xvals);
          fprintf(gdp ,"XDEF %d LEVELS ", xsize);
          j = 0;
          for ( i = 0; i < xsize; i++ )
            {
              fprintf(gdp, "%7.3f ", xvals[i]);
              j++;
              if ( j == 6 )
                {
                  fprintf(gdp, "\n");
                  j = 0;
                  if ( i != xsize-1 ) fprintf(gdp, "               ");
                }
            }
          if ( j ) fprintf(gdp, "\n");

          free(xvals);
        }
      else
        {
          if ( IS_EQUAL(xinc, 0) ) xinc = 360.0/xsize;
          fprintf(gdp, "XDEF %d LINEAR %f %f\n", xsize, xfirst, xinc);
        }
    }

  /* YDEF */

  if ( gridtype != GRID_LCC )
    {
      yfirst = gridInqYval(gridID, 0);
      yinc   = gridInqYinc(gridID);
      if ( gridtype == GRID_GAUSSIAN ) yinc = 0;

      if ( IS_EQUAL(yinc, 0) && gridInqYvals(gridID, NULL) )
        {
          yvals = (double*) malloc(ysize*sizeof(double));
          gridInqYvals(gridID, yvals);
          fprintf(gdp ,"YDEF %d LEVELS ", ysize);
          j = 0;
          if ( yvals[0] > yvals[ysize-1] )
            {
              *yrev = TRUE;
              for ( i = ysize-1; i >= 0; i-- )
                {
                  fprintf(gdp, "%7.3f ", yvals[i]);
                  j++;
                  if ( j == 6 )
                    {
                      fprintf(gdp, "\n");
                      j = 0;
                      if ( i != 0 ) fprintf(gdp, "               ");
                    }
                }
            }
          else
            {
              for ( i = 0; i < ysize; i++ )
                {
                  fprintf(gdp, "%7.3f ", yvals[i]);
                  j++;
                  if ( j == 6 )
                    {
                      fprintf(gdp, "\n");
                      j = 0;
                      if ( i != ysize-1 ) fprintf(gdp, "               ");
                    }
                }
            }

          if ( j ) fprintf(gdp, "\n");

          free(yvals);
        }
      else
        {
          if ( IS_EQUAL(yinc, 0) ) yinc = 180.0/ysize;
          if ( yinc < 0)
            {
              *yrev = TRUE;
              fprintf(gdp, "YDEF %d LINEAR %f %f\n", ysize, yfirst + yinc * (ysize-1 ), -yinc);
            }
          else
            fprintf(gdp, "YDEF %d LINEAR %f %f\n", ysize, yfirst, yinc);
        }
    }
}

static
void ctl_zdef(FILE *gdp, int vlistID, int *zrev)
{
  int i, j, index;
  int zaxisIDmax = -1, nlevmax;
  int nzaxis, zaxisID, nlev;
  int lplev = FALSE;
  double *levels, level0, levinc = 0;

  *zrev = FALSE;
  nzaxis  = vlistNzaxis(vlistID);

  nlevmax = 0;
  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID, index);
      nlev    = zaxisInqSize(zaxisID);
      if ( nlev > nlevmax )
        {
          nlevmax = nlev;
          zaxisIDmax = zaxisID;
        }
    }

  levels = (double*) malloc(nlevmax*sizeof(double));
  zaxisInqLevels(zaxisIDmax, levels);
  if ( zaxisInqType(zaxisIDmax) == ZAXIS_PRESSURE ) lplev = TRUE;
  level0 = levels[0];
  if ( nlevmax > 1 )
    {
      if ( levels[0] < levels[1] && zaxisInqType(zaxisIDmax) != ZAXIS_HYBRID )
        *zrev = TRUE;

      levinc = levels[1] - levels[0];

      if ( IS_EQUAL(levinc, 1) ) *zrev = FALSE;

      for ( i = 1; i < nlevmax; i++ )
        {
          if ( IS_NOT_EQUAL(levinc, (levels[i] - levels[i-1])) )
            {
              levinc = 0;
              break;
            }
        }
    }

  if ( IS_NOT_EQUAL(levinc, 0) )
    fprintf(gdp,"ZDEF %d LINEAR %g %g\n", nlevmax, level0, levinc);
  else
    {
      fprintf(gdp, "ZDEF %d LEVELS ", nlevmax);
      j  = 0;
      /* zrev not needed !!!
      if ( *zrev )
        {
          for ( i = nlevmax-1; i >=0 ; i-- )
            {
              if ( lplev ) fprintf(gdp, "%g ", levels[i]/100);
              else         fprintf(gdp, "%d ", (int) levels[i]);
              j++;
              if ( j == 10 )
                {
                  fprintf(gdp, "\n");
                  j = 0;
                  if ( i != 0 ) fprintf(gdp, "               ");
                }
            }
        }
      else
      */
        {
          for ( i = 0; i < nlevmax ; i++ )
            {
              if ( lplev ) fprintf(gdp, "%g ", levels[i]/100);
              else         fprintf(gdp, "%g ", levels[i]);
              j++;
              if ( j == 10 )
                {
                  fprintf(gdp, "\n");
                  j = 0;
                  if ( i != (nlevmax-1) ) fprintf(gdp, "               ");
                }
            }
        }
      if ( j ) fprintf(gdp, "\n");
    }

  free(levels);
}

static
void ctl_options(FILE *gdp, int yrev, int zrev, int sequential, int bigendian, int littleendian, int flt64, int cal365day)
{
  /* if ( filetype == FILETYPE_GRB ) zrev = FALSE; */

  if ( yrev || zrev || sequential || bigendian || littleendian || flt64 )
    {
      fprintf(gdp, "OPTIONS");
      if ( yrev )         fprintf(gdp, " yrev");
      if ( zrev )         fprintf(gdp, " zrev");
      if ( sequential )   fprintf(gdp, " sequential");
      if ( bigendian )    fprintf(gdp, " big_endian");
      if ( littleendian ) fprintf(gdp, " little_endian");
      if ( flt64 )        fprintf(gdp, " flt64");
      if ( cal365day )    fprintf(gdp, " 365_day_calendar");
      fprintf(gdp, "\n");
    }
}

static
void ctl_undef(FILE *gdp, int vlistID)
{
  double missval;

  missval = vlistInqVarMissval(vlistID, 0);
  fprintf(gdp, "UNDEF  %g\n", missval);
}

static
void ctl_vars(FILE *gdp, int filetype, int vlistID, int nvarsout, int *vars)
{
  int varID, nvars;
  int ltype, code;
  int zaxisID, nlev;
  int i, j;
  int len;
  char varname[CDI_MAX_NAME], varlongname[CDI_MAX_NAME], varunits[CDI_MAX_NAME];

  nvars   = vlistNvars(vlistID);

  fprintf(gdp, "VARS  %d\n", nvarsout);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vars[varID] == TRUE )
        {
          zaxisID = vlistInqVarZaxis(vlistID, varID);
          ltype   = zaxisInqLtype(zaxisID);
          nlev    = zaxisInqSize(zaxisID);
          vlistInqVarName(vlistID, varID, varname);

          len = (int) strlen(varname);
          for ( i = 0; i < len; i++ )
            if ( varname[i] == '-' ) break;

          if ( i < len )
            for ( j = i; j < len; j++ )
              varname[j] = varname[j+1];

          vlistInqVarLongname(vlistID, varID, varlongname);
          vlistInqVarUnits(vlistID, varID, varunits);
          fprintf(gdp, "%-15s", varname);

          if ( nlev == 1 ) nlev = 0;

          fprintf(gdp, "  %3d", nlev);

          if ( filetype == FILETYPE_GRB )
            {
              code = vlistInqVarCode(vlistID, varID);
              /*
              if      ( ltype == ZAXIS_SURFACE )  ltype = 1;
              else if ( ltype == ZAXIS_PRESSURE ) ltype = 99;
              else if ( nlev == 1 )  ltype = 1;
              else ltype = 99;
              */
              fprintf(gdp, "  %d,%d", code, ltype);
            }
          else if ( filetype == FILETYPE_NC )
            {
              int xyz = vlistInqVarXYZ(vlistID, varID);

              fprintf(gdp, "  ");
              if ( vlistInqVarTsteptype(vlistID, varID) != TSTEP_CONSTANT )
                fprintf(gdp, "t,");
              if ( xyz == 321 )
                {
                  if ( nlev > 0 ) fprintf(gdp, "z,");
                  fprintf(gdp, "y,x");
                }
              else if ( xyz == 312 )
                {
                  if ( nlev > 0 ) fprintf(gdp, "z,");
                  fprintf(gdp, "x,y");
                }
              else if ( xyz == 231 )
                {
                  fprintf(gdp, "y,");
                  if ( nlev > 0 ) fprintf(gdp, "z,");
                  fprintf(gdp, "x");
                }
              else if ( xyz == 132 )
                {
                  fprintf(gdp, "x,");
                  if ( nlev > 0 ) fprintf(gdp, "z,");
                  fprintf(gdp, "y");
                }
              else
                {
                  if ( nlev > 0 ) fprintf(gdp, "z,");
                  fprintf(gdp, "y,x");
                }
            }
          else
            fprintf(gdp, "  99");

          if ( varlongname[0] == 0 )
            fprintf(gdp, "  %s", varname);
          else
            fprintf(gdp, "  %s", varlongname);

          if ( varunits[0] != 0 )
            fprintf(gdp, "  [%s]", varunits);

          fprintf(gdp, "\n");
        }
    }

  fprintf(gdp, "ENDVARS\n");
}

static
void write_map_grib1(const char *ctlfile, int map_version, int nrecords, int *intnum, float *fltnum, off_t *bignum)
{
  int i;
  struct gaindx indx;
  struct gaindxb indxb;
  FILE *mapfp;
  int hinum[5];

  memset(&indx, 0, sizeof(struct gaindx));

  mapfp = fopen(ctlfile, "w");
  if ( mapfp == NULL ) cdoAbort("Open failed on %s", ctlfile);

  indx.type   = map_version;
  indx.hfnum  = 0;
  if ( map_version == 4 )
    {
      indx.hinum  = 5;
      indx.intnum = nrecords;
      indxb.bignum = 2 * nrecords;
    }
  else
    {
      indx.hinum  = 4;
      indx.intnum = 3 * nrecords;
      indxb.bignum = 0;
    }
  indx.fltnum = 3 * nrecords;

  indx.hipnt  = NULL;
  indx.hfpnt  = NULL;
  indx.intpnt = NULL;
  indx.fltpnt = NULL;
  indxb.bigpnt = NULL;

  hinum[0] = map_version;
  hinum[1] = 1;
  hinum[2] = nrecords;
  hinum[3] = 255;
  hinum[4] = indxb.bignum;

  if ( map_version == 2 )
    {
      int nb, bcnt, rc, j;
      float fdum;
      unsigned char *map;
      unsigned char ibmfloat[4];

      /* calculate the size of the ver==1 index file */

      nb = 2 + (indx.hinum*4) +  /* version in byte 2, then 4 ints with number of each data type */
        indx.hinum*sizeof(int) +
        indx.hfnum*sizeof(int) +
        indx.intnum*sizeof(int) +
        indx.fltnum*sizeof(float);

      /* add additional info */

      nb += 7;      /* base time (+ sec)  for compatibility with earlier version 2 maps */
      nb += 8*4;    /* grvals for time <-> grid conversion */

      map = (unsigned char*) malloc(nb);

      bcnt = 0;
      Put1Byte(map, bcnt, 0);
      Put1Byte(map, bcnt, map_version);

      Put4Byte(map, bcnt, indx.hinum);
      Put4Byte(map, bcnt, indx.hfnum);
      Put4Byte(map, bcnt, indx.intnum);
      Put4Byte(map, bcnt, indx.fltnum);

      Put2Byte(map, bcnt, 0);   /* initial year   */
      Put1Byte(map, bcnt, 0);   /* initial month  */
      Put1Byte(map, bcnt, 0);   /* initial day    */
      Put1Byte(map, bcnt, 0);   /* initial hour   */
      Put1Byte(map, bcnt, 0);   /* initial minute */
      Put1Byte(map, bcnt, 0);   /* initial second */

      if( indx.hinum )
        for ( i = 0; i < indx.hinum; i++ )
          Put4Byte(map, bcnt, hinum[i]);

      if( indx.hfnum ) {
        /* blank for now */
      }

      for ( i = 0; i < nrecords; i++ )
        {
          PutInt(map, bcnt, (int) bignum[i*2]);
          PutInt(map, bcnt, (int) bignum[i*2+1]);
          PutInt(map, bcnt, intnum[i]);
        }

      for ( i = 0; i < indx.fltnum; i++)
        {
          fdum= fltnum[i];
          rc = flt2ibm(fdum, ibmfloat);
          if ( rc < 0 ) cdoAbort("overflow in IBM float conversion");
          for ( j = 0; j < 4; j++ ) map[bcnt++] = ibmfloat[j];
        }

      /* write out the factors for converting from grid to absolute time */

      for ( i = 0; i < 8; i++)
        {
          fdum = 0;
          rc = flt2ibm(fdum, ibmfloat);
          if ( rc < 0 ) cdoAbort("overflow in IBM float conversion");
          for ( j = 0; j < 4; j++ ) map[bcnt++] = ibmfloat[j];
        }

      fwrite(map, 1, bcnt, mapfp);

      free(map);
    }
  else
    {
      fwrite(&indx, sizeof(struct gaindx), 1, mapfp);
      if ( indx.hinum > 0 )  fwrite(hinum, sizeof(int), indx.hinum, mapfp);
      if ( map_version == 1 )
        {
          int *intnumbuf;
          intnumbuf = (int*) malloc(indx.intnum*sizeof(int));
          for ( i = 0; i < nrecords; i++ )
            {
              intnumbuf[i*3+0] = (int) bignum[i*2];
              intnumbuf[i*3+1] = (int) bignum[i*2+1];
              intnumbuf[i*3+2] = intnum[i];
            }
          if ( indx.intnum > 0 ) fwrite(intnumbuf, sizeof(int), indx.intnum, mapfp);
          free(intnumbuf);
          if ( indx.fltnum > 0 ) fwrite(fltnum, sizeof(float), indx.fltnum, mapfp);
        }
      else
        {
          if ( indx.intnum  > 0 ) fwrite(intnum, sizeof(int), indx.intnum, mapfp);
          if ( indx.fltnum  > 0 ) fwrite(fltnum, sizeof(float), indx.fltnum, mapfp);
          if ( indxb.bignum > 0 ) fwrite(bignum, sizeof(off_t), indxb.bignum, mapfp);
        }
    }

  fclose(mapfp);
}

/*
static
void write_map_grib2(const char *ctlfile, int map_version, int nrecords, int *intnum, float *fltnum, off_t *bignum)
{
}
*/

void *Gradsdes(void *argument)
{
  int gridID = -1;
  int gridtype = -1;
  int index;
  int varID;
  int recID, levelID;
  int nrecs;
  int vdate, vtime;
  char *idxfile = NULL;
  char varname[CDI_MAX_NAME];
  int yrev = FALSE;
  int zrev = FALSE;
  int xyheader = 0;
  int nrecords = 0;
  int bigendian = FALSE, littleendian = FALSE;
  int flt64 = 0;
  int cal365day = 0;
  int sequential = FALSE;
  char Time[30], Incr[10] = {"1mn"}, *IncrKey[] = {"mn","hr","dy","mo","yr"};
  int isd, imn, ihh, iyy, imm, idd;
  int isds = 0, imns = 0, ihhs = 0, iyys = 0, imms = 0, idds = 0;
  int imn0 = 0, ihh0 = 0, iyy0 = 0, imm0 = 0, idd0 = 0;
  int idmn, idhh, idmm, idyy, iddd;
  int dt=1, iik=0, mdt = 0;
  int gridsize = 0;
  long checksize = 0;
  int nmiss;
  int prec;
  int map_version = 2;
  int maxrecs = 0;
  int monavg = -1;
  int *intnum = NULL;
  float *fltnum = NULL;
  off_t *bignum = NULL;
  double *array = NULL;
  const char *cmons[]={"jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"};

  cdoInitialize(argument);

  int GRADSDES = cdoOperatorAdd("gradsdes",  0, 0, NULL);
  int DUMPMAP  = cdoOperatorAdd("dumpmap",   0, 0, NULL);

  UNUSED(GRADSDES);

  int operatorID = cdoOperatorID();

  const char *datfile = cdoStreamName(0)->args;
  size_t len = strlen(datfile);
  char *ctlfile = (char *) malloc(len+10);
  strcpy(ctlfile, datfile);

  if ( cdoStreamName(0)->args[0] == '-' )
    cdoAbort("This operator does not work with pipes!");

  if ( operatorID == DUMPMAP )
    {
      dumpmap();
      goto END_LABEL;
    }

  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  if ( operatorArgc() == 1 )
    {
      map_version = parameter2int(operatorArgv()[0]);
      if ( map_version != 1 && map_version != 2 && map_version != 4 )
        cdoAbort("map_version=%d unsupported!", map_version);
    }
  else
    {
      if ( fileSize(cdoStreamName(0)->args) > 2147483647L ) map_version = 4;
    }

  if ( cdoVerbose ) cdoPrint("GrADS GRIB map version: %d", map_version);

  if ( map_version == 4 && sizeof(off_t) != 8 )
    cdoAbort("GrADS GRIB map version %d requires size of off_t to be 8! The size of off_t is %ld.",
             map_version, sizeof(off_t));

  int streamID = streamOpenRead(cdoStreamName(0));

  int vlistID = streamInqVlist(streamID);

  int nvars   = vlistNvars(vlistID);
  int ntsteps = vlistNtsteps(vlistID);
  int ngrids  = vlistNgrids(vlistID);

  int filetype  = streamInqFiletype(streamID);
  int byteorder = streamInqByteorder(streamID);

  if ( filetype == FILETYPE_NC2 || filetype == FILETYPE_NC4 ) filetype = FILETYPE_NC;

  if ( filetype != FILETYPE_SRV &&
       filetype != FILETYPE_EXT &&
       filetype != FILETYPE_IEG &&
       filetype != FILETYPE_GRB )
    {
      if ( filetype == FILETYPE_NC )
        //        cdoAbort("Unsupported file format: netCDF");
        ;
      else if ( filetype == FILETYPE_GRB2 )
        //cdoAbort("Unsupported file format: GRIB2");
        ;
      else
        cdoAbort("Unsupported file format!");
    }

  /* find the first lonlat or Gaussian grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID, index);
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_LONLAT   ||
           gridtype == GRID_GAUSSIAN ||
           gridtype == GRID_LCC  ) break;
    }

  if ( index == ngrids )
    cdoAbort("No Lon/Lat, Gaussian or Lambert grid found (%s data unsupported)!", gridNamePtr(gridtype));

  /* select all variables with used gridID */
  int *vars = (int*) malloc(nvars*sizeof(int));
  int *recoffset = (int*) malloc(nvars*sizeof(int));
  int nvarsout = 0;
  int nrecsout = 0;
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarGrid(vlistID, varID) == gridID )
        {
          if ( filetype == FILETYPE_SRV ||
               filetype == FILETYPE_EXT ||
               filetype == FILETYPE_IEG )
            {
              prec = vlistInqVarDatatype(vlistID, varID);
              if ( prec == DATATYPE_FLT64 ) flt64 = 1;
            }
          vars[varID] = TRUE;
          recoffset[varID] = nrecsout;
          nvarsout++;
          nrecsout += zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
          if ( ntsteps != 1 && ntsteps != 0 && vlistInqVarTsteptype(vlistID, varID) == TSTEP_CONSTANT )
            cdoAbort("Unsupported GrADS record structure! Variable %d has only 1 time step.",
                     vlistInqVarCode(vlistID, varID));
        }
      else
        {
          vlistInqVarName(vlistID, varID, varname);
          cdoPrint("Unsupported grid type >%s<, skipped variable %s!",
                   gridNamePtr(gridInqType(vlistInqVarGrid(vlistID, varID))), varname);
          vars[varID] = FALSE;
        }
    }

  if ( filetype != FILETYPE_GRB && nvars != nvarsout )
    cdoAbort("Too many different grids!");

  if ( filetype == FILETYPE_SRV )
    {
      xyheader = 40;
      if ( flt64 ) xyheader = 72;
      sequential = TRUE;
      if ( byteorder == CDI_BIGENDIAN )    bigendian = TRUE;
      if ( byteorder == CDI_LITTLEENDIAN ) littleendian = TRUE;
    }

  if ( filetype == FILETYPE_EXT )
    {
      xyheader = 24;
      if ( flt64 ) xyheader = 40;
      sequential = TRUE;
      if ( byteorder == CDI_BIGENDIAN )    bigendian = TRUE;
      if ( byteorder == CDI_LITTLEENDIAN ) littleendian = TRUE;
    }

  if ( filetype == FILETYPE_IEG )
    {
      xyheader = 644;
      if ( flt64 ) xyheader = 1048;
      sequential = TRUE;
      if ( byteorder == CDI_BIGENDIAN )    bigendian = TRUE;
      if ( byteorder == CDI_LITTLEENDIAN ) littleendian = TRUE;
    }

  /* ctl file name */
  repl_filetypeext(ctlfile, filetypeext(filetype), ".ctl");

  /* open ctl file*/
  FILE *gdp = fopen(ctlfile, "w");
  if ( gdp == NULL ) cdoAbort("Open failed on %s", ctlfile);

  /* VERSION */
  fprintf(gdp, "* Generated by CDO operator gradsdes\n");
  fprintf(gdp, "*\n");
  /*
#if defined(VERSION)
  fprintf(gdp, "* Generated by CDO version %s\n", VERSION);
  fprintf(gdp, "*\n");
#endif
  */
  /* DSET */
  if ( datfile[0] == '/' )
    fprintf(gdp, "DSET  %s\n", datfile);
  else
    {
      datfile = strrchr(datfile, '/');
      if ( datfile == 0 ) datfile = cdoStreamName(0)->args;
      else                datfile++;
      fprintf(gdp, "DSET  ^%s\n", datfile);
    }

  /*
   * DTYPE Print file type
   * INDEX Print filename of the control/index file .ctl/.idx
   */
  if ( filetype == FILETYPE_GRB ||  filetype == FILETYPE_GRB2 )
    {
      idxfile = strdup(ctlfile);
      char *pidxfile = idxfile;

      // print GRIB[12] file type
      // generate the index file
      if ( filetype == FILETYPE_GRB )
        {
          fprintf(gdp, "DTYPE  GRIB\n");
          repl_filetypeext(idxfile, ".ctl", ".gmp");
        }
      else if ( filetype == FILETYPE_GRB2 )
        {
          fprintf(gdp, "DTYPE  GRIB2\n");
          repl_filetypeext(pidxfile, ".ctl", ".idx");
        }

      // print file name of index file
      if ( datfile[0] == '/' )
        fprintf(gdp, "INDEX  %s\n", pidxfile);
      else
        {
          pidxfile = strrchr(pidxfile, '/');
          if ( pidxfile == 0 ) pidxfile = idxfile;
          else                 pidxfile++;
          fprintf(gdp, "INDEX  ^%s\n", pidxfile);
        }

      gridsize = vlistGridsizeMax(vlistID);
      array = (double*) malloc(gridsize*sizeof(double));
    }
  else if ( filetype == FILETYPE_NC )
    {
      fprintf(gdp, "DTYPE  netCDF\n");
    }

  /* XYHEADER */
  if ( xyheader ) fprintf(gdp, "XYHEADER  %d\n", xyheader);

  /* TIME */

  int taxisID = vlistInqTaxis(vlistID);

  if ( taxisInqCalendar(taxisID) == CALENDAR_365DAYS ) cal365day = 1;

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      if ( tsID == 0 )
        {
          cdiDecodeDate(vdate, &iyys, &imms, &idds);
          cdiDecodeTime(vtime, &ihhs, &imns, &isds);

          if ( imms < 1 || imms > 12 )  imms=1;

          ihh0 = ihhs;
          imn0 = imns;
          iyy0 = iyys;
          imm0 = imms;
          idd0 = idds;
        }

      if ( tsID == 1 )
        {
          cdiDecodeDate(vdate, &iyy, &imm, &idd);
          cdiDecodeTime(vtime, &ihh, &imn, &isd);

          idmn = imn - imns;
          idhh = ihh - ihhs;
          iddd = idd - idds;
          idmm = imm - imms;
          idyy = iyy - iyys;

          if ( idmn != 0 )
            {
              dt = idmn + (idhh + (iddd + (idmm*30 + idyy*12)*30)*24)*60;
            }
          else if ( idhh != 0 )
            {
              dt = idhh + (iddd + (idmm + idyy*12)*30)*24;
              iik = 1;
            }
          else if ( iddd != 0 )
            {
              dt = iddd + (idmm + idyy*12)*30;
              iik = 2;
            }
          else if ( idmm != 0 )
            {
              dt = idmm + idyy*12;
              iik = 3;
            }
          else if ( idyy != 0 )
            {
              dt = idyy;
              iik = 4;
            }

          if ( dt <= 0 ) dt = 1;
        }

      if ( tsID > 0 && tsID < 6 && iik != 3 && (monavg == TRUE || monavg == -1) )
        {
          cdiDecodeDate(vdate, &iyy, &imm, &idd);
          cdiDecodeTime(vtime, &ihh, &imn, &isd);

          idmn = imn - imns;
          idhh = ihh - ihhs;
          iddd = idd - idds;
          idmm = imm - imms;
          idyy = iyy - iyys;

          if ( iddd < 0 ) iddd *= -1;
          if ( idyy > 0 ) idmm += idyy*12;

          if ( /*idmn == 0 && idhh == 0 &&*/ (iddd == 0 || iddd == 1 || idd > 27 ) &&
               idmm > 0 && (mdt == 0 || idmm == mdt) )
            {
              mdt = idmm;
              monavg = TRUE;
            }
          else
            {
              monavg = FALSE;
            }
          /*          
          printf("monavg %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n",
                 tsID, monavg, mdt, imm , imms, idmm, iyy, iyys, idyy, idd, idds, iddd);
          */
          imns = imn;
          ihhs = ihh;
          idds = idd;
          imms = imm;
          iyys = iyy;
        }

      if ( filetype == FILETYPE_GRB )
        {
          nrecords += nrecsout;
          if ( nrecords >= maxrecs )
            {
              maxrecs = nrecords;
              intnum = (int*) realloc(intnum, 1*maxrecs*sizeof(int));
              fltnum = (float*) realloc(fltnum, 3*maxrecs*sizeof(float));
              bignum = (off_t*) realloc(bignum, 2*maxrecs*sizeof(off_t));
            }

          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(streamID, &varID, &levelID);
              if ( vars[varID] == TRUE )
                {
                  streamReadRecord(streamID, array, &nmiss);

                  index = (tsID*nrecsout + recoffset[varID] + levelID);

                  streamInqGRIBinfo(streamID, &intnum[index], &fltnum[index*3], &bignum[index*2]);

                  if ( map_version != 4 )
                    {
                      checksize = (long)bignum[index*2] + (long)gridsize*intnum[index]/8;
                      if ( checksize < 0L || checksize > 2147483647L )
                        {
                          nrecords -= nrecsout;
                          cdoWarning("File size limit reached for GrADS GRIB map_version=%d! Only the first %d time steps (2GB) are processed.", map_version, tsID);
                          goto LABEL_STOP;
                        }
                    }
                }
            }
        }

      tsID++;
    }

 LABEL_STOP:

  /* XYDEF */
  ctl_xydef(gdp, gridID, &yrev);

  /* ZDEF */
  ctl_zdef(gdp, vlistID, &zrev);

  /* TDEF */

  if ( monavg == TRUE )
    {
      dt = mdt;
      iik = 3;
      if ( idd0 > 28 )
        {
          /* int iddx = idd0; */
          idd0 = 1;
          cdoPrint("Reset start date to %02d:%02dZ%02d%s%04d",
                   ihh0, imn0, idd0, cmons[imm0-1], iyy0);
        }
    }

  sprintf(Time, "%02d:%02dZ%02d%s%04d", ihh0, imn0, idd0, cmons[imm0-1], iyy0);
  sprintf(Incr, "%d%s", dt, IncrKey[iik]);

  fprintf(gdp, "TDEF %d LINEAR %s %s\n", tsID, Time, Incr);

  /* TITLE */

  int xsize  = gridInqXsize(gridID);
  int ysize  = gridInqYsize(gridID);

  int res = 0;
  if ( gridtype == GRID_GAUSSIAN ) res = nlat2ntr(ysize);

  if ( res )
    fprintf(gdp, "TITLE  %s  T%d grid\n", datfile, res);
  else
    fprintf(gdp, "TITLE  %s  %dx%d grid\n", datfile, xsize, ysize);

  /* OPTIONS */
  ctl_options(gdp, yrev, zrev, sequential, bigendian, littleendian, flt64, cal365day);

  /* UNDEF */
  ctl_undef(gdp, vlistID);

  /* VARS */
  ctl_vars(gdp, filetype, vlistID, nvarsout, vars);


  /* INDEX file */
  if ( filetype == FILETYPE_GRB )
    {
      write_map_grib1(idxfile, map_version, nrecords, intnum, fltnum, bignum);
    }
  if ( filetype == FILETYPE_GRB2 )
    {
      cdoAbort("\nThe fileformat GRIB2 is not fully supported yet\nfor the gradsdes operator.\nThe .ctl file %s was generated.\nYou can add the necessary .idx file by running\n\tgribmap -i %s", ctlfile, ctlfile);
      // write_map_grib2(idxfile, map_version, nrecords, intnum, fltnum, bignum);
    }

  streamClose(streamID);

  if ( ctlfile ) free(ctlfile);
  if ( idxfile ) free(idxfile);

  if ( vars ) free(vars);
  if ( recoffset ) free(recoffset);
  if ( array ) free(array);
  if ( intnum ) free(intnum);
  if ( fltnum ) free(fltnum);

 END_LABEL:

  cdoFinish();

  return (0);
}
