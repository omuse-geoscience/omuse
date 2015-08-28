#if defined (HAVE_CONFIG_H)
#  include "../src/config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include <cdi.h>
int      vlistInqVarMissvalUsed(int vlistID, int varID);
#ifndef DBL_IS_NAN
#if  defined  (HAVE_DECL_ISNAN)
#  define DBL_IS_NAN(x)     (isnan(x))
#elif  defined  (FP_NAN)
#  define DBL_IS_NAN(x)     (fpclassify(x) == FP_NAN)
#else
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#endif

#ifndef DBL_IS_EQUAL
/*#define DBL_IS_EQUAL(x,y) (!(x < y || y < x)) */
#  define DBL_IS_EQUAL(x,y) (DBL_IS_NAN(x)||DBL_IS_NAN(y)?(DBL_IS_NAN(x)&&DBL_IS_NAN(y)?1:0):!(x < y || y < x))
#endif

#ifndef IS_EQUAL
#  define IS_NOT_EQUAL(x,y) (x < y || y < x)
#  define IS_EQUAL(x,y)     (!IS_NOT_EQUAL(x,y))
#endif


#include "printinfo.h"

void cdiDefTableID(int tableID);

int getopt(int argc, char *const argv[], const char *optstring);

extern char *optarg;
extern int optind, opterr, optopt;


char *Progname;

int DefaultFileType  = CDI_UNDEFID;
int DefaultDataType  = CDI_UNDEFID;
int DefaultByteorder = CDI_UNDEFID;

int comptype  = COMPRESS_NONE;  // Compression type
int complevel = 0;              // Compression level


static
void version(void)
{
  int   filetypes[] = {FILETYPE_SRV, FILETYPE_EXT, FILETYPE_IEG, FILETYPE_GRB, FILETYPE_GRB2, FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, FILETYPE_NC4C};
  char *typenames[] = {        "srv",        "ext",        "ieg",        "grb",        "grb2",        "nc",        "nc2",        "nc4",        "nc4c"};

  fprintf(stderr, "CDI version 1.8\n");
#if defined (COMPILER)
  fprintf(stderr, "Compiler: %s\n", COMPILER);
#endif
#if defined (COMP_VERSION)
  fprintf(stderr, " version: %s\n", COMP_VERSION);
#endif
#if defined (USER_NAME) && defined(HOST_NAME) && defined(SYSTEM_TYPE)
  fprintf(stderr, "Compiled: by %s on %s (%s) %s %s\n",
	  USER_NAME, HOST_NAME, SYSTEM_TYPE, __DATE__, __TIME__);
#endif

  fprintf(stderr, "filetype: ");
  for ( size_t i = 0; i < sizeof(filetypes)/sizeof(int); ++i )
    if ( cdiHaveFiletype(filetypes[i]) ) fprintf(stderr, "%s ", typenames[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "    with:");
#if defined (HAVE_LIBPTHREAD)
  fprintf(stderr, " PTHREADS");
#endif
#if defined (_OPENMP)
  fprintf(stderr, " OpenMP");
#endif
#if  defined  (HAVE_NETCDF4)
  fprintf(stderr, " NC4");
#endif
#if  defined  (HAVE_LIBNC_DAP)
  fprintf(stderr, " OPeNDAP");
#endif
#if defined (HAVE_LIBSZ)
  fprintf(stderr, " SZ");
#endif
#if defined (HAVE_LIBZ)
  fprintf(stderr, " Z");
#endif
#if defined (HAVE_LIBJASPER)
  fprintf(stderr, " JASPER");
#endif
#if defined (HAVE_LIBPROJ)
  fprintf(stderr, " PROJ.4");
#endif
#if defined (HAVE_LIBDRMAA)
  fprintf(stderr, " DRMAA");
#endif
#if defined (HAVE_LIBCURL)
  fprintf(stderr, " CURL");
#endif
  fprintf(stderr, "\n");
  cdiPrintVersion();
  fprintf(stderr, "\n");
/*
  1.0.0   6 Feb 2001 : initial version
  1.1.0  30 Jul 2003 : missing values implemented
  1.2.0   8 Aug 2003 : changes for CDI library version 0.7.0
  1.3.0  10 Feb 2004 : changes for CDI library version 0.7.9
  1.4.0   5 May 2004 : changes for CDI library version 0.8.1 (error handling)
  1.4.1  18 Sep 2004 : netCDF 2 support
  1.4.2  22 Mar 2005 : change level from int to double
  1.4.3  11 Apr 2005 : change date and time format to ISO
  1.5.0  22 Nov 2005 : IEG support
  1.5.1  21 Feb 2006 : added option -s for short info
  1.6.0   1 Aug 2006 : added option -z szip for SZIP compression of GRIB records
  1.6.1  27 Feb 2007 : short info with ltype for GENERIC zaxis
  1.6.2   3 Jan 2008 : changes for CDI library version 1.1.0 (compress)
  1.6.3  26 Mar 2008 : call streamDefTimestep also if ntsteps = 0 (buf fix)
  1.7.0  11 Apr 2008 : added option -z zip for deflate compression of netCDF4 variables
  1.7.1   1 Nov 2009 : added option -z jpeg for JPEG compression of GRIB2 records
  1.7.2  14 Nov 2012 : added optional compression level -z zip[_1-9]
*/
}

static
void usage(void)
{
  const char *name;
  int id;

  fprintf(stderr, "usage : %s  [Option]  [ifile]  [ofile]\n", Progname);

  fprintf(stderr, "\n");
  fprintf(stderr, "  Options:\n");
  fprintf(stderr, "    -d             Print debugging information\n");
  fprintf(stderr, "    -f <format>    Format of the output file. (grb, grb2, nc, nc2, nc4, nc4c, src, ext or ieg)\n");
  fprintf(stderr, "    -s             give short information if ofile is missing\n");
  fprintf(stderr, "    -t <table>     Parameter table name/file\n");
  fprintf(stderr, "                   Predefined tables: ");
  for ( id = 0; id < tableInqNumber(); id++ )
    if ( (name = tableInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");

  fprintf(stderr, "    -V             Print version number\n");
  fprintf(stderr, "    -z szip        SZIP compression of GRIB1 records\n");
  fprintf(stderr, "       jpeg        JPEG compression of GRIB2 records\n");
  fprintf(stderr, "        zip[_1-9]  Deflate compression of netCDF4 variables\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Report bugs to <http://code.zmaw.de/projects/cdi>\n");
}


static
void printInfo(int vdate, int vtime, char *varname, double level,
	       int datasize, int number, int nmiss, double missval, const double *data, int vardis)
{
  static int rec = 0;
  int i, ivals = 0, imiss = 0;
  double arrmean, arrmin, arrmax;
  char vdatestr[32], vtimestr[32];

  if ( ! rec )
  {
    if ( vardis )
      fprintf(stdout,
    "   Rec :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name\n");
/*   ----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+ */
    else
      fprintf(stdout,
    "   Rec :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter ID\n");
/*   ----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+ */
  }

  date2str(vdate, vdatestr, sizeof(vdatestr));
  time2str(vtime, vtimestr, sizeof(vtimestr));

  fprintf(stdout, "%6d :%s %s %7g ", ++rec, vdatestr, vtimestr, level);

  fprintf(stdout, "%8d ", datasize);

  fprintf(stdout, "%7d :", nmiss);

  if ( number == CDI_REAL )
    {
      if ( nmiss > 0 )
	{
	  arrmean = 0;
	  arrmin  =  1.e300;
	  arrmax  = -1.e300;
	  for ( i = 0; i < datasize; i++ )
	    {
	      if ( !DBL_IS_EQUAL(data[i], missval) )
		{
		  if ( data[i] < arrmin ) arrmin = data[i];
		  if ( data[i] > arrmax ) arrmax = data[i];
		  arrmean += data[i];
		  ivals++;
		}
	    }
	  imiss = datasize - ivals;
	  datasize = ivals;
	}
      else
	{
	  arrmean = data[0];
	  arrmin  = data[0];
	  arrmax  = data[0];
	  for ( i = 1; i < datasize; i++ )
	    {
	      if ( data[i] < arrmin ) arrmin = data[i];
	      if ( data[i] > arrmax ) arrmax = data[i];
	      arrmean += data[i];
	    }
	}

      if ( datasize > 0 ) arrmean /= datasize;

      fprintf(stdout, "%#12.5g%#12.5g%#12.5g", arrmin, arrmean, arrmax);
    }
  else
    {
      int nvals_r = 0, nvals_i = 0;
      double arrsum_r, arrsum_i, arrmean_r = 0, arrmean_i = 0;
      arrsum_r = 0;
      arrsum_i = 0;

      for ( i = 0; i < datasize; i++ )
	{
	  if ( !DBL_IS_EQUAL(data[i*2], missval) )
	    {
	      arrsum_r += data[i*2];
	      nvals_r++;
	    }
	  if ( !DBL_IS_EQUAL(data[i*2+1], missval) )
	    {
	      arrsum_i += data[i*2+1];
	      nvals_i++;
	    }
	}

      imiss = datasize - nvals_r;

      if ( nvals_r > 0 ) arrmean_r = arrsum_r / nvals_r;
      if ( nvals_i > 0 ) arrmean_i = arrsum_i / nvals_i;
      fprintf(stdout, "  -  (%#12.5g,%#12.5g)  -", arrmean_r, arrmean_i);
    }

  fprintf(stdout, " : %-11s\n", varname);

  if ( imiss != nmiss && nmiss > 0 )
    fprintf(stdout, "Found %d of %d missing values!\n", imiss, nmiss);
}

#define MAXCHARS 82

const char * tunit2str(int tunits)
{
  if      ( tunits == TUNIT_YEAR )       return ("years");
  else if ( tunits == TUNIT_MONTH )      return ("months");
  else if ( tunits == TUNIT_DAY )        return ("days");
  else if ( tunits == TUNIT_12HOURS )    return ("12hours");
  else if ( tunits == TUNIT_6HOURS )     return ("6hours");
  else if ( tunits == TUNIT_3HOURS )     return ("3hours");
  else if ( tunits == TUNIT_HOUR )       return ("hours");
  else if ( tunits == TUNIT_30MINUTES )  return ("30minutes");
  else if ( tunits == TUNIT_QUARTER )    return ("15minutes");
  else if ( tunits == TUNIT_MINUTE )     return ("minutes");
  else if ( tunits == TUNIT_SECOND )     return ("seconds");
  else                                   return ("unknown");
}


const char* calendar2str(int calendar)
{
  if      ( calendar == CALENDAR_STANDARD )  return ("standard");
  else if ( calendar == CALENDAR_PROLEPTIC ) return ("proleptic_gregorian");
  else if ( calendar == CALENDAR_360DAYS )   return ("360_day");
  else if ( calendar == CALENDAR_365DAYS )   return ("365_day");
  else if ( calendar == CALENDAR_366DAYS )   return ("366_day");
  else                                       return ("unknown");
}

static
void limit_string_length(char* string, size_t maxlen)
{
  string[maxlen-1] = 0;
  size_t len = strlen(string);

  if ( len > 10 )
    {
      for ( size_t i = 3; i < len; ++i )
	if ( string[i] == ' ' || string[i] == ',' || (i>10 && string[i] == '.') )
	  {
	    string[i] = 0;
	    break;
	  }
    }
}

static
void printShortinfo(int streamID, int vlistID, int vardis)
{
  int varID;
  int gridsize = 0;
  int gridID, zaxisID, param;
  int vdate, vtime;
  int ntsteps;
  int levelsize;
  int tsteptype, taxisID;
  char tmpname[CDI_MAX_NAME];
  char varname[CDI_MAX_NAME];
  const char *modelptr, *instptr;
  int datatype;
  int year, month, day, hour, minute, second;
  char pstr[4];
  char paramstr[32];

      fprintf(stdout, "   File format");
      fprintf(stdout, " : ");
      printFiletype(streamID, vlistID);

      //vlistPrint(vlistID);
      int nvars = vlistNvars(vlistID);
      int nsubtypes = vlistNsubtypes(vlistID);

      if ( nsubtypes > 0 )
        fprintf(stdout, "   Var : Institut Source   Ttype    Subtypes Num  Levels Num  Gridsize Num Dtype : ");
      else
        fprintf(stdout, "   Var : Institut Source   Ttype    Levels Num  Gridsize Num Dtype : ");

      if ( vardis )
	fprintf(stdout, "Parameter name\n");
      else
	fprintf(stdout, "Parameter ID\n");

      for ( varID = 0; varID < nvars; varID++ )
	{
	  param   = vlistInqVarParam(vlistID, varID);
	  gridID  = vlistInqVarGrid(vlistID, varID);
	  zaxisID = vlistInqVarZaxis(vlistID, varID);

	  fprintf(stdout, "%6d : ", varID + 1);

	  /* institute info */
	  instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
	  strcpy(tmpname, "unknown");
	  if ( instptr ) strncpy(tmpname, instptr, CDI_MAX_NAME);
	  limit_string_length(tmpname, CDI_MAX_NAME);
	  fprintf(stdout, "%-8s ", tmpname);

	  /* source info */
	  modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
	  strcpy(tmpname, "unknown");
	  if ( modelptr ) strncpy(tmpname, modelptr, CDI_MAX_NAME);
	  limit_string_length(tmpname, CDI_MAX_NAME);
	  fprintf(stdout, "%-8s ", tmpname);

	  /* tsteptype */
	  tsteptype = vlistInqVarTsteptype(vlistID, varID);
	  if      ( tsteptype == TSTEP_CONSTANT ) fprintf(stdout, "%-8s ", "constant");
	  else if ( tsteptype == TSTEP_INSTANT  ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_INSTANT2 ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_INSTANT3 ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_MIN      ) fprintf(stdout, "%-8s ", "min");
	  else if ( tsteptype == TSTEP_MAX      ) fprintf(stdout, "%-8s ", "max");
	  else if ( tsteptype == TSTEP_AVG      ) fprintf(stdout, "%-8s ", "avg");
	  else if ( tsteptype == TSTEP_ACCUM    ) fprintf(stdout, "%-8s ", "accum");
	  else if ( tsteptype == TSTEP_RANGE    ) fprintf(stdout, "%-8s ", "range");
	  else if ( tsteptype == TSTEP_DIFF     ) fprintf(stdout, "%-8s ", "diff");
	  else                                    fprintf(stdout, "%-8s ", "unknown");

          if ( nsubtypes > 0 )
            {
              int subtypeID = vlistInqVarSubtype(vlistID, varID);
              int subtypesize = subtypeInqSize(subtypeID);
              fprintf(stdout, " %6d  ", subtypesize);
              fprintf(stdout, "%3d ", vlistSubtypeIndex(vlistID, subtypeID) + 1);
            }

	  /* layer info */
	  levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, "%6d ", levelsize);
	  fprintf(stdout, "%3d ", vlistZaxisIndex(vlistID, zaxisID) + 1);

	  /* grid info */
	  gridsize = gridInqSize(gridID);
	  fprintf(stdout, "%9d ", gridsize);
	  fprintf(stdout, "%3d ", vlistGridIndex(vlistID, gridID) + 1);

	  /* datatype */
	  datatype = vlistInqVarDatatype(vlistID, varID);
	  if      ( datatype == DATATYPE_PACK   ) strcpy(pstr, "P0");
	  else if ( datatype > 0 && datatype <= 32  ) sprintf(pstr, "P%d", datatype);
	  else if ( datatype == DATATYPE_CPX32  ) strcpy(pstr, "C32");
	  else if ( datatype == DATATYPE_CPX64  ) strcpy(pstr, "C64");
	  else if ( datatype == DATATYPE_FLT32  ) strcpy(pstr, "F32");
	  else if ( datatype == DATATYPE_FLT64  ) strcpy(pstr, "F64");
	  else if ( datatype == DATATYPE_INT8   ) strcpy(pstr, "I8");
	  else if ( datatype == DATATYPE_INT16  ) strcpy(pstr, "I16");
	  else if ( datatype == DATATYPE_INT32  ) strcpy(pstr, "I32");
	  else if ( datatype == DATATYPE_UINT8  ) strcpy(pstr, "U8");
	  else if ( datatype == DATATYPE_UINT16 ) strcpy(pstr, "U16");
	  else if ( datatype == DATATYPE_UINT32 ) strcpy(pstr, "U32");
	  else                                    strcpy(pstr, "-1");

	  fprintf(stdout, " %-3s", pstr);

	  if ( vlistInqVarCompType(vlistID, varID) == COMPRESS_NONE )
	    fprintf(stdout, "  ");
	  else
	    fprintf(stdout, "z ");

	  /* parameter info */
	  fprintf(stdout, ": ");

	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  if (vardis)
            {
              vlistInqVarName(vlistID, varID, varname);
              fprintf(stdout, "%-11s", varname);
            }
	  else
	    fprintf(stdout, "%-11s", paramstr);

	  fprintf(stdout, "\n");
	}

      fprintf(stdout, "   Grid coordinates");
      fprintf(stdout, " :\n");

      printGridInfo(vlistID);

      fprintf(stdout, "   Vertical coordinates");
      fprintf(stdout, " :\n");

      printZaxisInfo(vlistID);

      if ( nsubtypes > 0 )
        {
          fprintf(stdout, "   Subtypes");
          fprintf(stdout, " :\n");

          printSubtypeInfo(vlistID);
        }

      taxisID = vlistInqTaxis(vlistID);
      ntsteps = vlistNtsteps(vlistID);

      if ( ntsteps != 0 )
	{
	  if ( ntsteps == CDI_UNDEFID )
	    fprintf(stdout, "   Time coordinate :  unlimited steps\n");
	  else
	    fprintf(stdout, "   Time coordinate :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

	  if ( taxisID != CDI_UNDEFID )
	    {
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
                  int calendar, tunits;

		  vdate = taxisInqRdate(taxisID);
		  vtime = taxisInqRtime(taxisID);

		  cdiDecodeDate(vdate, &year, &month, &day);
		  cdiDecodeTime(vtime, &hour, &minute, &second);

		  fprintf(stdout, "     RefTime = %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
			  year, month, day, hour, minute, second);

		  tunits = taxisInqTunit(taxisID);
		  if ( tunits != CDI_UNDEFID )  fprintf(stdout, "  Units = %s", tunit2str(tunits));

		  calendar = taxisInqCalendar(taxisID);
		  if ( calendar != CDI_UNDEFID )  fprintf(stdout, "  Calendar = %s", calendar2str(calendar));

		  if ( taxisHasBounds(taxisID) )
		    fprintf(stdout, "  Bounds = true");

		  fprintf(stdout, "\n");
		}
	    }

	  fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

	  printTimesteps(streamID, taxisID, 0);

	  fprintf(stdout, "\n");
	}
}


#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )


static
void setDefaultDataType(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int nbits = -1;
  enum {D_UINT, D_INT, D_FLT, D_CPX};
  int dtype = -1;

  if      ( *datatypestr == 'i' || *datatypestr == 'I' )
    {
      dtype = D_INT;
      datatypestr++;
    }
  else if ( *datatypestr == 'u' || *datatypestr == 'U' )
    {
      dtype = D_UINT;
      datatypestr++;
    }
  else if ( *datatypestr == 'f' || *datatypestr == 'F' )
    {
      dtype = D_FLT;
      datatypestr++;
    }
  else if ( *datatypestr == 'c' || *datatypestr == 'C' )
    {
      dtype = D_CPX;
      datatypestr++;
    }

  if ( isdigit((int) *datatypestr) )
    {
      nbits = atoi(datatypestr);
      if ( nbits < 10 )
	datatypestr += 1;
      else
	datatypestr += 2;

      if ( dtype == -1 )
	{
	  if      ( nbits > 0 && nbits < 32 ) DefaultDataType = nbits;
	  else if ( nbits == 32 )
	    {
	      if ( DefaultFileType == FILETYPE_GRB )
		DefaultDataType = DATATYPE_PACK32;
	      else
		DefaultDataType = DATATYPE_FLT32;
	    }
	  else if ( nbits == 64 ) DefaultDataType = DATATYPE_FLT64;
	  else
	    {
	      fprintf(stderr, "Unsupported number of bits %d!\n", nbits);
	      fprintf(stderr, "Use I8/I16/I32/F32/F64 for nc/nc2/nc4/nc4c; F32/F64 for grb2/srv/ext/ieg; P1 - P24 for grb/grb2.\n");
	      exit(EXIT_FAILURE);
	    }
	}
      else
	{
	  if ( dtype == D_INT )
	    {
	      if      ( nbits ==  8 ) DefaultDataType = DATATYPE_INT8;
	      else if ( nbits == 16 ) DefaultDataType = DATATYPE_INT16;
	      else if ( nbits == 32 ) DefaultDataType = DATATYPE_INT32;
	      else
		{
		  fprintf(stderr, "Unsupported number of bits = %d for datatype INT!\n", nbits);
		  exit(EXIT_FAILURE);
		}
	    }
	  /*
	  else if ( dtype == D_UINT )
	    {
	      if      ( nbits ==  8 ) DefaultDataType = DATATYPE_UINT8;
	      else
		{
		  fprintf(stderr, "Unsupported number of bits = %d for datatype UINT!\n", nbits);
		  exit(EXIT_FAILURE);
		}
	    }
	  */
	  else if ( dtype == D_FLT )
	    {
	      if      ( nbits == 32 ) DefaultDataType = DATATYPE_FLT32;
	      else if ( nbits == 64 ) DefaultDataType = DATATYPE_FLT64;
	      else
		{
		  fprintf(stderr, "Unsupported number of bits = %d for datatype FLT!\n", nbits);
		  exit(EXIT_FAILURE);
		}
	    }
	  else if ( dtype == D_CPX )
	    {
	      if      ( nbits == 32 ) DefaultDataType = DATATYPE_CPX32;
	      else if ( nbits == 64 ) DefaultDataType = DATATYPE_CPX64;
	      else
		{
		  fprintf(stderr, "Unsupported number of bits = %d for datatype CPX!\n", nbits);
		  exit(EXIT_FAILURE);
		}
	    }
	}
    }

  if ( *datatypestr != 0 )
    {
      if ( *datatypestr == 'l' || *datatypestr == 'L' )
	{
	  if ( IsBigendian() ) DefaultByteorder = CDI_LITTLEENDIAN;
	  datatypestr++;
	}
      else if ( *datatypestr == 'b' || *datatypestr == 'B' )
	{
	  if ( ! IsBigendian() ) DefaultByteorder = CDI_BIGENDIAN;
	  datatypestr++;
	}
      else
	{
	  fprintf(stderr, "Unsupported character in number of bytes: >%s< !\n", datatypestr);
	  exit(EXIT_FAILURE);
	}
    }
}

static
void setDefaultFileType(char *filetypestr)
{
  if ( filetypestr )
    {
      char *ftstr = filetypestr;

      if      ( memcmp(filetypestr, "grb2", 4)  == 0 ) { ftstr += 4; DefaultFileType = FILETYPE_GRB2;}
      else if ( memcmp(filetypestr, "grb1", 4)  == 0 ) { ftstr += 4; DefaultFileType = FILETYPE_GRB; }
      else if ( memcmp(filetypestr, "grb",  3)  == 0 ) { ftstr += 3; DefaultFileType = FILETYPE_GRB; }
      else if ( memcmp(filetypestr, "nc2",  3)  == 0 ) { ftstr += 3; DefaultFileType = FILETYPE_NC2; }
      else if ( memcmp(filetypestr, "nc4c", 4)  == 0 ) { ftstr += 4; DefaultFileType = FILETYPE_NC4C;}
      else if ( memcmp(filetypestr, "nc4",  3)  == 0 ) { ftstr += 3; DefaultFileType = FILETYPE_NC4; }
      else if ( memcmp(filetypestr, "nc",   2)  == 0 ) { ftstr += 2; DefaultFileType = FILETYPE_NC;  }
      else if ( memcmp(filetypestr, "srv",  3)  == 0 ) { ftstr += 3; DefaultFileType = FILETYPE_SRV; }
      else if ( memcmp(filetypestr, "ext",  3)  == 0 ) { ftstr += 3; DefaultFileType = FILETYPE_EXT; }
      else if ( memcmp(filetypestr, "ieg",  3)  == 0 ) { ftstr += 3; DefaultFileType = FILETYPE_IEG; }
      else
	{
	  fprintf(stderr, "Unsupported filetype %s!\n", filetypestr);
	  fprintf(stderr, "Available filetypes: grb, grb2, nc, nc2, nc4, nc4c, srv, ext and ieg\n");
	  exit(EXIT_FAILURE);
	}

      if ( DefaultFileType != CDI_UNDEFID && *ftstr != 0 )
	{
	  if ( *ftstr == '_' )
	    {
	      ftstr++;

	      setDefaultDataType(ftstr);
	    }
	  else
	    {
	      fprintf(stderr, "Unexpected character >%c< in file type >%s<!\n", *ftstr, filetypestr);
	      fprintf(stderr, "Use format[_nbits] with:\n");
	      fprintf(stderr, "    format = grb, grb2, nc, nc2, nc4, nc4c, srv, ext or ieg\n");
	      fprintf(stderr, "    nbits  = 32/64 for grb2/nc/nc2/nc4/nc4c/srv/ext/ieg; 1 - 24 for grb/grb2\n");
	      exit(EXIT_FAILURE);
	    }
	}
    }
}

static
int handle_error(int cdiErrno, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  printf("\n");
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  return (cdiErrno);
}

static
void defineCompress(const char *arg)
{
  size_t len = strlen(arg);

  if      ( strncmp(arg, "szip", len) == 0 )
    {
      comptype = COMPRESS_SZIP;
    }
  else if ( strncmp(arg, "jpeg", len) == 0 )
    {
      comptype = COMPRESS_JPEG;
    }
  else if ( strncmp(arg, "gzip", len) == 0 )
    {
      comptype = COMPRESS_GZIP;
      complevel = 6;
    }
  else if ( strncmp(arg, "zip", 3) == 0 )
    {
      comptype = COMPRESS_ZIP;
      if ( len == 5 && arg[3] == '_' && isdigit(arg[4]) )
	complevel = atoi(&arg[4]);
      else
        complevel = 1;
    }
  else
    fprintf(stderr, "%s compression unsupported!\n", arg);
}



int main(int argc, char *argv[])
{
  int c;
  char *fname1 = NULL;
  char *fname2 = NULL;
  char *rTable = NULL;
  char *wTable = NULL;
  int Move = 0;
  int Record = 0;
  int Debug = 0;
  /* int Quiet  = 0; */
  int Vardis = 0;
  int Version = 0;
  int Longinfo = 0;
  int Shortinfo = 0;
  int varID;
  int itableID = CDI_UNDEFID, otableID = CDI_UNDEFID;
  int Info = 1;
  char varname[CDI_MAX_NAME];
  char paramstr[32];

  Progname = strrchr(argv[0], '/');
  if (Progname == 0) Progname = argv[0];
  else               Progname++;

  while ( (c = getopt(argc, argv, "b:f:t:w:z:cdhlMmqRrsvVZ")) != EOF )
    {
      switch (c)
	{
	case 'b':
	  setDefaultDataType(optarg);
	  break;
	case 'd':
	  Debug = 1;
	  break;
	case 'f':
	  setDefaultFileType(optarg);
	  break;
	case 'h':
	  usage();
	  exit (0);
	case 'l':
	  Longinfo = 1;
	  break;
	case 'M':
	  cdiDefGlobal("HAVE_MISSVAL", 1);
	  break;
	case 'm':
	  Move = 1;
	  break;
	case 'q':
	  /* Quiet = 1; */
	  break;
	case 'R':
	  cdiDefGlobal("REGULARGRID", 1);
	  break;
	case 'r':
	  Record = 1;
	  break;
	case 's':
	  Shortinfo = 1;
	  break;
	case 't':
	  rTable = optarg;
	  break;
	case 'v':
	  Vardis = 1;
	  break;
	case 'V':
	  Version = 1;
	  break;
	case 'w':
	  wTable = optarg;
	  break;
	case 'z':
	  defineCompress(optarg);
	  break;
	}
    }

  if ( optind < argc ) fname1 = argv[optind++];
  if ( optind < argc ) fname2 = argv[optind++];
  if ( optind < argc ) fprintf(stderr, "optind: %d argc: %d\n", optind, argc);

  if ( Debug || Version ) version();

  if ( Debug ) cdiDebug(Debug);

  if ( rTable )
    {
      itableID = tableInq(-1, 0, rTable);
      if ( itableID != CDI_UNDEFID ) cdiDefTableID(itableID);
      otableID = itableID;
    }

  if ( fname1 == NULL && ! (Debug || Version) )
    {
      usage();
      exit (0);
    }

  if ( fname1 )
    {
      double *data = NULL;
      double missval;
      double level;
      int nmiss;
      int number;
      int datasize = 0;
      int streamID1 = CDI_UNDEFID;
      int streamID2 = CDI_UNDEFID;
      int filetype;
      int gridID, zaxisID;
      int param;
      int vdate, vtime;
      int nrecs, nvars;
      int levelID, levelsize;
      int nts = 0;
      int gridsize = 0;
      int recID;
      int tsID;
      int ntsteps = 0;
      int taxisID1, taxisID2 = CDI_UNDEFID;
      int vlistID1, vlistID2 = CDI_UNDEFID;

      streamID1 = streamOpenRead(fname1);
      if ( streamID1 < 0 )
	return (handle_error(streamID1, "Open failed on %s", fname1));

      vlistID1 = streamInqVlist(streamID1);

      if ( Longinfo )
	{
	  int ngrids, nzaxis;
	  vlistPrint(vlistID1);
	  ngrids = vlistNgrids(vlistID1);
	  nzaxis = vlistNzaxis(vlistID1);
	  for ( gridID = 0; gridID < ngrids; gridID++ ) gridPrint(gridID, gridID, 1);
	  for ( zaxisID = 0; zaxisID < nzaxis; zaxisID++ ) zaxisPrint(zaxisID, zaxisID);
	}

      nvars   = vlistNvars(vlistID1);
      taxisID1 = vlistInqTaxis(vlistID1);
      ntsteps = vlistNtsteps(vlistID1);

      if (Debug)
        fprintf(stderr, "nvars   = %d\n"
                "ntsteps = %d\n", nvars, ntsteps);

      if ( fname2 )
        {
          vlistID2 = vlistDuplicate(vlistID1);
          taxisID2 = taxisDuplicate(taxisID1);
          vlistDefTaxis(vlistID2, taxisID2);
        }

      for ( varID = 0; varID < nvars; varID++)
	{
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  if ( gridsize > datasize ) datasize = gridsize;
	  if ( fname2 )
	    {
	      if ( DefaultDataType != CDI_UNDEFID )
		vlistDefVarDatatype(vlistID2, varID, DefaultDataType);
	    }
	}

      if ( fname2 )
	{
	  Info = 0;
	  filetype = streamInqFiletype(streamID1);

	  if ( DefaultFileType != CDI_UNDEFID )
	    filetype = DefaultFileType;

	  streamID2 = streamOpenWrite(fname2, filetype);
	  if ( streamID2 < 0 )
	    return (handle_error(streamID2, "Open failed on %s", fname2));

	  if ( DefaultByteorder != CDI_UNDEFID )
	    streamDefByteorder(streamID2, DefaultByteorder);

	  if ( comptype != COMPRESS_NONE )
	    {
	      streamDefCompType(streamID2, comptype);
	      streamDefCompLevel(streamID2, complevel);
	    }

	  streamDefVlist(streamID2, vlistID2);

	  if ( otableID == CDI_UNDEFID ) otableID = itableID;
	}

      if ( vlistNumber(vlistID1) != CDI_REAL ) datasize *= 2;
      data = (double *) malloc((size_t)datasize * sizeof (double));

      /*
	nts = cdiInqTimeSize(streamID1);
      */
      if (Debug)
	printf("nts = %d\n"
               "streamID1 = %d, streamID2 = %d\n",
               nts, streamID1, streamID2);

      if ( Shortinfo )
	{
	  Info = 0;
	  printShortinfo(streamID1, vlistID1, Vardis);
	}

      tsID = 0;
      if ( Info || fname2 )
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) > 0 )
	{
	  if ( fname2 /* && ntsteps != 0*/ )
            {
              taxisCopyTimestep(taxisID2, taxisID1);
              streamDefTimestep(streamID2, tsID);
            }
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  if ( Debug )
	    fprintf(stdout, "tsID = %d nrecs = %d date = %d time = %d\n", tsID, nrecs, vdate, vtime);

	  if ( Record )
	    {
	      for ( recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  streamReadRecord(streamID1, data, &nmiss);

		  number   = vlistInqVarNumber(vlistID1, varID);
		  gridID   = vlistInqVarGrid(vlistID1, varID);
		  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
		  param    = vlistInqVarParam(vlistID1, varID);

		  cdiParamToString(param, paramstr, sizeof(paramstr));

		  if ( Vardis ) vlistInqVarName(vlistID1, varID, varname);
		  else          strcpy(varname, paramstr);
		  /*
		  printf("varID=%d, param=%d, gridID=%d, zaxisID=%d levelID=%d\n",
			 varID, param, gridID, zaxisID, levelID);
		  */
		  gridsize = gridInqSize(gridID);
		  level    = zaxisInqLevel(zaxisID, levelID);
		  missval  = vlistInqVarMissval(vlistID1, varID);

		  if ( Info )
		    printInfo(vdate, vtime, varname, level, gridsize, number, nmiss, missval, data, Vardis);

		  if ( fname2 )
		    {
		      streamDefRecord(streamID2, varID, levelID);
		      if ( Move )
			streamCopyRecord(streamID2, streamID1);
		      else
			streamWriteRecord(streamID2, data, nmiss);
		    }
	      	}
	    }
	  else
	    {
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT && tsID > 0 ) continue;

		  number   = vlistInqVarNumber(vlistID1, varID);
		  gridID   = vlistInqVarGrid(vlistID1, varID);
		  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
		  param    = vlistInqVarParam(vlistID1, varID);

		  cdiParamToString(param, paramstr, sizeof(paramstr));

		  if ( Vardis ) vlistInqVarName(vlistID1, varID, varname);
		  else          strcpy(varname, paramstr);

		  if ( Debug )
		    fprintf(stdout, "varID = %d param = %d gridID = %d zaxisID = %d\n",
			    varID, param, gridID, zaxisID);

		  gridsize = gridInqSize(gridID);
		  missval  = vlistInqVarMissval(vlistID1, varID);

		  levelsize = zaxisInqSize(zaxisID);
		  for ( levelID = 0; levelID < levelsize; levelID++ )
		    {
		      level = zaxisInqLevel(zaxisID, levelID);
		      streamReadVarSlice(streamID1, varID, levelID, data, &nmiss);

		      if ( Info )
			printInfo(vdate, vtime, varname, level, gridsize, number, nmiss, missval, data, Vardis);

		      if ( fname2 )
			streamWriteVarSlice(streamID2, varID, levelID, data, nmiss);
		    }
		}
	    }
	  tsID++;
        }

      free(data);

      /* fprintf(stderr, "%ld\n", (long) streamNvals(streamID1)); */

      if ( fname2 )
	{
	  streamClose(streamID2);
	  vlistDestroy(vlistID2);
	  taxisDestroy(taxisID2);
	}
      streamClose(streamID1);
    }

  if ( wTable )
    tableWrite(wTable, itableID);

  return (0);
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
