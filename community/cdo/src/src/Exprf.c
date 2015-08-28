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

      Exprf      expr            Evaluate expressions
      Exprf      exprf           Evaluate expressions from script file
      Exprf      aexpr           Append evaluated expressions
      Exprf      aexprf          Append evaluated expressions from script file
*/
/*
Operatoren: +, -, *, \, ^
Functions: sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
Functions: min, max, avg, std, var
Constansts: M_PI, M_E
*/

#include <sys/types.h> /* stat */
#include <sys/stat.h>  /* stat */
#include <unistd.h>    /* stat */

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "expr.h"


void *Expr(void *argument)
{
  int operatorID;
  char *exprs = NULL;
  const char *exprf = NULL;
  int streamID1, streamID2 = CDI_UNDEFID;
  int offset;
  int nrecs, nvars, nvars1, nvars2;
  int gridID, zaxisID;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int gridsize, nlevel;
  int nmiss;
  int taxisID1, taxisID2;
  //int lwarn = TRUE;
  double missval;
  double *array = NULL;
  double *single1, *single2;
  parse_parm_t parse_arg;
  void *scanner;
  int yy_scan_string(const char *str, void *scanner);

  cdoInitialize(argument);

  yylex_init(&scanner);
  yyset_extra(&parse_arg, scanner);


# define REPLACES_VARIABLES(id) cdoOperatorF1(id)
# define READS_COMMAND_LINE(id) cdoOperatorF2(id)

  cdoOperatorAdd("expr",   1, 1, "expressions");
  cdoOperatorAdd("exprf",  1, 0, "expr script filename");
  cdoOperatorAdd("aexpr",  0, 1, "expressions");
  cdoOperatorAdd("aexprf", 0, 0, "expr script filename");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( READS_COMMAND_LINE(operatorID) )
    {
      size_t slen;

      slen = strlen(operatorArgv()[0]);
      exprs = (char*) malloc(slen+2);
      strcpy(exprs, operatorArgv()[0]);
      if ( exprs[slen-1] != ';' )
	{
	  exprs[slen]   = ';';
	  exprs[slen+1] = 0;
	}
    }
  else
    {
      int ichar, ipos = 0;
      FILE *fp;
      size_t fsize;
      struct stat filestat;

      exprf = operatorArgv()[0];

      /* Open script file for reading */
      if( (fp = fopen(exprf, "r")) == NULL ) cdoAbort("Open failed on %s", exprf);

      if ( stat(exprf, &filestat) != 0 ) cdoAbort("Stat failed on %s", exprf);

      fsize = (size_t) filestat.st_size;
      exprs = (char*) malloc(fsize+1);

      while ( (ichar = fgetc(fp)) != EOF ) exprs[ipos++] = ichar;

      exprs[ipos] = 0;

      if ( ipos == 0 ) cdoAbort("%s is empty!", exprf);

      fclose(fp);
    }

  if ( cdoVerbose ) cdoPrint(exprs);


  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars1 = vlistNvars(vlistID1);

  if ( REPLACES_VARIABLES(operatorID) )
    {
      vlistID2 = vlistCreate();
      nvars = 0;
    }
  else
    {
      vlistID2 = vlistDuplicate(vlistID1);
      nvars = nvars1;
    }

  parse_arg.init = 1;
  parse_arg.vlistID1 = vlistID1;
  parse_arg.vlistID2 = vlistID2;
  parse_arg.nvars1   = 0;
  parse_arg.debug    = 0;
  if ( cdoVerbose ) parse_arg.debug    = 1;
  parse_arg.gridID2  = -1;
  parse_arg.zaxisID2 = -1;
  parse_arg.tsteptype2  = -1;
   
  /* Set all input variables to 'needed' if replacing is switched off */
  for ( varID = 0; varID < nvars1; varID++ )
    parse_arg.var_needed[varID] = ! REPLACES_VARIABLES(operatorID);

  yy_scan_string(exprs, scanner);
  yyparse(&parse_arg, scanner);

  parse_arg.init = 0;

  nvars2 = vlistNvars(vlistID2);
  if ( nvars2 == 0 ) cdoAbort("No output variable found!");

  if ( cdoVerbose ) vlistPrint(vlistID2);

  if ( cdoVerbose )
    for ( varID = 0; varID < nvars1; varID++ )
      if ( parse_arg.var_needed[varID] )
	printf("Needed var: %d %s\n", varID, parse_arg.var[varID]);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  parse_arg.vardata1 = (double**) malloc(nvars1*sizeof(double*));
  parse_arg.vardata2 = (double**) malloc(nvars2*sizeof(double*));

  for ( varID = 0; varID < nvars1; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      /* parse_arg.missval[varID] = vlistInqVarMissval(vlistID1, varID); */

      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      if ( parse_arg.var_needed[varID] )
	parse_arg.vardata1[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
      else
	parse_arg.vardata1[varID] = NULL;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      parse_arg.vardata2[varID] = parse_arg.vardata1[varID];
    }

  for ( varID = nvars; varID < nvars2; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID2, varID);
      zaxisID = vlistInqVarZaxis(vlistID2, varID);

      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      parse_arg.vardata2[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
    }

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( varID = 0; varID < nvars; varID++ ) parse_arg.nmiss[varID] = 0;

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( parse_arg.var_needed[varID] )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single1  = parse_arg.vardata1[varID] + offset;
	      streamReadRecord(streamID1, single1, &nmiss);
	      parse_arg.nmiss[varID] += nmiss;
	      /*
	      if ( nmiss && lwarn )
		{
		  cdoWarning("Missing values unsupported for this operator!");
		  lwarn = FALSE;
		}
	      */
	    }
	}

      for ( varID = nvars; varID < nvars2; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID2, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID2, varID);
	  gridsize = gridInqSize(gridID);
	  nlevel   = zaxisInqSize(zaxisID);

	  memset(parse_arg.vardata2[varID], 0, gridsize*nlevel*sizeof(double));
	}

      yy_scan_string(exprs, scanner);
      yyparse(&parse_arg, scanner);

      for ( varID = 0; varID < nvars2; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID2, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID2, varID);
	  missval  = vlistInqVarMissval(vlistID2, varID);

	  gridsize = gridInqSize(gridID);
	  nlevel   = zaxisInqSize(zaxisID);
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      long i;
	      offset   = gridsize*levelID;
	      single2  = parse_arg.vardata2[varID] + offset;

	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(single2[i], missval) ) nmiss++;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  yylex_destroy(scanner);

  if ( array ) free(array);
  if ( exprs ) free(exprs);

  cdoFinish();

  return (0);
}
