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

#include <string.h>

static int    gargc = 0;
static char **gargv;

static char CommandLine[1024];

void initCommandLine(void)
{
  int iarg;
  char *pargv;
  size_t len, offset = 0;
  /*
  time_t tp;
  struct tm *ltime;

  tp = time(NULL);

  if ( tp != -1 )
    {
      ltime = localtime(&tp);
      offset = strftime(CommandLine, 1024, "%d %b %Y : ", ltime);
    }
    */
  for ( iarg = 0; iarg < gargc; iarg++ )
    {
      if ( iarg == 0 )
	{
	  pargv = strrchr(gargv[iarg], '/');
	  if ( pargv == 0 ) pargv = gargv[0];
	  else              pargv++;
	}
      else
	pargv = gargv[iarg];

      len = strlen(pargv);
      if ( offset+len+1 > 1024 ) break;
      memcpy(CommandLine+offset, pargv, len);
      offset += len;
      CommandLine[offset] = ' ';
      offset++;
    }

  CommandLine[offset-1] = '\0';
}

char *commandLine(void)
{
  static int init = 0;

  if ( init == 0 )
    {
      initCommandLine();
      init = 1;
    }

  return (CommandLine);
}

void setCommandLine(int argc, char **argv)
{
  gargc = argc;
  gargv = argv;
}
