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

int readline(FILE *fp, char *line, int len)
{
  int ichar, ipos = 0;

  while ( (ichar = fgetc(fp)) != EOF )
    {
      if ( ichar == '\r' ) break;
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
	{
	  fprintf(stderr, "readline Warning: end of line not found (maxlen = %d)!\n", len);
	  break;
	}
    }
  line[ipos] = 0;

  if ( feof(fp) && ipos == 0 ) return (0);

  return (1);
}
