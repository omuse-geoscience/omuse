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

#ifndef _CDO_GETOPT_H
#define _CDO_GETOPT_H

struct cdo_option {
  /* name of long option */
  const char *name;
  /*
   * one of no_argument, required_argument, and optional_argument:
   * whether option takes an argument
   */
  int has_arg;
  /* if not NULL, set *flag to val when option found */
  int *flag;
  /* if flag not NULL, value to set *flag to; else return value */
  int val;
};


#define  no_argument        1   // no argument to the option is expect
#define  required_argument  2   // an argument to the option is required
#define  optional_argument  3   // an argument to the option may be presented.

int cdo_getopt(int argc, char * const *argv, const char *optstring);
int cdo_getopt_long(int argc, char * const *argv, const char *optstring, const struct cdo_option *longopts, int *longindex);

#endif  /* _CDO_GETOPT_H */
