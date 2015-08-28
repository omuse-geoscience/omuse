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
#include <string.h>
#include <assert.h>
#include "cdo_getopt.h"
#include "util.h" // refactor: necessary for accessing global vars

static int CDO_optopt = 0;
static int CDO_optreset = 0; // TODO refactor: var without semantic effect

int cdo_getopt(int argc, char * const *argv, const char *optstring)
{
  static int optpos = 0;
  int optval = -1, value;
  int opthasarg = 0;
  int optstrlen = strlen(optstring);
  int iargc;

  assert(argv != NULL);
  assert(optstring != NULL);

  CDO_optarg = NULL;

  while ( optpos < optstrlen && CDO_optind < argc )
    {
      value = optstring[optpos];
      optpos++;
      if ( optstring[optpos] == ':' )
	{
	  opthasarg = 1;
	  optpos++;
	}
      else
	opthasarg = 0;

      for ( iargc = 1; iargc < argc; iargc++ )
	{
	  if ( *argv[iargc] == '-' && argv[iargc][2] == 0 )
	    {
	      if ( (argv[iargc][1]) == value )
		{
		  optval = value;
		  CDO_optind++;
		  if ( opthasarg )
		    {
		      CDO_optarg = argv[iargc+1];
		      CDO_optind++;
		    }
		  break;
		}
	    }
	}

      if ( iargc < argc ) break;
    }

  if ( opthasarg && CDO_optarg == NULL ) optval = ':';

  return (optval);
}

/* http://www.opensource.apple.com/source/Kerberos/Kerberos-47/KerberosFramework/Kerberos5/Sources/util/windows/getopt_long.c */

static
char *__progname(char *nargv0)
{
  char * tmp;

  assert(nargv0 != NULL);

  tmp = strrchr(nargv0, '/');
  if (tmp)
    tmp++;
  else
    tmp = nargv0;
  
  return(tmp);
}

#define	BADCH	(int)'?'
#define	BADARG	(int)':'
#define	EMSG	""

static
int cdo_getopt_internal(int nargc, char * const *nargv, const char *ostr)
{
  static char *place = EMSG;		/* option letter processing */
  char *oli;				/* option letter list index */

  assert(nargv != NULL);
  assert(ostr != NULL);

  if ( CDO_optreset || !*place )
    {		/* update scanning pointer */
      CDO_optreset = 0;
      if ( CDO_optind >= nargc || *(place = nargv[CDO_optind]) != '-' )
	{
	  place = EMSG;
	  return (-1);
	}
      
      if ( place[1] && *++place == '-' )
	{	/* found "--" */
	  /* ++CDO_optind; */
	  place = EMSG;
	  return (-2);
	}

      if ( place[1] ) return (-1);
    } 

  /* option letter okay? */
  if ( (CDO_optopt = (int)*place++) == (int)':' || !(oli = strchr(ostr, CDO_optopt)) )
    {
      /* if the user didn't specify '-' as an option, assume it means -1. */
      if ( CDO_optopt == (int)'-' ) return (-1);

      if ( !*place ) ++CDO_optind;

      if ( CDO_opterr && *ostr != ':' )
	(void)fprintf(stderr, "%s: illegal option -- %c\n", __progname(nargv[0]), CDO_optopt);

      return (BADCH);
    }

  if ( *++oli != ':' )
    {		        	/* don't need argument */
      CDO_optarg = NULL;
      if (!*place)
	++CDO_optind;
    }
  else
    {				/* need an argument */
      if ( *place )		/* no white space */
	CDO_optarg = place;
      else if (nargc <= ++CDO_optind)
	{               	/* no arg */
	  place = EMSG;
	  if ( (CDO_opterr) && (*ostr != ':') )
	    fprintf(stderr, "%s: option requires an argument -- %c\n", __progname(nargv[0]), CDO_optopt);
	  return (BADARG);
	}
      else  /* white space */
	CDO_optarg = nargv[CDO_optind];

      place = EMSG;
      ++CDO_optind;
    }
 
  return (CDO_optopt);		/* dump back option letter */
}

int cdo_getopt_long(int argc, char * const *argv, const char *optstring, const struct cdo_option *longopts, int *longindex)
{
  int retval;

  assert(argv != NULL);
  assert(optstring != NULL);
  assert(longopts != NULL);
  /* index may be NULL */

  if ( (retval = cdo_getopt_internal(argc, argv, optstring)) == -2 )
    {
      char *current_argv = argv[CDO_optind++] + 2, *has_equal;
      int i, current_argv_len, match = -1;

      if ( *current_argv == '\0' ) return(-1);

      if ( (has_equal = strchr(current_argv, '=')) != NULL )
	{
	  current_argv_len = has_equal - current_argv;
	  has_equal++;
	}
      else
	current_argv_len = strlen(current_argv);

      for ( i = 0; longopts[i].name; i++ )
	{ 
	  if ( strncmp(current_argv, longopts[i].name, current_argv_len) )
	    continue;

	  if ( strlen(longopts[i].name) == (unsigned)current_argv_len )
	    { 
	      match = i;
	      break;
	    }

	  if ( match == -1 ) match = i;
	}
	
      if ( match != -1 )
	{
	  if ( longopts[match].has_arg == required_argument ||
	       longopts[match].has_arg == optional_argument )
	    {
	      if ( has_equal )
		CDO_optarg = has_equal;
	      else
		CDO_optarg = argv[CDO_optind++];
	    }
	  
	  if ( (longopts[match].has_arg == required_argument) && (CDO_optarg == NULL) )
	    {
	      /* Missing argument, leading: indicates no error should be generated */
	      if ( (CDO_opterr) && (*optstring != ':') )
		fprintf(stderr, "%s: option requires an argument -- %s\n", __progname(argv[0]), current_argv);
	      return (BADARG);
	    }
	}
      else
	{ /* No matching argument */
	  if ( (CDO_opterr) && (*optstring != ':') )
	    fprintf(stderr, "%s: illegal option -- %s\n", __progname(argv[0]), current_argv);

	  return (BADCH);
	}
	
      if ( longopts[match].flag )
	{
	  *longopts[match].flag = longopts[match].val;
	  retval = 0;
	}
      else 
	retval = longopts[match].val;
  
      if ( longindex ) *longindex = match;
    }
  
  return (retval);
}
