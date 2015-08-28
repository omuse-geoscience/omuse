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

/* key value list */

/*
  list type 1:
  ============

  &parameter
    name:              ps
    standard_name:     surface_air_pressure
    units:             Pa
    cell_methods:      "time: mean"
    cell_measures:     "area: areacella"
    long_name:         "Surface Air Pressure"
    comment:           "not, in general, the same as mean sea-level pressure"
    valid_min:         4.791e+04
    valid_max:         1.119e+05

  list type 2: one entry on each line
  ============

  variable_entry:    ps
    standard_name:     surface_air_pressure
    units:             Pa
    cell_methods:      time: mean
    cell_measures:     area: areacella
    long_name:         Surface Air Pressure
    comment:           not, in general, the same as mean sea-level pressure
    valid_min:         4.791e+04
    valid_max:         1.119e+05
*/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* strdup */
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>


#define MAX_KVLISTS    4096
#define MAX_KVELEMENTS 1024

typedef struct {
  char name[128];
  char *value;
} kvelement_t;

typedef struct {
  char name[128];
  int num_elements;
  kvelement_t elements[MAX_KVELEMENTS];
} kvlist_t;

typedef struct {
  char *filename;
  char *buffer;
  char *bufferp;
  size_t buffersize;
  int num_lists;
  kvlist_t lists[MAX_KVLISTS];
} kvl_t;


static
int kvlNewList(kvl_t *kvl)
{
  assert(kvl != NULL);

  kvl->num_lists++;
  assert(kvl->num_lists < MAX_KVLISTS);

  return (kvl->num_lists-1);
}

static
int kvlNewListElement(kvl_t *kvl, int listID)
{
  assert(kvl != NULL);

  kvl->lists[listID].num_elements++;
  assert(kvl->lists[listID].num_elements < MAX_KVELEMENTS);

  return (kvl->lists[listID].num_elements-1);
}

static
int kvlAddList(kvl_t *kvl, const char *name)
{
  int listID = kvlNewList(kvl);

  strcpy(kvl->lists[listID].name, name);

  return (listID);
}

static
void kvlAddListElement(kvl_t *kvl, int listID, const char *name, const char *value)
{
  int elemID = kvlNewListElement(kvl, listID);

  strcpy(kvl->lists[listID].elements[elemID].name, name);
  kvl->lists[listID].elements[elemID].value = strdup(value);
}

static
char *readLineFromBuffer(char *buffer, size_t *buffersize, char *line, size_t len)
{
  int ichar;
  size_t ipos = 0;

  while ( *buffersize )
    {
      ichar = *buffer;
      (*buffersize)--;
      buffer++;
      if ( ichar == '\r' ) break;
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
        {
          fprintf(stderr, "readLineFromBuffer: end of line not found (maxlen = %ld)!\n", len);
          break;
        }
    }
  line[ipos] = 0;

  if ( *buffersize == 0 && ipos == 0 ) buffer = NULL;

  return (buffer);
}

static
void pfree(void *ptr)
{
  if ( ptr ) free(ptr);
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return (pline);
}

static
char *getElementName(char *pline, char *name)
{
  int pos = 0, len;

  while ( isspace((int) *pline) ) pline++;
  len = strlen(pline);
  while ( pos < len && !isspace((int) *(pline+pos)) && *(pline+pos) != '=' && *(pline+pos) != ':' ) pos++;

  strncpy(name, pline, pos);
  name[pos] = 0;

  pline += pos;
  return (pline);
}

static
char *getElementValue(char *pline)
{
  int len;

  while ( isspace((int) *pline) ) pline++;
  len = strlen(pline);
  while ( isspace((int) *(pline+len-1)) && len ) { *(pline+len-1) = 0; len--;}

  return (pline);
}

static
void kvlParseBuffer(kvl_t *kvl)
{
  char line[4096];
  char name[256];
  char *pline;
  char *buffer = kvl->buffer;
  size_t buffersize = kvl->buffersize;
  int linenumber = 0;
  char listkey1[] = "axis_entry:";
  char listkey2[] = "variable_entry:";
  int listtype = 0;
  int listID = -1;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '#' || *pline == '!' || *pline == '\0' ) continue;
      //  len = (int) strlen(pline);
      if ( listtype == 0 && *pline == '&' )
	{
	  listtype = 1;
	}
      
      if ( strncmp(pline, listkey1, strlen(listkey1)) == 0 )
	{
	  pline += strlen(listkey1);

	  listtype = 2;

	  listID = kvlAddList(kvl, "axis");

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlAddListElement(kvl, listID, "name", pline);
	}
      else if ( strncmp(pline, listkey2, strlen(listkey2)) == 0 )
	{
	  pline += strlen(listkey2);

	  listtype = 2;

	  listID = kvlAddList(kvl, "variable");

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlAddListElement(kvl, listID, "name", pline);
	}
      else
	{
	  pline = getElementName(pline, name);
	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( listID == -1 ) listID = kvlAddList(kvl, "global");

	  if ( *pline ) kvlAddListElement(kvl, listID, name, pline);

	    {
	      //fprintf(stderr, "%d skip line %3d: %s\n", newlist, linenumber, pline);
	    }
	}

      //   printf("%s\n", pline);
    }
}


void *kvlParseFile(const char *filename)
{
  kvl_t *kvl = NULL;
  FILE *fp;
  char *buffer;
  size_t filesize;
  size_t nitems;

  assert(filename != NULL);

  filesize = fileSize(filename);

  fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return (kvl);
    }

  buffer = (char*) malloc(filesize);
  nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return (kvl);
    }
 
  kvl = (kvl_t*) calloc(1, sizeof(kvl_t));
  kvl->buffer = buffer;
  kvl->buffersize = filesize;
  kvl->filename = strdup(filename);

  kvlParseBuffer(kvl);
  
  return ((void *) kvl);
}


void kvlDelete(void *kvlist)
{
  kvl_t *kvl = (kvl_t *) kvlist;

  assert(kvl != NULL);

  pfree(kvl->filename);
  pfree(kvl->buffer);

  free(kvl);
}


int kvlGetNumLists(void *kvlist)
{
  kvl_t *kvl = (kvl_t *) kvlist;

  assert(kvl != NULL);

  return(kvl->num_lists);
}


const char *kvlGetListName(void *kvlist, int listID)
{
  kvl_t *kvl = (kvl_t *) kvlist;
  char *listname = NULL;
  
  assert(listID < kvl->num_lists);
  
  listname = kvl->lists[listID].name;

  return (listname);
}

int kvlGetListNumElements(void *kvlist, int listID)
{
  kvl_t *kvl = (kvl_t *) kvlist;
  int nelements = 0;

  assert(listID < kvl->num_lists);

  nelements = kvl->lists[listID].num_elements;

  return (nelements);
}


const char *kvlGetListElementName(void *kvlist, int listID, int elemID)
{
  kvl_t *kvl = (kvl_t *) kvlist;
  char *ename = NULL;

  assert(listID < kvl->num_lists);
  assert(elemID < kvl->lists[listID].num_elements);

  ename = kvl->lists[listID].elements[elemID].name;

  return (ename);
}


const char *kvlGetListElementValue(void *kvlist, int listID, int elemID)
{
  kvl_t *kvl = (kvl_t *) kvlist;
  char *evalue = NULL;

  assert(listID < kvl->num_lists);
  assert(elemID < kvl->lists[listID].num_elements);

  evalue = kvl->lists[listID].elements[elemID].value;

  return (evalue);
}

/*
int main(int argc, char *argv[])
{
  char *filename;
  void *kvlist;
  int nlists, listID;
  int nelements, elemID;
  const char *listname;
  const char *ename;
  const char *evalue;

  if ( argc != 2 ) 
    {
      fprintf(stderr, "usage: kvlist filename\n");
      return (1);
    }

  filename = argv[1];

  printf("Parse file: %s\n", filename);

  kvlist = kvlParseFile(filename);
  nlists = kvlGetNumLists(kvlist);
  printf("# Number of lists: %d\n", nlists);
  for ( listID = 0; listID < nlists; ++listID )
    {
      listname = kvlGetListName(kvlist, listID);
      nelements = kvlGetListNumElements(kvlist, listID);
      printf("# list ID: %d;   Number of elements: %d\n", listID, nelements);
      printf("&%s\n", listname);
      for ( elemID = 0; elemID < nelements; ++elemID )
	{
	  ename  = kvlGetListElementName(kvlist, listID, elemID);
	  evalue = kvlGetListElementValue(kvlist, listID, elemID);
	  printf("  %s = %s\n", ename, evalue);
	}
      printf("/\n");
    }

  kvlDelete(kvlist);

  return (0);
}
*/
