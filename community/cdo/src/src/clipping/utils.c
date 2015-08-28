/**
 * @file utils.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://redmine.dkrz.de/doc/YAC/html/index.html
 *
 * This file is part of YAC.
 *
 * YAC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with YAC.  If not, see <http://www.gnu.org/licenses/gpl.txt>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static void ** pointer_lookup_table = NULL;
static unsigned pointer_lookup_table_size = 0;

unsigned yac_pointer_to_unique_id(void * pointer) {

   pointer_lookup_table = realloc (pointer_lookup_table,
      ++pointer_lookup_table_size * sizeof(pointer_lookup_table[0]));

   pointer_lookup_table[pointer_lookup_table_size-1] = pointer;

   return pointer_lookup_table_size - 1;
}

void * yac_unique_id_to_pointer(unsigned id) {

   if (id < pointer_lookup_table_size)
      return pointer_lookup_table[id];
   else
      return NULL;
}

void yac_free_pointer_unique_lookup() {

  free(pointer_lookup_table);
  pointer_lookup_table = NULL;
  pointer_lookup_table_size = 0;
}

void yac_internal_abort_message ( char * text, char * file, int line )
{
  fprintf(stderr, "%s \n", text); 
  fprintf(stderr, "Aborting in file %s, line %i ...\n", file, line );
  exit(EXIT_FAILURE);
}

void yac_abort_message ( char * text, char * file, int line ) {
	yac_internal_abort_message ( text, file, line );
}

/* ------------------------------------

   Source: http://www.cse.yorku.ca/~oz/hash.html 

unsigned long hash(const char *str) {

  unsigned long hash = 5381;
  int c;

  while ((c = *str++))
    hash = ((hash << 5) + hash) + c; 

  return hash;
}

  ------------------------------------*/

/*  Source http://snipplr.com/view/9022/string-hash-table/ */

#define NHASH 29989 //Use a prime number!
#define MULT 31

unsigned int yac_hash(const char *str) {
  unsigned int h = 0;
  for(; *str; str++)
    h = MULT * h + *str;
  return h % NHASH;
}
