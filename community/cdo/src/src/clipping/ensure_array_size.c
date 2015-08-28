/**
 * @file ensure_array_size.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include <stdio.h>
#include <stdlib.h>

#include "ensure_array_size.h"

void
yac_realloc_array(void **array, size_t elem_size, size_t *curr_array_size,
                  size_t requested_size)
{
  const size_t array_inc_size = (1024 + elem_size - 1)/ elem_size;
  *curr_array_size = array_inc_size
    * ((requested_size + array_inc_size) / array_inc_size);
  *array = realloc(*array, *curr_array_size * elem_size);
  if (*array == NULL && (*curr_array_size * elem_size)) {
     fprintf(stderr, "error in realloc\n");
     exit(EXIT_FAILURE);
  }
}
