/**
 * @file grid_cell.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 *             Thomas Jahns <jahns@dkrz.de>
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
#include <string.h>
#include <stdio.h>

#include "grid_cell.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "geometry.h"

void yac_init_grid_cell(struct grid_cell * cell) {

   cell->coordinates_x = NULL;
   cell->coordinates_y = NULL;
   cell->coordinates_xyz = NULL;
   cell->edge_type = NULL;
   cell->num_corners = 0;
   cell->array_size = 0;
}

void yac_copy_grid_cell(struct grid_cell in_cell, struct grid_cell * out_cell) {

   if (in_cell.num_corners > out_cell->array_size) {

      free(out_cell->coordinates_x);
      free(out_cell->coordinates_y);
      free(out_cell->coordinates_xyz);
      free(out_cell->edge_type);
      out_cell->coordinates_x = malloc(in_cell.num_corners *
                                       sizeof(*(out_cell->coordinates_x)));
      out_cell->coordinates_y = malloc(in_cell.num_corners *
                                       sizeof(*(out_cell->coordinates_y)));
      out_cell->coordinates_xyz = malloc(3 * in_cell.num_corners *
                                         sizeof(*(out_cell->coordinates_xyz)));
      out_cell->edge_type = malloc(in_cell.num_corners *
                                   sizeof(*(out_cell->edge_type)));
      out_cell->array_size = in_cell.num_corners;
   }

   memcpy(out_cell->coordinates_x, in_cell.coordinates_x,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_x)));
   memcpy(out_cell->coordinates_y, in_cell.coordinates_y,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_y)));
   memcpy(out_cell->coordinates_xyz, in_cell.coordinates_xyz,
          3 * in_cell.num_corners * sizeof(*(out_cell->coordinates_xyz)));
   memcpy(out_cell->edge_type, in_cell.edge_type,
          in_cell.num_corners * sizeof(*(out_cell->edge_type)));
   out_cell->num_corners = in_cell.num_corners;
}

void yac_free_grid_cell(struct grid_cell * cell) {

   if (cell->coordinates_x != NULL) free(cell->coordinates_x);
   if (cell->coordinates_y != NULL) free(cell->coordinates_y);
   if (cell->coordinates_xyz != NULL) free(cell->coordinates_xyz);
   if (cell->edge_type != NULL) free(cell->edge_type);

   yac_init_grid_cell(cell);
}
