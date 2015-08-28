/**
 * @file points.h
 * @brief Structs and interfaces for points
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

#ifndef POINTS_H
#define POINTS_H

/** \example test_points.c
 * This shows how to work with struct points.
 */

#include "grid.h"

enum yac_location {

   CELL =   0,
   CORNER = 1,
   EDGE =   2,
};

struct points {

   enum yac_location location;

   double * coordinates_x, * coordinates_y;

   struct grid * base_grid;

   struct grid * point_grid;
};

enum yac_location yac_get_location(int const location);

/**
 * initialises a struct points
 * @param[in,out] points        points to be initialised
 * @param[in]     base_grid     grid on which the points situated
 * @param[in]     location      location of the points
 * @param[in]     coordinates_x x coordinates of the points
 * @param[in]     coordinates_y y coordinates of the points
 *
 * \remarks - for CORNER and CELL points the coordinates need to be provided in the same
 *            way that is used to initialise the base_grid
 * \remarks - for EDGE points the coordinates need to have an entry for each edge (in
 *            the order of the local ids)
 */
void yac_init_points(struct points * points, struct grid * base_grid, enum yac_location location,
                     double * coordinates_x, double * coordinates_y);

/**
 * returns a grid that has the points defined by the struct points as corners
 * @param[in] points points for which the point grid is to be generated
 * @return a point grid
 *
 * \remarks does not work for EDGE
 */
struct grid * yac_get_point_grid(struct points * points);
struct grid * yac_get_base_grid(struct points * points);

void yac_get_coordinate_array_sizes (struct grid * grid, enum yac_location location, unsigned * sizes);

/**
 * returns the number of points in the struct points
 * @param[in] points
 * @return number of points
 */
unsigned yac_get_data_size(struct points points);

void yac_get_point_coordinates (struct points * points, unsigned local_point_id,
                                double * coordinates);

void yac_free_points(struct points * points);

#endif // POINTS_H
