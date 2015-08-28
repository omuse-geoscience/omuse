/**
 * @file clipping.h
 * @brief Structs and interfaces for cell clipping
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

#ifndef CLIPPING_H
#define CLIPPING_H

#include "points.h"

/** \example test_clipping.c
 * This contains some examples on how to use the cell_clipping routine.
 */


/**
  * \brief cell clipping to get the cells describing the intersections
  *
  * The routine takes (a list of) source cells and a target cell. It sets the
  * target cell data and does some further initialisation. Thus it needs to be
  * called for each new target cell intersection calculation
  *
  * The vertices of source and target cells can be either provided in a clockwise
  * or anticlockwise sense. However, the same sense must be used for source and
  * target cells.
  *
  * @param[in] N              number of source cells
  * @param[in] source_cell    list of source cells
  * @param[in] target_cell    target cell
  * @param[in] overlap_buffer buffer for the overlaps between the target and
  *                           the source cells
  *
  * \remark source and target cells have to be convex
  * \remark cells in overlap_buffer can be concave
  * \remark overlap_buffer must contain valid grid_cells (have to be initialised
  *         using \ref yac_init_grid_cell; initialisation have to be done only once,
  *         in consecutive calls, the cells can be reused with have to be
  *         reinitialised)
  *
 **/
void yac_cell_clipping (unsigned N,
                        struct grid_cell * source_cell,
                        struct grid_cell target_cell,
                        struct grid_cell * overlap_buffer);

/** \example test_partial_areas.c
 * This contains examples on how to use \ref yac_compute_overlap_areas.
 */

/** \example test_compute_overlap_area.c
 * This contains examples on how to use \ref yac_compute_overlap_areas.
 */

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with triangular target cells. This is required for
  *        conservative remapping
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a convex target cell. As
  * a triangle is always convex we recommend to use this routine only for
  * triangular target cells. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N             number of source cells
  * @param[in]  source_cell   list of source cells
  * @param[in]  target_cell   target cell
  * @param[out] partial_areas list of N partial weights, one weight for each
  *                           source-target intersection
  *
  * \remark source and target cell have to be convex
  *
 **/
void yac_compute_overlap_areas (unsigned N,
                                struct grid_cell * source_cell,
                                struct grid_cell target_cell,
                                double * partial_areas);

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with arbitrary target cells, this is required for conservative
  *        remapping.
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a target cell. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N             number of source cells
  * @param[in]  source_cell   list of source cells
  * @param[in]  target_cell   target cell
  * @param[in]  target_node_x x-coordinate of target cell node or center point
  * @param[in]  target_node_y y-coordinate of target cell node or center point
  * @param[out] partial_areas list of N partial weights, one weight for each
  *                           source-target intersection
  *
  * \remark source and target cell have to be convex
  *
 **/
void yac_compute_concave_overlap_areas (unsigned N,
                                        struct grid_cell * source_cell,
                                        struct grid_cell target_cell,
                                        double * target_node_x,
                                        double * target_node_y,
                                        double * partial_areas);
/**
  * \brief correct interpolation weights
  *
  * Returns weights with a sum close to 1.0
  *
  * @param[in]  N                 number of source cells
  * @param[out] weight            list of N partial weights
  *
 **/
void yac_correct_weights (unsigned N, double * weight);

#endif // CLIPPING_H
