/**
 * @file clipping.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "geometry.h"
#include "clipping.h"
#include "area.h"
#include "ensure_array_size.h"
#include "utils.h"

//#define VERBOSE

static double const tol = 1.0e-12;

enum cell_type {
  LON_LAT_CELL,
  LAT_CELL,
  GREAT_CIRCLE_CELL,
  MIXED_CELL
};

struct point_list_element {

  double vec_coords[3];
  enum yac_edge_type edge_type; // type of edge with next corner
  int to_be_removed;
  struct point_list_element * next;
};

struct point_list {

  struct point_list_element * first;
  struct point_list_element * last;
  struct point_list_element * free_elements;
};

/* internal helper routines for working with linked lists of points */

static void init_point_list(struct point_list * list);

static void reset_point_list(struct point_list * list);

static void generate_point_list(struct point_list * list, struct grid_cell cell);

static struct point_list_element *
get_free_point_list_element(struct point_list * list);

static void remove_points(struct point_list * list);
static unsigned remove_zero_length_edges(struct point_list * list);

static void free_point_list(struct point_list * list);

static unsigned get_cell_points_ordering(struct point_list * cell);

static void generate_overlap_cell(struct point_list * list,
                                  struct grid_cell * cell);

static enum cell_type get_cell_type(struct grid_cell target_cell);

/* ------------------------- */

void yac_compute_overlap_areas (unsigned N,
                                struct grid_cell * source_cell,
                                struct grid_cell target_cell,
                                double * partial_areas) {

  static struct grid_cell * overlap_buffer = NULL;
  static unsigned overlap_buffer_size = 0;

  // ensure that there are enough buffer cells

  if (overlap_buffer_size < N) {

    unsigned old_overlap_buffer_size = overlap_buffer_size;

    ENSURE_ARRAY_SIZE(overlap_buffer, overlap_buffer_size, N);

    for (; old_overlap_buffer_size < overlap_buffer_size;
         ++old_overlap_buffer_size)
      yac_init_grid_cell(overlap_buffer + old_overlap_buffer_size);
  }

  /* Do the clipping and get the cell for the overlapping area */

  yac_cell_clipping ( N, source_cell, target_cell, overlap_buffer);

  /* Get the partial areas for the overlapping regions */

  for (unsigned n = 0; n < N; n++) {
    partial_areas[n] = yac_huiliers_area (overlap_buffer[n]);
    // we cannot use pole_area because it is rather inaccurate for great circle
    // edges that are nearly circles of longitude
    //partial_areas[n] = pole_area (overlap_buffer[n]);
  }

#ifdef VERBOSE
  for (unsigned n = 0; n < N; n++)
    printf("overlap area : %lf\n", partial_areas[n]);
#endif
}

/* ------------------------- */

void yac_compute_concave_overlap_areas (unsigned N,
                                        struct grid_cell * source_cell,
                                        struct grid_cell target_cell,
                                        double * target_node_x,
                                        double * target_node_y,
                                        double * partial_areas) {
  enum cell_type target_cell_type;

  if ( target_cell.num_corners > 3 )
    target_cell_type = get_cell_type (target_cell);

  if ( target_cell.num_corners < 4 || target_cell_type == LON_LAT_CELL ) {
    yac_compute_overlap_areas (N, source_cell, target_cell, partial_areas);
    return;
  }

  if ( target_node_x == NULL || target_node_y == NULL )
    yac_internal_abort_message("ERROR: missing target point coordinates "
                               "(x_coordinates == NULL || y_coordinates == NULL)",
                               __FILE__ , __LINE__);

  struct grid_cell target_partial_cell =
    {.coordinates_x   = (double[3]){-1},
     .coordinates_y   = (double[3]){-1},
     .coordinates_xyz = (double[3*3]){-1},
     .edge_type       = (enum yac_edge_type[3]) {GREAT_CIRCLE},
     .num_corners     = 3};

  static struct grid_cell * overlap_buffer = NULL;
  static unsigned overlap_buffer_size = 0;

  // ensure that there are enough buffer cells

  if (overlap_buffer_size < N) {

    unsigned old_overlap_buffer_size = overlap_buffer_size;

    ENSURE_ARRAY_SIZE(overlap_buffer, overlap_buffer_size, N);

    for (; old_overlap_buffer_size < overlap_buffer_size;
         ++old_overlap_buffer_size)
      yac_init_grid_cell(overlap_buffer + old_overlap_buffer_size);
  }

  /* Do the clipping and get the cell for the overlapping area */

  for ( unsigned n = 0; n < N; n++) partial_areas[n] = 0.0;

  // common node point to all partial target cells
  target_partial_cell.coordinates_x[0] = *target_node_x;
  target_partial_cell.coordinates_y[0] = *target_node_y;

  LLtoXYZ ( *target_node_x, *target_node_y, target_partial_cell.coordinates_xyz );

  for ( unsigned num_corners = 0; num_corners < target_cell.num_corners; ++num_corners ) {

    unsigned corner_a = num_corners;
    unsigned corner_b = (num_corners+1)%target_cell.num_corners;

    // skip clipping and area calculation for degenerated triangles
    //
    // If this is not sufficient, instead we can try something like:
    //
    //     struct point_list target_list
    //     init_point_list(&target_list);
    //     generate_point_list(&target_list, target_cell);
    //     struct grid_cell temp_target_cell;
    //     generate_overlap_cell(target_list, temp_target_cell);
    //     free_point_list(&target_list);
    //
    // and use temp_target_cell for triangulation.
    //
    // Compared to the if statement below the alternative seems
    // to be quite costly.

    if ( ( ( fabs(target_cell.coordinates_xyz[0+3*corner_a]-target_cell.coordinates_xyz[0+3*corner_b]) < tol ) &&
           ( fabs(target_cell.coordinates_xyz[1+3*corner_a]-target_cell.coordinates_xyz[1+3*corner_b]) < tol ) &&
           ( fabs(target_cell.coordinates_xyz[2+3*corner_a]-target_cell.coordinates_xyz[2+3*corner_b]) < tol ) ) ||
    	 ( ( fabs(target_cell.coordinates_xyz[0+3*corner_a]-target_partial_cell.coordinates_xyz[0]) < tol    ) &&
           ( fabs(target_cell.coordinates_xyz[1+3*corner_a]-target_partial_cell.coordinates_xyz[1]) < tol    ) &&
           ( fabs(target_cell.coordinates_xyz[2+3*corner_a]-target_partial_cell.coordinates_xyz[2]) < tol    ) ) ||
         ( ( fabs(target_cell.coordinates_xyz[0+3*corner_b]-target_partial_cell.coordinates_xyz[0]) < tol    ) &&
           ( fabs(target_cell.coordinates_xyz[1+3*corner_b]-target_partial_cell.coordinates_xyz[1]) < tol    ) &&
           ( fabs(target_cell.coordinates_xyz[2+3*corner_b]-target_partial_cell.coordinates_xyz[2]) < tol    ) ) )
       continue;

    target_partial_cell.coordinates_x[1] = target_cell.coordinates_x[corner_a];
    target_partial_cell.coordinates_y[1] = target_cell.coordinates_y[corner_a];
    target_partial_cell.coordinates_x[2] = target_cell.coordinates_x[corner_b];
    target_partial_cell.coordinates_y[2] = target_cell.coordinates_y[corner_b];

    target_partial_cell.coordinates_xyz[0+3*1] = target_cell.coordinates_xyz[0+3*corner_a];
    target_partial_cell.coordinates_xyz[1+3*1] = target_cell.coordinates_xyz[1+3*corner_a];
    target_partial_cell.coordinates_xyz[2+3*1] = target_cell.coordinates_xyz[2+3*corner_a];
    target_partial_cell.coordinates_xyz[0+3*2] = target_cell.coordinates_xyz[0+3*corner_b];
    target_partial_cell.coordinates_xyz[1+3*2] = target_cell.coordinates_xyz[1+3*corner_b];
    target_partial_cell.coordinates_xyz[2+3*2] = target_cell.coordinates_xyz[2+3*corner_b];

    yac_cell_clipping ( N, source_cell, target_partial_cell, overlap_buffer);

    /* Get the partial areas for the overlapping regions as sum over the partial target cells. */

    for (unsigned n = 0; n < N; n++) {
      partial_areas[n] += yac_huiliers_area (overlap_buffer[n]);
      // we cannot use pole_area because it is rather inaccurate for great circle
      // edges that are nearly circles of longitude
      //partial_areas[n] = pole_area (overlap_buffer[n]);
    }
  }

#ifdef VERBOSE
  for (unsigned n = 0; n < N; n++)
	    printf("overlap area %i: %lf \n", n, partial_areas[n]);
#endif
}

/* ------------------------- */

static double dotproduct(double a[], double b[]) {

  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void compute_norm_vector(double a[], double b[], double norm[]) {

  crossproduct_ld(a, b, norm);

  if ((fabs(norm[0]) < tol) &&
      (fabs(norm[1]) < tol) &&
      (fabs(norm[2]) < tol))
    yac_internal_abort_message("ERROR: a and b are identical -> no norm vector\n",
                               __FILE__, __LINE__);

  double scale = 1.0 / sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);

  norm[0] *= scale;
  norm[1] *= scale;
  norm[2] *= scale;
}

static void compute_lat_circle_z_value(double a[], double b[], double z[]) {

  double temp[3];

  crossproduct_ld(a, b, temp);

  z[0] = 0;
  z[1] = 0;

  if (temp[2] > 0)
    z[2] = 1.0 - a[2];
  else
    z[2] = -1.0 - a[2];
}

/**
 * Determines whether a given point is on a hemisphere that is defined by a plane
 * through the middle of the sphere.\n
 * The plane is defined by its norm vector.
 * @param[in] point point to be checked
 * @param[in] norm_vec norm vector of the plane dividing the sphere
 * @returns  0 if the point is not inside the hemisphere
 *           1 if the point is inside the hemisphere
 *           2 if the point is in the plane
 */
static unsigned is_inside_gc(double point[], double norm_vec[]) {

  double dot;

  // the product is defined as follows
  // a * b = |a| * |b| * cos(alpha)
  // where alpha is the angle between a and b

  dot = dotproduct(point, norm_vec);

  // if the point is on the line
  if (fabs(dot) <= yac_angle_tol * 1e-3)
    return 2;

  return dot < 0;
}

static unsigned is_inside_latc(double point[], double z) {

  double temp = fabs(point[2] + z);

  if (fabs(1.0 - temp) <= yac_angle_tol * 1e-3) return 2;
  else return temp < 1.0;
}

static unsigned is_inside(double point[], double help_vec[],
                          enum yac_edge_type edge_type,
                          unsigned cell_points_ordering) {

  unsigned ret_val = 0;

  switch (edge_type) {

    case (LON_CIRCLE) :
    case (GREAT_CIRCLE) :
      ret_val = is_inside_gc(point, help_vec);
      break;
    case (LAT_CIRCLE) :
      ret_val = is_inside_latc(point, help_vec[2]);
      break;
    default:
      yac_internal_abort_message("invalid edge type\n", __FILE__, __LINE__);
  };

  if (ret_val == 2) return 2;
  else return ret_val ^ cell_points_ordering;
}

static enum cell_type get_cell_type(struct grid_cell target_cell) {

  int count_lat_edges = 0, count_great_circle_edges = 0;

   if ((target_cell.num_corners == 4) &&
       ((target_cell.edge_type[0] == LAT_CIRCLE &&
         target_cell.edge_type[1] == LON_CIRCLE &&
         target_cell.edge_type[2] == LAT_CIRCLE &&
         target_cell.edge_type[3] == LON_CIRCLE) ||
        (target_cell.edge_type[0] == LON_CIRCLE &&
         target_cell.edge_type[1] == LAT_CIRCLE &&
         target_cell.edge_type[2] == LON_CIRCLE &&
         target_cell.edge_type[3] == LAT_CIRCLE)))
      return LON_LAT_CELL;
   else
      for (unsigned i = 0; i < target_cell.num_corners; ++i)
         if (target_cell.edge_type[i] == LON_CIRCLE ||
             target_cell.edge_type[i] == GREAT_CIRCLE)
            count_great_circle_edges++;
         else
            count_lat_edges++;

   if (count_lat_edges && count_great_circle_edges)
      return MIXED_CELL;
   else if (count_lat_edges)
      return LAT_CELL;
   else
      return GREAT_CIRCLE_CELL;
}

static void get_edge_middle_point_lat(double a[3], double b[3],
                                      double middle[3]) {

  middle[0] = a[0] + b[0];
  middle[1] = a[1] + b[1];
  middle[2] = a[2];

  double length = sqrt(middle[0] * middle[0] + middle[1] * middle[1]);

  if (length < tol) {

    length = 1;
    middle[0] = 1;
    middle[1] = 0;
  }

  double scale = sqrt(1.0 - middle[2] * middle[2]) / length;

  middle[0] *= scale;
  middle[1] *= scale;
}

static void get_edge_middle_point_gc(double a[3], double b[3],
                                     double middle[3]) {

  middle[0] = a[0] + b[0];
  middle[1] = a[1] + b[1];
  middle[2] = a[2] + b[2];

  double length = sqrt(middle[0] * middle[0] +
                       middle[1] * middle[1] +
                       middle[2] * middle[2]);

  if (length < tol) {

    middle[0] = 1;
    middle[1] = 0;
    middle[2] = 0;
    return;
  }

  double scale = 1.0 / length;

  middle[0] *= scale;
  middle[1] *= scale;
  middle[2] *= scale;
}

static void get_edge_middle_point (double a[3], double b[3],
                                   enum yac_edge_type edge_type, double middle[3]) {

  switch (edge_type) {

    case (LAT_CIRCLE):

      get_edge_middle_point_lat(a, b, middle);
      break;

    case (LON_CIRCLE):
    case (GREAT_CIRCLE):

      get_edge_middle_point_gc(a, b, middle);
      break;

    default:
      yac_internal_abort_message("ERROR: invalid edge type\n", __FILE__, __LINE__);
      middle[0] = -1; // program should never reach this point
      middle[1] = -1;
      middle[2] = -1;
  };
}

/**
 * cell clipping using Sutherland-Hodgman algorithm;
 */
static void point_list_clipping (struct point_list * source_list, int source_ordering,
                                 struct point_list target_list, int target_ordering,
                                 unsigned nct, double * tgt_edge_norm_vec) {

  struct {
    double * edge_norm_vec;
    struct point_list_element * point;
  } tgt_points[nct];

  // to avoid some problems that can occur close the the pole, we process the
  // target lat-circle edges at the end
  {

    int count = 0;

    struct point_list_element * curr_tgt_point = target_list.first;

    for (int i = 0; i < nct; ++i, curr_tgt_point = curr_tgt_point->next) {

      if (curr_tgt_point->edge_type == LAT_CIRCLE) continue;

      tgt_points[count].edge_norm_vec = tgt_edge_norm_vec + i * 3;
      tgt_points[count++].point = curr_tgt_point;
    }

    if (count != nct) {
      for (int i = 0; i < nct; ++i, curr_tgt_point = curr_tgt_point->next) {

        if (curr_tgt_point->edge_type != LAT_CIRCLE) continue;

        tgt_points[count].edge_norm_vec = tgt_edge_norm_vec + i * 3;
        tgt_points[count++].point = curr_tgt_point;
      }
    }
  }

  for (int i = 0; i < nct; ++i) {

    struct point_list_element * curr_src_point = source_list->first;
    struct point_list_element * prev_src_point = source_list->last;

    unsigned prev_is_inside, curr_is_inside;

    prev_is_inside = is_inside(prev_src_point->vec_coords,
                               tgt_points[i].edge_norm_vec,
                               tgt_points[i].point->edge_type, target_ordering);

    // for all edges of the target cell
    do {

      curr_is_inside = is_inside(curr_src_point->vec_coords,
                                 tgt_points[i].edge_norm_vec,
                                 tgt_points[i].point->edge_type,
                                 target_ordering);

      double p[3], q[3];
      int intersect = -1;

      if ((curr_is_inside + prev_is_inside == 1) ||
          (((prev_src_point->edge_type == LAT_CIRCLE) ^
            (tgt_points[i].point->edge_type == LAT_CIRCLE)) &&
           (prev_is_inside + curr_is_inside < 4))) {

        // get intersection points
        intersect = yac_intersect_vec(prev_src_point->edge_type,
                                      prev_src_point->vec_coords,
                                      curr_src_point->vec_coords,
                                      tgt_points[i].point->edge_type,
                                      tgt_points[i].point->vec_coords,
                                      tgt_points[i].point->next->vec_coords,
                                      p, q);

        // if both edges are on an identical great circle
        if ((intersect != -1) && (intersect & (1 << 4))) {

          prev_is_inside = 2;
          curr_is_inside = 2;
        }
      }

      // if the current edges change from inside/outside to outside/inside
      if (curr_is_inside + prev_is_inside == 1) {

        // if there is an no intersection
        if ((intersect == -1) || ((intersect & 3) == 0)) {

          if (intersect == -1) intersect = 0;

          intersect |= 1;

          // which of the two end points of the source edge is closer to the
          // target edge

          struct point_list_element * temp_src_point;

          switch (tgt_points[i].point->edge_type) {

            case (LON_CIRCLE) :
            case (GREAT_CIRCLE) :
              temp_src_point =
                (fabs(dotproduct(prev_src_point->vec_coords,
                                 tgt_points[i].edge_norm_vec)) <
                 fabs(dotproduct(curr_src_point->vec_coords,
                                 tgt_points[i].edge_norm_vec)))?
                prev_src_point:curr_src_point;
              break;
            case (LAT_CIRCLE) :
              temp_src_point =
                (fabs(1.0 - fabs(prev_src_point->vec_coords[2] +
                                 tgt_points[i].edge_norm_vec[2])) <
                 fabs(1.0 - fabs(curr_src_point->vec_coords[2] +
                                 tgt_points[i].edge_norm_vec[2])))?
                prev_src_point:curr_src_point;
              break;
            default:
              yac_internal_abort_message("invalid edge type\n", __FILE__, __LINE__);
              return;
          };

          p[0] = temp_src_point->vec_coords[0];
          p[1] = temp_src_point->vec_coords[1];
          p[2] = temp_src_point->vec_coords[2];
        }

        // if there are two intersection points with the source edge
        if ((intersect & ((1 << 0) | (1 << 1))) == ((1 << 0) | (1 << 1))) {

          if (!((prev_src_point->edge_type == LAT_CIRCLE) ^
                (tgt_points[i].point->edge_type == LAT_CIRCLE))) {

            yac_internal_abort_message("ERROR: ...this should not have happened...\n",
                                       __FILE__, __LINE__);
          }

          if (((get_vector_angle(prev_src_point->vec_coords, p) < yac_angle_tol) &&
               (get_vector_angle(curr_src_point->vec_coords, q) < yac_angle_tol)) ||
              ((get_vector_angle(prev_src_point->vec_coords, q) < yac_angle_tol) &&
               (get_vector_angle(curr_src_point->vec_coords, p) < yac_angle_tol))) {

            prev_is_inside = 2;
            curr_is_inside = 2;

          } else {

            // which of the two end points of the source edge is closer to the
            // target edge

            int prev_is_closer;

            switch (tgt_points[i].point->edge_type) {

              case (LON_CIRCLE) :
              case (GREAT_CIRCLE) :
                prev_is_closer =
                  fabs(dotproduct(prev_src_point->vec_coords,
                                  tgt_points[i].edge_norm_vec)) <
                  fabs(dotproduct(curr_src_point->vec_coords,
                                  tgt_points[i].edge_norm_vec));
                break;
              case (LAT_CIRCLE) :
                prev_is_closer =
                  fabs(1.0 - fabs(prev_src_point->vec_coords[2] +
                                  tgt_points[i].edge_norm_vec[2])) <
                  fabs(1.0 - fabs(curr_src_point->vec_coords[2] +
                                  tgt_points[i].edge_norm_vec[2]));
                break;
              default:
                yac_internal_abort_message("invalid edge type\n", __FILE__, __LINE__);
                return;
            };

            if (prev_is_closer)
              prev_is_inside = curr_is_inside;
            else
              curr_is_inside = prev_is_inside;
          }

        // if p or q is on the source edge
        } else {

          struct point_list_element * intersect_point;

          // if the previous point was inside or current edge is the last one
          if (prev_is_inside ||
              (curr_is_inside && (prev_src_point == source_list->last))) {

            intersect_point = get_free_point_list_element(source_list);
            prev_src_point->next = intersect_point;
            intersect_point->next = curr_src_point;

            if (prev_src_point == source_list->last)
              source_list->last = intersect_point;

          } else
            intersect_point = prev_src_point;

          if (prev_is_inside)
            intersect_point->edge_type = tgt_points[i].point->edge_type;
          else
            intersect_point->edge_type = prev_src_point->edge_type;

          if (intersect & (1 << 0)) {

            intersect_point->vec_coords[0] = p[0];
            intersect_point->vec_coords[1] = p[1];
            intersect_point->vec_coords[2] = p[2];

          // if q is on the source edge
          } else if (intersect & (1 << 1)) {

            intersect_point->vec_coords[0] = q[0];
            intersect_point->vec_coords[1] = q[1];
            intersect_point->vec_coords[2] = q[2];

          }

          if (intersect_point == prev_src_point)
            prev_is_inside = 1;
        }
      }

      // if the one edge is a circle of latitude while the other is not
      // and both corners are not directly on the edge
      if (((prev_src_point->edge_type == LAT_CIRCLE) ^
           (tgt_points[i].point->edge_type == LAT_CIRCLE)) &&
          ((prev_is_inside + curr_is_inside == 0) ||
           (prev_is_inside + curr_is_inside == 2) ||
           (prev_is_inside + curr_is_inside == 3))) {

        // if there is an intersection possible
        if ((intersect != -1) && (intersect != 0)) {

          // if there are two intersection points with the source edge
          if ((intersect & ((1 << 0) | (1 << 1))) == ((1 << 0) | (1 << 1))) {

            struct point_list_element * intersect_points[2];

            // if the previous point was inside or current edge is the last one
            if ((prev_is_inside || prev_src_point == source_list->last) &&
                (prev_is_inside != 2)) {

              intersect_points[0] = get_free_point_list_element(source_list);
              prev_src_point->next = intersect_points[0];
              intersect_points[0]->next = curr_src_point;

              if (prev_src_point == source_list->last)
                source_list->last = intersect_points[0];

            } else {
              intersect_points[0] = prev_src_point;
              intersect_points[0]->to_be_removed = 0;
            }

            // second intersection point
            intersect_points[1] = get_free_point_list_element(source_list);

            if (intersect_points[0] == source_list->last)
              source_list->first = intersect_points[1];

            intersect_points[1]->next = intersect_points[0]->next;
            intersect_points[0]->next = intersect_points[1];

            int p_is_first = get_vector_angle(prev_src_point->vec_coords, p) <
                             get_vector_angle(prev_src_point->vec_coords, q);
            enum yac_edge_type prev_src_point_edge_type =
              prev_src_point->edge_type;

            intersect_points[!p_is_first]->vec_coords[0] = p[0];
            intersect_points[!p_is_first]->vec_coords[1] = p[1];
            intersect_points[!p_is_first]->vec_coords[2] = p[2];
            intersect_points[p_is_first]->vec_coords[0] = q[0];
            intersect_points[p_is_first]->vec_coords[1] = q[1];
            intersect_points[p_is_first]->vec_coords[2] = q[2];
            intersect_points[(prev_is_inside != 0) &&
                             (curr_is_inside != 0)]->edge_type =
              prev_src_point->edge_type;
            intersect_points[(prev_is_inside == 0) ||
                             (curr_is_inside == 0)]->edge_type =
              tgt_points[i].point->edge_type;

            if (prev_is_inside == 2 || curr_is_inside == 2) {

              int tgt_edge_inside_src;

              {
                double edge_middle[3];

                get_edge_middle_point(p, q, tgt_points[i].point->edge_type,
                                      edge_middle);

                double norm_vec[3];

                switch (prev_src_point_edge_type) {

                  case (LON_CIRCLE) :
                  case (GREAT_CIRCLE) :
                    compute_norm_vector(prev_src_point->vec_coords,
                                        curr_src_point->vec_coords,
                                        norm_vec);
                    break;
                  case (LAT_CIRCLE):
                    compute_lat_circle_z_value(prev_src_point->vec_coords,
                                               curr_src_point->vec_coords,
                                               norm_vec);
                    break;
                  default:
                    norm_vec[0] = 0.0, norm_vec[1] = 0.0, norm_vec[2] = 0.0;
                    yac_internal_abort_message("invalid edge type\n", __FILE__, __LINE__);
                };

                tgt_edge_inside_src = is_inside(edge_middle, norm_vec,
                                                prev_src_point_edge_type,
                                                source_ordering) != 1;
              }

              // if one source point is on the target edge and the other is inside
              if (prev_is_inside + curr_is_inside == 3) {

                // if the current source point is on the target edge, then the
                // second intersection point is just a dummy that is identical
                // to the current source point but might have the wrong edge
                // type
                if (curr_is_inside == 2)
                  intersect_points[1]->to_be_removed = 1;

                if (curr_is_inside == 2)
                  curr_src_point->to_be_removed = tgt_edge_inside_src;
                if (prev_is_inside == 2)
                  prev_src_point->to_be_removed = tgt_edge_inside_src;
              }
              if (prev_is_inside == 0 && curr_is_inside == 2) {
                if (tgt_edge_inside_src)
                  intersect_points[1]->to_be_removed = 1;
                else
                  curr_src_point->to_be_removed = 1;
              }
            }

            // if the previous point has been reused for an intersection
            if (intersect_points[0] == prev_src_point)
              prev_is_inside = 1;
          }
        }

      // if the one edge is a circle of latitude while the other is not
      // and both corners are directly on the edge
      } else if (((prev_src_point->edge_type == LAT_CIRCLE) ^
                  (tgt_points[i].point->edge_type == LAT_CIRCLE)) &&
                 (prev_is_inside == 2) && (curr_is_inside == 2)) {

        double cross_src_z, cross_tgt_z;

        cross_src_z = (long double)prev_src_point->vec_coords[0] *
                      (long double)curr_src_point->vec_coords[1] -
                      (long double)prev_src_point->vec_coords[1] *
                      (long double)curr_src_point->vec_coords[0];
        cross_tgt_z = (long double)tgt_points[i].point->vec_coords[0] *
                      (long double)tgt_points[i].point->next->vec_coords[1] -
                      (long double)tgt_points[i].point->vec_coords[1] *
                      (long double)tgt_points[i].point->next->vec_coords[0];

        int same_ordering = source_ordering == target_ordering;
        int same_direction = (cross_src_z > 0) == (cross_tgt_z > 0);

        // if source and target cell have the same ordering and both
        // edges have the same direction or if both cells have different
        // ordering and the edges have different directions, then we might
        // have to change the edge type of the source edge
        if (same_ordering == same_direction) {

          // well...it works...do not ask  ;-)
          // ((edge is on south hemisphere) XOR (direction of source edge) XOR
          //  (ordering of source cell))
          if ((curr_src_point->vec_coords[2] > 0) ^
              (cross_src_z < 0) ^ source_ordering)
            prev_src_point->edge_type = LAT_CIRCLE;
          else
            prev_src_point->edge_type = GREAT_CIRCLE;
        }
      }

      // if the previous points was on the target edge and the current
      // one is outside
      if (prev_is_inside == 2 && curr_is_inside == 0)
        prev_src_point->edge_type = tgt_points[i].point->edge_type;

      if (!prev_is_inside)
        prev_src_point->to_be_removed = 1;

      prev_src_point = curr_src_point;
      curr_src_point = curr_src_point->next;
      prev_is_inside = curr_is_inside;

    } while ((prev_src_point != source_list->last) &&
             (source_list->first != NULL));

    // remove all points that are to be deleted
    remove_points(source_list);

    // if there are no more corners in the source cell
    if (source_list->first == NULL) break;
  }
}

static void copy_point_list(struct point_list in, struct point_list * out) {

  reset_point_list(out);

  struct point_list_element * curr = in.first;

  if (curr == NULL) return;

  struct point_list_element * new = get_free_point_list_element(out);
  out->first = new;
  *new = *curr;
  curr = curr->next;

  do {

    new->next = get_free_point_list_element(out);
    new = new->next;
    *new = *curr;
    curr = curr->next;

  } while (curr != in.first);

  new->next = out->first;
  out->last = new;
}

void yac_cell_clipping (unsigned N,
                        struct grid_cell * source_cell,
                        struct grid_cell target_cell,
                        struct grid_cell * overlap_buffer) {

  unsigned ncs;               /* number of vertices of source cell */
  unsigned nct;               /* number of vertices of target cell */

  struct point_list target_list, source_list, temp_list;

  unsigned target_ordering; /* ordering of target cell corners */
  unsigned source_ordering; /* ordering of source cell corners */

  double * norm_vec; /* norm vector for temporary target edge plane */

  enum cell_type tgt_cell_type = get_cell_type(target_cell);

  if (tgt_cell_type == MIXED_CELL)
    yac_internal_abort_message("invalid target cell type (cell contains edges consisting "
                               "of great circles and circles of latitude)\n",
			       __FILE__, __LINE__);

  init_point_list(&temp_list);

  // generate point list for target cell (clip cell)
  init_point_list(&target_list);
  generate_point_list(&target_list, target_cell);

  // if there is no target cell (e.g. if all edges of target cell have a length
  // of zero)
  if (target_list.first == NULL) {
    free_point_list(&target_list);
    for (unsigned i = 0; i < N; ++i)
      overlap_buffer[i].num_corners = 0;
    return;
  }

  struct point_list_element * prev_tgt_point = target_list.first;
  struct point_list_element * curr_tgt_point = target_list.first->next;

  nct = 0;
  do {
    nct++;
    prev_tgt_point = prev_tgt_point->next;
  } while (prev_tgt_point != target_list.first);

  norm_vec = malloc(3 * nct * sizeof(*norm_vec));

  // compute norm vectors for all edges
  // or for lat circle edges a special z value
  for (unsigned i = 0; i < nct; ++i) {

    switch (prev_tgt_point->edge_type) {

      case (LON_CIRCLE) :
      case (GREAT_CIRCLE) :
        compute_norm_vector(prev_tgt_point->vec_coords, curr_tgt_point->vec_coords,
                            norm_vec + 3 * i);
        break;
      case (LAT_CIRCLE):
        compute_lat_circle_z_value(prev_tgt_point->vec_coords, curr_tgt_point->vec_coords,
                            norm_vec + 3 * i);
        break;
      default:
        yac_internal_abort_message("invalid edge type\n", __FILE__, __LINE__);
    };
    prev_tgt_point = curr_tgt_point;
    curr_tgt_point = curr_tgt_point->next;
  }

  // compute target direction
  target_ordering = get_cell_points_ordering(&target_list);

  init_point_list(&source_list);

  // for all source cells
  for (unsigned n = 0; n < N; n++ ) {

    overlap_buffer[n].num_corners = 0;

    enum cell_type src_cell_type = get_cell_type(source_cell[n]);

    if (src_cell_type == MIXED_CELL)
      yac_internal_abort_message("invalid source cell type (cell contains edges consisting "
                                 "of great circles and circles of latitude)\n",
                                 __FILE__, __LINE__);

    if (source_cell[n].num_corners < 2)
      continue;

    // generate point list for current source list
    generate_point_list(&source_list, source_cell[n]);

    {
      struct point_list_element * src_point = source_list.first;

      ncs = 0;
      do {
        ncs++;
        src_point = src_point->next;
      } while (src_point != source_list.first);
    }

    // compute source direction
    source_ordering = get_cell_points_ordering(&source_list);

    struct point_list * overlap;

    // in this case the source an target cell are both LAT_CELL's than the
    // bigger one has to be the target cell
    // a similar problem occurs when the target cell is a LAT_CELL and the
    // source is a GREAT_CIRCLE_CELL which overlaps with the pole that is also
    // include in the target cell
    if (((tgt_cell_type == LAT_CELL) && (src_cell_type == GREAT_CIRCLE_CELL)) ||
        ((tgt_cell_type == LAT_CELL) && (src_cell_type == LAT_CELL) &&
         (fabs(target_cell.coordinates_y[0]) >
          fabs(source_cell[n].coordinates_y[0])))) {

      copy_point_list(target_list, &temp_list);

      double temp_norm_vec[3*ncs];
      struct point_list_element * src_point = source_list.first;


      for (unsigned i = 0; i < ncs; ++i) {
        switch (src_point->edge_type) {

          case (LON_CIRCLE) :
          case (GREAT_CIRCLE) :
            compute_norm_vector(src_point->vec_coords,
                                src_point->next->vec_coords,
                                temp_norm_vec + 3 * i);
            break;
          case (LAT_CIRCLE):
            compute_lat_circle_z_value(src_point->vec_coords,
                                       src_point->next->vec_coords,
                                       temp_norm_vec + 3 * i);
            break;
          default:
            yac_internal_abort_message("invalid edge type\n", __FILE__, __LINE__);
        };
        src_point = src_point->next;
      }

      point_list_clipping(&temp_list, target_ordering,
                          source_list, source_ordering, ncs, temp_norm_vec);

      overlap = &temp_list;

    } else {

      point_list_clipping(&source_list, source_ordering,
                          target_list, target_ordering, nct, norm_vec);

      overlap = &source_list;
    }

    if (overlap->first != NULL)
      generate_overlap_cell(overlap, overlap_buffer + n);
  }

  free(norm_vec);
  free_point_list(&source_list);
  free_point_list(&target_list);
  free_point_list(&temp_list);
}

/* ---------------------------------------------------- */

void yac_correct_weights ( unsigned nSourceCells, double * weight ) {

  static unsigned maxIter = 10; // number of iterations to get better accuracy of the weights
  static double const tol = 1.0e-15;

  unsigned n;
  unsigned iter;

  double weight_diff;

  for ( iter = 1; iter < maxIter; iter++ ) {

    weight_diff = 1.0;

    for ( n = 0; n < nSourceCells; n++ )
      weight_diff -= weight[n];

#ifdef VERBOSE
    printf ("weight sum is %.15f \n", weight_sum);
    printf ("weights are  ");
    for (unsigned i = 0; i < nSourceCells; ++i)
      printf (" %.15f", weight[i]);
    printf("\n");
#endif

    if ( fabs(weight_diff) < tol ) break;

    for ( n = 0; n < nSourceCells; n++ )
      weight[n] += weight[n] * weight_diff;
  }
#ifdef VERBOSE
  if ( fabs(weight_diff) > tol )
    printf ("weight sum is %.15f \n", weight_sum);
#endif
}

/* ---------------------------------------------------- */

static unsigned get_cell_points_ordering(struct point_list * cell) {

  if ((cell->first == NULL) || (cell->first == cell->last))
    yac_internal_abort_message("ERROR: invalid cell\n", __FILE__, __LINE__);

  double norm_vec[3];
  struct point_list_element * curr = cell->first;

  compute_norm_vector(curr->vec_coords, curr->next->vec_coords, norm_vec);

  curr = curr->next;

  if (curr->next == cell->first)
    yac_internal_abort_message("ERROR: invalid cell\n", __FILE__, __LINE__);

  do {
    curr = curr->next;

    double dot = dotproduct(curr->vec_coords, norm_vec);

    if (fabs(dot) > tol)
      return dot > 0;

  } while (curr != cell->first);

  yac_internal_abort_message("ERROR: could not determine order of points in cell\n",
                             __FILE__, __LINE__);

  return -1;
}

static void init_point_list(struct point_list * list) {

  list->first = NULL;
  list->last = NULL;
  list->free_elements = NULL;
}

static void reset_point_list(struct point_list * list) {

  if (list->first != NULL) {

    list->last->next = list->free_elements;
    list->free_elements = list->first;

    list->first = NULL;
    list->last = NULL;
  } else
    list->last = NULL;
}

static void remove_points(struct point_list * list) {

  struct point_list_element * curr = list->first;
  struct point_list_element * prev;
  struct point_list_element * remaining_points = NULL;

  if (curr == NULL) return;

  unsigned num_edges = 1;

  while (curr != list->last) {

    curr = curr->next;
    num_edges++;
  }

  for (unsigned i = 0; i < num_edges; ++i) {

    prev = curr;
    curr = curr->next;

    if (curr->to_be_removed) {
      prev->next = curr->next;
      curr->next = list->free_elements;
      list->free_elements = curr;
      curr = prev;
    } else
      remaining_points = curr;
  }

  if ((remaining_points != NULL) &&
      (remaining_points == remaining_points->next)) {

    remaining_points->next = list->free_elements;
    list->free_elements = remaining_points;
    remaining_points = NULL;
  }

  list->last = remaining_points;
  if (remaining_points != NULL)
    list->first = remaining_points->next;
  else
    list->first = NULL;
}

//! returns number of edges/corners
static unsigned remove_zero_length_edges(struct point_list * list) {

  struct point_list_element * curr = list->first;
  double const tol = 1e-10;

  if (curr == NULL) return 0;

  unsigned num_edges = 1;
  unsigned temp_num_edges;

  while (curr != list->last) {

    curr = curr->next;
    num_edges++;
  }

  temp_num_edges = num_edges;

  for (unsigned i = 0; i < num_edges; ++i) {

    // if both points are nearly identical (angle between them is very small)
    if (!curr->to_be_removed &&
        (get_vector_angle(curr->vec_coords, curr->next->vec_coords) < tol)) {
      curr->to_be_removed = 1;
      temp_num_edges--;
    } else if (curr->edge_type == LAT_CIRCLE &&
               curr->next->edge_type == LAT_CIRCLE &&
               (fabs(get_angle(atan2(curr->vec_coords[1] ,
                                     curr->vec_coords[0]),
                               atan2(curr->next->next->vec_coords[1] ,
                                     curr->next->next->vec_coords[0]))) <
                M_PI_2)) {
      curr->next->to_be_removed = 1;
      temp_num_edges--;
    }

    curr = curr->next;
  }

  remove_points(list);

  return temp_num_edges;
}

static void generate_point_list(struct point_list * list,
                                struct grid_cell cell) {

  reset_point_list(list);

  if (cell.num_corners == 0) return;

  struct point_list_element * curr = get_free_point_list_element(list);

  list->first = curr;
  curr->vec_coords[0] = cell.coordinates_xyz[0+0*3];
  curr->vec_coords[1] = cell.coordinates_xyz[1+0*3];
  curr->vec_coords[2] = cell.coordinates_xyz[2+0*3];

  for (unsigned i = 1; i < cell.num_corners; ++i) {

    curr->next = get_free_point_list_element(list);
    curr->edge_type = cell.edge_type[i - 1];
    curr = curr->next;

    curr->vec_coords[0] = cell.coordinates_xyz[0+i*3];
    curr->vec_coords[1] = cell.coordinates_xyz[1+i*3];
    curr->vec_coords[2] = cell.coordinates_xyz[2+i*3];
    curr->edge_type = cell.edge_type[i];
  }

  curr->next = list->first;
  list->last = curr;

  remove_zero_length_edges(list);
}

static struct point_list_element *
get_free_point_list_element(struct point_list * list) {

  struct point_list_element * element;

  if (list->free_elements == NULL) {

    for (int i = 0; i < 7; ++i) {

      element = malloc(1 * sizeof(*element));

      element->next = list->free_elements;
      list->free_elements = element;
    }

    element = malloc(1 * sizeof(*element));

  } else {

    element = list->free_elements;
    list->free_elements = list->free_elements->next;
  }

  element->next = NULL;
  element->to_be_removed = 0;

  return element;
}

static void free_point_list(struct point_list * list) {

  struct point_list_element * element;

  if (list->first != NULL) {

    list->last->next = NULL;

    while (list->first != NULL) {

      element = list->first;
      list->first = element->next;
      free(element);
    }
  }

  while (list->free_elements != NULL) {

    element = list->free_elements;
    list->free_elements = element->next;
    free(element);
  }

  list->first = NULL;
  list->last = NULL;
  list->free_elements = NULL;
}

static int is_empty_gc_cell(struct point_list * list, unsigned num_edges) {

  double const tol = 1e-6;

  if ((num_edges == 2) &&
      (list->first->edge_type != LAT_CIRCLE) &&
      (list->last->edge_type != LAT_CIRCLE))
    return 1;

  struct point_list_element * curr = list->first;

  for (unsigned i = 0; i < num_edges; ++i) {

    if (curr->edge_type == LAT_CIRCLE) return 0;
    curr = curr->next;
  }

  double ref_norm[3];

  compute_norm_vector(curr->vec_coords, curr->next->vec_coords, ref_norm);

  for (unsigned i = 0; i < num_edges-1; ++i) {

    curr = curr->next;

    double norm[3];

    compute_norm_vector(curr->vec_coords, curr->next->vec_coords, norm);

    if (((fabs(ref_norm[0] - norm[0]) > tol) ||
         (fabs(ref_norm[1] - norm[1]) > tol) ||
         (fabs(ref_norm[2] - norm[2]) > tol)) &&
        ((fabs(ref_norm[0] + norm[0]) > tol) ||
         (fabs(ref_norm[1] + norm[1]) > tol) ||
         (fabs(ref_norm[2] + norm[2]) > tol)))
      return 0;
  }

  return 1;
}

static void generate_overlap_cell(struct point_list * list,
                                  struct grid_cell * cell) {

  //! \todo test whether all points of the cell are on a single
  //!       great circle --> empty cell

  unsigned num_edges = remove_zero_length_edges(list);

  if ((num_edges < 2) ||
      is_empty_gc_cell(list, num_edges)){

    reset_point_list(list);
    cell->num_corners = 0;
    return;
  }

  if (num_edges > cell->array_size) {
    free(cell->coordinates_x);
    free(cell->coordinates_y);
    free(cell->coordinates_xyz);
    free(cell->edge_type);
    cell->coordinates_x = malloc(num_edges * sizeof(*cell->coordinates_x));
    cell->coordinates_y = malloc(num_edges * sizeof(*cell->coordinates_y));
    cell->coordinates_xyz = malloc(3 * num_edges * sizeof(*cell->coordinates_xyz));
    cell->edge_type = malloc(num_edges * sizeof(*cell->edge_type));
    cell->array_size = num_edges;
  }
  cell->num_corners = num_edges;

  struct point_list_element * curr = list->first;

  for (unsigned i = 0; i < num_edges; ++i) {


    XYZtoLL(curr->vec_coords, cell->coordinates_x+i, cell->coordinates_y+i);
    cell->coordinates_xyz[0+i*3] = curr->vec_coords[0];
    cell->coordinates_xyz[1+i*3] = curr->vec_coords[1];
    cell->coordinates_xyz[2+i*3] = curr->vec_coords[2];
    cell->edge_type[i] = curr->edge_type;

    curr = curr->next;
  }
}
