/**
 * @file grid.h
 * @brief Structs and interfaces to handle grids
 *
 * struct grid is used as an abstraction for all types of grids. It contains longitude
 * and latitude data of the corners and the mapping between cells and their corners and
 * edges. \n
 * In addition to the struct grid definition this file contains the declaration of a
 * set of routines that are valid for all grid types and do certain operations on grids.
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

#ifndef GRID_H
#define GRID_H

/** \example test_grid.c
 * This shows how to work with struct grid.
 */

#include "dep_list.h"
#include "math.h"
#include "grid_cell.h"

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif

#ifndef  M_PI_2
#define  M_PI_2      1.57079632679489661923132169163975144  /* pi/2 */
#endif

#define EARTH_RADIUS 6371.2290

static double const EarthRadius  = EARTH_RADIUS;
static double const rad          = M_PI / 180.0;
static double const deg          = 180.0 / M_PI;
static double const EarthRadius2 = EARTH_RADIUS * EARTH_RADIUS / 2.0 ;

// forward declaration required by grid_vtable
struct grid;
struct bounding_circle;

struct grid_vtable {

   struct grid * (*copy)(struct grid *);
   void (*get_2d_extent)(struct grid *, double (* extent)[2]);
   void (*get_grid_cell) (struct grid *, unsigned, struct grid_cell *);
   void (*get_grid_cell2) (struct grid *, unsigned, struct grid_cell *,
                           struct bounding_circle * bnd_circle);
   unsigned (*get_size_x_coords)(struct grid *);
   unsigned (*get_size_y_coords)(struct grid *);
   double const * (*get_x_coords)(struct grid *);
   double const * (*get_y_coords)(struct grid *);
   void (*set_x_coords)(struct grid *, double *);
   void (*set_y_coords)(struct grid *, double *);
   unsigned (*get_size_cell_grid_x_coords)(struct grid *);
   unsigned (*get_size_cell_grid_y_coords)(struct grid *);
   unsigned (*get_num_grid_cells)(struct grid *);
   unsigned (*get_num_grid_corners)(struct grid *);
   unsigned (*get_num_cell_corners)(struct grid *, unsigned);
   unsigned (*get_num_corner_cells)(struct grid *, unsigned);
   unsigned (*get_num_grid_edges)(struct grid *);
   unsigned (*get_num_corner_edges)(struct grid *, unsigned);
   unsigned (*get_num_cell_edges) (struct grid *, unsigned);
   unsigned const * (*get_corner_edges)(struct grid *, unsigned);
   unsigned const * (*get_cell_edge_indices)(struct grid *, unsigned);
   enum yac_edge_type (*get_edge_type)(struct grid *, unsigned);
   unsigned const * (*get_cell_corner_indices)(struct grid *, unsigned);
   unsigned const * (*get_corner_cell_indices)(struct grid *, unsigned);
   unsigned const * (*get_cell_x_coord_indices) (struct grid *, unsigned);
   unsigned const * (*get_cell_y_coord_indices) (struct grid *, unsigned);
   double (*get_corner_x_coord) (struct grid *, unsigned);
   double (*get_corner_y_coord) (struct grid *, unsigned);
   unsigned (*get_corner_x_coord_index) (struct grid *, unsigned);
   unsigned (*get_corner_y_coord_index) (struct grid *, unsigned);
   int (*get_aux_grid_cell)(struct grid *, unsigned, unsigned *, enum yac_edge_type *);
   struct dep_list (*get_cell_neigh_dep_list)(struct grid *);
   void (*get_boundary_corners) (struct grid *, unsigned *, unsigned *);
   struct grid * (*generate_cell_grid)(struct grid *, double *, double *);
   struct grid * (*generate_subgrid)(struct grid *, unsigned *, unsigned,
                                     unsigned **, unsigned **, unsigned **);
   void (*pack_grid)(struct grid *, double **, unsigned, unsigned *, unsigned *,
                     unsigned **, unsigned, unsigned *, unsigned *);
   struct grid_search * (*get_grid_search)(struct grid * grid);
   void (*delete)(struct grid *);
};

struct grid {

   struct grid_vtable *vtable;
};

/**
 * computes the 2d %grid extent of a given %grid
 * @param[in] grid is a pointer to a declared grid object
 * @param[out] extent is the computed 2d extent
 */
void yac_get_2d_grid_extent(struct grid * grid, double (* extent)[2]);

/**
 * gets a %grid cell from a grid object
 * @param[in] grid
 * @param[in] cell_index local cell index of the requested cell
 * @param[out] cell requested cell (grid_cell object has be initialised once before)
 * @see init_grid_cell
 */
void yac_get_grid_cell (struct grid * grid, unsigned cell_index,
                        struct grid_cell * cell);

/**
 * gets a %grid cell from a grid object\n
 * in addition the bounding circle for this grid cell is also returned
 * @param[in] grid
 * @param[in] cell_index local cell index of the requested cell
 * @param[out] cell requested cell (grid_cell object has be initialised once before)
 * @param[out] bnd_circle bounding circle of input cell
 * @see init_grid_cell
 */
void yac_get_grid_cell2 (struct grid * grid, unsigned cell_index,
                         struct grid_cell * cell,
                         struct bounding_circle * bnd_circle);

/**
 * @param[in] grid
 * @return size of the array returned by \ref yac_get_x_coords
 */
unsigned yac_get_size_x_coords(struct grid * grid);

/**
 * @param[in] grid
 * @return size of the array returned by \ref yac_get_y_coords
 */
unsigned yac_get_size_y_coords(struct grid * grid);

/**
 * returns an array that contains the x coordinates of all corners of the grid
 * @param[in] grid
 * @return x coordinates of all corners of the grid
 */
double const * yac_get_x_coords(struct grid * grid);

/**
 * returns an array that contains the y coordinates of all corners of the grid
 * @param[in] grid
 * @return y coordinates of all corners of the grid
 */
double const * yac_get_y_coords(struct grid * grid);

/**
 * sets the array that contains the x coordinates of all corners of the grid
 * @param[in] grid
 * @param[in] x_coords
 */
void yac_set_x_coords(struct grid * grid, double * x_coords);

/**
 * set the array that contains the y coordinates of all corners of the grid
 * @param[in] grid
 * @param[in] y_coords
 */
void yac_set_y_coords(struct grid * grid, double * y_coords);

/**
 * gets the size of the x coordinate array of a cell %grid generated by \ref yac_generate_cell_grid for the given %grid object
 * @param[in] grid
 * @return size of x coordinate for a potential cell grid
 * @see generate_cell_grid
 */
unsigned yac_get_size_cell_grid_x_coords(struct grid * grid);

/**
 * gets the size of the y coordinate array of a cell %grid generated by \ref yac_generate_cell_grid for the given %grid object
 * @param[in] grid
 * @return size of y coordinate array for a potential cell grid
 * @see generate_cell_grid
 */
unsigned yac_get_size_cell_grid_y_coords(struct grid * grid);

/**
 * gets the number of cells in the given %grid
 * @param[in] grid
 * @return number of cells in the %grid
 */
unsigned yac_get_num_grid_cells (struct grid * grid);

/**
 * gets the number of corners in the given %grid
 * @param[in] grid
 * @return number of corners in the %grid
 */
unsigned yac_get_num_grid_corners (struct grid * grid); 

/**
 * gets the number of corners for a cell of a given %grid
 * @param[in] grid
 * @param[in] cell_index local id of the respective cell
 * @return number of corners for the specified cell in %grid
 */
unsigned yac_get_num_cell_corners (struct grid * grid, unsigned cell_index);

/**
 * gets the number of cells associated to a corner of a given %grid
 * @param[in] grid
 * @param[in] corner_index local id of the respective corner
 * @return number of cells for the specified corner in %grid
 */
unsigned yac_get_num_corner_cells (struct grid * grid, unsigned corner_index);

/**
 * gets the number of edges in the given %grid
 * @param[in] grid
 * @return number of edges in %grid
 */
unsigned yac_get_num_grid_edges (struct grid * grid);

/**
 * gets the number of edges that are associated with a corner of a given %grid
 * this is identical to the number of corners directly connected to the given corner
 * @param[in] grid
 * @param[in] corner_index local index of the respective corner
 * @return number of edges associated with given corner
 */
unsigned yac_get_num_corner_edges (struct grid * grid, unsigned corner_index);

/**
 * gets the number of edges for a cell of a given %grid (identical to \ref yac_get_num_cell_corners)
 * @param[in] grid
 * @param[in] cell_index local id of the respective cell
 * @return number of edges for the specified cell in %grid
 * @see get_num_cell_corners
 */
unsigned yac_get_num_cell_edges (struct grid * grid, unsigned cell_index);

/**
 * gets neighbour corners of a corner of a given %grid
 * @param[in] grid
 * @param[in] corner_index local id of a corner in %grid
 * @return array containing the local ids of all corners directly linked to the given corner
 * @see get_num_corner_edges
 */
unsigned const * yac_get_corner_edges (struct grid * grid, unsigned corner_index);

/**
 * gets the %edges of a cell for a given %grid
 * @param[in] grid
 * @param[in] cell_index local id of the respective cell
 * @return array that contains the local ids of all edges for the given cell
 * @see get_num_cell_edges
 */
unsigned const * yac_get_cell_edge_indices (struct grid * grid, unsigned cell_index);

/**
 * gets the type of an %edge for a given %grid
 * @param[in] grid
 * @param[in] edge_index local id of the respective edge
 * @return type of the specified %edge
 */
enum yac_edge_type yac_get_edge_type(struct grid * grid, unsigned edge_index);

/**
 * gets the corners of a cell for a given grid
 * @param[in] grid
 * @param[in] cell_index local id of the respective cell
 * @return array that contains the local ids of all corners for the given cell
 * @see get_num_cell_corners
 */
unsigned const * yac_get_cell_corner_indices (struct grid * grid, unsigned cell_index);

/**
 * gets the cells associated with a corner for a given grid
 * @param[in] grid
 * @param[in] corner_index local id of the respective corner
 * @return array that contains the local ids of all cells for the given cell
 * @see get_num_corner_cells
 */
unsigned const * yac_get_corner_cell_indices (struct grid * grid, unsigned corner_index);


/**
 * gets indices to access x coordinate array for all corners of a cell for a given %grid
 * @param[in] grid
 * @param[in] cell_index local id of the respective cell
 * @return array containing an index for each corner for the respective cell; these indices can be used to access x coordinate array (these indices are not the local corner ids)
 * @see get_num_cell_corners
 */
unsigned const * yac_get_cell_x_coord_indices (struct grid * grid, unsigned cell_index);

/**
 * gets indices to access y coordinate array for all corners of a cell for a given %grid
 * @param[in] grid
 * @param[in] cell_index local id of the respective cell
 * @return array containing an index for each corner for the respective cell; these indices can be used to access y coordinate array (these indices are not the local corner ids)
 * @see get_num_cell_corners
 */
unsigned const * yac_get_cell_y_coord_indices (struct grid * grid, unsigned cell_index);

/**
 * gets the x coordinate for a corner of a given %grid
 * @param[in] grid
 * @param[in] corner_index local id of a corner in %grid
 * @return x coordinate of the given corner
 * @see get_num_grid_corners
 */
double yac_get_corner_x_coord (struct grid * grid, unsigned corner_index);

/**
 * gets the y coordinate for a corner of a given %grid
 * @param[in] grid
 * @param[in] corner_index local id of a corner in %grid
 * @return y coordinate of the given corner
 * @see get_num_grid_corners
 */
double yac_get_corner_y_coord (struct grid * grid, unsigned corner_index);

/**
 * gets an index that can used to access x coordinate array for a corner of a given %grid
 * @param[in] grid
 * @param[in] corner_index local id of a corner in %grid
 * @return index that can be used to get the actual coordinate from x coordinate array (this is not the local corner id)
 * @see get_num_grid_corners
 */
unsigned yac_get_corner_x_coord_index (struct grid * grid, unsigned corner_index);

/**
 * gets an index that can used to access y coordinate array for a corner of a given %grid
 * @param[in] grid
 * @param[in] corner_index local id of a corner in %grid
 * @return index that can be used to get the actual coordinate from y coordinate array (this is not the local corner id)
 * @see get_num_grid_corners
 */
unsigned yac_get_corner_y_coord_index (struct grid * grid, unsigned corner_index);

/**
 * computes the cell, whose corners are the cells, which have the given corner
 * in common
 *
 * @param[in]  grid
 * @param[in]  corner_index         index of the corner around which the auxiliary cell
 *                                  is to be built
 * @param[out] cell_indices cell indices of the auxiliary cell
 * @param[out] edge_type type of the edges of the auxiliary cell
 * @return return 0 in case there is no auxiliary cell
 * @see get_num_corner_cells
 */
int yac_get_aux_grid_cell(struct grid * grid, unsigned corner_index,
                          unsigned * cell_indices, enum yac_edge_type * edge_type);

/**
 * gets a dependency list containing all cell neighbourhood information for the given %grid
 * @param[in] grid
 * @return dependency list that contains for each local cell id the local ids of all neighbour cells
 */
struct dep_list yac_get_cell_neigh_dep_list(struct grid * grid);

/**
 * makes a complete copy of a grid
 * @param[in] grid grid to be copied
 * @returns copy of the given grid
 */
struct grid * yac_copy_grid(struct grid * grid);

/** \example test_polecover.c
 * This contains examples for cell_covers_pole.
 */

/**
 * checks whether a given cell covers a pole on the globe
 * @param[in] num_corners number of corners for the given cell
 * @param[in] corners_lon longitude coordinates of the given cell
 * @param[in] corners_lat latitude coordinates of the given cell
 * @return 0 if the given cell does not cover a pole
 */
unsigned yac_cell_covers_pole (unsigned num_corners, double * const corners_lon, 
                                                     double * const corners_lat);

/**
 * get local ids of all boundary corners of the given grid
 * @param[in] grid
 * @param[in,out] bnd_corners array with the boundary corners (user must provide this array, it must be big enough to be able to hole all boundary cells)
 * @param[out] num_bnd_corners number of boundary corners written to bnd_corners
 */
void yac_get_boundary_corners (struct grid * grid, unsigned * bnd_corners,
                               unsigned * num_bnd_corners);

/**
 * generates a grid that has a corner for each cell of the input %grid
 * (the corners of the new %grid are in the same order as the cells of the input %grid)
 * @param[in] grid input %grid
 * @param[in] coordinates_x longitude data for the corners of the cell %grid (the layout of this array depends on the %grid type)
 * @param[in] coordinates_y latitude data for the corners of the cell %grid (the layout of this array depends on the %grid type)
 * @returns the generated cell grid
 */
struct grid * yac_generate_cell_grid(struct grid * grid, double * coordinates_x, 
                                     double * coordinates_y);

/**
 * generates a grid the contains the selected cells from the input grid (it may actually contain more cells)
 * @param[in] grid input %grid
 * @param[in] selected_local_cell_ids local ids of all cells that are to be included in the subgrid
 * @param[in] num_local_cells num of local ids in selected_local_cell_ids
 * @param[out] local_cell_ids pointer to an array containing the local ids of the input %grid for all cells of subgrid (the user is responsible for freeing the memory of the array)
 * @param[out] local_corner_ids pointer to an array containing the local ids of the input %grid for all corners of subgrid (the user is responsible for freeing the memory of the array)
 * @param[out] local_edge_ids pointer to an array containing the local ids of the input %grid for all edges of subgrid (the user is responsible for freeing the memory of the array)
 * @returns the generate subgrid
 */
struct grid * yac_generate_subgrid(struct grid * grid, unsigned * selected_local_cell_ids,
                                   unsigned num_local_cells, unsigned ** local_cell_ids,
                                   unsigned ** local_corner_ids, unsigned ** local_edge_ids);

/**
 * packs given grid into a double and unsigned buffer
 * @param[in] grid
 * @param[in,out] dble_buf pointer to double buffer (buffer is reallocated by this routine if necessary)
 * @param[in] dble_buf_offset number of elements already in the buffer (this data is not overwritten by this routine)
 * @param[out] dble_buf_data_size number of elements added to the buffer by this routine
 * @param[in,out] dble_buf_size size of dble_buf (this is adjusted in case this routine reallocates dble_buf)
 * @param[in,out] uint_buf pointer to double buffer (buffer is reallocated by this routine if necessary)
 * @param[in] uint_buf_offset number of elements already in the buffer (this data is not overwritten by this routine)
 * @param[out] uint_buf_data_size number of elements added to the buffer by this routine
 * @param[in,out] uint_buf_size size of uint_buf (this is adjusted in case this routine reallocates uint_buf)
 * @see \ref yac_unpack_grid
 */
void yac_pack_grid(struct grid * grid, double ** dble_buf,
                   unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                   unsigned * dble_buf_size, unsigned ** uint_buf,
                   unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                   unsigned * uint_buf_size);

/**
 * unpacks a grid from the provided double and unsigned buffers; the unpacked grid is identical
 * to the one packed into these buffers by \ref yac_pack_grid
 * @param[in] dble_buf
 * @param[out] dble_buf_data_size number of elements extracted from dble_buf
 * @param[in] uint_buf
 * @param[out] uint_buf_data_size number of elements extracted from uint_buf
 * @returns unpacked grid
 * @see pack_grid
 */
struct grid * yac_unpack_grid(double * dble_buf, unsigned * dble_buf_data_size,
                              unsigned * uint_buf, unsigned * uint_buf_data_size);

/**
 * generates a grid search object for the given grid
 * @param[in] grid
 * @return grid_search object for the provided grid
 * @remark the grid_search object returned by this routine is deleted by
 *         \ref yac_delete_grid
 */
struct grid_search * yac_get_grid_search(struct grid * grid);

/**
 * frees all memory associated with the given %grid (memory objects provided
 * by the user must be taken care of by him; e.g. coordinates_x and coordinates_y
 * provided to \ref yac_reg2d_grid_new)
 * @param[in] grid
 */
void yac_delete_grid(struct grid * grid);

#undef EARTH_RADIUS

#endif // GRID_H
